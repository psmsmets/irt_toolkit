!*****************************************************************************80
!
!                           I E R F 3 D _ R A Y T R A C E
!
!  Module:       / (include)
!
!  Programmer:   Pieter S. M. Smets
!                R&D depart. of Seismology and Acoustics - Koninklijk Nederlands
!                Meteorologisch Instituut (KNMI)
!                De Bilt, The Netherlands
!
!  Date:         August 2, 2017
!
!  Language:     Fortran-90
!
!  Description:  Include file with ierf3d raytrace procedure.
!
!
!*****************************************************************************80
!
!.....print header
!
  if (verb) print "('> ',a)", 'Ray tracing'
  if (verb) print "(4x,a,i6)"   , '#rays                  = ', rays_count
!
  azimuth_vec   = (/ ( azimuth_first  + i * azimuth_step , i=0, azimuth_count  - 1 ) /)
  elevation_vec = (/ ( elevation_first + i * elevation_step, i=0, elevation_count-1 ) /)
!
!.....initialize
!
  f0       = 0._double
  progress = 0
  cnt      = 0
!
  if ( rays_count .lt. 10 ) then
    dprogress = 50
  else if ( rays_count .lt. 20 ) then
    dprogress = 25
  else
    dprogress = 10
  end if
!
! initialize geometrical spreading
!
  call initialize_spreading ( tloss_re_unit, A0 )
!    
!.....Set stop conditions
!
  phMin      = 0.d0
  phMax      = sc%lonrange * D2R
  thMin      = ( 90.d0 - sc%latmax ) * D2R
  thMax      = ( 90.d0 - sc%latmin ) * D2R
  rMin       = RE + atmo%hmin
  rMax       = RE + atmo%hmax
  refl_hmin  = 0.3d0 / frq
  refl_hmin  = 0.3d0 / frq
  if ( refl_hmin .lt. .5d0 ) refl_hmin = .5d0
!     
!.....Source parameters, from lonlat to spherical
!
  if ( ltopo ) then
    call topo_get ( lon0, lat0, h, hdLon, hdLat )
    if ( h.gt.h0 ) h0=h
  end if
  call llh2sph( lat0, lon0, h0, f0(1), f0(2), f0(3), atmo%lonshift )
!
!.....Parallel loop
!
!$OMP PARALLEL do &
!$OMP& PRIVATE ( i, ia, ie, nofCaustics, f, d, s, ud, absorption, alpha, tloss, jacobian_det, &
!$OMP&           refl_hmax, altmin, altmax, altmean, nofs, nofrefl, elev1, azi1, atmo_specs, error, &
!$OMP&           k1, k2, k3, k4, df, lat, lon, h, rmin, hdLon, hdLat, dd, A, px, py, pz, rk_stepvar, &
!$OMP&           inclPrev, jacobian_det_is_zero, factor, rho, c, rho0, c0 ) &
!$OMP& SHARED  ( f0, azimuth_vec, elevation_vec, store_lon, store_lat, store_alt, store_capp, store_bazi, &
!$OMP&           store_incl, store_tt, store_tloss, store_d, store_hmin, store_hmax, &
!$OMP&           store_hmean, store_refl, cnt, progress )
!
!
! Start loop for all rays
!
  do i = 1, rays_count
!
!   get azimuth and elevation index
!
    ia = ceiling( 1.d0 * i / elevation_count )
    ie = i - ( ia - 1 ) * elevation_count
!
!........reset counters
!
    nofCaustics  = 0_int32
!
!........reset values
!
    f            = f0
    d            = 0.d0
    s            = 0.d0
    ud           = 0.d0
    absorption   = 0.d0
    tloss        = 0.d0
    jacobian_det = 0.d0
    refl_hmax    = 0.d0
    altmin       = 0.d0
    altmax       = 0.d0
    altmean      = 0.d0
    nofs         = 0
    nofrefl      = 0
!
!........get elevation and azimuth and convert to radians
!
    elev1        = (  90.d0 - elevation_vec( ie )  ) * d2r
    azi1         = (  azimuth_vec( ia ) - 90.d0 )  * d2r
!
!........initialize EOM for ray equations and transmission loss
!
    call initialize_eom( f, azi1, elev1, tloss_geometrical, atmo_specs, error )
!
    if ( error .eq. 0 ) then
!
!...........init Jacobian / geometrical spreading
!
      if (tloss_geometrical) then
        c0 = re_unit_conversion / sqrt( ( f(4)/f(3)/dsin(f(2)) )**2 &
          + ( f(5)/f(3) )** 2 + f(6)** 2 ) 
        rho0 = 1.d0
        rho  = 1.d0
        if (tloss_density) rho0 = atmo_specs(6)
      end if
!          
!...........show progress
!
      cnt = cnt + 1
!$OMP CRITICAL
      if ( verb ) then
        if ( anint( cnt * 100.d0 / rays_count ) .ge. progress ) then
          print &
            & "(3x,'-->',1x,'azimuth = ',F7.2"// &
            & ",' deg (',I3,1x,'%)')", &
            & azimuth_vec(ia), progress
          progress = progress + dprogress
        end if
      end if
!$OMP END CRITICAL 
!          
!.....Loop to integrate a ray (rk4)
!
      do
!
!.....adapt stepsize based on elevation angle (max resolution vertically)
!       elevation =  0 deg --> rk_stepvar = rk_stepsize
!       elevation = 30 deg --> rk_stepvar =  2 * rk_stepsize
!       elevation = 60 deg --> rk_stepvar =  5 * rk_stepsize
!       elevation = 90 deg --> rk_stepvar = 10 * rk_stepsize
!
        px = f(4)/f(3)/dsin(f(2))
        py = f(5)/f(3)
        pz = f(6)
        inclPrev   = dabs( dasin( dabs(pz) / dsqrt( px ** 2 + py ** 2 + pz ** 2 ) ) * r2d )
!
        rk_stepvar = rk_stepsize * ( 1.d0 + 9.d0 * ( inclPrev / 90.d0 ) ** 2 )
!
!.......time increase
!
        s = s + rk_stepvar
!
!.......count steps
!
        nofs = nofs + 1
!
!.......4 rk approximations
!
        call GET_EOM_DF ( f, df, tloss_geometrical, error )
        if ( error .eq. 0 ) then
          k1 = rk_stepvar * df
          call GET_EOM_DF ( f + .5d0 * k1, df, tloss_geometrical, error )
          if ( error .eq. 0 ) then
            k2 = rk_stepvar * df
            call GET_EOM_DF ( f + .5d0 * k2, df, tloss_geometrical, error )
            if ( error .eq. 0 ) then
              k3 = rk_stepvar * df
              call GET_EOM_DF ( f + k3, df, tloss_geometrical, error )
              if ( error .eq. 0 ) then
                k4 = rk_stepvar * df
!
!...............accumulated value
!
                f = f + ( k1 + 2.d0 * k2 + 2.d0 * k3 + k4 ) * A6TH
!
!...............correct for continuous longitude, and North and South Pole
!
                if ( atmo%l360 ) call twopirange( f(1) )
!
!...............get spherical distance
!
!!                call sphdistance(  f0(1), f0(2), f(1), f(2), d  )
!
!...............get euclidean distance
!
                call euclidean_distance_sph( &
                  & f0(1), f0(2), f0(3), f(1), f(2), f(3), ud &
                  & )
!
!...............Get attenuation (transmission loss)
!
                if ( tloss_attenuation ) then
!
!.................geometrical spreading (jacobian)
!
                  if (tloss_geometrical) then
                    call get_eom_df ( f, df, .true., error, atmo_specs )
                    call get_jacobian ( f, jacobian_det, nofCaustics, jacobian_det_is_zero )
                    if ( .not. jacobian_det_is_zero ) then
                      jacobian_det=jacobian_det
                      c = re_unit_conversion / sqrt( ( f(4)/f(3)/dsin(f(2)) )**2 &
                        + ( f(5)/f(3) )** 2 + f(6)** 2 )
                      if (tloss_density) rho = atmo_specs(6)
                      call get_spreading ( elev1, rho0, c0, rho, c, nofcaustics, jacobian_det, a )
                    end if
                  end if
!
!.................absorption
!
                  if (tloss_absorption) then
                    call atmo_absorption( (f(3) - RE), atmo_specs(1)*1000, atmo_specs(5), &
                      & atmo_specs(6), tloss_frequency, .true., alpha )
                    absorption = absorption + alpha
                  end if
!
!.................total transmission loss
!
                  call GET_TLOSS ( A0, A, absorption, tloss_frequency, f(19), tloss )
!
                end if
!
!...............Get earth radius at end point
!
                if ( ltopo ) then
                  call TOPO_GET ( lonrangef( f(1) * r2d - atmo%lonshift ), &
                    & 90.d0 - f(2) * r2d, h, hdLon, hdLat )
                  if ( h .lt. atmo%hMin ) h = atmo%hMin
                  rMin = RE + h
                end if
!
              end if
            end if
          end if
        end if
!
!.......Stop ray iteration conditions
!
        if (  (f(1).gt.phMax-1.d-8) .or. (f(1).lt.phMin+1.d-8) .or. &
          &  (f(2).gt.thMax-1.d-8) .or. (f(2).lt.thMin+1.d-8) .or. &
          &  (f(3).gt.rMax-1.d-8 ) .or. ud.gt.sc%d-1.d-8 .or. &
          &  error.eq.1 ) then
          exit
        else
!
!.........Reflection candidate
!
          if ( f(3).lt.rMin+1.d-8 .or. error.eq.-1 ) then
!
!...........Surface wave?
!
            if ( refl_hmax.gt.refl_hmin ) then
!
!.............Reflection count up
!
              nofrefl = nofrefl + 1
!
!.............Topography ?
!
              if ( ltopo ) then
!
!...............Topography to spherical
!
                hdLon =  hdLon/f(3)/dsin(f(2))
                hdLat = -hdLat/f(3)
!
!...............New slowness (reversed p_r and incorporating topography angle)
!
                factor = (  f(4) * hdLon * (  f(3) * dsin( f(2) )  )**( -2 )    &
                  &         + f(5) * hdLat * f(3)**( -2 ) + f(6)             ) &
                  &   / (  hdLon * (  f(3) * dsin( f(2) )  )**( -2 )           &
                  &         + hdLat * f(3)**( -2 ) + 1.d0                    ) &
                  &   * ( -2.d0 )
                f(4) = f(4) + factor * hdLon
                f(5) = f(5) + factor * hdLat
                f(6) = f(6) + factor
!
!.............No topography: simply reverse p_r
!
              else
!
                f(6) = -f(6)
!
              end if
!
!.............Set height to minimal height.
!
              f(3) = rMin
!
!.............Store minimal reflection height to the surface 
!
              if ( nofrefl .eq. 1 ) then
                altmin = refl_hmax
              else if ( refl_hmax .lt. altmin ) then
                altmin = refl_hmax
              end if
!
!.............Reset maximum height to the surface for each reflection
!
              refl_hmax = 0.d0
!
            else
!
!.............Surface wave, stop ray iteration
!
              exit
!
            end if ! surface wave?
!
          end if ! reflection candidate?
!
        end if ! end if exit conditions?
!
!.......Get maximum height in reflection above the surface 
!
        h = f(3) - rMin
        if ( h .gt. refl_hmax ) refl_hmax = h
        altmean = altmean + h
!
!.......Get maximum height above sea level 
!
        h = f(3) - RE
        if ( h .gt. altmax ) altmax = h
!
!.....End rk4 loop
!
      end do
!
!.....Write store f to output
!
      call sph2llh( f(1), f(2), f(3), lat, lon, h, atmo%lonshift )
!
      px = f(4)/f(3)/dsin(f(2))
      py = f(5)/f(3)
      pz = f(6)
!
      call euclidean_distance_sph( ph1, th1, r1, f(1), f(2), f(3), dd )
!            call euclidean_distance_llh( lat1, lon1, h1, lat, lon, h, dd )
!
      store_d(ia,ie) = dd
!
      if (clip.and.dd.gt.clip_max_offset) cycle
!
      store_lat(ia,ie)   = lat
      store_lon(ia,ie)   = lon
      store_alt(ia,ie)   = h
      store_tt(ia,ie)    = f(19)
      store_tloss(ia,ie) = tloss
      if (altmax.gt.0.d0) store_hmax(ia,ie) = altmax
      if (altmin.gt.0.d0) store_hmin(ia,ie) = altmin
      store_hmean(ia,ie) = altmean/nofs
      store_refl(ia,ie)  = nofrefl
      store_capp(ia,ie)  = 1000.d0 / dsqrt( px**2 + py**2 )
      store_bazi(ia,ie)  = lonrangef( datan2( px, -py )*r2d + 180.d0 )
      store_incl(ia,ie)  = dasin( dabs(pz) / dsqrt( px**2 + py**2 + pz**2 ) ) * r2d
!
!...End if okay initial ray conditions ( error == 0 )
!
    end if
!
! End ray loop
!
  end do
!
! End parallel loop
!
!$OMP END PARALLEL DO
!
  if (verb) print "(A)", ''
!
