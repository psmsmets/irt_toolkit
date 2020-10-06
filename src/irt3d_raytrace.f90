! --------------------------------------------------------------------------------------------------
!
! RAYTRACE (rk4 integration)
!
! --------------------------------------------------------------------------------------------------
!
      IF ( verbose ) PRINT "(A)"         , 'RAYTRACING'
      IF ( verbose ) PRINT "(3x,A,1x,I6)", '#RAYS        = ', nofRays
!
      progress = 0
!
      if ( nofRays .lt. 10 ) then
        dprogress = 50
      else if ( nofRays .lt. 20 ) then
        dprogress = 25
      else
        dprogress = 10
      end if
!
!.....make 2d array of angle combinations
!
      ALLOCATE (  angles( nofRays, 2 )  )
!
      angles( :, 1 ) = (/   ( ( elevMin + i * elevStep ) - &
        & CEILING (  DBLE( i / nofElev )  ) * nofElev * elevStep, &
        & i = 0, nofRays - 1, 1 )  /)
      angles( :, 2 ) = (/  &
        & ( aziMin + CEILING (  DBLE( i / nofElev )  ) * aziStep, &
        & i = 0, nofRays - 1, 1 )  /)
!     
!.....Write angles to file
!
      OPEN ( UNIT = 10, STATUS = 'replace', FILE = out_angle, &
        & FORM = 'formatted' ) 
      WRITE ( 10, "(F10.5,1x,F10.5)" ) TRANSPOSE ( angles )
      CLOSE ( 10 )   
!     
!.....Prepare output file
!
      if ( lsaverefl ) open ( unit = 11, status = 'replace', &
        & file = out_refl, form = 'unformatted', access = 'stream' ) 
      if ( lsaverays ) open ( unit = 12, status = 'replace', &
        & file = out_rays, form = 'unformatted', access = 'stream' ) 
! 
!.....initialize
!
      f0         = 0.d0
      store_rays = 0.d0
      store_refl = 0.d0
!
!.....initialize spreading
!
      IF ( jacobian ) CALL INITIALIZE_SPREADING ( re_unit, A0 )
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
      if ( refl_hmin .lt. .5d0 ) refl_hmin = .5d0
!     
!.....Source parameters, from lonlat to spherical
!
      call llh2sph( lat0, lon0, h0, f0(1), f0(2), f0(3), atmo%lonshift )
!
!.....start loop for all rays
!
      do i = 1, nofRays
!
!.......reset counters
!
        cnt          = 1
        nofRefl      = 0
        nofCaustics  = 0
!
!.......reset values
!
        f            = f0
        d            = 0.d0
        s            = 0.d0
        absorption   = 0.d0
        tloss        = 0.d0
        jacobian_det = 0.d0
        altMax       = 0.d0
!
!.......reset logical
!
        lreflection  = .false.
!
!.......get elevation and azimuth and convert to radians
!
        elev1        = (  90.d0 - angles( i, 1 )  ) * D2R
        azi1         = (  angles( i, 2 ) - 90.d0 )  * D2R
!
!.......initialize EOM for ray equations and transmission loss
!
        call initialize_eom ( f, azi1, elev1, jacobian, atmo_specs, ierr )
!
        if ( ierr .eq. 0 ) then
!
!.........init Jacobian / geometrical spreading
!
          if (jacobian) then
            c0 = re_unit_conversion / sqrt( ( f(4)/f(3)/dsin(f(2)) )**2 &
              + ( f(5)/f(3) )** 2 + f(6)** 2 )
            rho0 = 1.d0
            rho  = 1.d0
            if (jacobian_density) rho0 = atmo_specs(6)
          end if
!
!.........store rays
!
          IF ( lsaveRays ) store_rays( 1, : ) = (/ DBLE ( i ), f( 19 ), lon0, lat0, h0, &
!!               & f(4) * f(3) * DSIN ( f(2) ), f(5) * f(3), f(6), 0.d0 /)
            & 0.d0, 0.d0 /)
!          
!.........show progress
!
          if ( verbose ) then
            if ( anint( i * 100.d0 / nofrays ) .ge. progress ) then
              print &
                & "(3x,'-->',1x,'azimuth = ',F7.2"// &
                & ",' deg (',I3,1x,'%)')", &
                & angles(i,2), progress
              progress = progress + dprogress
            end if
          end if
!          
!.........Loop to integrate a ray (rk4)
!
          DO
!
!...........adapt stepsize based on elevation angle (max resolution vertically)
!                 elevation =  0 deg --> rk_stepvar = rk_stepsize
!                 elevation = 30 deg --> rk_stepvar =  2 * rk_stepsize
!                 elevation = 60 deg --> rk_stepvar =  5 * rk_stepsize
!                 elevation = 90 deg --> rk_stepvar = 10 * rk_stepsize
!
            px = f(4)/f(3)/dsin(f(2))
            py = f(5)/f(3)
            pz = f(6)
            inclPrev   = dabs( dasin( dabs( pz ) / dsqrt( px ** 2 + py ** 2 + pz ** 2 ) ) * r2d )
!
            rk_stepvar = rk_stepsize * ( 1.d0 + 9.d0 * ( inclPrev / 90.d0 ) ** 2 )
!
!...........time increase
!
            s = s + rk_stepvar
!
!...........4 rk approximations
!
            CALL GET_EOM_DF ( f, df, jacobian, ierr )
            IF ( ierr .EQ. 0 ) THEN
              k1 = rk_stepvar * df
              CALL GET_EOM_DF ( f + .5d0 * k1, df, jacobian, ierr )
              IF ( ierr .EQ. 0 ) THEN
                k2 = rk_stepvar * df
                CALL GET_EOM_DF ( f + .5d0 * k2, df, jacobian, ierr )
                IF ( ierr .EQ. 0 ) THEN
                  k3 = rk_stepvar * df
                  CALL GET_EOM_DF ( f + k3, df, jacobian, ierr )
                  IF ( ierr .EQ. 0 ) THEN
                    k4 = rk_stepvar * df
!
!...................accumulated value
!
                    f = f + ( k1 + 2.d0 * k2 + 2.d0 * k3 + k4 ) * A6TH
!
!...................correct for continuous longitude, and North and South Pole
!
                    IF ( atmo%l360 ) CALL TWOPIRANGE ( f(1) )
!
!...................get distance (in km!)
!
                    CALL SPHDISTANCE (  f0(1), f0(2), f(1), f(2), d  )
!
!...................get attenuation (transmission loss): spreading ( jacobian ) and absorption
!
                    IF ( jacobian ) THEN
!
                      CALL GET_EOM_DF ( f, df, .TRUE., ierr, atmo_specs )
                      CALL GET_JACOBIAN ( f, jacobian_det, nofCaustics, jacobian_det_is_zero )
!
                      if ( .not. jacobian_det_is_zero ) then
                        jacobian_det=jacobian_det
                        c = re_unit_conversion / sqrt( ( f(4)/f(3)/dsin(f(2)) )**2 &
                          + ( f(5)/f(3) )** 2 + f(6)** 2 )
                        if (jacobian_density) rho =  atmo_specs(6)
                        call get_spreading ( elev1, rho0, c0, rho, c, nofcaustics, jacobian_det, a )
                      end if
!
!.....................absorption
!
!!                      CALL ATMO_ABSORPTION ( (f(3) - RE), atmo_specs(1)*1000, atmo_specs(5), &
!!                        & atmo_specs(6), frq, .TRUE., alpha )
!!                      absorption = absorption + alpha
!
!.....................transmission loss
!
                      CALL GET_TLOSS ( A0, A, absorption, frq, f( 19 ), tloss )
!
                    END IF
!
!...................get earth radius at end point
!
                    IF ( ltopo ) THEN
                      call TOPO_GET ( lonrangef( f(1) * r2d - atmo%lonshift ), &
                        & 90.d0 - f(2) * r2d, h, hdLon, hdLat )
                      IF ( h .LT. atmo%hMin ) h = atmo%hMin
                      rMin = RE + h
                    END IF
!
                  END IF
                END IF
              END IF
            END IF
!
!...........stop ray iteration conditions
!
            if (  (f(1).gt.phMax-1.d-8) .or. (f(1).lt.phMin+1.d-8) .or. &
               &  (f(2).gt.thMax-1.d-8) .or. (f(2).lt.thMin+1.d-8) .or. &
               &  (f(3).gt.rMax-1.d-8 ) .or. d.gt.sc%d-1.d-8 .or. &
               &  ierr.eq.1 ) then
                  exit
            else
!
!.............Reflection candidate
!
              if ( f(3).lt.rmin+1.d-8 .or. ierr.eq.-1 ) then
!
!...............Reflection (on or just below surface of the earth).
!
                IF ( altMax .GT. refl_hmin ) THEN
!
                  lreflection = .TRUE.
!
!.................topography ?
!
                  IF ( ltopo ) THEN
!
!...................to spherical
!
                    hdLon =   hdLon / (  f(3) * DSIN ( f(2) )  )
                    hdLat = - hdLat / f(3)
!
!...................New slowness (reversed p_r and incorporating topography angle)
!
                    factor = (  f(4) * hdLon * (  f(3) * DSIN ( f(2) )  )**( -2 )    &
                      &         + f(5) * hdLat * f(3)**( -2 ) + f(6)              ) &
                      &   / (  hdLon * (  f(3) * DSIN ( f(2) )  )**( -2 )             &
                      &         + hdLat * f(3)**( -2 ) + 1.d0                         ) &
                      &   * ( -2.d0 )
                    f(4) = f(4) + factor * hdLon
                    f(5) = f(5) + factor * hdLat
                    f(6) = f(6) + factor
!
!.................No topography: simply reverse p_r
!
                  ELSE
!
                    f(6) = - f(6)
!
                  END IF
!
!.................Set height to minimal height.
!
                  f(3) = rMin
!
!...............Surface wave, does not count as reflection
!
                ELSE
!
                  EXIT ! stop ray iteration
!
                END IF
!
!.............Catch near reflection points!
!
              ELSE IF (  f(3) .LT. rMin + refl_hmin * .99d0 .AND. cnt .GT. 10 &
                & .AND. altMax .GT. refl_hmin ) THEN
                lreflection = store_rays( cnt, 5 ) - store_rays( cnt - 1, 5 ) .LT. 0.d0 &
                  & .AND. f(3) - RE - store_rays( cnt, 5 ) .GT. 0.d0
              END IF
!
!.............ray step counter up
!
              cnt = cnt + 1
!
!.............store ray position
!                 (/ ray number,time,lon,lat,h,p_lon,p_lat,p_h,0 /)
!
              IF ( lsaveRays .AND. lreflection ) store_rays( cnt, : ) = &
                & (/  DBLE ( i ), f( 19 ), LONRANGEF (  f(1) * R2D - atmo%lonshift ), &
                &     90.d0 - f(2) * R2D, f(3) - RE, tloss, altMax &
!!                     &     f(4) / f(3) / DSIN ( f(2) ), f(5) / f(3), f(6) &
                & /)
              IF ( lsaveRays .AND. .NOT. lreflection ) store_rays( cnt, : ) = &
                & (/  DBLE ( i ), f( 19 ), LONRANGEF (  f(1) * R2D - atmo%lonshift ), &
                &     90.d0 - f(2) * R2D, f(3) - RE, tloss, 0.d0 &
!!                     &     f(4) / f(3) / DSIN ( f(2) ), f(5) / f(3), f(6) &
                & /)
!
!.............Max height
!
              IF ( f(3) - RE .GT. altMax ) altMax = f(3) - RE
!
!.............Reset altitude reflection height
!
              IF ( lreflection ) THEN
!
                lreflection = .FALSE.
!
!...............increase reflection counter
!
                nofRefl = nofRefl + 1
!
!...............store reflection
!                    (/ ray number,time,lon,lat,h,p_lon,p_lat,p_h,refl altitude /)
!
                store_refl( nofRefl, : ) = &
                  & (/  DBLE ( i ), f( 19 ), LONRANGEF (  f(1) * R2D - atmo%lonshift  ), &
                  &     90.d0 - f(2) * R2D, h, tloss, altMax &
!!                        &     f(4) / f(3) / DSIN( f(2) ), f(5) / f(3), f(6) &
                  &  /)
! 
!...............Reset height properties
!
                altMax  = 0.d0
!
              END IF
!
            END IF ! ENDIF exit conditions
!
!.........end rk4 loop
!
          END DO
!          
!.........write ray data
!
          IF ( lsaveRefl ) WRITE ( 11 ) TRANSPOSE (  store_refl( 1 : nofRefl, : )  )
          IF ( lsaveRays ) WRITE ( 12 ) TRANSPOSE (  store_rays( 1 : cnt, : )  )
!
!........end if okay initial ray conditions ( ierr == 0 )
!
        END IF
!
!.....end ray loop
!
      END DO
!     
!.....Close output files
!
      IF ( lsaveRays ) CLOSE ( 11 )
      IF ( lsaveRefl ) CLOSE ( 12 )
!
      IF ( verbose ) PRINT "(A)", ''
!
