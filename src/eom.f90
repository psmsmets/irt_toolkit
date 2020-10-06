! =========================================================================================
!
! MODULE EOM
!
! =========================================================================================
!
MODULE EOM
!
CONTAINS
! =========================================================================================
!
!.......SUBROUTINE GET_DF 
!
!       Get derivatives of equation system for ray equations: f( 1:6 ) and df( 1:6 ),
!       and the jacobian if set: f( 7:18 ) and df( 7:18 )
!
!****************************
   SUBROUTINE GET_EOM_DF ( f, df, jacobian, ierr, atmo_specs )
!****************************
!
      USE TOOLS_ATMO
      USE MATH, ONLY: SPH2LLH;
      USE MATH, ONLY: TWOPIRANGEF;
      USE MATH, ONLY: CSC;
      USE MATH, ONLY: COT;
      IMPLICIT NONE
!    
!.....dummy variables
!
      DOUBLE PRECISION, DIMENSION ( 19 ), INTENT ( IN  )            :: f
      DOUBLE PRECISION, DIMENSION ( 19 ), INTENT ( OUT )            :: df
      LOGICAL                           , INTENT ( IN  )            :: jacobian
      integer(kind=4)                       , INTENT ( OUT )            :: ierr
      DOUBLE PRECISION, DIMENSION ( 6 ) , INTENT ( OUT ), OPTIONAL  :: atmo_specs
!
!.....local variables
!
      DOUBLE PRECISION                              :: lat1, dlon1, h1
      DOUBLE PRECISION, DIMENSION ( 0:2, 0:2, 0:2 ) :: C, Wph, Wth, Wr, P, Rho
!
! ---
!
!.....get atmospheric values
!          
      CALL SPH2LLH (  TWOPIRANGEF ( f(1) ), f(2), f(3), lat1, dlon1, h1 )      ! spherical to lonlat
      CALL ATMO_GET ( dlon1, lat1, h1, .TRUE., ierr, C, Wph, Wth, Wr, P, Rho ) ! get atmospheric specifiations
!
      IF ( ierr .NE. 0 ) RETURN
!
!.....if requested, return atmo specs (no derivs!)
!
      IF (  PRESENT ( atmo_specs )  ) atmo_specs = (/ C(0,0,0), Wph(0,0,0), Wth(0,0,0), Wr(0,0,0), P(0,0,0), Rho(0,0,0)  /)
!          
!.....wind conversion to spherical
!
       Wph = Wph / f(3) / dsin( f(2) )  ! wind in azimuth direction (rad/s)
       Wth = Wth / f(3)
!
!......Derivative structure
!
!         df(  1 ) = d X_ph / dS
!         df(  2 ) = d X_th / dS
!         df(  3 ) = d X_r  / dS
!
!         df(  4 ) = d P_ph / dS
!         df(  5 ) = d P_th / dS
!         df(  6 ) = d P_r  / dS
!
!         df(  7 ) = d X_ph / dA
!         df(  8 ) = d X_th / dA
!         df(  9 ) = d X_r  / dA
!
!         df( 10 ) = d X_ph / dE
!         df( 11 ) = d X_th / dE
!         df( 12 ) = d X_r  / dE
!
!         df( 13 ) = d P_ph / dA
!         df( 14 ) = d P_th / dA
!         df( 15 ) = d P_r  / dA
!
!         df( 16 ) = d P_ph / dE
!         df( 17 ) = d P_th / dE
!         df( 18 ) = d P_r  / dE
!
! 
!.....d X_ph / dS
!
      df( 1 ) = &
         & Wph(0,0,0) * (  1 - f(4) * Wph(0,0,0) - f(5) * Wth(0,0,0) - f(6) * Wr(0,0,0)  ) &
         & + (  C(0,0,0) ** 2  ) * f(4) / (   (  f(3) * dsin( f(2) )  ) ** 2   )
!
!.....d X_th / dS
!
      df( 2 ) = &
         & Wth(0,0,0) * (  1 - f(4) * Wph(0,0,0) - f(5) * Wth(0,0,0) - f(6) * Wr(0,0,0)  ) &
         & + (  C(0,0,0) ** 2  ) * f(5) / (  f(3) ** 2  )
!
!.....d X_r / dS
!
      df( 3 ) = &
         & Wr(0,0,0) * (  1 - f(4) * Wph(0,0,0) - f(5) * Wth(0,0,0) - f(6) * Wr(0,0,0)  ) &
         & + (  C(0,0,0) ** 2  ) * f(6)
!
!.....d P_ph / dS
!
      df( 4 ) = &
         & - C(0,0,0) * (   (  f(4) / f(3) / dsin( f(2) )  ) ** 2 + (  f(5) / f(3)  ) ** 2 + f(6) ** 2   ) * C(1,0,0) &
         & + (  1 - f(4) * Wph(0,0,0) - f(5) * Wth(0,0,0) - f(6) * Wr(0,0,0)  ) &
         & * (  - f(4) * Wph(1,0,0) - f(5) * Wth(1,0,0) - f(6) * Wr(1,0,0)  ) 
!
!.....d P_th / dS
!
      df( 5 ) = &
         & C(0,0,0) ** 2 * f(4) ** 2 / (  f(3)**2  ) / (  dsin( f(2) ) ** 2  ) / dtan( f(2) ) &
         &  - C(0,0,0) * (   (  f(4) / f(3) / dsin( f(2) )  ) ** 2 + (  f(5) / f(3)  ) ** 2 + f(6) ** 2   ) * C(0,1,0) &
         & + (  1 - f(4) * Wph(0,0,0) - f(5) * Wth(0,0,0) - f(6) * Wr(0,0,0)  ) &
         & * (  - f(4) * Wph(0,1,0) - f(5) * Wth(0,1,0) - f(6) * Wr(0,1,0)  )
!
!.....d P_r / dS
!
      df( 6 ) = &
         & C(0,0,0) ** 2 * (   f(5) ** 2 / (  f(3) ** 3  ) + f(4) ** 2 / (  f(3) ** 3 ) / (  dsin( f(2) ) ** 2  )   ) &
         & - C(0,0,0) * (   (  f(4) / f(3) / dsin( f(2) )  ) ** 2 + (  f(5) / f(3)  ) ** 2 + f(6) ** 2   ) * C(0,0,1) &
         & + (  1 - f(4) * Wph(0,0,0) - f(5) * Wth(0,0,0) - f(6) * Wr(0,0,0)  ) &
         & * (  - f(4) * Wph(0,0,1) - f(5) * Wth(0,0,1) - f(6) * Wr(0,0,1)  )
!
      IF ( jacobian ) THEN
!
!........d X_ph / dA
!
         df( 7 ) = &
            & 0.5*(C(0,0,0)**2*((4*CSC(f(2))**2*f(9)*f(4))/f(3)**3 + (2*CSC(f(2))**2*f(13))/f(3)**2 &
            & - (2*CSC(f(2))**4*f(4)*(2*f(9)*dsin(f(2))**2*f(3)**3 + f(3)**2*(2*f(9)*dsin(f(2))**2*f(3) &
            & + 2*dcos(f(2))*f(8)*dsin(f(2))*f(3)**2)))/f(3)**6) + (4*C(0,0,0)*CSC(f(2))**2*f(4)*(f(7)*C(1,0,0) &
            & + f(8)*C(0,1,0) + f(9)*C(0,0,1)))/f(3)**2 + 2*(1 - f(6)*Wr(0,0,0) - f(5)*Wth(0,0,0) &
            & - f(4)*Wph(0,0,0))*(f(7)*Wph(1,0,0) + f(8)*Wph(0,1,0) + f(9)*Wph(0,0,1)) &
            & + 2*Wph(0,0,0)*(-(f(15)*Wr(0,0,0)) - f(14)*Wth(0,0,0) - f(13)*Wph(0,0,0) - f(6)*(f(7)*Wr(1,0,0) &
            & + f(8)*Wr(0,1,0) + f(9)*Wr(0,0,1)) - f(5)*(f(7)*Wth(1,0,0) + f(8)*Wth(0,1,0) + f(9)*Wth(0,0,1)) &
            & - f(4)*(f(7)*Wph(1,0,0) + f(8)*Wph(0,1,0) + f(9)*Wph(0,0,1))))
!
!........d X_th / dA
!
         df( 8 ) = &
            & 0.5*(C(0,0,0)**2*((2*f(14))/f(3)**2 + (2*CSC(f(2))**2*f(5)*(2*f(9)*dsin(f(2))**2*f(3) &
            & + 2*dcos(f(2))*f(8)*dsin(f(2))*f(3)**2))/f(3)**4 - (2*CSC(f(2))**2*f(5)*(2*f(9)*dsin(f(2))**2*f(3)**3 &
            & + f(3)**2*(2*f(9)*dsin(f(2))**2*f(3) + 2*dcos(f(2))*f(8)*dsin(f(2))*f(3)**2)))/f(3)**6) &
            & + (4*C(0,0,0)*f(5)*(f(7)*C(1,0,0) + f(8)*C(0,1,0) + f(9)*C(0,0,1)))/f(3)**2 + 2*(1 - f(6)*Wr(0,0,0) &
            & - f(5)*Wth(0,0,0) - f(4)*Wph(0,0,0))*(f(7)*Wth(1,0,0) + f(8)*Wth(0,1,0) + f(9)*Wth(0,0,1)) &
            & + 2*Wth(0,0,0)*(-(f(15)*Wr(0,0,0)) - f(14)*Wth(0,0,0) - f(13)*Wph(0,0,0) - f(6)*(f(7)*Wr(1,0,0) &
            & + f(8)*Wr(0,1,0) + f(9)*Wr(0,0,1)) - f(5)*(f(7)*Wth(1,0,0) + f(8)*Wth(0,1,0) + f(9)*Wth(0,0,1)) &
            & - f(4)*(f(7)*Wph(1,0,0) + f(8)*Wph(0,1,0) + f(9)*Wph(0,0,1))))
!
!........d X_r / dA
!
         df( 9 ) = &
            & 0.5*(2*C(0,0,0)**2*f(15) + 4*C(0,0,0)*f(6)*(f(7)*C(1,0,0) + f(8)*C(0,1,0) + f(9)*C(0,0,1)) &
            & + 2*(1 - f(6)*Wr(0,0,0) - f(5)*Wth(0,0,0) - f(4)*Wph(0,0,0))*(f(7)*Wr(1,0,0) + f(8)*Wr(0,1,0) &
            & + f(9)*Wr(0,0,1)) + 2*Wr(0,0,0)*(-(f(15)*Wr(0,0,0)) - f(14)*Wth(0,0,0) - f(13)*Wph(0,0,0) &
            & - f(6)*(f(7)*Wr(1,0,0) + f(8)*Wr(0,1,0) + f(9)*Wr(0,0,1)) - f(5)*(f(7)*Wth(1,0,0) + f(8)*Wth(0,1,0) &
            & + f(9)*Wth(0,0,1)) - f(4)*(f(7)*Wph(1,0,0) + f(8)*Wph(0,1,0) + f(9)*Wph(0,0,1))))
!
!........d X_ph / dE
!
         df( 10 ) = &
            & 0.5*(C(0,0,0)**2*((4*CSC(f(2))**2*f(12)*f(4))/f(3)**3 + (2*CSC(f(2))**2*f(16))/f(3)**2 &
            & - (2*CSC(f(2))**4*f(4)*(2*f(12)*dsin(f(2))**2*f(3)**3 + f(3)**2*(2*f(12)*dsin(f(2))**2*f(3) &
            & + 2*dcos(f(2))*f(11)*dsin(f(2))*f(3)**2)))/f(3)**6) + (4*C(0,0,0)*CSC(f(2))**2*f(4)*(f(10)*C(1,0,0) &
            & + f(11)*C(0,1,0) + f(12)*C(0,0,1)))/f(3)**2 + 2*(1 - f(6)*Wr(0,0,0) - f(5)*Wth(0,0,0) &
            & - f(4)*Wph(0,0,0))*(f(10)*Wph(1,0,0) + f(11)*Wph(0,1,0) + f(12)*Wph(0,0,1)) &
            & + 2*Wph(0,0,0)*(-(f(18)*Wr(0,0,0)) - f(17)*Wth(0,0,0) - f(16)*Wph(0,0,0) - f(6)*(f(10)*Wr(1,0,0) &
            & + f(11)*Wr(0,1,0) + f(12)*Wr(0,0,1)) - f(5)*(f(10)*Wth(1,0,0) + f(11)*Wth(0,1,0) + f(12)*Wth(0,0,1)) &
            & - f(4)*(f(10)*Wph(1,0,0) + f(11)*Wph(0,1,0) + f(12)*Wph(0,0,1))))
!
!........d X_th / dE
!
         df( 11 ) = &
            & 0.5*(C(0,0,0)**2*((2*f(17))/f(3)**2 + (2*CSC(f(2))**2*f(5)*(2*f(12)*dsin(f(2))**2*f(3) &
            & + 2*dcos(f(2))*f(11)*dsin(f(2))*f(3)**2))/f(3)**4 - (2*CSC(f(2))**2*f(5)*(2*f(12)*dsin(f(2))**2*f(3)**3 &
            & + f(3)**2*(2*f(12)*dsin(f(2))**2*f(3) + 2*dcos(f(2))*f(11)*dsin(f(2))*f(3)**2)))/f(3)**6) &
            & + (4*C(0,0,0)*f(5)*(f(10)*C(1,0,0) + f(11)*C(0,1,0) + f(12)*C(0,0,1)))/f(3)**2 &
            & + 2*(1 - f(6)*Wr(0,0,0) - f(5)*Wth(0,0,0) - f(4)*Wph(0,0,0))*(f(10)*Wth(1,0,0) &
            & + f(11)*Wth(0,1,0) + f(12)*Wth(0,0,1)) + 2*Wth(0,0,0)*(-(f(18)*Wr(0,0,0)) - f(17)*Wth(0,0,0) &
            & - f(16)*Wph(0,0,0) - f(6)*(f(10)*Wr(1,0,0) + f(11)*Wr(0,1,0) + f(12)*Wr(0,0,1)) &
            & - f(5)*(f(10)*Wth(1,0,0) + f(11)*Wth(0,1,0) + f(12)*Wth(0,0,1)) - f(4)*(f(10)*Wph(1,0,0) &
            & + f(11)*Wph(0,1,0) + f(12)*Wph(0,0,1))))
!
!........d X_r / dE
!
         df( 12 ) = &
            & 0.5*(2*C(0,0,0)**2*f(18) + 4*C(0,0,0)*f(6)*(f(10)*C(1,0,0) + f(11)*C(0,1,0) + f(12)*C(0,0,1)) &
            & + 2*(1 - f(6)*Wr(0,0,0) - f(5)*Wth(0,0,0) - f(4)*Wph(0,0,0))*(f(10)*Wr(1,0,0) + f(11)*Wr(0,1,0) &
            & + f(12)*Wr(0,0,1)) + 2*Wr(0,0,0)*(-(f(18)*Wr(0,0,0)) - f(17)*Wth(0,0,0) - f(16)*Wph(0,0,0) &
            & - f(6)*(f(10)*Wr(1,0,0) + f(11)*Wr(0,1,0) + f(12)*Wr(0,0,1)) - f(5)*(f(10)*Wth(1,0,0) &
            & + f(11)*Wth(0,1,0) + f(12)*Wth(0,0,1)) - f(4)*(f(10)*Wph(1,0,0) + f(11)*Wph(0,1,0) + f(12)*Wph(0,0,1))))
!
!........d P_ph / dA
!
         df( 13 ) = &
            & -0.5*(2*C(0,0,0)*(2*f(15)*f(6) + (f(14)*f(5))/f(3)**2 + (CSC(f(2))**2*f(13)*f(4))/f(3)**2 &
            & + f(5)*(f(14)/f(3)**2 + (CSC(f(2))**2*f(5)*(2*f(9)*dsin(f(2))**2*f(3) &
            & + 2*dcos(f(2))*f(8)*dsin(f(2))*f(3)**2))/f(3)**4 - (CSC(f(2))**2*f(5)*(2*f(9)*dsin(f(2))**2*f(3)**3 &
            & + f(3)**2*(2*f(9)*dsin(f(2))**2*f(3) + 2*dcos(f(2))*f(8)*dsin(f(2))*f(3)**2)))/f(3)**6) &
            & + f(4)*((2*CSC(f(2))**2*f(9)*f(4))/f(3)**3 + (CSC(f(2))**2*f(13))/f(3)**2 &
            & - (CSC(f(2))**4*f(4)*(2*f(9)*dsin(f(2))**2*f(3)**3 + f(3)**2*(2*f(9)*dsin(f(2))**2*f(3) &
            & + 2*dcos(f(2))*f(8)*dsin(f(2))*f(3)**2)))/f(3)**6))*C(1,0,0) + 2*(f(6)**2 + f(5)**2/f(3)**2 &
            & + (CSC(f(2))**2*f(4)**2)/f(3)**2)*C(1,0,0)*(f(7)*C(1,0,0) + f(8)*C(0,1,0) + f(9)*C(0,0,1)) &
            & - 2*(-(f(6)*Wr(1,0,0)) - f(5)*Wth(1,0,0) - f(4)*Wph(1,0,0))*(-(f(15)*Wr(0,0,0)) - f(14)*Wth(0,0,0) &
            & - f(13)*Wph(0,0,0) - f(6)*(f(7)*Wr(1,0,0) + f(8)*Wr(0,1,0) + f(9)*Wr(0,0,1)) - f(5)*(f(7)*Wth(1,0,0) &
            & + f(8)*Wth(0,1,0) + f(9)*Wth(0,0,1)) - f(4)*(f(7)*Wph(1,0,0) + f(8)*Wph(0,1,0) + f(9)*Wph(0,0,1))) &
            & + 2*C(0,0,0)*(f(6)**2 + f(5)**2/f(3)**2 + (CSC(f(2))**2*f(4)**2)/f(3)**2)*(f(7)*C(2,0,0) + f(8)*C(1,1,0) &
            & + f(9)*C(1,0,1)) - 2*(1 - f(6)*Wr(0,0,0) - f(5)*Wth(0,0,0) - f(4)*Wph(0,0,0))*(-(f(15)*Wr(1,0,0)) &
            & - f(14)*Wth(1,0,0) - f(13)*Wph(1,0,0) - f(6)*(f(7)*Wr(2,0,0) + f(8)*Wr(1,1,0) + f(9)*Wr(1,0,1)) &
            & - f(5)*(f(7)*Wth(2,0,0) + f(8)*Wth(1,1,0) + f(9)*Wth(1,0,1)) - f(4)*(f(7)*Wph(2,0,0) + f(8)*Wph(1,1,0) &
            & + f(9)*Wph(1,0,1))))
!
!........d P_th / dA
!
         df( 14 ) = &
            & -0.5*(C(0,0,0)**2*((-2*Cot(f(2))*CSC(f(2))**2*f(13)*f(4))/f(3)**2 &
            & + f(5)*((-2*Cot(f(2))*CSC(f(2))**2*f(5)*(2*f(9)*dsin(f(2))**2*f(3) &
            & + 2*dcos(f(2))*f(8)*dsin(f(2))*f(3)**2))/f(3)**4 + (2*Cot(f(2))*CSC(f(2))**2*f(5)*(2*f(9)*dsin(f(2))**2*f(3)**3 &
            & + f(3)**2*(2*f(9)*dsin(f(2))**2*f(3) + 2*dcos(f(2))*f(8)*dsin(f(2))*f(3)**2)))/f(3)**6 &
            & + (CSC(f(2))**2*f(5)*(4*dcos(f(2))*f(9)*dsin(f(2))*f(3) + f(8)*(2*dcos(f(2))**2*f(3)**2 &
            & - 2*dsin(f(2))**2*f(3)**2)))/f(3)**4 - (CSC(f(2))**2*f(5)*(4*dcos(f(2))*f(9)*dsin(f(2))*f(3)**3 &
            & + f(3)**2*(4*dcos(f(2))*f(9)*dsin(f(2))*f(3) + f(8)*(2*dcos(f(2))**2*f(3)**2 - 2*dsin(f(2))**2*f(3)**2))))/f(3)**6) &
            & + f(4)*((-4*Cot(f(2))*CSC(f(2))**2*f(9)*f(4))/f(3)**3 - (2*Cot(f(2))*CSC(f(2))**2*f(13))/f(3)**2 &
            & + (4*Cot(f(2))*CSC(f(2))**4*f(4)*(2*f(9)*dsin(f(2))**2*f(3)**3 + f(3)**2*(2*f(9)*dsin(f(2))**2*f(3) &
            & + 2*dcos(f(2))*f(8)*dsin(f(2))*f(3)**2)))/f(3)**6 - (CSC(f(2))**4*f(4)*(4*dcos(f(2))*f(9)*dsin(f(2))*f(3)**3 &
            & + f(3)**2*(4*dcos(f(2))*f(9)*dsin(f(2))*f(3) + f(8)*(2*dcos(f(2))**2*f(3)**2 - 2*dsin(f(2))**2*f(3)**2))))/f(3)**6)) &
            & + 2*C(0,0,0)*(2*f(15)*f(6) + (f(14)*f(5))/f(3)**2 + (CSC(f(2))**2*f(13)*f(4))/f(3)**2 + f(5)*(f(14)/f(3)**2 &
            & + (CSC(f(2))**2*f(5)*(2*f(9)*dsin(f(2))**2*f(3) + 2*dcos(f(2))*f(8)*dsin(f(2))*f(3)**2))/f(3)**4 &
            & - (CSC(f(2))**2*f(5)*(2*f(9)*dsin(f(2))**2*f(3)**3 + f(3)**2*(2*f(9)*dsin(f(2))**2*f(3) &
            & + 2*dcos(f(2))*f(8)*dsin(f(2))*f(3)**2)))/f(3)**6) + f(4)*((2*CSC(f(2))**2*f(9)*f(4))/f(3)**3 &
            & + (CSC(f(2))**2*f(13))/f(3)**2 - (CSC(f(2))**4*f(4)*(2*f(9)*dsin(f(2))**2*f(3)**3 &
            & + f(3)**2*(2*f(9)*dsin(f(2))**2*f(3) + 2*dcos(f(2))*f(8)*dsin(f(2))*f(3)**2)))/f(3)**6))*C(0,1,0) &
            & - (4*C(0,0,0)*Cot(f(2))*CSC(f(2))**2*f(4)**2*(f(7)*C(1,0,0) + f(8)*C(0,1,0) + f(9)*C(0,0,1)))/f(3)**2 &
            & + 2*(f(6)**2 + f(5)**2/f(3)**2 + (CSC(f(2))**2*f(4)**2)/f(3)**2)*C(0,1,0)*(f(7)*C(1,0,0) + f(8)*C(0,1,0) &
            & + f(9)*C(0,0,1)) - 2*(-(f(6)*Wr(0,1,0)) - f(5)*Wth(0,1,0) - f(4)*Wph(0,1,0))*(-(f(15)*Wr(0,0,0)) &
            & - f(14)*Wth(0,0,0) - f(13)*Wph(0,0,0) - f(6)*(f(7)*Wr(1,0,0) + f(8)*Wr(0,1,0) + f(9)*Wr(0,0,1)) &
            & - f(5)*(f(7)*Wth(1,0,0) + f(8)*Wth(0,1,0) + f(9)*Wth(0,0,1)) - f(4)*(f(7)*Wph(1,0,0) + f(8)*Wph(0,1,0) &
            & + f(9)*Wph(0,0,1))) + 2*C(0,0,0)*(f(6)**2 + f(5)**2/f(3)**2 + (CSC(f(2))**2*f(4)**2)/f(3)**2)*(f(7)*C(1,1,0) &
            & + f(8)*C(0,2,0) + f(9)*C(0,1,1)) - 2*(1 - f(6)*Wr(0,0,0) - f(5)*Wth(0,0,0) - f(4)*Wph(0,0,0))*(-(f(15)*Wr(0,1,0)) &
            & - f(14)*Wth(0,1,0) - f(13)*Wph(0,1,0) - f(6)*(f(7)*Wr(1,1,0) + f(8)*Wr(0,2,0) + f(9)*Wr(0,1,1)) &
            & - f(5)*(f(7)*Wth(1,1,0) + f(8)*Wth(0,2,0) + f(9)*Wth(0,1,1)) - f(4)*(f(7)*Wph(1,1,0) + f(8)*Wph(0,2,0) &
            & + f(9)*Wph(0,1,1))))
!
!........d P_r / dA
!
         df( 15 ) = &
            & -0.5*(C(0,0,0)**2*((-2*f(14)*f(5))/f(3)**3 - (2*CSC(f(2))**2*f(13)*f(4))/f(3)**3 + f(5)*((-2*f(14))/f(3)**3 &
            & + (CSC(f(2))**2*f(5)*(2*f(9)*dsin(f(2))**2 + 4*dcos(f(2))*f(8)*dsin(f(2))*f(3)))/f(3)**4 &
            & - (4*CSC(f(2))**2*f(5)*(2*f(9)*dsin(f(2))**2*f(3) + 2*dcos(f(2))*f(8)*dsin(f(2))*f(3)**2))/f(3)**5 &
            & - (CSC(f(2))**2*f(5)*(6*f(9)*dsin(f(2))**2*f(3)**2 + f(3)**2*(2*f(9)*dsin(f(2))**2 &
            & + 4*dcos(f(2))*f(8)*dsin(f(2))*f(3)) + 2*f(3)*(2*f(9)*dsin(f(2))**2*f(3) &
            & + 2*dcos(f(2))*f(8)*dsin(f(2))*f(3)**2)))/f(3)**6 + (6*CSC(f(2))**2*f(5)*(2*f(9)*dsin(f(2))**2*f(3)**3 &
            & + f(3)**2*(2*f(9)*dsin(f(2))**2*f(3) + 2*dcos(f(2))*f(8)*dsin(f(2))*f(3)**2)))/f(3)**7) &
            & + f(4)*((-6*CSC(f(2))**2*f(9)*f(4))/f(3)**4 - (2*CSC(f(2))**2*f(13))/f(3)**3 &
            & - (CSC(f(2))**4*f(4)*(6*f(9)*dsin(f(2))**2*f(3)**2 + f(3)**2*(2*f(9)*dsin(f(2))**2 &
            & + 4*dcos(f(2))*f(8)*dsin(f(2))*f(3)) + 2*f(3)*(2*f(9)*dsin(f(2))**2*f(3) &
            & + 2*dcos(f(2))*f(8)*dsin(f(2))*f(3)**2)))/f(3)**6 + (6*CSC(f(2))**4*f(4)*(2*f(9)*dsin(f(2))**2*f(3)**3 &
            & + f(3)**2*(2*f(9)*dsin(f(2))**2*f(3) + 2*dcos(f(2))*f(8)*dsin(f(2))*f(3)**2)))/f(3)**7)) + 2*C(0,0,0)*(2*f(15)*f(6) &
            & + (f(14)*f(5))/f(3)**2 + (CSC(f(2))**2*f(13)*f(4))/f(3)**2 + f(5)*(f(14)/f(3)**2 &
            & + (CSC(f(2))**2*f(5)*(2*f(9)*dsin(f(2))**2*f(3) + 2*dcos(f(2))*f(8)*dsin(f(2))*f(3)**2))/f(3)**4 &
            & - (CSC(f(2))**2*f(5)*(2*f(9)*dsin(f(2))**2*f(3)**3 + f(3)**2*(2*f(9)*dsin(f(2))**2*f(3) &
            & + 2*dcos(f(2))*f(8)*dsin(f(2))*f(3)**2)))/f(3)**6) + f(4)*((2*CSC(f(2))**2*f(9)*f(4))/f(3)**3 &
            & + (CSC(f(2))**2*f(13))/f(3)**2 - (CSC(f(2))**4*f(4)*(2*f(9)*dsin(f(2))**2*f(3)**3 &
            & + f(3)**2*(2*f(9)*dsin(f(2))**2*f(3) + 2*dcos(f(2))*f(8)*dsin(f(2))*f(3)**2)))/f(3)**6))*C(0,0,1) &
            & + 2*C(0,0,0)*((-2*f(5)**2)/f(3)**3 - (2*CSC(f(2))**2*f(4)**2)/f(3)**3)*(f(7)*C(1,0,0) + f(8)*C(0,1,0) &
            & + f(9)*C(0,0,1)) + 2*(f(6)**2 + f(5)**2/f(3)**2 + (CSC(f(2))**2*f(4)**2)/f(3)**2)*C(0,0,1)*(f(7)*C(1,0,0) &
            & + f(8)*C(0,1,0) + f(9)*C(0,0,1)) - 2*(-(f(6)*Wr(0,0,1)) - f(5)*Wth(0,0,1) - f(4)*Wph(0,0,1))*(-(f(15)*Wr(0,0,0)) &
            & - f(14)*Wth(0,0,0) - f(13)*Wph(0,0,0) - f(6)*(f(7)*Wr(1,0,0) + f(8)*Wr(0,1,0) + f(9)*Wr(0,0,1)) &
            & - f(5)*(f(7)*Wth(1,0,0) + f(8)*Wth(0,1,0) + f(9)*Wth(0,0,1)) - f(4)*(f(7)*Wph(1,0,0) + f(8)*Wph(0,1,0) &
            & + f(9)*Wph(0,0,1))) + 2*C(0,0,0)*(f(6)**2 + f(5)**2/f(3)**2 + (CSC(f(2))**2*f(4)**2)/f(3)**2)*(f(7)*C(1,0,1) &
            & + f(8)*C(0,1,1) + f(9)*C(0,0,2)) - 2*(1 - f(6)*Wr(0,0,0) - f(5)*Wth(0,0,0) - f(4)*Wph(0,0,0))*(-(f(15)*Wr(0,0,1)) &
            & - f(14)*Wth(0,0,1) - f(13)*Wph(0,0,1) - f(6)*(f(7)*Wr(1,0,1) + f(8)*Wr(0,1,1) + f(9)*Wr(0,0,2)) &
            & - f(5)*(f(7)*Wth(1,0,1) + f(8)*Wth(0,1,1) + f(9)*Wth(0,0,2)) - f(4)*(f(7)*Wph(1,0,1) + f(8)*Wph(0,1,1) &
            & + f(9)*Wph(0,0,2))))
!
!........d P_ph / dE
!
         df( 16 ) = &
            & -0.5*(2*C(0,0,0)*(2*f(18)*f(6) + (f(17)*f(5))/f(3)**2 + (CSC(f(2))**2*f(16)*f(4))/f(3)**2 + f(5)*(f(17)/f(3)**2 &
            & + (CSC(f(2))**2*f(5)*(2*f(12)*dsin(f(2))**2*f(3) + 2*dcos(f(2))*f(11)*dsin(f(2))*f(3)**2))/f(3)**4 &
            & - (CSC(f(2))**2*f(5)*(2*f(12)*dsin(f(2))**2*f(3)**3 + f(3)**2*(2*f(12)*dsin(f(2))**2*f(3) &
            & + 2*dcos(f(2))*f(11)*dsin(f(2))*f(3)**2)))/f(3)**6) + f(4)*((2*CSC(f(2))**2*f(12)*f(4))/f(3)**3 &
            & + (CSC(f(2))**2*f(16))/f(3)**2 - (CSC(f(2))**4*f(4)*(2*f(12)*dsin(f(2))**2*f(3)**3 &
            & + f(3)**2*(2*f(12)*dsin(f(2))**2*f(3) + 2*dcos(f(2))*f(11)*dsin(f(2))*f(3)**2)))/f(3)**6))*C(1,0,0) &
            & + 2*(f(6)**2 + f(5)**2/f(3)**2 + (CSC(f(2))**2*f(4)**2)/f(3)**2)*C(1,0,0)*(f(10)*C(1,0,0) + f(11)*C(0,1,0) &
            & + f(12)*C(0,0,1)) - 2*(-(f(6)*Wr(1,0,0)) - f(5)*Wth(1,0,0) - f(4)*Wph(1,0,0))*(-(f(18)*Wr(0,0,0)) &
            & - f(17)*Wth(0,0,0) - f(16)*Wph(0,0,0) - f(6)*(f(10)*Wr(1,0,0) + f(11)*Wr(0,1,0) + f(12)*Wr(0,0,1)) &
            & - f(5)*(f(10)*Wth(1,0,0) + f(11)*Wth(0,1,0) + f(12)*Wth(0,0,1)) - f(4)*(f(10)*Wph(1,0,0) + f(11)*Wph(0,1,0) &
            & + f(12)*Wph(0,0,1))) + 2*C(0,0,0)*(f(6)**2 + f(5)**2/f(3)**2 + (CSC(f(2))**2*f(4)**2)/f(3)**2)*(f(10)*C(2,0,0) &
            & + f(11)*C(1,1,0) + f(12)*C(1,0,1)) - 2*(1 - f(6)*Wr(0,0,0) - f(5)*Wth(0,0,0) &
            & - f(4)*Wph(0,0,0))*(-(f(18)*Wr(1,0,0)) - f(17)*Wth(1,0,0) - f(16)*Wph(1,0,0) - f(6)*(f(10)*Wr(2,0,0) &
            & + f(11)*Wr(1,1,0) + f(12)*Wr(1,0,1)) - f(5)*(f(10)*Wth(2,0,0) + f(11)*Wth(1,1,0) + f(12)*Wth(1,0,1)) &
            & - f(4)*(f(10)*Wph(2,0,0) + f(11)*Wph(1,1,0) + f(12)*Wph(1,0,1))))
!
!........d P_th / dE
!
         df( 17 ) = &
            & -0.5*(C(0,0,0)**2*((-2*Cot(f(2))*CSC(f(2))**2*f(16)*f(4))/f(3)**2 &
            & + f(5)*((-2*Cot(f(2))*CSC(f(2))**2*f(5)*(2*f(12)*dsin(f(2))**2*f(3) &
            & + 2*dcos(f(2))*f(11)*dsin(f(2))*f(3)**2))/f(3)**4 + (2*Cot(f(2))*CSC(f(2))**2*f(5)*(2*f(12)*dsin(f(2))**2*f(3)**3 &
            & + f(3)**2*(2*f(12)*dsin(f(2))**2*f(3) + 2*dcos(f(2))*f(11)*dsin(f(2))*f(3)**2)))/f(3)**6 &
            & + (CSC(f(2))**2*f(5)*(4*dcos(f(2))*f(12)*dsin(f(2))*f(3) + f(11)*(2*dcos(f(2))**2*f(3)**2 &
            & - 2*dsin(f(2))**2*f(3)**2)))/f(3)**4 - (CSC(f(2))**2*f(5)*(4*dcos(f(2))*f(12)*dsin(f(2))*f(3)**3 &
            & + f(3)**2*(4*dcos(f(2))*f(12)*dsin(f(2))*f(3) + f(11)*(2*dcos(f(2))**2*f(3)**2 &
            & - 2*dsin(f(2))**2*f(3)**2))))/f(3)**6) + f(4)*((-4*Cot(f(2))*CSC(f(2))**2*f(12)*f(4))/f(3)**3 &
            & - (2*Cot(f(2))*CSC(f(2))**2*f(16))/f(3)**2 + (4*Cot(f(2))*CSC(f(2))**4*f(4)*(2*f(12)*dsin(f(2))**2*f(3)**3 &
            & + f(3)**2*(2*f(12)*dsin(f(2))**2*f(3) + 2*dcos(f(2))*f(11)*dsin(f(2))*f(3)**2)))/f(3)**6 &
            & - (CSC(f(2))**4*f(4)*(4*dcos(f(2))*f(12)*dsin(f(2))*f(3)**3 + f(3)**2*(4*dcos(f(2))*f(12)*dsin(f(2))*f(3) &
            & + f(11)*(2*dcos(f(2))**2*f(3)**2 - 2*dsin(f(2))**2*f(3)**2))))/f(3)**6)) + 2*C(0,0,0)*(2*f(18)*f(6) &
            & + (f(17)*f(5))/f(3)**2 + (CSC(f(2))**2*f(16)*f(4))/f(3)**2 + f(5)*(f(17)/f(3)**2 &
            & + (CSC(f(2))**2*f(5)*(2*f(12)*dsin(f(2))**2*f(3) + 2*dcos(f(2))*f(11)*dsin(f(2))*f(3)**2))/f(3)**4 &
            & - (CSC(f(2))**2*f(5)*(2*f(12)*dsin(f(2))**2*f(3)**3 + f(3)**2*(2*f(12)*dsin(f(2))**2*f(3) &
            & + 2*dcos(f(2))*f(11)*dsin(f(2))*f(3)**2)))/f(3)**6) + f(4)*((2*CSC(f(2))**2*f(12)*f(4))/f(3)**3 &
            & + (CSC(f(2))**2*f(16))/f(3)**2 - (CSC(f(2))**4*f(4)*(2*f(12)*dsin(f(2))**2*f(3)**3 &
            & + f(3)**2*(2*f(12)*dsin(f(2))**2*f(3) + 2*dcos(f(2))*f(11)*dsin(f(2))*f(3)**2)))/f(3)**6))*C(0,1,0) &
            & - (4*C(0,0,0)*Cot(f(2))*CSC(f(2))**2*f(4)**2*(f(10)*C(1,0,0) + f(11)*C(0,1,0) + f(12)*C(0,0,1)))/f(3)**2 &
            & + 2*(f(6)**2 + f(5)**2/f(3)**2 + (CSC(f(2))**2*f(4)**2)/f(3)**2)*C(0,1,0)*(f(10)*C(1,0,0) + f(11)*C(0,1,0) &
            & + f(12)*C(0,0,1)) - 2*(-(f(6)*Wr(0,1,0)) - f(5)*Wth(0,1,0) - f(4)*Wph(0,1,0))*(-(f(18)*Wr(0,0,0)) &
            & - f(17)*Wth(0,0,0) - f(16)*Wph(0,0,0) - f(6)*(f(10)*Wr(1,0,0) + f(11)*Wr(0,1,0) + f(12)*Wr(0,0,1)) &
            & - f(5)*(f(10)*Wth(1,0,0) + f(11)*Wth(0,1,0) + f(12)*Wth(0,0,1)) - f(4)*(f(10)*Wph(1,0,0) + f(11)*Wph(0,1,0) &
            & + f(12)*Wph(0,0,1))) + 2*C(0,0,0)*(f(6)**2 + f(5)**2/f(3)**2 + (CSC(f(2))**2*f(4)**2)/f(3)**2)*(f(10)*C(1,1,0) &
            & + f(11)*C(0,2,0) + f(12)*C(0,1,1)) - 2*(1 - f(6)*Wr(0,0,0) - f(5)*Wth(0,0,0) &
            & - f(4)*Wph(0,0,0))*(-(f(18)*Wr(0,1,0)) - f(17)*Wth(0,1,0) - f(16)*Wph(0,1,0) - f(6)*(f(10)*Wr(1,1,0) &
            & + f(11)*Wr(0,2,0) + f(12)*Wr(0,1,1)) - f(5)*(f(10)*Wth(1,1,0) + f(11)*Wth(0,2,0) + f(12)*Wth(0,1,1)) &
            & - f(4)*(f(10)*Wph(1,1,0) + f(11)*Wph(0,2,0) + f(12)*Wph(0,1,1))))
!
!........d P_r / dE
!
         df( 18 ) = &
            & -0.5*(C(0,0,0)**2*((-2*f(17)*f(5))/f(3)**3 - (2*CSC(f(2))**2*f(16)*f(4))/f(3)**3 + f(5)*((-2*f(17))/f(3)**3 &
            & + (CSC(f(2))**2*f(5)*(2*f(12)*dsin(f(2))**2 + 4*dcos(f(2))*f(11)*dsin(f(2))*f(3)))/f(3)**4 &
            & - (4*CSC(f(2))**2*f(5)*(2*f(12)*dsin(f(2))**2*f(3) + 2*dcos(f(2))*f(11)*dsin(f(2))*f(3)**2))/f(3)**5 &
            & - (CSC(f(2))**2*f(5)*(6*f(12)*dsin(f(2))**2*f(3)**2 + f(3)**2*(2*f(12)*dsin(f(2))**2 &
            & + 4*dcos(f(2))*f(11)*dsin(f(2))*f(3)) + 2*f(3)*(2*f(12)*dsin(f(2))**2*f(3) &
            & + 2*dcos(f(2))*f(11)*dsin(f(2))*f(3)**2)))/f(3)**6 + (6*CSC(f(2))**2*f(5)*(2*f(12)*dsin(f(2))**2*f(3)**3 &
            & + f(3)**2*(2*f(12)*dsin(f(2))**2*f(3) + 2*dcos(f(2))*f(11)*dsin(f(2))*f(3)**2)))/f(3)**7) &
            & + f(4)*((-6*CSC(f(2))**2*f(12)*f(4))/f(3)**4 - (2*CSC(f(2))**2*f(16))/f(3)**3 &
            & - (CSC(f(2))**4*f(4)*(6*f(12)*dsin(f(2))**2*f(3)**2 + f(3)**2*(2*f(12)*dsin(f(2))**2 &
            & + 4*dcos(f(2))*f(11)*dsin(f(2))*f(3)) + 2*f(3)*(2*f(12)*dsin(f(2))**2*f(3) &
            & + 2*dcos(f(2))*f(11)*dsin(f(2))*f(3)**2)))/f(3)**6 + (6*CSC(f(2))**4*f(4)*(2*f(12)*dsin(f(2))**2*f(3)**3 &
            & + f(3)**2*(2*f(12)*dsin(f(2))**2*f(3) + 2*dcos(f(2))*f(11)*dsin(f(2))*f(3)**2)))/f(3)**7)) &
            & + 2*C(0,0,0)*(2*f(18)*f(6) + (f(17)*f(5))/f(3)**2 + (CSC(f(2))**2*f(16)*f(4))/f(3)**2 + f(5)*(f(17)/f(3)**2 &
            & + (CSC(f(2))**2*f(5)*(2*f(12)*dsin(f(2))**2*f(3) + 2*dcos(f(2))*f(11)*dsin(f(2))*f(3)**2))/f(3)**4 &
            & - (CSC(f(2))**2*f(5)*(2*f(12)*dsin(f(2))**2*f(3)**3 + f(3)**2*(2*f(12)*dsin(f(2))**2*f(3) &
            & + 2*dcos(f(2))*f(11)*dsin(f(2))*f(3)**2)))/f(3)**6) + f(4)*((2*CSC(f(2))**2*f(12)*f(4))/f(3)**3 &
            & + (CSC(f(2))**2*f(16))/f(3)**2 - (CSC(f(2))**4*f(4)*(2*f(12)*dsin(f(2))**2*f(3)**3 &
            & + f(3)**2*(2*f(12)*dsin(f(2))**2*f(3) + 2*dcos(f(2))*f(11)*dsin(f(2))*f(3)**2)))/f(3)**6))*C(0,0,1) &
            & + 2*C(0,0,0)*((-2*f(5)**2)/f(3)**3 - (2*CSC(f(2))**2*f(4)**2)/f(3)**3)*(f(10)*C(1,0,0) + f(11)*C(0,1,0) &
            & + f(12)*C(0,0,1)) + 2*(f(6)**2 + f(5)**2/f(3)**2 + (CSC(f(2))**2*f(4)**2)/f(3)**2)*C(0,0,1)*(f(10)*C(1,0,0) &
            & + f(11)*C(0,1,0) + f(12)*C(0,0,1)) - 2*(-(f(6)*Wr(0,0,1)) - f(5)*Wth(0,0,1) &
            & - f(4)*Wph(0,0,1))*(-(f(18)*Wr(0,0,0)) - f(17)*Wth(0,0,0) - f(16)*Wph(0,0,0) - f(6)*(f(10)*Wr(1,0,0) &
            & + f(11)*Wr(0,1,0) + f(12)*Wr(0,0,1)) - f(5)*(f(10)*Wth(1,0,0) + f(11)*Wth(0,1,0) + f(12)*Wth(0,0,1)) &
            & - f(4)*(f(10)*Wph(1,0,0) + f(11)*Wph(0,1,0) + f(12)*Wph(0,0,1))) + 2*C(0,0,0)*(f(6)**2 + f(5)**2/f(3)**2 &
            & + (CSC(f(2))**2*f(4)**2)/f(3)**2)*(f(10)*C(1,0,1) + f(11)*C(0,1,1) + f(12)*C(0,0,2)) - 2*(1 - f(6)*Wr(0,0,0) &
            & - f(5)*Wth(0,0,0) - f(4)*Wph(0,0,0))*(-(f(18)*Wr(0,0,1)) - f(17)*Wth(0,0,1) - f(16)*Wph(0,0,1) &
            & - f(6)*(f(10)*Wr(1,0,1) + f(11)*Wr(0,1,1) + f(12)*Wr(0,0,2)) - f(5)*(f(10)*Wth(1,0,1) + f(11)*Wth(0,1,1) &
            & + f(12)*Wth(0,0,2)) - f(4)*(f(10)*Wph(1,0,1) + f(11)*Wph(0,1,1) + f(12)*Wph(0,0,2))))
!
      ELSE
!
          df( 7:18 ) = 0.d0
!
      END IF
!
!.....d t / d S = dX / dS * P 
!
!!    df( 19 ) = SUM (  f( 4:6 ) * df( 1:3 )  )
      df( 19 ) = df(1)*f(4) + df(2)*f(5) + df(3)*f(6)
!!        df( 19 ) = dsqrt( (df(1)*f(4))**2 + (df(2)*f(5))**2 + (df(3)*f(6))**2 )
!
      RETURN
!
   END SUBROUTINE GET_EOM_DF
!
! =========================================================================================


! =========================================================================================
!
!.......SUBROUTINE INITIALIZE_EOM 
!
!       Get derivatives of equation system for ray equations: f0( 1:6 ) and df0( 1:6 ),
!       and the jacobian if set: f0( 7:18 ) and df0( 7:18 )
!
!****************************
   SUBROUTINE INITIALIZE_EOM ( f0, azi, elev, jacobian, atmo_specs, ierr )
!****************************
!
      USE TOOLS_ATMO
      USE MATH, ONLY: SPH2LLH;
      USE MATH, ONLY: TWOPIRANGEF;
      USE MATH, ONLY: CSC;
      USE MATH, ONLY: COT;
      IMPLICIT NONE
!    
!.....dummy variables
!
      DOUBLE PRECISION, DIMENSION ( 19 ), INTENT ( INOUT )            :: f0
      DOUBLE PRECISION                  , INTENT ( IN    )            :: azi, elev
      LOGICAL                           , INTENT ( IN    )            :: jacobian
      integer(kind=4)                       , INTENT ( OUT   )            :: ierr
      DOUBLE PRECISION, DIMENSION ( 6 ) , INTENT ( OUT   ), OPTIONAL  :: atmo_specs
!
!.....local variables
!
      DOUBLE PRECISION                              :: lat1, dlon1, h1, p_den 
      DOUBLE PRECISION, DIMENSION ( 0:2, 0:2, 0:2 ) :: C, Wph, Wth, Wr, P, Rho
!
! ---
!
!.....get atmospheric values
!          
      CALL SPH2LLH (  TWOPIRANGEF ( f0(1) ), f0(2), f0(3), lat1, dlon1, h1 ) ! spherical to lonlat
      CALL ATMO_GET ( dlon1, lat1, h1, .TRUE., ierr, C, Wph, Wth, Wr, P, Rho ) ! get atmospheric specifiations
!
      IF ( ierr .NE. 0 ) RETURN
!
!.....if requested, return atmo specs (no derivs!)
!
      IF (  PRESENT ( atmo_specs )  ) atmo_specs = (/ C(0,0,0), Wph(0,0,0), Wth(0,0,0), Wr(0,0,0), P(0,0,0), Rho(0,0,0)  /)
!          
!.....wind conversion to spherical
!
       Wph = Wph / f0(3) / dsin(f0(2))  ! wind in azimuth direction (rad/s)
       Wth = Wth / f0(3)
!
!.....slowness denominator
!
      p_den  = 1 / abs(  Wph(0,0,0) * dsin(elev) * dcos(azi) * ( f0(3) * dsin(f0(2)) ) &
         & + Wth(0,0,0) * dsin(elev) * dsin(azi) * f0(3) + Wr(0,0,0) * dcos(elev) + C(0,0,0)  )
!!      p_den = 1 / C(0,0,0)
!
!.....d P_ph / dS
!
      f0( 4 ) = dsin(elev) * dcos(azi) * p_den * f0(3) * dsin(f0(2))
!
!.....d P_th / dS
!
      f0( 5 ) = dsin(elev) * dsin(azi) * p_den * f0(3)
!
!.....d P_r / dS
!
      f0( 6 ) = dcos(elev) * p_den
!
!.....jacobian initial conditions for position are equal to zero!
!
      f0 ( 7:12 ) = 0.d0
!
      IF ( jacobian ) THEN
!
!........d P_ph / dA
!
         f0( 13 ) = &
            & -((dsin(azi)*dsqrt(f0(3)**2*dsin(f0(2))**2)*dsin(elev)*(-(dsqrt(f0(3)**2)*dsin(azi)*dsin(elev)*Wth(0,0,0)) &
            & + dcos(azi)*dsqrt(f0(3)**2*dsin(f0(2))**2)*dsin(elev)*Wph(0,0,0)))/(C(0,0,0) &
            & + dsqrt(f0(3)**2)*dcos(azi)*dsin(elev)*Wth(0,0,0) &
            & + dsin(azi)*dsqrt(f0(3)**2*dsin(f0(2))**2)*dsin(elev)*Wph(0,0,0))**2) &
            & + (dcos(azi)*dsqrt(f0(3)**2*dsin(f0(2))**2)*dsin(elev))/(C(0,0,0) &
            & + dsqrt(f0(3)**2)*dcos(azi)*dsin(elev)*Wth(0,0,0) &
            & + dsin(azi)*dsqrt(f0(3)**2*dsin(f0(2))**2)*dsin(elev)*Wph(0,0,0))
!
!........d P_th / dA
!
         f0( 14 ) = &
            & -((dsqrt(f0(3)**2)*dcos(azi)*dsin(elev)*(-(dsqrt(f0(3)**2)*dsin(azi)*dsin(elev)*Wth(0,0,0)) &
            & + dcos(azi)*dsqrt(f0(3)**2*dsin(f0(2))**2)*dsin(elev)*Wph(0,0,0)))/(C(0,0,0) &
            & + dsqrt(f0(3)**2)*dcos(azi)*dsin(elev)*Wth(0,0,0) &
            & + dsin(azi)*dsqrt(f0(3)**2*dsin(f0(2))**2)*dsin(elev)*Wph(0,0,0))**2) &
            & - (dsqrt(f0(3)**2)*dsin(azi)*dsin(elev))/(C(0,0,0) &
            & + dsqrt(f0(3)**2)*dcos(azi)*dsin(elev)*Wth(0,0,0) &
            & + dsin(azi)*dsqrt(f0(3)**2*dsin(f0(2))**2)*dsin(elev)*Wph(0,0,0))
!
!........d P_r / dA
!
         f0( 15 ) = &
            & ((dcos(elev)*(-(dsqrt(f0(3)**2)*dsin(azi)*dsin(elev)*Wth(0,0,0)) &
            & + dcos(azi)*dsqrt(f0(3)**2*dsin(f0(2))**2)*dsin(elev)*Wph(0,0,0)))/(C(0,0,0) &
            & + dsqrt(f0(3)**2)*dcos(azi)*dsin(elev)*Wth(0,0,0) &
            & + dsin(azi)*dsqrt(f0(3)**2*dsin(f0(2))**2)*dsin(elev)*Wph(0,0,0))**2)
!
!........d P_ph / dE
!
         f0( 16 ) = &
            & -((dsin(azi)*dsqrt(f0(3)**2*dsin(f0(2))**2)*dsin(elev)*(dsqrt(f0(3)**2)*dcos(azi)*dcos(elev)*Wth(0,0,0) &
            & + dcos(elev)*dsin(azi)*dsqrt(f0(3)**2*dsin(f0(2))**2)*Wph(0,0,0)))/(C(0,0,0) &
            & + dsqrt(f0(3)**2)*dcos(azi)*dsin(elev)*Wth(0,0,0) &
            & + dsin(azi)*dsqrt(f0(3)**2*dsin(f0(2))**2)*dsin(elev)*Wph(0,0,0))**2) &
            & + (dcos(elev)*dsin(azi)*dsqrt(f0(3)**2*dsin(f0(2))**2))/(C(0,0,0) &
            & + dsqrt(f0(3)**2)*dcos(azi)*dsin(elev)*Wth(0,0,0) &
            & + dsin(azi)*dsqrt(f0(3)**2*dsin(f0(2))**2)*dsin(elev)*Wph(0,0,0))
!
!........d P_th / dE
!
         f0( 17 ) = &
            & -((dsqrt(f0(3)**2)*dcos(azi)*dsin(elev)*(dsqrt(f0(3)**2)*dcos(azi)*dcos(elev)*Wth(0,0,0) &
            & + dcos(elev)*dsin(azi)*dsqrt(f0(3)**2*dsin(f0(2))**2)*Wph(0,0,0)))/(C(0,0,0) &
            & + dsqrt(f0(3)**2)*dcos(azi)*dsin(elev)*Wth(0,0,0) &
            & + dsin(azi)*dsqrt(f0(3)**2*dsin(f0(2))**2)*dsin(elev)*Wph(0,0,0))**2) &
            & + (dsqrt(f0(3)**2)*dcos(azi)*dcos(elev))/(C(0,0,0) &
            & + dsqrt(f0(3)**2)*dcos(azi)*dsin(elev)*Wth(0,0,0) &
            & + dsin(azi)*dsqrt(f0(3)**2*dsin(f0(2))**2)*dsin(elev)*Wph(0,0,0))
!
!........d P_r / dE
!
         f0( 18 ) = &
            & -((dcos(elev)*(dsqrt(f0(3)**2)*dcos(azi)*dcos(elev)*Wth(0,0,0) &
            & + dcos(elev)*dsin(azi)*dsqrt(f0(3)**2*dsin(f0(2))**2)*Wph(0,0,0)))/(C(0,0,0) &
            & + dsqrt(f0(3)**2)*dcos(azi)*dsin(elev)*Wth(0,0,0) &
            & + dsin(azi)*dsqrt(f0(3)**2*dsin(f0(2))**2)*dsin(elev)*Wph(0,0,0))**2) &
            & - dsin(elev)/(C(0,0,0) + dsqrt(f0(3)**2)*dcos(azi)*dsin(elev)*Wth(0,0,0) &
            & + dsin(azi)*dsqrt(f0(3)**2*dsin(f0(2))**2)*dsin(elev)*Wph(0,0,0))
!
      ELSE
!
          f0( 13:18 ) = 0.d0
!
      END IF
!
!....time starts at zero!
!
     f0( 19 ) = 0.d0
!
      RETURN
!
   END SUBROUTINE INITIALIZE_EOM
!
! =========================================================================================


! =========================================================================================
!
!.......SUBROUTINE INITIALIZE_SPREADING 
!
!       Get spreading at step S, a sphere at 1m or 1km, respectively, depending on the RE 1 unit.
!
!****************************
   SUBROUTINE INITIALIZE_SPREADING ( re_unit, A0 )
!****************************
!
      USE SYMS
      IMPLICIT NONE
!    
!.....dummy variables
!
      CHARACTER ( LEN = * )             , INTENT ( IN  )  :: re_unit
      DOUBLE PRECISION                  , INTENT ( OUT )  :: A0
!
! ---
!
!.....Transmission loss start unit
!
      SELECT CASE ( re_unit )
         CASE ('km')
            A0 = 1.d0 / ( 4 * PI * 1.d-3 )
         CASE DEFAULT
            A0 = 1.d0 / ( 4 * PI * 1.d-6 )
      END SELECT
!
      RETURN
!
   END SUBROUTINE INITIALIZE_SPREADING
!
! =========================================================================================


! =========================================================================================
!
!.......SUBROUTINE GET_JACOBIAN_DETERMINANT
!
!       Get spreading in terms of the Jacobian determinant at time step t.
!
!****************************
   SUBROUTINE GET_JACOBIAN_DETERMINANT ( f, J )
!****************************
!
      IMPLICIT NONE
!    
!.....dummy variables
!
      DOUBLE PRECISION, DIMENSION ( 19 ), INTENT ( IN  )  :: f
      DOUBLE PRECISION                  , INTENT ( OUT )  :: J
!
!.....Derivative structure
!
!        f(  1 ) = d X_ph / dS
!        f(  2 ) = d X_th / dS
!        f(  3 ) = d X_r  / dS
!
!        f(  4 ) = d P_ph / dS
!        f(  5 ) = d P_th / dS
!        f(  6 ) = d P_r  / dS
!
!        f(  7 ) = d X_ph / dA
!        f(  8 ) = d X_th / dA
!        f(  9 ) = d X_r  / dA
!
!        f( 10 ) = d X_ph / dE
!        f( 11 ) = d X_th / dE
!        f( 12 ) = d X_r  / dE
!
!        f( 13 ) = d P_ph / dA
!        f( 14 ) = d P_th / dA
!        f( 15 ) = d P_r  / dA
!
!        f( 16 ) = d P_ph / dE
!        f( 17 ) = d P_th / dE
!        f( 18 ) = d P_r  / dE
!
! ---
!
!.....Jacobian (determinant) from df. See Jensen et. al 2011, p. 166. 
!
      J = f( 1 ) * f( 11 ) * f ( 9 ) + f( 10 ) * f( 8 ) * f( 3 ) + f( 2 ) * f( 12 ) * f( 7 ) &
         & - f( 3 ) * f( 11 ) * f( 7 ) - f( 2 ) * f( 10 ) * f( 9 ) - f( 12 ) * f( 8 ) * f( 1 )
!
      RETURN
!
   END SUBROUTINE GET_JACOBIAN_DETERMINANT
!
! =========================================================================================


! =========================================================================================
!
!.......SUBROUTINE GET_JACOBIAN 
!
!       Get spreading in terms of the Jacobian determinant at time step t.
!
!****************************
   SUBROUTINE GET_JACOBIAN ( f, J, nofCaustics, lzero )
!****************************
!
      IMPLICIT NONE
!    
!.....dummy variables
!
      DOUBLE PRECISION, DIMENSION ( 19 ), INTENT ( IN    )  :: f
      integer(kind=4)                   , INTENT ( INOUT )  :: nofCaustics
      DOUBLE PRECISION                  , INTENT ( INOUT )  :: J
      LOGICAL                           , INTENT ( OUT   )  :: lzero
!
!.....local variables
!
      DOUBLE PRECISION :: J0
!
! ---
!
!.....copy previous determinant
!
      J0 = J
!
!.....Jacobian (determinant) from df. See Jensen et. al 2011, p. 166. 
!
      call get_jacobian_determinant( f, J )
!
!.....Zero determinant --> caustic!
!
      IF ( abs(J).lt.1.d-8 ) THEN
         nofCaustics = nofCaustics + 1
         lzero = .TRUE.
      ELSE IF ( J0.LT.0.d0 .AND. J.GT.0.d0 &
         & .OR. J0.GT.0.d0 .AND. J.LT.0.d0 ) THEN
         nofCaustics = nofCaustics + 1
         lzero = .FALSE.
      ELSE
         lzero = .FALSE.
      END IF
!
      RETURN
!
   END SUBROUTINE GET_JACOBIAN
!
! =========================================================================================

! =========================================================================================
!
!.......SUBROUTINE GET_JACOBIAN_OLD 
!
!       Get spreading in terms of the Jacobian determinant at time step t.
!
!****************************
   SUBROUTINE GET_JACOBIAN_OLD ( f, d, azi, J, nofCaustics, lzero )
!****************************
!
      USE SYMS
      IMPLICIT NONE
!    
!.....dummy variables
!
      DOUBLE PRECISION, DIMENSION ( 19 ), INTENT ( IN    )  :: f
      DOUBLE PRECISION                  , INTENT ( IN    )  :: d, azi
      integer(kind=4)                   , INTENT ( INOUT )  :: nofCaustics
      DOUBLE PRECISION                  , INTENT ( INOUT )  :: J
      LOGICAL                           , INTENT ( OUT   )  :: lzero
!
!.....local variables
!
      DOUBLE PRECISION :: J0
!
!.....Derivative structure
!
!        f(  1 ) = d X_ph / dS
!        f(  2 ) = d X_th / dS
!        f(  3 ) = d X_r  / dS
!
!        f(  4 ) = d P_ph / dS
!        f(  5 ) = d P_th / dS
!        f(  6 ) = d P_r  / dS
!
!        f(  7 ) = d X_ph / dA
!        f(  8 ) = d X_th / dA
!        f(  9 ) = d X_r  / dA
!
!        f( 10 ) = d X_ph / dE
!        f( 11 ) = d X_th / dE
!        f( 12 ) = d X_r  / dE
!
!        f( 13 ) = d P_ph / dA
!        f( 14 ) = d P_th / dA
!        f( 15 ) = d P_r  / dA
!
!        f( 16 ) = d P_ph / dE
!        f( 17 ) = d P_th / dE
!        f( 18 ) = d P_r  / dE
!
! ---
!
!.....copy previous determinant
!
      J0 = J
!
!.....Weir implimentation from past 
!
      J = - d * dsin( azi ) * (  f( 11 ) * f( 3 ) - f( 12 ) * f( 2 )  ) &
         &  - d * dcos( azi ) * (  f( 10 ) * f( 3 ) - f( 12 ) * f( 1 )  )
!
!.....Zero determinant --> caustic!
!
      IF ( abs(J).lt.1.d-8 ) THEN
         nofCaustics = nofCaustics + 1
         lzero = .TRUE.
      ELSE IF ( J0.LT.0.d0 .AND. J.GT.0.d0 &
         & .OR. J0.GT.0.d0 .AND. J.LT.0.d0 ) THEN
         nofCaustics = nofCaustics + 1
         lzero = .FALSE.
      ELSE
         lzero = .FALSE.
      END IF
!
      RETURN
!
   END SUBROUTINE GET_JACOBIAN_OLD
!
! =========================================================================================


! =========================================================================================
!
!.......SUBROUTINE GET_SPREADING 
!
!       Get amplitude due to spreading spreading at step S, based on the Jacobian determinant,
!       density times speed of sound and the elevation angle.
!
!****************************
   SUBROUTINE GET_SPREADING ( elev, rho0, c0, rho, c, nofCaustics, det, A )
!****************************
!
      USE SYMS
      IMPLICIT NONE
!    
!.....dummy variables
!
      DOUBLE PRECISION                 , INTENT ( IN  )  :: elev, rho0, c0, rho, c, det
      integer(kind=4)                  , INTENT ( IN  )  :: nofCaustics
      DOUBLE PRECISION                 , INTENT ( OUT )  :: A
!
! ---
!
!.....Spreading per spherical area ( devided by 1 / 4 / PI )
!
      A = ( - 1.d0 ) ** nofCaustics * &
         & dsqrt( abs( ( dsin(elev)*rho*c )/( rho0*c0*det ) ) )/(4*pi)
!
      RETURN
!
   END SUBROUTINE GET_SPREADING
!
! =========================================================================================


! =========================================================================================
!
!.......SUBROUTINE GET_TLOSS 
!
!       Get spreading at step S, based on the Jacobian determinant, density times speed of
!       sound and the elevation angle.
!
!****************************
   SUBROUTINE GET_TLOSS ( A0, A, absorption, frq, t, TL )
!****************************
!
      USE TOOLS_ATMO
      USE SYMS
      IMPLICIT NONE
!    
!.....dummy variables
!
      DOUBLE PRECISION, INTENT ( IN  )  :: A0, A, absorption, frq, t
      DOUBLE PRECISION, INTENT ( OUT )  :: TL
!
!.....local variables
!
      complex(kind=8)  :: p
!
! ---
!
!.....Transmission Loss
!
      p  = A / A0 * exp( cmplx( 0.d0, twopi*frq*t, kind=8 ) )
      TL = - 20.d0 * log10( abs(p) ) !+ absorption
!
      RETURN
!
   END SUBROUTINE GET_TLOSS
!
! =========================================================================================
!
END MODULE EOM
