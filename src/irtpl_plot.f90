! =========================================================================================
!
! MODULE IRTPL_PLOT
!
! --> subroutines to read plot IRT3D output data.
!
! =========================================================================================
!
MODULE IRTPL_PLOT
!
   USE syms
!
CONTAINS
!
! =========================================================================================
!
!.......Subroutine plot_vertical 
!
!
!****************************
  SUBROUTINE MKGMT_VERTICAL ( prefix, dmax, hmax, tmax, omin, omax, resd, resh, &
    marker, marker_sym, re_unit, tstr )
!****************************
!
    IMPLICIT NONE
!          
!.....dummy variables
!
    CHARACTER ( LEN = * )                     , INTENT ( IN  )  :: prefix, re_unit, tstr
    CHARACTER ( LEN = * ), DIMENSION ( : )    , INTENT ( IN  )  :: marker_sym
    DOUBLE PRECISION                          , INTENT ( IN  )  :: dmax, hmax, tmax, omin, omax, resd, resh
    DOUBLE PRECISION     , DIMENSION ( : , : ), INTENT ( IN  )  :: marker
!          
!.....local variables
!
    integer(kind=4)       :: marker_cnt, i, munit = 20
    CHARACTER         :: str_tmax*8, str_omin*8, str_omax*8, str_hmax*8, str_dmax*8, str_resd*8
    CHARACTER         :: str_resh*8, dlabel*10, tlabel*10, olabel*10
    DOUBLE PRECISION  :: orange
!
! ----
!
!.....write parameters to string
!
    WRITE ( str_dmax, "(I8)" )   floor( dmax / 10 ) * 10
    WRITE ( str_hmax, "(F8.2)" ) hmax 
    WRITE ( str_tmax, "(I8)" )   CEILING ( tmax / 100 ) * 100
    WRITE ( str_omin, "(I8)" )   FLOOR ( omin / 10 ) * 10
    WRITE ( str_omax, "(I8)" )   CEILING ( omax / 10 ) * 10
    WRITE ( str_resd, "(F8.2)" ) resd
    WRITE ( str_resh, "(F8.2)" ) resh
!
!.....get distance label
!
    IF ( dmax .LE.   500 ) dlabel = 'a50f5    '
    IF ( dmax .GT.   500 ) dlabel = 'a100f10  '
    IF ( dmax .GT.  1000 ) dlabel = 'a250f25  '
    IF ( dmax .GT.  2000 ) dlabel = 'a500f50  '
    IF ( dmax .GT.  4000 ) dlabel = 'a1000f100'
    IF ( dmax .GT.  8000 ) dlabel = 'a2000f200'
!
!.....get time label
!
    IF ( tmax .LE.    500 ) tlabel = 'a50f10    '
    IF ( tmax .GT.    500 ) tlabel = 'a100f20   '
    IF ( tmax .GT.   1000 ) tlabel = 'a250f50   '
    IF ( tmax .GT.   2000 ) tlabel = 'a500f100  '
    IF ( tmax .GT.   4000 ) tlabel = 'a1000f200 '
    IF ( tmax .GT.   8000 ) tlabel = 'a2000f400 '
    IF ( tmax .GT.  16000 ) tlabel = 'a3000f500 '
    IF ( tmax .GT.  24000 ) tlabel = 'a4000f800 '
    IF ( tmax .GT.  32000 ) tlabel = 'a5000f1000'
!
!.....get offset label
!
    orange = omax - omin
    IF ( orange .LE.    50 ) olabel = 'a10f2  '
    IF ( orange .GT.    50 ) olabel = 'a20f4  '
    IF ( orange .GT.   100 ) olabel = 'a50f10 '
    IF ( orange .GT.   200 ) olabel = 'a50f10 '
    IF ( orange .GT.   400 ) olabel = 'a100f20'
    IF ( orange .GT.   700 ) olabel = 'a150f30'
    IF ( orange .GT.  1050 ) olabel = 'a200f40'
!
!.....get markers
!
    marker_cnt = SIZE ( marker, 1 )
!
    OPEN ( UNIT = 19, STATUS = 'replace', FILE = prefix // '_markersel.xy', FORM = 'formatted' )
!
    DO i = 1, marker_cnt
!
       IF ( abs(marker(i,6)-1.d0).lt.1.d-8 ) &
        & WRITE ( 19, "(F9.2,1x,F6.2,1x,A)" ) marker(i,3), 0.d0, TRIM (  marker_sym( i )  )
!
    ENDDO
!
    CLOSE ( 19 )
!
!.....open gmt script file
!
    OPEN ( UNIT = munit, STATUS = 'replace', FILE = prefix // '_vertical.gmt', FORM='formatted' )
!
!.....write gmtset parameters
!
    write(munit,"(a)") '#!/bin/bash'
    write(munit,"(a)") '# ***************************************************************************'
    write(munit,"(a)") '#'
    write(munit,"(a)") '#                              IRTPL GMT SCRIPT FILE'
    write(munit,"(a)") '#'
    write(munit,"(a)") '# ***************************************************************************'
    write(munit,"(a)") '#'
    write(munit,"(a)") '# GMT scripts to plot the vertical cross-section of the IRT3D raytrace data.'
    write(munit,"(a)") '#'
    write(munit,"(a)") '# v1.0 [Jan 2013, KNMI]'
    write(munit,"(a)") '#'
    write(munit,"(a)") '# ***'
    write(munit,"(a)") '#'
    write(munit,"(a)") '# GMT DEFAULTS '
    write(munit,"(a)") 'gmtset PAGE_ORIENTATION=landscape'
    write(munit,"(a)") 'gmtset ANNOT_FONT_SIZE_PRIMARY=10p ANNOT_FONT_SIZE_SECONDARY=10p'
    write(munit,"(a)") 'gmtset LABEL_FONT_SIZE=10p HEADER_FONT_SIZE=10p'
    write(munit,"(a)") 'gmtset COLOR_BACKGROUND=white COLOR_FOREGROUND=brown COLOR_NAN=gray'
!
!.....general
!
    write(munit,"(a)") '# NAMES'
    write(munit,"(a)") 'prefix="'// prefix // '"'
    write(munit,"(a)") 'ps="${prefix}_vert.ps"'
    write(munit,"(a)") 'cEF="${prefix}_ceff.cpt"'
    write(munit,"(a)") 'cTL="${prefix}_tloss.cpt"'
    write(munit,"(a)") ''
!
    write(munit,"(a)") '# GMT VARIABLES'
    write(munit,"(3a)") 'rTLoss=0/', trim(adjustl(str_dmax)), '/-30/180'
    write(munit,"(4a)") 'rTime=0/', trim(adjustl(str_dmax)), '/0/', trim(adjustl(str_tmax))
    write(munit,"(6a)") 'rOffset=0/', trim(adjustl(str_dmax)), '/', trim(adjustl(str_omin)), &
       '/', trim(adjustl(str_omax))
    write(munit,"(4a)") 'rCross=0/', trim(adjustl(str_dmax)), '/0/', trim(adjustl(str_hmax))
    write(munit,"(a)")  'scale=7.5c/-1.5c/10c/0.3ch'
    write(munit,"(2a)") 'bDist=', trim(dlabel)
    write(munit,"(2a)") 'bTime=', trim(tlabel)
    write(munit,"(2a)") 'bOffset=', trim(olabel)
!    write(munit,"(a)") ''
!
!.....make cpt file
!
    write(munit,"(a)") '# CPT FILE'
    write(munit,"(a)") 'makecpt -C$GMT_CEFF -T240/410/5 -D -M > $cEF'
    SELECT CASE ( re_unit )
       CASE ('km')
        write(munit,"(a)") 'makecpt -Cgray -T30/150/10 -D -M -Z > $cTL'
       CASE DEFAULT
        write(munit,"(a)") 'makecpt -Cgray -T90/210/10 -D -M -Z > $cTL'
    END SELECT
    write(munit,"(a)") ''
!
!.....plot 1.1 : tloss vs distance 
!
    write(munit,"(a)") '# TLOSS vs Distance'
    SELECT CASE ( re_unit )
       CASE ('km')
        write(munit,"(2a)") 'psbasemap -R${rTLoss} -JX6.5i/-1.8i -B:."":${bDist}', &
         & ':"Distance (km)":/50f10:"TL (dB Re 1km)":neWS -K -X3.0i -Y.75i > $ps' 
       CASE DEFAULT
        write(munit,"(2a)") 'psbasemap -R${rTLoss} -JX6.5i/-1.8i -B:."":${bDist}', &
         & ':"Distance (km)":/30f6:"TL (dB Re 1m)":neWS -K -X3.0i -Y.75i > $ps' 
    END SELECT
    write(munit,"(a)") 'psxy ${prefix}_rays_v_tloss.xy -R -JX -m -W0.5p,black -O -K >> $ps'
!
!.....plot 1.2 : time and offset
!
    write(munit,"(a)") '# TIME'
    write(munit,"(a)") 'gmtset BASEMAP_FRAME_RGB=red'
    write(munit,"(2a)") 'psbasemap -R$rTime -JX6.5i/1.8i -Y2.1i ', &
       '-B:."":${bDist}:"":/${bTime}:"Time (sec)":E -O -K >> $ps'
    write(munit,"(a)") 'psxy ${prefix}_time.xy -R -JX -Sc1.0p -Wred -Gred -O -K >> $ps'
!
    write(munit,"(a)") 'rm .gmtdefaults4'
    write(munit,"(a)") 'gmtset ANNOT_FONT_SIZE_PRIMARY=10p ANNOT_FONT_SIZE_SECONDARY=10p'
    write(munit,"(a)") 'gmtset LABEL_FONT_SIZE=10p HEADER_FONT_SIZE=10p'
    write(munit,"(a)") '# OFFSET'
    write(munit,"(a)") 'gmtset BASEMAP_FRAME_RGB=black'
    write(munit,"(a)") 'psbasemap -R$rOffset -JX -B:."":${bDist}:"Distance (km)":/${bOffset}:"Offset (km)":nWs -O -K >> $ps'
    write(munit,"(a)") 'psxy ${prefix}_offset.xy -R -JX -Sc1.0p -Wblack -Gblack -O -K >> $ps'
    write(munit,"(a)") ''
!
!.....plot 1.3 : cross_section
!
    write(munit,"(a)") '# CROSS-SECTION'
    write(munit,"(5a)") 'xyz2grd ${prefix}_ceff.xyz -G${prefix}_ceff.grd -R${rCross} -I', &
       trim(adjustl(str_resd)), '/', trim(adjustl(str_resh)), &
       ' -D:"distance"/"altitude"/"c_eff"/1/0/"Effective speed of sound"/='
    write(munit,"(2a)") 'grdimage ${prefix}_ceff.grd -R${rCross} -JX6.5i/2.5i', &
       ' -C$cEF -E300 -Sn/0 -Y2.1i -O -K >> $ps'
    write(munit,"(a)") 'psxy ${prefix}_rays_v.xy -R -JX -m -W+0.5p -Sp1.0p -C$cTL -O -K >> $ps'
    write(munit,"(a)") 'psxy ${prefix}_markersel.xy -R -JX -S16.0p -W0.5p,black -Gwhite -O -K >> $ps'
    write(munit,"(4a)") 'psbasemap -R -JX -B:."', trim(tstr), '":${bDist}', &
       ':"Distance (km)":/10f2:"Altitude (km)":news -O -K >> $ps'
    write(munit,"(2a)") 'psscale -D3.25i/3.15i/4.0i/.12ih -C$cEF -E', &
       ' -Ba20f4:"":/:"c@-eff@- (m/s)":neSW -O -K >> $ps'
    SELECT CASE ( re_unit )
       CASE ('km')
        write(munit,"(2a)") 'psscale -D6.8i/1.25i/-2.5i/.12i -C$cTL -E', &
           ' -Ba20f4:"TL (dB Re 1km)":/:"":neSW -O -K >> $ps'
       CASE DEFAULT
        write(munit,"(2a)") 'psscale -D6.8i/1.25i/-2.5i/.12i -C$cTL -E', &
           ' -Ba20f4:"TL (dB Re 1m)":/:"":neSW -O -K >> $ps'
    END SELECT
    write(munit,"(a)") ''
!
!.....plot 1.4 : profiles source
!
    write(munit,"(a)") '# PROFILES SOURCE'
    write(munit,"(4a)") 'psbasemap -R240/410/0/', trim(adjustl(str_hmax)), ' -JX1.3i/2.5i -X-1.7i', &
       ' -B:."":a50f10:"c@-T@- and c@-eff@- (m/s)":/10f2:"Altitude (km)":neWS -O -K >> $ps'
    write(munit,"(a)") 'psxy ${prefix}_source_c.xy -R -JX -W0.5p,gray -O -K >> $ps'
    write(munit,"(a)") 'psxy ${prefix}_source_ceff.xy -R -JX -W.5p,black -O >> $ps'
    write(munit,"(a)") ''
!
!.....plot 2 : tloss vs time 
!
!!      write(munit,"(a)") '# TLOSS vs time'
!!      write(munit,"(a)") 'psbasemap -R0/' // TRIM (  ADJUSTL (  str_tmax  )  ) // '/0/130' &
!!         & //' -JX6.5i/-3.0i -B:."":' // TRIM ( tlabel ) // ':"Time (sec)":/a20f4' &
!!         & // ':"TLOSS (dB Re 1m)":neWS -K -X2.5i > ${prefix}_vert_tloss.ps' 
!!      write(munit,"(a)") 'psxy ${prefix}_rays_v_tloss.xy -R -JX -m -W0.5p,black -O >> ' &
!!         & // prefix // '_vert_tloss.ps'
!!      write(munit,"(a)") ''
!
!.....cleanup
!
    write(munit,"(a)") '# CLEANUP FILES'
    write(munit,"(a)") 'rm *.cpt *.grd .gmtdefaults4'
!
!.....close gmt script file
!
    CLOSE ( munit )
!
!.....execute gmt script
!
    CALL SYSTEM( 'chmod 755 ' // prefix // '_vertical.gmt' )
    CALL SYSTEM( './' // prefix // '_vertical.gmt' )
!
    RETURN
!
   END SUBROUTINE MKGMT_VERTICAL
!
! =========================================================================================


! =========================================================================================
!
!.......Subroutine MKGMT_REFL
!
!
!****************************
   SUBROUTINE MKGMT_REFL ( prefix, projection, lat0, lon0, re_unit, tstr )
!****************************
!
    USE TOOLS_ATMO 
    USE IRTPL_READ
!
    IMPLICIT NONE
!          
!.....dummy variables
!
    CHARACTER ( LEN = * ), INTENT ( IN  )  :: prefix, projection, re_unit, tstr
    DOUBLE PRECISION     , INTENT ( IN  )  :: lat0, lon0
!          
!.....local variables
!
    integer(kind=4) :: munit, area, tmax
    CHARACTER   :: str_hmax*8, str_tmax*8, str_xmin*8, str_xmax*8, str_ymin*8, str_ymax*8, &
       & str_lat0*8, str_lon0*8, b_x*10, b_y*10, tlabel*10, coast*100, &
       & str_proj*30, str_range*50
!
! init
!
  munit=22
  tmax=1000
  area=0
!
  coast=''
!
!.....write parameters to string
!
    IF (sc%lonmin .LE. sc%lonmax ) WRITE ( str_xmin, "(F8.2)" ) sc%lonmin
    IF (sc%lonmin .GT. sc%lonmax ) WRITE ( str_xmin, "(F8.2)" ) sc%lonmin - 360.d0
    WRITE ( str_xmax, "(F8.2)" ) sc%lonmax
    IF ( sc%latmin .GT. -85.d0 ) WRITE ( str_ymin, "(F8.2)" ) sc%latmin
    IF ( sc%latmin .LE. -85.d0 ) WRITE ( str_ymin, "(F8.2)" ) -85.d0
    IF ( sc%latmax .LT.  85.d0 ) WRITE ( str_ymax, "(F8.2)" ) sc%latmax
    IF ( sc%latmax .GE.  85.d0 ) WRITE ( str_ymax, "(F8.2)" )  85.d0
    WRITE ( str_hmax, "(F8.2)" ) atmo%hmax
    IF (  ALLOCATED( store_refl )  ) tmax = CEILING(  MAXVAL( store_refl( :, 2 ), 1 )  )
    WRITE ( str_tmax, "(I8)" ) tmax
    WRITE ( str_lat0, "(F8.2)" ) lat0
    WRITE ( str_lon0, "(F8.2)" ) lon0
!
!.....get lat label
!
    IF ( sc%latrange .LE.   7 ) b_y = 'a1g.5 '
    IF ( sc%latrange .GT.   7 ) b_y = 'a2g1  '
    IF ( sc%latrange .GT.  15 ) b_y = 'a5g2.5'
    IF ( sc%latrange .GT.  30 ) b_y = 'a10g5 '
    IF ( sc%latrange .GT.  60 ) b_y = 'a20g10'
    IF ( sc%latrange .GT. 120 ) b_y = 'a30g15'
!
!.....get lon label
!
    IF ( sc%lonrange .LE.   7 ) b_x = 'a1g.5 '
    IF ( sc%lonrange .GT.   7 ) b_x = 'a2g1  '
    IF ( sc%lonrange .GT.  15 ) b_x = 'a5g2.5'
    IF ( sc%lonrange .GT.  30 ) b_x = 'a10g5 '
    IF ( sc%lonrange .GT.  60 ) b_x = 'a20g10'
    IF ( sc%lonrange .GT. 120 ) b_x = 'a30g15'
    IF ( sc%lonrange .GT. 180 ) b_x = 'a40g20'
    IF ( sc%lonrange .GT. 240 ) b_x = 'a60g30'
!
!.....set projection
!
    SELECT CASE ( projection )
       CASE ( 'sn' )
        str_proj  = 'S0/90/15c'
        IF ( sc%latmin .LT. 0.d0 ) WRITE ( str_ymin, "(F8.2)" ) 0.d0
        str_range = '0/360/' // TRIM ( ADJUSTL ( str_ymin ) ) // '/90'
        b_x       = 'a30g30'
       CASE ( 'ss' )
        str_proj  = 'S0/-90/15c'
        IF ( sc%latmax .GT. 0.d0 ) WRITE ( str_ymax, "(F8.2)" ) 0.d0
        str_range = '0/360/-90' // TRIM ( ADJUSTL ( str_ymax ) )
        b_x       = 'a30g30'
       CASE ( 'h' )
        str_range = 'd'
        str_proj  = 'W15c'
       CASE ( 'm' )
        str_proj  = 'M15c'
        str_range = TRIM (  ADJUSTL ( str_xmin )  ) // '/' // TRIM (  ADJUSTL ( str_xmax )  ) &
         & // '/' // TRIM (  ADJUSTL ( str_ymin )  ) // '/' // TRIM (  ADJUSTL ( str_ymax )  )
       CASE DEFAULT
        str_proj  = 'M15c'
        str_range = TRIM (  ADJUSTL ( str_xmin )  ) // '/' // TRIM (  ADJUSTL ( str_xmax )  ) &
         & // '/' // TRIM (  ADJUSTL ( str_ymin )  ) // '/' // TRIM (  ADJUSTL ( str_ymax )  )
    END SELECT
!
!.....get time label
!
    IF ( tmax .LE.    500 ) tlabel = 'a50f10    '
    IF ( tmax .GT.    500 ) tlabel = 'a100f20   '
    IF ( tmax .GT.   1000 ) tlabel = 'a250f50   '
    IF ( tmax .GT.   2000 ) tlabel = 'a500f100  '
    IF ( tmax .GT.   4000 ) tlabel = 'a1000f200 '
    IF ( tmax .GT.   8000 ) tlabel = 'a2000f400 '
    IF ( tmax .GT.  16000 ) tlabel = 'a3000f500 '
    IF ( tmax .GT.  24000 ) tlabel = 'a4000f800 '
    IF ( tmax .GT.  32000 ) tlabel = 'a5000f1000'
!
!.....set level of detail of coast, rivers and borders
!
    area = INT(  ( sc%lonmax - sc%lonmin ) * ( sc%latmax - sc%latmin ), KIND = 4  )
!
    IF ( area .LE.  500 ) coast = '-Df -Ia/thinner,white -Na/thinner,gray50'
    IF ( area .GT.  500 ) coast = '-Dh -Ir/thinner,white -A250 -N1/thin,gray50 -N2/thinner,gray75 -A250'
    IF ( area .GT. 10000 ) coast = '-Di -I1/thinner,white -I2/thinner,white -A500 -N1/thinner,gray50 -A500'
    IF ( area .GT. 10000 ) coast = '-Dl -A1000'
!
!.....open gmt script file
!
    OPEN ( UNIT = munit, STATUS = 'replace', FILE = prefix // '_refl.gmt', FORM = 'formatted' )
!
!.....write gmtset parameters
!
    write(munit,"(a)") '#!/bin/bash'
    write(munit,"(a)") '# ***************************************************************************'
    write(munit,"(a)") '#'
    write(munit,"(a)") '#                              IRTPL GMT SCRIPT FILE'
    write(munit,"(a)") '#'
    write(munit,"(a)") '# ***************************************************************************'
    write(munit,"(a)") '#'
    write(munit,"(a)") '# GMT scripts to plot the reflection height and time of the IRT3D raytrace data.'
    write(munit,"(a)") '#'
    write(munit,"(a)") '# v1.0 [Jan 2013, KNMI]'
    write(munit,"(a)") '#'
    write(munit,"(a)") '# ***'
    write(munit,"(a)") '#'
    write(munit,"(a)") '# GMT DEFAULTS '
    write(munit,"(a)") 'gmtset PAGE_ORIENTATION=portrait'
    write(munit,"(a)") 'gmtset ANNOT_FONT_SIZE_PRIMARY=10p ANNOT_FONT_SIZE_SECONDARY=10p'
    write(munit,"(a)") 'gmtset LABEL_FONT_SIZE=10p HEADER_FONT_SIZE=10p'
    write(munit,"(a)") 'gmtset BASEMAP_TYPE=FANCY PLOT_DEGREE_FORMAT=ddd:mm:ssF GRID_CROSS_SIZE_PRIMARY=0.05i'
    write(munit,"(a)") 'gmtset COLOR_NAN=gray'
    write(munit,"(a)") ''
!
!.....general
!
    write(munit,"(a)") '# NAMES'
    write(munit,"(a)") 'prefix="'// prefix // '"'
    write(munit,"(a)") 'cH="${prefix}_refl_h.cpt"'
    write(munit,"(a)") 'cT="${prefix}_refl_t.cpt"'
    write(munit,"(a)") 'cTL="${prefix}_refl_tloss.cpt"'
    write(munit,"(a)") ''
!
    write(munit,"(a)") '# GMT VARIABLES'
    write(munit,"(2a)") 'range=',trim( str_range )
    write(munit,"(2a)") 'proj=',trim( str_proj )
    write(munit,"(3a)") 'coast="',trim( coast ),'"'
    write(munit,"(a)")  'scale=7.5c/-1.5c/10c/0.3ch'
    write(munit,"(a)") ''
!
!.....make cpt file
!
    write(munit,"(a)") '# CPT FILE'
    write(munit,"(a)") 'makecpt -Crainbow -T0/' //trim( adjustl(str_hmax) )// '/5 -D -Z > $cH'
    write(munit,"(a)") 'makecpt -Crainbow -T0/' //trim( adjustl(str_tmax) )// '/250 -D -Z > $cT'
    SELECT CASE ( re_unit )
       CASE ('km')
        write(munit,"(a)") 'makecpt -C$GMT_WAXLER -T30/180/10 -I -D -M -Z > $cTL'
       CASE DEFAULT
        write(munit,"(a)") 'makecpt -C$GMT_WAXLER -T30/180/10 -I -D -M -Z > $cTL'
    END SELECT
    write(munit,"(a)") ''
!
!.....plot 1 : reflection height
!
    write(munit,"(a)") '# REFLECTION HEIGHT'
    write(munit,"(a)") 'ps="${prefix}_refl_h.ps"'
    write(munit,"(a)") 'pscoast -R$range -J$proj $coast -Sgray85 -Wthinnest,gray60 -Gwhite -Xc -Yc -K > $ps'
    write(munit,"(a)") 'psxy ${prefix}_rays_h.xy -R -J -m -W0.5p,black -O -K >> $ps'
    write(munit,"(a)") 'psxy ${prefix}_refl_h.xy -R -J -Sc2.5p -C$cH -O -K >> $ps'
    write(munit,"(a)") 'echo "' // trim(adjustl(str_lon0)) // ' ' // trim(adjustl(str_lat0)) &
       & // '" | psxy -R -J -Sa10p -Wgray30 -Ggray70 -O -K >> $ps'
    write(munit,"(a)") 'psxy ${prefix}_markers.xy -R -J -S8p -Wgray30 -Ggray70 -O -K'&
       & // ' -B:."' //trim(tstr)// ' - Refraction height":' // trim(b_x) // ':"":/' // trim(b_y) // ':"":neSW >> $ps'
    write(munit,"(a)") 'psscale -D$scale -C$cH -Ba10f2:"":/:"h (km)":neSW -O >> $ps'
    write(munit,"(a)") ''
!
!.....plot 2 : travel time
!
    write(munit,"(a)") '# TRAVEL TIME'
    write(munit,"(a)") 'ps="${prefix}_refl_t.ps"'
    write(munit,"(a)") 'pscoast -R$range -J$proj $coast -Sgray85 -Wthinnest,gray60 -Gwhite -Xc -Yc -K > $ps'
    write(munit,"(a)") 'psxy ${prefix}_rays_h.xy -R -J -m -W0.5p,black -O -K >> $ps'
    write(munit,"(a)") 'psxy ${prefix}_refl_t.xy -R -J -Sc2.5p -C$cT -O -K >> $ps'
    write(munit,"(a)") 'echo "' // TRIM (  ADJUSTL ( str_lon0 )  ) // ' ' // TRIM (  ADJUSTL ( str_lat0 )  ) &
       & // '" | psxy -R -J -Sa10p -Wgray30 -Ggray70 -O -K >> $ps'
    write(munit,"(a)") 'psxy ${prefix}_markers.xy -R -J -S8p -Wgray30 -Ggray70 -O -K' &
       & // ' -B:."' //trim(tstr)// ' - Travel time":' // TRIM(b_x) // ':"":/' // TRIM(b_y) // ':"":neSW >> $ps'
    write(munit,"(a)") 'psscale -D$scale -C$cT -B' // trim(tlabel) // ':"":/:"t (s)":neSW -O >> $ps'
    write(munit,"(a)") ''
!
!.....plot 3 : tloss
!
    write(munit,"(a)") '# Transmission Loss'
    write(munit,"(a)") 'ps="${prefix}_refl_tloss.ps"'
    write(munit,"(a)") 'pscoast -R$range -J$proj $coast -Sgray85 -Wthinnest,gray60 -Gwhite -Xc -Yc -K > $ps'
    write(munit,"(a)") 'psxy ${prefix}_rays_h.xy -R -J -m -W0.5p,black -O -K >> $ps'
    write(munit,"(a)") 'psxy ${prefix}_refl_tloss.xy -R -J -Sc2.5p -C$cTL -O -K >> $ps'
    write(munit,"(a)") 'echo "' // TRIM (  ADJUSTL ( str_lon0 )  ) // ' ' // TRIM (  ADJUSTL ( str_lat0 )  ) &
       & // '" | psxy -R -J -Sa10p -Wgray30 -Ggray70 -O -K >> $ps'
    write(munit,"(a)") 'psxy ${prefix}_markers.xy -R -J -S8p -Wgray30 -Ggray70 -O -K' &
       & // ' -B:."' //trim(tstr)// ' - Transmission Loss":' // trim(b_x) // ':"":/' // trim(b_y) // ':"":neSW >> $ps'
    SELECT CASE ( re_unit )
       CASE ('km')
         write(munit,"(a)") 'psscale -D$scale -C$cTL ' &
          & // '-Ba10f2:"":/:"TL (dB Re 1km)":neSW -E -O >> $ps'
       CASE DEFAULT
        write(munit,"(a)") 'psscale -D$scale -C$cTL ' &
         & // '-Ba10f2:"":/:"TL (dB Re 1m)":neSW -E -O >> $ps'
    END SELECT
    write(munit,"(a)") ''
!
!.....cleanup
!
    write(munit,"(a)") '# CLEANUP FILES'
    write(munit,"(a)") 'rm *.cpt'
!
!.....close gmt script file
!
    CLOSE ( munit )
!
!.....execute gmt script
!
    CALL SYSTEM( 'chmod 755 ' // prefix // '_refl.gmt' )
    CALL SYSTEM( './' // prefix // '_refl.gmt' )
!
    RETURN
!
   END SUBROUTINE MKGMT_REFL
!
! =========================================================================================
!
END MODULE IRTPL_PLOT
