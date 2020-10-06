!
!.....read additional markers (lat, lon, symbol)
!
      marker_cnt = 0
      marker     = 0.d0
!
      DO
!
!........read record
!
         READ ( 1,'(A)',IOSTAT = io ) dummy
!
!........error reading file
!
         IF (io .GT. 0)  THEN
!
           PRINT "(2A)", 'ERROR: Something wrong with the file ', &
              & TRIM( in_inputfile )
           STOP
!
!........end of file reached
!
         ELSE IF ( io .LT. 0 ) THEN
!
            EXIT
!
!........record ok
!
         ELSE
!
!...........marker count up
!
            marker_cnt = marker_cnt + 1
!
            IF ( marker_cnt .GT. 99 ) THEN
               PRINT "(2A)", 'WARNING: maximum of 999 markers reached. ', &
                  & 'Remaining markers will be ignored.'
               EXIT
            ENDIF
!
!...........lat
!
            i = SCAN( dummy( 1 : 256 ), ' ', .FALSE. )
            READ ( dummy( 1 : i ), * ) marker( marker_cnt, 1 )
            dummy = dummy(i+1:256)
!
!...........lon
!
            i = SCAN( dummy( 1 : 256 ), ' ', .FALSE. )
            READ ( dummy( 1 : i ), * ) marker( marker_cnt, 2 )
            dummy = dummy(i+1:256)
!
            CALL LONRANGE( marker( marker_cnt,2 ) )
!
!...........symbol
!
            READ ( dummy, * ) marker_sym( marker_cnt )
!
         ENDIF
!
      ENDDO
!
