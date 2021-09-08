!*==TIMESTAMP.spg  processed by SPAG 6.72Dc at 06:35 on  8 Sep 2021
      SUBROUTINE TIMESTAMP()
 
!*********************************************************************72
!
!c TIMESTAMP prints out the current YMDHMS date as a timestamp.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
      IMPLICIT NONE
!*--TIMESTAMP25
 
      CHARACTER*(8) ampm
      INTEGER d
      CHARACTER*(8) date
      INTEGER h
      INTEGER m
      INTEGER mm
      CHARACTER*(9) month(12)
      INTEGER n
      INTEGER s
      CHARACTER*(10) time
      INTEGER y
 
      SAVE month
 
      DATA month/'January  ' , 'February ' , 'March    ' , 'April    ' ,&
          &'May      ' , 'June     ' , 'July     ' , 'August   ' ,      &
          &'September' , 'October  ' , 'November ' , 'December '/
 
      CALL DATE_AND_TIME(date,time)
 
      READ (date,'(i4,i2,i2)') y , m , d
      READ (time,'(i2,i2,i2,1x,i3)') h , n , s , mm
 
      IF ( h.LT.12 ) THEN
         ampm = 'AM'
      ELSEIF ( h.EQ.12 ) THEN
         IF ( n.EQ.0 .AND. s.EQ.0 ) THEN
            ampm = 'Noon'
         ELSE
            ampm = 'PM'
         ENDIF
      ELSE
         h = h - 12
         IF ( h.LT.12 ) THEN
            ampm = 'PM'
         ELSEIF ( h.EQ.12 ) THEN
            IF ( n.EQ.0 .AND. s.EQ.0 ) THEN
               ampm = 'Midnight'
            ELSE
               ampm = 'AM'
            ENDIF
         ENDIF
      ENDIF
 
      WRITE (*,'(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)') d ,&
           & month(m) , y , h , ':' , n , ':' , s , '.' , mm , ampm
 
      END
