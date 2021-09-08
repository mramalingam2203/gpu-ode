!*==DSCAL.spg  processed by SPAG 6.72Dc at 05:14 on  8 Sep 2021
      SUBROUTINE DSCAL(N,Da,Dx,Incx)
 
!*********************************************************************72
!
!c DSCAL scales a vector by a constant.
!
!  Discussion:
!
!    This routine uses double precision real arithmetic.
!
!  Modified:
!
!    07 July 2007
!
!  Author:
!
!    Jack Dongarra
!
!  Reference:
!
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Basic Linear Algebra Subprograms for FORTRAN usage,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, pages 308-323, 1979.
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input, double precision SA, the multiplier.
!
!    Input/output, double precision X(*), the vector to be scaled.
!
!    Input, integer INCX, the increment between successive entries of X.
!
      IMPLICIT NONE
!*--DSCAL44
 
      DOUBLE PRECISION Da
      DOUBLE PRECISION Dx(*)
      INTEGER i , Incx , m , N , nincx
 
      IF ( N.LE.0 .OR. Incx.LE.0 ) RETURN
      IF ( Incx.EQ.1 ) THEN
!
!  code for increment equal to 1
!
!
!  clean-up loop
!
         m = MOD(N,5)
         IF ( m.NE.0 ) THEN
            DO i = 1 , m
               Dx(i) = Da*Dx(i)
            ENDDO
            IF ( N.LT.5 ) RETURN
         ENDIF
         DO i = m + 1 , N , 5
            Dx(i) = Da*Dx(i)
            Dx(i+1) = Da*Dx(i+1)
            Dx(i+2) = Da*Dx(i+2)
            Dx(i+3) = Da*Dx(i+3)
            Dx(i+4) = Da*Dx(i+4)
         ENDDO
      ELSE
!
!  code for increment not equal to 1
!
         nincx = N*Incx
         DO i = 1 , nincx , Incx
            Dx(i) = Da*Dx(i)
         ENDDO
         RETURN
      ENDIF
 
      END
