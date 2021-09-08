!*==DDOT.spg  processed by SPAG 6.72Dc at 05:11 on  8 Sep 2021
      FUNCTION DDOT(N,Dx,Incx,Dy,Incy)
 
!*********************************************************************72
!
!c DDOT forms the dot product of two vectors.
!
!  Discussion:
!
!    This routine uses double precision real arithmetic.
!
!    This routine uses unrolled loops for increments equal to one.
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
!    Input, integer N, the number of entries in the vectors.
!
!    Input, double precision DX(*), the first vector.
!
!    Input, integer INCX, the increment between successive entries in DX.
!
!    Input, double precision DY(*), the second vector.
!
!    Input, integer INCY, the increment between successive entries in DY.
!
!    Output, double precision DDOT, the sum of the product of the
!    corresponding entries of DX and DY.
!
      IMPLICIT NONE
!*--DDOT51
 
      DOUBLE PRECISION DDOT
      DOUBLE PRECISION Dx(*)
      DOUBLE PRECISION Dy(*)
      DOUBLE PRECISION dtemp
      INTEGER i , Incx , Incy , ix , iy , m , N
 
      DDOT = 0.0D0
      dtemp = 0.0D0
      IF ( N.LE.0 ) RETURN
      IF ( Incx.EQ.1 .AND. Incy.EQ.1 ) THEN
!
!  code for both increments equal to 1
!
!
!  clean-up loop
!
         m = MOD(N,5)
         IF ( m.NE.0 ) THEN
            DO i = 1 , m
               dtemp = dtemp + Dx(i)*Dy(i)
            ENDDO
            IF ( N.LT.5 ) THEN
 
               DDOT = dtemp
               GOTO 99999
            ENDIF
         ENDIF
         DO i = m + 1 , N , 5
            dtemp = dtemp + Dx(i)*Dy(i) + Dx(i+1)*Dy(i+1) + Dx(i+2)     &
     &              *Dy(i+2) + Dx(i+3)*Dy(i+3) + Dx(i+4)*Dy(i+4)
         ENDDO
         DDOT = dtemp
      ELSE
!
!  code for unequal increments or equal increments not equal to 1
!
         ix = 1
         iy = 1
         IF ( Incx.LT.0 ) ix = (-N+1)*Incx + 1
         IF ( Incy.LT.0 ) iy = (-N+1)*Incy + 1
         DO i = 1 , N
            dtemp = dtemp + Dx(ix)*Dy(iy)
            ix = ix + Incx
            iy = iy + Incy
         ENDDO
         DDOT = dtemp
         RETURN
      ENDIF
 
99999 END
