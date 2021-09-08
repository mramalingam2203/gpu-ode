!*==DAXPY.spg  processed by SPAG 6.72Dc at 05:10 on  8 Sep 2021
      SUBROUTINE DAXPY(N,Da,Dx,Incx,Dy,Incy)
 
!*********************************************************************72
!
!c DAXPY computes constant times a vector plus a vector.
!
!  Discussion:
!
!    This routine uses double precision real arithmetic.
!
!    This routine uses unrolled loops for increments equal to one.
!
!  Modified:
!
!    18 December 2008
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
!    Input, integer N, the number of elements in DX and DY.
!
!    Input, double precision DA, the multiplier of DX.
!
!    Input, double precision DX(*), the first vector.
!
!    Input, integer INCX, the increment between successive entries of DX.
!
!    Input/output, double precision DY(*), the second vector.
!    On output, DY(*) has been replaced by DY(*) + DA * DX(*).
!
!    Input, integer INCY, the increment between successive entries of DY.
!
      IMPLICIT NONE
!*--DAXPY51
 
      DOUBLE PRECISION Da
      DOUBLE PRECISION Dx(*)
      DOUBLE PRECISION Dy(*)
      INTEGER i
      INTEGER Incx
      INTEGER Incy
      INTEGER ix
      INTEGER iy
      INTEGER m
      INTEGER N
 
      IF ( N.LE.0 ) RETURN
 
      IF ( Da.EQ.0.0D0 ) RETURN
 
      IF ( Incx.NE.1 .OR. Incy.NE.1 ) THEN
 
         ix = 1
         iy = 1
         IF ( Incx.LT.0 ) ix = (-N+1)*Incx + 1
         IF ( Incy.LT.0 ) iy = (-N+1)*Incy + 1
 
         DO i = 1 , N
            Dy(iy) = Dy(iy) + Da*Dx(ix)
            ix = ix + Incx
            iy = iy + Incy
         ENDDO
 
      ELSE
 
         m = MOD(N,4)
 
         DO i = 1 , m
            Dy(i) = Dy(i) + Da*Dx(i)
         ENDDO
 
         DO i = m + 1 , N , 4
            Dy(i) = Dy(i) + Da*Dx(i)
            Dy(i+1) = Dy(i+1) + Da*Dx(i+1)
            Dy(i+2) = Dy(i+2) + Da*Dx(i+2)
            Dy(i+3) = Dy(i+3) + Da*Dx(i+3)
         ENDDO
 
      ENDIF
 
      END
