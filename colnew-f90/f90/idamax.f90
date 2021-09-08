!*==IDAMAX.spg  processed by SPAG 6.72Dc at 06:43 on  8 Sep 2021
      FUNCTION IDAMAX(N,Dx,Incx)
 
!*********************************************************************72
!
!c IDAMAX finds the index of element having maximum absolute value.
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
!    Input, double precision X(*), the vector to be examined.
!
!    Input, integer INCX, the increment between successive entries of SX.
!
!    Output, integer IDAMAX, the index of the element of SX of maximum
!    absolute value.
!
      IMPLICIT NONE
!*--IDAMAX45
 
      DOUBLE PRECISION Dx(*) , dmax
      INTEGER IDAMAX
      INTEGER i , Incx , ix , N
 
      IDAMAX = 0
      IF ( N.LT.1 .OR. Incx.LE.0 ) RETURN
      IDAMAX = 1
      IF ( N.EQ.1 ) RETURN
      IF ( Incx.EQ.1 ) THEN
!
!  code for increment equal to 1
!
         dmax = DABS(Dx(1))
         DO i = 2 , N
            IF ( dmax.LT.DABS(Dx(i)) ) THEN
               IDAMAX = i
               dmax = DABS(Dx(i))
            ENDIF
         ENDDO
         GOTO 99999
      ENDIF
!
!  code for increment not equal to 1
!
      ix = 1
      dmax = DABS(Dx(1))
      ix = ix + Incx
      DO i = 2 , N
         IF ( dmax.LT.DABS(Dx(ix)) ) THEN
            IDAMAX = i
            dmax = DABS(Dx(ix))
         ENDIF
         ix = ix + Incx
      ENDDO
      RETURN
 
99999 END
