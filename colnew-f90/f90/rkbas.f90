!*==RKBAS.spg  processed by SPAG 6.72Dc at 06:31 on  8 Sep 2021
      SUBROUTINE RKBAS(S,Coef,K,M,Rkb,Dm,Mode)
!
!**********************************************************************
!
!c RKBAS evaluates a mesh-independent Runge-Kutta basis.
!
!   purpose
!           evaluate mesh independent runge-kutta basis for given s
!
!   variables
!     s      - argument, i.e. the relative position for which
!              the basis is to be evaluated ( 0. .le. s .le. 1. ).
!     coef   - precomputed derivatives of the basis
!     k      - number of collocatin points per subinterval
!     m      - maximal order of the differential equation
!     rkb    - the runge-kutta basis (0-th to (m-1)-th derivatives )
!     dm     - basis elements for m-th derivative
!
!**********************************************************************
!
      IMPLICIT NONE
!*--RKBAS23
!*** Start of declarations inserted by SPAG
      REAL*8 Coef , Dm , p , Rkb , S , t
      INTEGER i , j , K , kpm1 , l , lb , M , Mode
!*** End of declarations inserted by SPAG
      DIMENSION Coef(K,1) , Rkb(7,1) , Dm(1) , t(10)
!
      IF ( K.EQ.1 ) THEN
         Rkb(1,1) = 1.0D0
         Dm(1) = 1.0D0
         GOTO 99999
      ENDIF
      kpm1 = K + M - 1
      DO i = 1 , kpm1
         t(i) = S/DFLOAT(i)
      ENDDO
      DO l = 1 , M
         lb = K + l + 1
         DO i = 1 , K
            p = Coef(1,i)
            DO j = 2 , K
               p = p*t(lb-j) + Coef(j,i)
            ENDDO
            Rkb(i,l) = p
         ENDDO
      ENDDO
      IF ( Mode.EQ.0 ) RETURN
      DO i = 1 , K
         p = Coef(1,i)
         DO j = 2 , K
            p = p*t(K+1-j) + Coef(j,i)
         ENDDO
         Dm(i) = p
      ENDDO
      RETURN
99999 END
