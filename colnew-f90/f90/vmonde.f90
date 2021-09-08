!*==VMONDE.spg  processed by SPAG 6.72Dc at 06:35 on  8 Sep 2021
      SUBROUTINE VMONDE(Rho,Coef,K)
      IMPLICIT NONE
!*--VMONDE4
!
!**********************************************************************
!
!c VMONDE solves a Vandermonde linear system.
!
!   purpose
!          solve vandermonde system v * x = e
!          with  v(i,j) = rho(j)**(i-1)/(i-1)! .
!
!**********************************************************************
!
      INTEGER K , i , ifac , j , km1 , kmi
      DOUBLE PRECISION Rho(K) , Coef(K)
!
      IF ( K.EQ.1 ) RETURN
      km1 = K - 1
      DO i = 1 , km1
         kmi = K - i
         DO j = 1 , kmi
            Coef(j) = (Coef(j+1)-Coef(j))/(Rho(j+i)-Rho(j))
         ENDDO
      ENDDO
!
      ifac = 1
      DO i = 1 , km1
         kmi = K + 1 - i
         DO j = 2 , kmi
            Coef(j) = Coef(j) - Rho(j+i-1)*Coef(j-1)
         ENDDO
         Coef(kmi) = DFLOAT(ifac)*Coef(kmi)
         ifac = ifac*i
      ENDDO
      Coef(1) = DFLOAT(ifac)*Coef(1)
      END
