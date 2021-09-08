!*==DMZSOL.spg  processed by SPAG 6.72Dc at 05:13 on  8 Sep 2021
      SUBROUTINE DMZSOL(Kd,Mstar,N,V,Z,Dmz)
!
!**********************************************************************
!
!c DMZSOL computes DMZ in a blockwise manner.
!
!   purpose
!          compute dmz in a blockwise manner
!          dmz(i) = dmz(i)  +  v(i) * z(i), i = 1,...,n
!
!**********************************************************************
!
      IMPLICIT NONE
!*--DMZSOL15
!*** Start of declarations inserted by SPAG
      REAL*8 Dmz , fact , V , Z
      INTEGER i , j , jz , Kd , l , Mstar , N
!*** End of declarations inserted by SPAG
      DIMENSION V(Kd,1) , Dmz(Kd,1) , Z(1)
!
      jz = 1
      DO i = 1 , N
         DO j = 1 , Mstar
            fact = Z(jz)
            DO l = 1 , Kd
               Dmz(l,i) = Dmz(l,i) + fact*V(l,jz)
            ENDDO
            jz = jz + 1
         ENDDO
      ENDDO
      END
