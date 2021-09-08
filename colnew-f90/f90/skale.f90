!*==SKALE.spg  processed by SPAG 6.72Dc at 06:33 on  8 Sep 2021
      SUBROUTINE SKALE(N,Mstar,Kd,Z,Xi,Scale,Dscale)
!
!**********************************************************************
!
!c SKALE provides a scaling for the state variables.
!
!   purpose
!            provide a proper scaling of the state variables, used
!            to control the damping factor for a newton iteration [2].
!
!   variables
!
!            n      = number of mesh subintervals
!            mstar  = number of unknomns in z(u(x))
!            kd     = number of unknowns in dmz
!            z      = the global unknown vector
!            xi     = the current mesh
!            scale  = scaling vector for z
!            dscale = scaling vector for dmz
!
!**********************************************************************
!
      IMPLICIT NONE
!*--SKALE25
!*** Start of declarations inserted by SPAG
      REAL*8 basm , Dscale , h , scal , Scale , Xi , Z
      INTEGER icomp , ID1 , ID2 , idmz , iz , j , K , Kd , l , M , mj , &
            & MMAx , Mstar , N , NCOmp , np1
!*** End of declarations inserted by SPAG
      DIMENSION Z(Mstar,1) , Scale(Mstar,1) , Dscale(Kd,1)
      DIMENSION Xi(1) , basm(5)
!
      COMMON /COLORD/ K , NCOmp , ID1 , ID2 , MMAx , M(20)
!
      basm(1) = 1.D0
      DO j = 1 , N
         iz = 1
         h = Xi(j+1) - Xi(j)
         DO l = 1 , MMAx
            basm(l+1) = basm(l)*h/DFLOAT(l)
         ENDDO
         DO icomp = 1 , NCOmp
            scal = (DABS(Z(iz,j))+DABS(Z(iz,j+1)))*.5D0 + 1.D0
            mj = M(icomp)
            DO l = 1 , mj
               Scale(iz,j) = basm(l)/scal
               iz = iz + 1
            ENDDO
            scal = basm(mj+1)/scal
            DO idmz = icomp , Kd , NCOmp
               Dscale(idmz,j) = scal
            ENDDO
         ENDDO
      ENDDO
      np1 = N + 1
      DO iz = 1 , Mstar
         Scale(iz,np1) = Scale(iz,N)
      ENDDO
      END
