!*==HORDER.spg  processed by SPAG 6.72Dc at 06:42 on  8 Sep 2021
      SUBROUTINE HORDER(I,Uhigh,Hi,Dmz,Ncomp,K)
!
!**********************************************************************
!
!c HORDER determines the highest order derivatives of the current solution.
!
!   purpose
!           determine highest order (piecewise constant) derivatives
!           of the current collocation solution
!
!   variables
!     hi     - the stepsize, hi = xi(i+1) - xi(i)
!     dmz    - vector of mj-th derivative of the solution
!     uhigh  - the array of highest order (piecewise constant)
!              derivatives of the approximate solution on
!              (xi(i),xi(i+1)), viz,
!                          (k+mj-1)
!              uhigh(j) = u   (x)    on (xi(i),xi(i+1))
!                          j
!
!**********************************************************************
!
      IMPLICIT NONE
!*--HORDER25
!*** Start of declarations inserted by SPAG
      REAL*8 COEf , Dmz , dn , fact , Hi , RHO , Uhigh
      INTEGER I , id , idmz , j , K , kin , Ncomp
!*** End of declarations inserted by SPAG
      DIMENSION Uhigh(1) , Dmz(1)
!
      COMMON /COLLOC/ RHO(7) , COEf(49)
!
      dn = 1.D0/Hi**(K-1)
!
!...  loop over the ncomp solution components
!
      DO id = 1 , Ncomp
         Uhigh(id) = 0.D0
      ENDDO
      kin = 1
      idmz = (I-1)*K*Ncomp + 1
      DO j = 1 , K
         fact = dn*COEf(kin)
         DO id = 1 , Ncomp
            Uhigh(id) = Uhigh(id) + fact*Dmz(idmz)
            idmz = idmz + 1
         ENDDO
         kin = kin + K
      ENDDO
      END
