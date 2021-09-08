!*==SBBLOK.spg  processed by SPAG 6.72Dc at 06:31 on  8 Sep 2021
      SUBROUTINE SBBLOK(Bloks,Integs,Nbloks,Ipivot,X)
      IMPLICIT NONE
!*--SBBLOK4
!*** Start of declarations inserted by SPAG
      INTEGER Nbloks
!*** End of declarations inserted by SPAG
!
!**********************************************************************
!
!c SBBLOK calls SUBFOR and SUBBAK so solve a block linear system.
!
!     calls subroutines  subfor  and  subbak .
!
!     supervises the solution (by forward and backward substitution) of
!     the linear system  a*x = b  for x, with the plu factorization of
!     a  already generated in  fcblok .  individual blocks of
!     equations are solved via  subfor  and  subbak .
!
!    parameters
!       bloks, integs, nbloks, ipivot    are as on return from fcblok.
!       x       on input: the right hand side, in dense storage
!               on output: the solution vector
!
!*********************************************************************
!
      INTEGER Integs(3,Nbloks) , Ipivot(1) , i , index , indexx , j ,   &
            & last , nbp1 , ncol , nrow
      DOUBLE PRECISION Bloks(1) , X(1)
!
!...  forward substitution pass
!
      index = 1
      indexx = 1
      DO i = 1 , Nbloks
         nrow = Integs(1,i)
         last = Integs(3,i)
         CALL SUBFOR(Bloks(index),Ipivot(indexx),nrow,last,X(indexx))
         index = nrow*Integs(2,i) + index
         indexx = indexx + last
      ENDDO
!
!...  back substitution pass
!
      nbp1 = Nbloks + 1
      DO j = 1 , Nbloks
         i = nbp1 - j
         nrow = Integs(1,i)
         ncol = Integs(2,i)
         last = Integs(3,i)
         index = index - nrow*ncol
         indexx = indexx - last
         CALL SUBBAK(Bloks(index),nrow,ncol,last,X(indexx))
      ENDDO
      END
