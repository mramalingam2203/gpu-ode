!*==FCBLOK.spg  processed by SPAG 6.72Dc at 07:07 on  8 Sep 2021
!----------------------------------------------------------------------
!                            p a r t  5
!          we list here a modified (column oriented, faster)
!          version of the package solveblok of de boor - weiss [5].
!          we also give a listing of the linpack
!          routines dgefa und dgesl used by colnew.
!----------------------------------------------------------------------
!
      SUBROUTINE FCBLOK(Bloks,Integs,Nbloks,Ipivot,Scrtch,Info)
      IMPLICIT NONE
!*--FCBLOK12
!*** Start of declarations inserted by SPAG
      INTEGER indexx , Nbloks
!*** End of declarations inserted by SPAG
 
!**********************************************************************
!
!c FCBLOK calls subroutines FACTRB and SHIFTB.
!
!     fcblok  supervises the plu factorization with pivoting of
!     scaled rows of the almost block diagonal matrix stored in the
!     arrays  bloks  and  integs .
!
!     factrb = subprogram which carries out steps 1,...,last of gauss
!            elimination (with pivoting) for an individual block.
!     shiftb = subprogram which shifts the remaining rows to the top of
!            the next block
!
!     parameters
!      bloks   an array that initially contains the almost block diago-
!            nal matrix  a  to be factored, and on return contains the
!            computed factorization of  a .
!      integs  an integer array describing the block structure of  a .
!      nbloks  the number of blocks in  a .
!      ipivot  an integer array of dimension   sum (integs(3,n) ; n=1,
!            ...,nbloks) which, on return, contains the pivoting stra-
!            tegy used.
!      scrtch  work area required, of length  max (integs(1,n) ; n=1,
!            ...,nbloks).
!      info    output parameter;
!            = 0  in case matrix was found to be nonsingular.
!            otherwise,
!            = n if the pivot element in the nth gauss step is zero.
!
!**********************************************************************
!
      INTEGER Integs(3,Nbloks) , Ipivot(1) , Info , i , index , indexn ,&
            & last , ncol , nrow
      DOUBLE PRECISION Bloks(1) , Scrtch(1)
      Info = 0
      indexx = 1
      indexn = 1
      i = 1
!
!...  loop over the blocks.  i  is loop index
!
 100  index = indexn
      nrow = Integs(1,i)
      ncol = Integs(2,i)
      last = Integs(3,i)
!
!...       carry out elimination on the i-th block until next block
!...       enters, i.e., for columns 1,...,last  of i-th block.
!
      CALL FACTRB(Bloks(index),Ipivot(indexx),Scrtch,nrow,ncol,last,    &
                & Info)
!
!...       check for having reached a singular block or the last block
!
      IF ( Info.NE.0 ) THEN
         Info = Info + indexx - 1
      ELSE
         IF ( i.EQ.Nbloks ) RETURN
         i = i + 1
         indexn = nrow*ncol + index
         indexx = indexx + last
!
!...       put the rest of the i-th block onto the next block
!
         CALL SHIFTB(Bloks(index),nrow,ncol,last,Bloks(indexn),         &
                   & Integs(1,i),Integs(2,i))
         GOTO 100
      ENDIF
      END
