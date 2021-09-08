!*==SHIFTB.spg  processed by SPAG 6.72Dc at 06:32 on  8 Sep 2021
      SUBROUTINE SHIFTB(Ai,Nrowi,Ncoli,Last,Ai1,Nrowi1,Ncoli1)
      IMPLICIT NONE
!*--SHIFTB4
!*** Start of declarations inserted by SPAG
      INTEGER Ncoli , Ncoli1 , Nrowi , Nrowi1
!*** End of declarations inserted by SPAG
!
!*********************************************************************
!
!c SHIFTB shifts rows in the current block.
!
!     shifts the rows in current block, ai, not used as pivot rows, if
!     any, i.e., rows  (last+1),..., (nrowi), onto the first mmax =
!      = nrow-last  rows of the next block, ai1, with column last+j of
!      ai  going to column j , j=1,...,jmax=ncoli-last. the remaining
!     columns of these rows of ai1 are zeroed out.
!
!                                picture
!
!          original situation after         results in a new block i+1
!          last = 2 columns have been       created and ready to be
!          done in factrb (assuming no      factored by next factrb
!          interchanges of rows)            call.
!                      1
!                 x  x 1x  x  x           x  x  x  x  x
!                      1
!                 0  x 1x  x  x           0  x  x  x  x
!     block i          1                       ---------------
!     nrowi = 4   0  0 1x  x  x           0  0 1x  x  x  0  01
!     ncoli = 5        1                       1             1
!     last = 2    0  0 1x  x  x           0  0 1x  x  x  0  01
!     -------------------------------          1             1   new
!                      1x  x  x  x  x          1x  x  x  x  x1  block
!                      1                       1             1   i+1
!     block i+1        1x  x  x  x  x          1x  x  x  x  x1
!     nrowi1= 5        1                       1             1
!     ncoli1= 5        1x  x  x  x  x          1x  x  x  x  x1
!     -------------------------------          1-------------1
!                      1
!
!*********************************************************************
!
      INTEGER Last , j , jmax , jmaxp1 , m , mmax
      DOUBLE PRECISION Ai(Nrowi,Ncoli) , Ai1(Nrowi1,Ncoli1)
      mmax = Nrowi - Last
      jmax = Ncoli - Last
      IF ( mmax.LT.1 .OR. jmax.LT.1 ) RETURN
!
!...  put the remainder of block i into ai1
!
      DO j = 1 , jmax
         DO m = 1 , mmax
            Ai1(m,j) = Ai(Last+m,Last+j)
         ENDDO
      ENDDO
      IF ( jmax.EQ.Ncoli1 ) RETURN
!
!...  zero out the upper right corner of ai1
!
      jmaxp1 = jmax + 1
      DO j = jmaxp1 , Ncoli1
         DO m = 1 , mmax
            Ai1(m,j) = 0.D0
         ENDDO
      ENDDO
      END
