!*==FACTRB.spg  processed by SPAG 6.72Dc at 07:06 on  8 Sep 2021
      SUBROUTINE FACTRB(W,Ipivot,D,Nrow,Ncol,Last,Info)
      IMPLICIT NONE
!*--FACTRB4
!*** Start of declarations inserted by SPAG
      INTEGER Nrow
!*** End of declarations inserted by SPAG
!
!********************************************************************
!
!c FACTRB constructs a partial PLU factorization of a matrix.
!
!     adapted from p.132 of  element.numer.analysis  by conte-de boor
!
!     constructs a partial plu factorization, corresponding to steps
!      1,..., last   in gauss elimination, for the matrix  w  of
!      order ( nrow ,  ncol ), using pivoting of scaled rows.
!
!     parameters
!       w       contains the (nrow,ncol) matrix to be partially factored
!               on input, and the partial factorization on output.
!       ipivot  an integer array of length last containing a record of
!               the pivoting strategy used; explicit interchanges
!               are used for pivoting.
!       d       a work array of length nrow used to store row sizes
!               temporarily.
!       nrow    number of rows of w.
!       ncol    number of columns of w.
!       last    number of elimination steps to be carried out.
!       info    on output, zero if the matrix is found to be non-
!               singular, in case a zero pivot was encountered in row
!               n,  info = n on output.
!
!**********************************************************************
!
      INTEGER Ipivot(Nrow) , Ncol , Last , Info , i , j , k , l , kp1
      DOUBLE PRECISION W(Nrow,Ncol) , D(Nrow) , colmax , t , s
      DOUBLE PRECISION DABS , DMAX1
!
!...  initialize  d
!
      DO i = 1 , Nrow
         D(i) = 0.D0
      ENDDO
      DO j = 1 , Ncol
         DO i = 1 , Nrow
            D(i) = DMAX1(D(i),DABS(W(i,j)))
         ENDDO
      ENDDO
!
!...  gauss elimination with pivoting of scaled rows, loop over
!...  k=1,.,last
!
      k = 1
!
!...  as pivot row for k-th step, pick among the rows not yet used,
!...  i.e., from rows  k ,..., nrow , the one whose k-th entry
!...  (compared to the row size) is largest. then, if this row
!...  does not turn out to be row k, interchange row k with this
!...  particular row and redefine ipivot(k).
!
 100  IF ( D(k).EQ.0.D0 ) THEN
!
!...  singularity flag set
!
         Info = k
      ELSEIF ( k.EQ.Nrow ) THEN
!
!...  if  last  .eq. nrow , check now that pivot element in last row
!...  is nonzero.
!
         IF ( DABS(W(Nrow,Nrow))+D(Nrow).GT.D(Nrow) ) RETURN
         Info = k
      ELSE
         l = k
         kp1 = k + 1
         colmax = DABS(W(k,k))/D(k)
!
!...       find the (relatively) largest pivot
!
         DO i = kp1 , Nrow
            IF ( DABS(W(i,k)).GT.colmax*D(i) ) THEN
               colmax = DABS(W(i,k))/D(i)
               l = i
            ENDIF
         ENDDO
         Ipivot(k) = l
         t = W(l,k)
         s = D(l)
         IF ( l.NE.k ) THEN
            W(l,k) = W(k,k)
            W(k,k) = t
            D(l) = D(k)
            D(k) = s
         ENDIF
!
!...       if pivot element is too small in absolute value, declare
!...       matrix to be noninvertible and quit.
!
         IF ( DABS(t)+D(k).LE.D(k) ) THEN
            Info = k
         ELSE
!
!...       otherwise, subtract the appropriate multiple of the pivot
!...       row from remaining rows, i.e., the rows (k+1),..., (nrow)
!...       to make k-th entry zero. save the multiplier in its place.
!...       for high performance do this operations column oriented.
!
            t = -1.0D0/t
            DO i = kp1 , Nrow
               W(i,k) = W(i,k)*t
            ENDDO
            DO j = kp1 , Ncol
               t = W(l,j)
               IF ( l.NE.k ) THEN
                  W(l,j) = W(k,j)
                  W(k,j) = t
               ENDIF
               IF ( t.NE.0.D0 ) THEN
                  DO i = kp1 , Nrow
                     W(i,j) = W(i,j) + W(i,k)*t
                  ENDDO
               ENDIF
            ENDDO
            k = kp1
!
!...       check for having reached the next block.
!
            IF ( k.LE.Last ) GOTO 100
            RETURN
         ENDIF
      ENDIF
      END
