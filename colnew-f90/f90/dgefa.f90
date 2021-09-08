!*==DGEFA.spg  processed by SPAG 6.72Dc at 05:12 on  8 Sep 2021
      SUBROUTINE DGEFA(A,Lda,N,Ipvt,Info)
      IMPLICIT NONE
!*--DGEFA4
 
!*********************************************************************72
!
!c DGEFA factors a double precision matrix by gaussian elimination.
!
!     dgefa is usually called by dgeco, but it can be called
!     directly with a saving in time if  rcond  is not needed.
!     (time for dgeco) = (1 + 9/n)*(time for dgefa) .
!
!     on entry
!
!        a       double precision(lda, n)
!                the matrix to be factored.
!
!        lda     integer
!                the leading dimension of the array  a .
!
!        n       integer
!                the order of the matrix  a .
!
!     on return
!
!        a       an upper triangular matrix and the multipliers
!                which were used to obtain it.
!                the factorization can be written  a = l*u  where
!                l  is a product of permutation and unit lower
!                triangular matrices and  u  is upper triangular.
!
!        ipvt    integer(n)
!                an integer vector of pivot indices.
!
!        info    integer
!                = 0  normal value.
!                = k  if  u(k,k) .eq. 0.0 .  this is not an error
!                     condition for this subroutine, but it does
!                     indicate that dgesl or dgedi will divide by zero
!                     if called.  use  rcond  in dgeco for a reliable
!                     indication of singularity.
!
!     linpack. this version dated 08/14/78 .
!     cleve moler, university of new mexico, argonne national lab.
!
      INTEGER Lda , N , Ipvt(1) , Info
      DOUBLE PRECISION A(Lda,1)
      DOUBLE PRECISION t
      INTEGER IDAMAX , j , k , kp1 , l , nm1
!
!     gaussian elimination with partial pivoting
!
      Info = 0
      nm1 = N - 1
      IF ( nm1.GE.1 ) THEN
         DO k = 1 , nm1
            kp1 = k + 1
!
!        find l = pivot index
!
            l = IDAMAX(N-k+1,A(k,k),1) + k - 1
            Ipvt(k) = l
!
!        zero pivot implies this column already triangularized
!
            IF ( A(l,k).EQ.0.0D0 ) THEN
               Info = k
            ELSE
!
!           interchange if necessary
!
               IF ( l.NE.k ) THEN
                  t = A(l,k)
                  A(l,k) = A(k,k)
                  A(k,k) = t
               ENDIF
!
!           compute multipliers
!
               t = -1.0D0/A(k,k)
               CALL DSCAL(N-k,t,A(k+1,k),1)
!
!           row elimination with column indexing
!
               DO j = kp1 , N
                  t = A(l,j)
                  IF ( l.NE.k ) THEN
                     A(l,j) = A(k,j)
                     A(k,j) = t
                  ENDIF
                  CALL DAXPY(N-k,t,A(k+1,k),1,A(k+1,j),1)
               ENDDO
            ENDIF
         ENDDO
      ENDIF
      Ipvt(N) = N
      IF ( A(N,N).EQ.0.0D0 ) Info = N
      END
