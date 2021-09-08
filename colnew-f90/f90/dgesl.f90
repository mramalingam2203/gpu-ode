!*==DGESL.spg  processed by SPAG 6.72Dc at 05:13 on  8 Sep 2021
      SUBROUTINE DGESL(A,Lda,N,Ipvt,B,Job)
      IMPLICIT NONE
!*--DGESL4
 
!*********************************************************************72
!
!c DGESL solves a linear system factored by DGEFA.
!
!  DGESL solves the double precision system
!     a * x = b  or  trans(a) * x = b
!     using the factors computed by dgeco or dgefa.
!
!     on entry
!
!        a       double precision(lda, n)
!                the output from dgeco or dgefa.
!
!        lda     integer
!                the leading dimension of the array  a .
!
!        n       integer
!                the order of the matrix  a .
!
!        ipvt    integer(n)
!                the pivot vector from dgeco or dgefa.
!
!        b       double precision(n)
!                the right hand side vector.
!
!        job     integer
!                = 0         to solve  a*x = b ,
!                = nonzero   to solve  trans(a)*x = b  where
!                            trans(a)  is the transpose.
!
!     on return
!
!        b       the solution vector  x .
!
!     error condition
!
!        a division by zero will occur if the input factor contains a
!        zero on the diagonal.  technically this indicates singularity
!        but it is often caused by improper arguments or improper
!        setting of lda .  it will not occur if the subroutines are
!        called correctly and if dgeco has set rcond .gt. 0.0
!        or dgefa has set info .eq. 0 .
!
!     to compute  inverse(a) * c  where  c  is a matrix
!     with  p  columns
!           call dgeco(a,lda,n,ipvt,rcond,z)
!           if (rcond is too small) go to ...
!           do 10 j = 1, p
!              call dgesl(a,lda,n,ipvt,c(1,j),0)
!        10 continue
!
!     linpack. this version dated 08/14/78 .
!     cleve moler, university of new mexico, argonne national lab.
!
      INTEGER Lda , N , Ipvt(1) , Job
      DOUBLE PRECISION A(Lda,1) , B(1)
      DOUBLE PRECISION DDOT , t
      INTEGER k , kb , l , nm1
!
      nm1 = N - 1
      IF ( Job.NE.0 ) THEN
!
!        job = nonzero, solve  trans(a) * x = b
!        first solve  trans(u)*y = b
!
         DO k = 1 , N
            t = DDOT(k-1,A(1,k),1,B(1),1)
            B(k) = (B(k)-t)/A(k,k)
         ENDDO
!
!        now solve trans(l)*x = y
!
         IF ( nm1.GE.1 ) THEN
            DO kb = 1 , nm1
               k = N - kb
               B(k) = B(k) + DDOT(N-k,A(k+1,k),1,B(k+1),1)
               l = Ipvt(k)
               IF ( l.NE.k ) THEN
                  t = B(l)
                  B(l) = B(k)
                  B(k) = t
               ENDIF
            ENDDO
         ENDIF
      ELSE
!
!        job = 0 , solve  a * x = b
!        first solve  l*y = b
!
         IF ( nm1.GE.1 ) THEN
            DO k = 1 , nm1
               l = Ipvt(k)
               t = B(l)
               IF ( l.NE.k ) THEN
                  B(l) = B(k)
                  B(k) = t
               ENDIF
               CALL DAXPY(N-k,t,A(k+1,k),1,B(k+1),1)
            ENDDO
         ENDIF
!
!        now solve  u*x = y
!
         DO kb = 1 , N
            k = N + 1 - kb
            B(k) = B(k)/A(k,k)
            t = -B(k)
            CALL DAXPY(k-1,t,A(1,k),1,B(1),1)
         ENDDO
      ENDIF
      END
