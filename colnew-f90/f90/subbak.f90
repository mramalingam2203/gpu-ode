!*==SUBBAK.spg  processed by SPAG 6.72Dc at 06:33 on  8 Sep 2021
      SUBROUTINE SUBBAK(W,Nrow,Ncol,Last,X)
      IMPLICIT NONE
!*--SUBBAK4
!*** Start of declarations inserted by SPAG
      INTEGER kb , Ncol , Nrow
!*** End of declarations inserted by SPAG
!
!*********************************************************************
!
!c SUBBAK carries out backsubstitution for the current block.
!
!    parameters
!       w, ipivot, nrow, ncol, last  are as on return from factrb.
!       x(1),...,x(ncol)  contains, on input, the right side for the
!               equations in this block after backsubstitution has been
!               carried up to but not including equation (last).
!               means that x(j) contains the right side of equation (j)
!               as modified during elimination, j=1,...,last, while
!               for j .gt. last, x(j) is already a component of the
!               solution vector.
!       x(1),...,x(ncol) contains, on output, the components of the
!               solution corresponding to the present block.
!
!*********************************************************************
!
      INTEGER Last , i , j , k , km1 , lm1 , lp1
      DOUBLE PRECISION W(Nrow,Ncol) , X(Ncol) , t
!
      lp1 = Last + 1
      IF ( lp1.LE.Ncol ) THEN
         DO j = lp1 , Ncol
            t = -X(j)
            IF ( t.NE.0.D0 ) THEN
               DO i = 1 , Last
                  X(i) = X(i) + W(i,j)*t
               ENDDO
            ENDIF
         ENDDO
      ENDIF
      IF ( Last.NE.1 ) THEN
         lm1 = Last - 1
         DO kb = 1 , lm1
            km1 = Last - kb
            k = km1 + 1
            X(k) = X(k)/W(k,k)
            t = -X(k)
            IF ( t.NE.0.D0 ) THEN
               DO i = 1 , km1
                  X(i) = X(i) + W(i,k)*t
               ENDDO
            ENDIF
         ENDDO
      ENDIF
      X(1) = X(1)/W(1,1)
      END
