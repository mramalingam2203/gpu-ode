!*==SUBFOR.spg  processed by SPAG 6.72Dc at 06:34 on  8 Sep 2021
      SUBROUTINE SUBFOR(W,Ipivot,Nrow,Last,X)
      IMPLICIT NONE
!*--SUBFOR4
!*** Start of declarations inserted by SPAG
      INTEGER i , Last , Nrow
!*** End of declarations inserted by SPAG
!
!**********************************************************************
!
!c SUBFOR carries out the forward pass of substitution for the current block.
!
!     carries out the forward pass of substitution for the current
!     block, i.e., the action on the right side corresponding to the
!     elimination carried out in  factrb  for this block.
!
!    parameters
!       w, ipivot, nrow, last  are as on return from factrb.
!       x(j)  is expected to contain, on input, the right side of j-th
!             equation for this block, j=1,...,nrow.
!       x(j)  contains, on output, the appropriately modified right
!             side of equation (j) in this block, j=1,...,last and
!             for j=last+1,...,nrow.
!
!*********************************************************************
!
      INTEGER Ipivot(Last) , ip , k , kp1 , lstep
      DOUBLE PRECISION W(Nrow,Last) , X(Nrow) , t
!
      IF ( Nrow.EQ.1 ) RETURN
      lstep = MIN0(Nrow-1,Last)
      DO k = 1 , lstep
         kp1 = k + 1
         ip = Ipivot(k)
         t = X(ip)
         X(ip) = X(k)
         X(k) = t
         IF ( t.NE.0.D0 ) THEN
            DO i = kp1 , Nrow
               X(i) = X(i) + W(i,k)*t
            ENDDO
         ENDIF
      ENDDO
      END
