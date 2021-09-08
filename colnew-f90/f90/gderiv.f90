!*==GDERIV.spg  processed by SPAG 6.72Dc at 06:41 on  8 Sep 2021
      SUBROUTINE GDERIV(Gi,Nrow,Irow,Zval,Dgz,Mode,DGSUB)
!
!**********************************************************************
!
!c GDERIV constructs a row of the collocation matrix.
!
!   purpose:
!
!      construct a collocation matrix row according to mode:
!      mode = 1  -  a row corresponding to a initial condition
!                   (i.e. at the left end of the subinterval).
!      mode = 2  -  a row corresponding to a final condition.
!
!   variables:
!
!      gi     - the sub-block of the global bvp matrix in
!               which the equations are to be formed.
!      nrow   - no. of rows in gi.
!      irow   - the row in gi to be used for equations.
!      zval   - z(xi)
!      dg     - the derivatives of the side condition.
!
!**********************************************************************
      IMPLICIT NONE
!*--GDERIV26
!*** Start of declarations inserted by SPAG
      REAL*8 ALEft , ARIght , dg , Dgz , dot , Gi , ZETa , Zval
      INTEGER ICAre , IDUm , IGUess , Irow , ITEr , IZEta , j , KD ,    &
            & KDUm , LIMit , M , MMAx , Mode , MSTar , NDUm , NONlin ,  &
            & Nrow
!*** End of declarations inserted by SPAG
      DIMENSION Gi(Nrow,1) , Zval(1) , Dgz(1) , dg(40)
!
      COMMON /COLORD/ KDUm , NDUm , MSTar , KD , MMAx , M(20)
      COMMON /COLSID/ ZETa(40) , ALEft , ARIght , IZEta , IDUm
      COMMON /COLNLN/ NONlin , ITEr , LIMit , ICAre , IGUess
!
!...  zero jacobian dg
!
      DO j = 1 , MSTar
         dg(j) = 0.D0
      ENDDO
!
!...  evaluate jacobian dg
!
      CALL DGSUB(IZEta,Zval,dg)
!
!...  evaluate  dgz = dg * zval  once for a new mesh
!
      IF ( NONlin.NE.0 .AND. ITEr.LE.0 ) THEN
         dot = 0.D0
         DO j = 1 , MSTar
            dot = dot + dg(j)*Zval(j)
         ENDDO
         Dgz(IZEta) = dot
      ENDIF
!
!...  branch according to  m o d e
!
      IF ( Mode.EQ.2 ) THEN
!
!...  handle a final condition
!
         DO j = 1 , MSTar
            Gi(Irow,j) = 0.D0
            Gi(Irow,MSTar+j) = dg(j)
         ENDDO
         GOTO 99999
      ENDIF
!
!...  provide coefficients of the j-th linearized side condition.
!...  specifically, at x=zeta(j) the j-th side condition reads
!...  dg(1)*z(1) + ... +dg(mstar)*z(mstar) + g = 0
!
!
!...  handle an initial condition
!
      DO j = 1 , MSTar
         Gi(Irow,j) = dg(j)
         Gi(Irow,MSTar+j) = 0.D0
      ENDDO
      RETURN
99999 END
