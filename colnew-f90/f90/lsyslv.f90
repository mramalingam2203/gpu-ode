!*==LSYSLV.spg  processed by SPAG 6.72Dc at 06:44 on  8 Sep 2021
!---------------------------------------------------------------------
!                            p a r t  3
!          collocation system setup routines
!---------------------------------------------------------------------
!
      SUBROUTINE LSYSLV(Msing,Xi,Xiold,Z,Dmz,Delz,Deldmz,G,W,V,Rhs,Dmzo,&
                      & Integs,Ipvtg,Ipvtw,Rnorm,Mode,FSUB,DFSUB,GSUB,  &
                      & DGSUB,GUESS)
!*********************************************************************
!
!c LSYSLV controls the solution of linear systems.
!
!   purpose
!         this routine controls the set up and solution of a linear
!      system of collocation equations.
!         the matrix  g  is cast into an almost block diagonal
!      form by an appropriate ordering of the columns and solved
!      using the package of de boor-weiss [5]. the matrix is composed
!      of n blocks. the i-th block has the size
!                  integs(1,i) * integs(2,i).
!      it contains in its last rows the linearized collocation
!      equations, condensed as described in [2],
!      and the linearized side conditions corresponding to
!      the i-th subinterval.  integs(3,i)  steps of gaussian
!      elimination are applied to it to achieve a  partial plu
!      decomposition.  the right hand side vector is put into  rhs
!      and the solution vector is returned in  delz and deldmz.
!
!         lsyslv operates according to one of 5 modes:
!      mode = 0 - set up the collocation matrices  v , w , g
!                 and the right hand side  rhs ,  and solve.
!                 (for linear problems only.)
!      mode = 1 - set up the collocation matrices  v , w , g
!                 and the right hand sides  rhs  and  dmzo ,
!                 and solve. also set up  integs .
!                 (first iteration of nonlinear problems only).
!      mode = 2 - set up  rhs  only and compute its norm.
!      mode = 3 - set up  v, w, g  only and solve system.
!      mode = 4 - perform forward and backward substitution only
!                 (do not set up the matrices nor form the rhs).
!
!   variables
!
!      ig,izeta  - pointers to g,zeta respectively
!                       (necessary to keep track of blocks of g
!                       during matrix manipulations)
!      idmz,irhs,iv,iw - pointers to  rhs,v,w rspectively
!      df    - partial derivatives of f from dfsub
!      rnorm - euclidean norm of rhs
!      lside - number of side conditions in current and previous blocks
!      iguess = 1 when current soln is user specified via  guess
!             = 0 otherwise
!
!*********************************************************************
      IMPLICIT NONE
!*--LSYSLV57
!*** Start of declarations inserted by SPAG
      REAL*8 ACOl , ALEft , ARIght , ASAve , at , B , COEf , Deldmz ,   &
           & Delz , df , DFSUB , DGSUB , dgz , dmval , Dmz , Dmzo ,     &
           & dummy , f , G , gval
      REAL*8 h , hrho , paramsh , paramsl , PREcis , RHO , Rhs , Rnorm ,&
           & V , value , W , xcol , Xi , xii , Xiold , Z , ZETa , zval
      INTEGER i , ICAre , idmz , idmzo , ig , IGUess , Integs , iold ,  &
            & IOUt , IPRint , Ipvtg , Ipvtw , irhs , ITEr , iv , iw ,   &
            & iz , izet , IZEta , IZSave
      INTEGER j , jj , K , KD , l , LIMit , lside , lw , M , m1 , MMAx ,&
            & Mode , Msing , MSTar , N , ncol , NCOmp , NDMz , NMAx ,   &
            & NOLd
      INTEGER NONlin , nrow , NZ
!*** End of declarations inserted by SPAG
      DIMENSION Z(1) , Dmz(1) , Delz(1) , Deldmz(1) , Xi(1) , Xiold(1)
      DIMENSION G(1) , W(1) , V(1) , Rhs(1) , Dmzo(1) , dummy(1)
      DIMENSION Integs(3,1) , Ipvtg(1) , Ipvtw(1)
      DIMENSION zval(40) , f(40) , dgz(40) , dmval(20) , df(800) ,      &
              & at(28)
!
      COMMON /COLOUT/ PREcis , IOUt , IPRint
      COMMON /COLLOC/ RHO(7) , COEf(49)
      COMMON /COLORD/ K , NCOmp , MSTar , KD , MMAx , M(20)
      COMMON /COLSID/ ZETa(40) , ALEft , ARIght , IZEta , IZSave
      COMMON /COLAPR/ N , NOLd , NMAx , NZ , NDMz
      COMMON /COLNLN/ NONlin , ITEr , LIMit , ICAre , IGUess
      COMMON /COLBAS/ B(28) , ACOl(28,7) , ASAve(28,4)
!
      EXTERNAL DFSUB , DGSUB
!
      m1 = Mode + 1
      IF ( m1.EQ.2 .OR. m1.EQ.3 .OR. m1.EQ.4 ) THEN
      ELSEIF ( m1.EQ.5 ) THEN
         GOTO 400
      ELSE
!
!...  linear problem initialization
!
         DO i = 1 , MSTar
            zval(i) = 0.D0
         ENDDO
      ENDIF
!
!...  initialization
!
      idmz = 1
      idmzo = 1
      irhs = 1
      ig = 1
      iw = 1
      iv = 1
      IZEta = 1
      lside = 0
      iold = 1
      ncol = 2*MSTar
      Rnorm = 0.D0
      IF ( Mode.LE.1 ) THEN
!
!...  build integs (describing block structure of matrix)
!
         DO i = 1 , N
            Integs(2,i) = ncol
            IF ( i.LT.N ) THEN
               Integs(3,i) = MSTar
 10            IF ( lside.NE.MSTar ) THEN
                  IF ( ZETa(lside+1).LT.Xi(i)+PREcis ) THEN
                     lside = lside + 1
                     GOTO 10
                  ENDIF
               ENDIF
            ELSE
               Integs(3,N) = ncol
               lside = MSTar
            ENDIF
            nrow = MSTar + lside
            Integs(1,i) = nrow
         ENDDO
      ENDIF
      IF ( Mode.NE.2 ) THEN
!
!...  zero the matrices to be computed
!
         lw = KD*KD*N
         DO l = 1 , lw
            W(l) = 0.D0
         ENDDO
      ENDIF
!
!...  the do loop 290 sets up the linear system of equations.
!
      DO i = 1 , N
!
!...       construct a block of  a  and a corresponding piece of  rhs.
!
         xii = Xi(i)
         h = Xi(i+1) - Xi(i)
         nrow = Integs(1,i)
!
!...       go thru the ncomp collocation equations and side conditions
!...       in the i-th subinterval
!
 50      IF ( IZEta.LE.MSTar ) THEN
            IF ( ZETa(IZEta).LE.xii+PREcis ) THEN
!
!...       build equation for a side condition.
!
               IF ( Mode.NE.0 ) THEN
                  IF ( IGUess.EQ.1 ) THEN
!
!...       case where user provided current approximation
!
                     CALL GUESS(xii,zval,dmval)
!
!...       other nonlinear case
!
                  ELSEIF ( Mode.NE.1 ) THEN
                     CALL APPROX(i,xii,zval,at,dummy,Xi,N,Z,Dmz,K,NCOmp,&
                               & MMAx,M,MSTar,1,dummy,0)
                     IF ( Mode.EQ.3 ) THEN
!
!...       build a row of  a  corresponding to a boundary point
!
                        CALL GDERIV(G(ig),nrow,IZEta,zval,dgz,1,DGSUB)
                        IZEta = IZEta + 1
                        GOTO 50
                     ENDIF
                  ELSE
                     CALL APPROX(iold,xii,zval,at,COEf,Xiold,NOLd,Z,Dmz,&
                               & K,NCOmp,MMAx,M,MSTar,2,dummy,0)
                  ENDIF
               ENDIF
!
!...       find  rhs  boundary value.
!
               CALL GSUB(IZEta,zval,gval)
               Rhs(NDMz+IZEta) = -gval
               Rnorm = Rnorm + gval**2
               IF ( Mode.EQ.2 ) THEN
                  IZEta = IZEta + 1
               ELSE
                  CALL GDERIV(G(ig),nrow,IZEta,zval,dgz,1,DGSUB)
                  IZEta = IZEta + 1
               ENDIF
               GOTO 50
            ENDIF
         ENDIF
!
!...       assemble collocation equations
!
         DO j = 1 , K
            hrho = h*RHO(j)
            xcol = xii + hrho
!
!...         this value corresponds to a collocation (interior)
!...         point. build the corresponding  ncomp  equations.
!
            IF ( Mode.EQ.0 ) THEN
!
!...         the linear case
!
               CALL FSUB(paramsl,paramsh,xcol,zval,Rhs(irhs))
!  			 type (params) :: paramsL,paramsH
 
               irhs = irhs + NCOmp
            ELSE
               IF ( IGUess.EQ.1 ) THEN
!
!...         use initial approximation provided by the user.
!
                  CALL GUESS(xcol,zval,Dmzo(irhs))
!
!...         find  rhs  values
!
               ELSEIF ( Mode.NE.1 ) THEN
!
!...         evaluate former collocation solution
!
                  CALL APPROX(i,xcol,zval,ACOl(1,j),COEf,Xi,N,Z,Dmz,K,  &
                            & NCOmp,MMAx,M,MSTar,4,dummy,0)
                  IF ( Mode.EQ.3 ) GOTO 60
!
!...         fill in  rhs  values (and accumulate its norm).
!
                  CALL FSUB(paramsl,paramsh,xcol,zval,f)
!  			 type (params) :: paramsL,paramsH
 
                  DO jj = 1 , NCOmp
                     value = Dmz(irhs) - f(jj)
                     Rhs(irhs) = -value
                     Rnorm = Rnorm + value**2
                     irhs = irhs + 1
                  ENDDO
                  GOTO 100
               ELSE
                  CALL APPROX(iold,xcol,zval,at,COEf,Xiold,NOLd,Z,Dmz,K,&
                            & NCOmp,MMAx,M,MSTar,2,Dmzo(irhs),1)
               ENDIF
!
               CALL FSUB(paramsl,paramsh,xcol,zval,f)
!  			 type (params) :: paramsL,paramsH
 
               DO jj = 1 , NCOmp
                  value = Dmzo(irhs) - f(jj)
                  Rhs(irhs) = -value
                  Rnorm = Rnorm + value**2
                  irhs = irhs + 1
               ENDDO
            ENDIF
!
!...         fill in ncomp rows of  w and v
!
 60         CALL VWBLOK(xcol,hrho,j,W(iw),V(iv),Ipvtw(idmz),KD,zval,df, &
                      & ACOl(1,j),Dmzo(idmzo),NCOmp,DFSUB,Msing)
            IF ( Msing.NE.0 ) RETURN
 100     ENDDO
!
!...       build global bvp matrix  g
!
         IF ( Mode.NE.2 ) CALL GBLOCK(h,G(ig),nrow,IZEta,W(iw),V(iv),KD,&
                                    & dummy,Deldmz(idmz),Ipvtw(idmz),1)
         IF ( i.LT.N ) THEN
!
!...       update counters -- i-th block completed
!
            ig = ig + nrow*ncol
            iv = iv + KD*MSTar
            iw = iw + KD*KD
            idmz = idmz + KD
            IF ( Mode.EQ.1 ) idmzo = idmzo + KD
            GOTO 300
         ELSE
            IZSave = IZEta
         ENDIF
 150     IF ( IZEta.GT.MSTar ) GOTO 300
!
!...       build equation for a side condition.
!
         IF ( Mode.NE.0 ) THEN
            IF ( IGUess.EQ.1 ) THEN
!
!...       case where user provided current approximation
!
               CALL GUESS(ARIght,zval,dmval)
!
!...       other nonlinear case
!
            ELSEIF ( Mode.NE.1 ) THEN
               CALL APPROX(N+1,ARIght,zval,at,COEf,Xi,N,Z,Dmz,K,NCOmp,  &
                         & MMAx,M,MSTar,1,dummy,0)
               IF ( Mode.EQ.3 ) THEN
!
!...       build a row of  a  corresponding to a boundary point
!
                  CALL GDERIV(G(ig),nrow,IZEta+MSTar,zval,dgz,2,DGSUB)
                  GOTO 200
               ENDIF
            ELSE
               CALL APPROX(NOLd+1,ARIght,zval,at,COEf,Xiold,NOLd,Z,Dmz, &
                         & K,NCOmp,MMAx,M,MSTar,1,dummy,0)
            ENDIF
         ENDIF
!
!...       find  rhs  boundary value.
!
         CALL GSUB(IZEta,zval,gval)
         Rhs(NDMz+IZEta) = -gval
         Rnorm = Rnorm + gval**2
         IF ( Mode.NE.2 ) CALL GDERIV(G(ig),nrow,IZEta+MSTar,zval,dgz,2,&
                                    & DGSUB)
 200     IZEta = IZEta + 1
         GOTO 150
 300  ENDDO
!
!...       assembly process completed
!
      IF ( Mode.NE.0 .AND. Mode.NE.3 ) THEN
         Rnorm = DSQRT(Rnorm/DFLOAT(NZ+NDMz))
         IF ( Mode.EQ.2 ) RETURN
      ENDIF
!
!...  solve the linear system.
!
!...  matrix decomposition
!
      CALL FCBLOK(G,Integs,N,Ipvtg,df,Msing)
!
!...  check for singular matrix
!
      Msing = -Msing
      IF ( Msing.NE.0 ) RETURN
!
!...  perform forward and backward substitution for mode=4 only.
!
 400  DO l = 1 , NDMz
         Deldmz(l) = Rhs(l)
      ENDDO
      iz = 1
      idmz = 1
      iw = 1
      izet = 1
      DO i = 1 , N
         nrow = Integs(1,i)
         IZEta = nrow + 1 - MSTar
         IF ( i.EQ.N ) IZEta = IZSave
 450     IF ( izet.EQ.IZEta ) THEN
            h = Xi(i+1) - Xi(i)
            CALL GBLOCK(h,G(1),nrow,IZEta,W(iw),V(1),KD,Delz(iz),       &
                      & Deldmz(idmz),Ipvtw(idmz),2)
            iz = iz + MSTar
            idmz = idmz + KD
            iw = iw + KD*KD
            IF ( i.GE.N ) THEN
 460           IF ( izet.LE.MSTar ) THEN
                  Delz(iz-1+izet) = Rhs(NDMz+izet)
                  izet = izet + 1
                  GOTO 460
               ENDIF
            ENDIF
         ELSE
            Delz(iz-1+izet) = Rhs(NDMz+izet)
            izet = izet + 1
            GOTO 450
         ENDIF
      ENDDO
!
!...  perform forward and backward substitution for mode=0,2, or 3.
!
      CALL SBBLOK(G,Integs,N,Ipvtg,Delz)
!
!...  finaly find deldmz
!
      CALL DMZSOL(KD,MSTar,N,V,Delz,Deldmz)
!
      IF ( Mode.NE.1 ) RETURN
      DO l = 1 , NDMz
         Dmz(l) = Dmzo(l)
      ENDDO
      iz = 1
      idmz = 1
      iw = 1
      izet = 1
      DO i = 1 , N
         nrow = Integs(1,i)
         IZEta = nrow + 1 - MSTar
         IF ( i.EQ.N ) IZEta = IZSave
 500     IF ( izet.EQ.IZEta ) THEN
            h = Xi(i+1) - Xi(i)
            CALL GBLOCK(h,G(1),nrow,IZEta,W(iw),df,KD,Z(iz),Dmz(idmz),  &
                      & Ipvtw(idmz),2)
            iz = iz + MSTar
            idmz = idmz + KD
            iw = iw + KD*KD
            IF ( i.GE.N ) THEN
 510           IF ( izet.LE.MSTar ) THEN
                  Z(iz-1+izet) = dgz(izet)
                  izet = izet + 1
                  GOTO 510
               ENDIF
            ENDIF
         ELSE
            Z(iz-1+izet) = dgz(izet)
            izet = izet + 1
            GOTO 500
         ENDIF
      ENDDO
      CALL SBBLOK(G,Integs,N,Ipvtg,Z)
!
!...  finaly find dmz
!
      CALL DMZSOL(KD,MSTar,N,V,Z,Dmz)
!
      END
