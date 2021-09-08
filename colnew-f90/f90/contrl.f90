!*==CONTRL.spg  processed by SPAG 6.72Dc at 06:53 on  8 Sep 2021
      SUBROUTINE CONTRL(Xi,Xiold,Z,Dmz,Rhs,Delz,Deldmz,Dqz,Dqdmz,G,W,V, &
                      & Valstr,Slope,Scale,Dscale,Accum,Ipvtg,Integs,   &
                      & Ipvtw,Nfxpnt,Fixpnt,Iflag,FSUB,DFSUB,GSUB,DGSUB,&
                      & GUESS)
 
!*********************************************************************72
!
!c CONTRL is the driver for COLNEW.
!
!   purpose
!     this subroutine is the actual driver.  the nonlinear iteration
!     strategy is controlled here ( see [4] ). upon convergence, errchk
!     is called to test for satisfaction of the requested tolerances.
!
!   variables
!
!     check  - maximum tolerance value, used as part of criteria for
!              checking for nonlinear iteration convergence
!     relax  - the relaxation factor for damped newton iteration
!     relmin - minimum allowable value for relax  (otherwise the
!              jacobian is considered singular).
!     rlxold - previous relax
!     rstart - initial value for relax when problem is sensitive
!     ifrz   - number of fixed jacobian iterations
!     lmtfrz - maximum value for ifrz before performing a reinversion
!     iter   - number of iterations (counted only when jacobian
!              reinversions are performed).
!     xi     - current mesh
!     xiold  - previous mesh
!     ipred  = 0  if relax is determined by a correction
!            = 1  if relax is determined by a prediction
!     ifreez = 0  if the jacobian is to be updated
!            = 1  if the jacobian is currently fixed (frozen)
!     iconv  = 0  if no previous convergence has been obtained
!            = 1  if convergence on a previous mesh has been obtained
!     icare  =-1  no convergence occurred (used for regular problems)
!            = 0  a regular problem
!            = 1  a sensitive problem
!            = 2  used for continuation (see description of ipar(10)
!                 in colnew).
!     rnorm  - norm of rhs (right hand side) for current iteration
!     rnold  - norm of rhs for previous iteration
!     anscl  - scaled norm of newton correction
!     anfix  - scaled norm of newton correction at next step
!     anorm  - scaled norm of a correction obtained with jacobian fixed
!     nz     - number of components of  z  (see subroutine approx)
!     ndmz   - number of components of  dmz  (see subroutine approx)
!     imesh  - a control variable for subroutines newmsh and errchk
!            = 1  the current mesh resulted from mesh selection
!                 or is the initial mesh.
!            = 2  the current mesh resulted from doubling the
!                 previous mesh
!
!**********************************************************************
!
      IMPLICIT NONE
!*--CONTRL58
!*** Start of declarations inserted by SPAG
      REAL*8 Accum , ALEft , andif , anfix , anorm , anscl , arg ,      &
           & ARIght , check , Deldmz , Delz , DFSUB , DGSUB , Dmz ,     &
           & Dqdmz , Dqz , Dscale , dummy , fact , factor
      REAL*8 Fixpnt , FSUB , G , GSUB , GUESS , PREcis , relax ,        &
           & relmin , Rhs , rlxold , rnold , rnorm , ROOt , rstart ,    &
           & Scale , Slope , TOL , TOLin , V , Valstr
      REAL*8 W , WGTerr , WGTmsh , Xi , Xiold , Z , ZETa
      INTEGER i , ICAre , iconv , icor , IDUm , ifin , Iflag , ifreez , &
            & ifrz , IGUess , imesh , Integs , inz , IOUt , ipred ,     &
            & IPRint , Ipvtg , Ipvtw , it , ITEr
      INTEGER iz , IZEta , j , JTOl , K , KD , LIMit , lj , lmtfrz ,    &
            & LTOl , M , MMAx , MSHalt , MSHflg , MSHlmt , MSHnum ,     &
            & msing , MSTar , N , NCOmp
      INTEGER NDMz , Nfxpnt , NMAx , noconv , NOLd , NONlin , np1 ,     &
            & NTOl , NZ
!*** End of declarations inserted by SPAG
      DIMENSION Xi(1) , Xiold(1) , Z(1) , Dmz(1) , Rhs(1)
      DIMENSION G(1) , W(1) , V(1) , Valstr(1) , Slope(1) , Accum(1)
      DIMENSION Delz(1) , Deldmz(1) , Dqz(1) , Dqdmz(1) , Fixpnt(1)
      DIMENSION dummy(1) , Scale(1) , Dscale(1)
      DIMENSION Integs(1) , Ipvtg(1) , Ipvtw(1)
!
      COMMON /COLOUT/ PREcis , IOUt , IPRint
      COMMON /COLORD/ K , NCOmp , MSTar , KD , MMAx , M(20)
      COMMON /COLAPR/ N , NOLd , NMAx , NZ , NDMz
      COMMON /COLMSH/ MSHflg , MSHnum , MSHlmt , MSHalt
      COMMON /COLSID/ ZETa(40) , ALEft , ARIght , IZEta , IDUm
      COMMON /COLNLN/ NONlin , ITEr , LIMit , ICAre , IGUess
      COMMON /COLEST/ TOL(40) , WGTmsh(40) , WGTerr(40) , TOLin(40) ,   &
                    & ROOt(40) , JTOl(40) , LTOl(40) , NTOl
!
      EXTERNAL FSUB , DFSUB , GSUB , DGSUB , GUESS
!
!...  constants for control of nonlinear iteration
!
      relmin = 1.D-3
      rstart = 1.D-2
      lmtfrz = 4
!
!...  compute the maximum tolerance
!
      check = 0.D0
      DO i = 1 , NTOl
         check = DMAX1(TOLin(i),check)
      ENDDO
      imesh = 1
      iconv = 0
      IF ( NONlin.EQ.0 ) iconv = 1
      icor = 0
      noconv = 0
      msing = 0
!
!...  the main iteration begins here .
!...  loop 20 is executed until error tolerances are satisfied or
!...  the code fails (due to a singular matrix or storage limitations)
!
!
!...       initialization for a new mesh
!
 100  ITEr = 0
      IF ( NONlin.GT.0 ) THEN
!
!...       iteration loop for nonlinear case
!...       define the initial relaxation parameter (= relax)
!
         relax = 1.D0
!
!...       check for previous convergence and problem sensitivity
!
         IF ( ICAre.EQ.1 .OR. ICAre.EQ.(-1) ) relax = rstart
         IF ( iconv.EQ.0 ) GOTO 500
!
!...       convergence on a previous mesh has been obtained.    thus
!...       we have a very good initial approximation for the newton
!...       process.    proceed with one full newton and then iterate
!...       with a fixed jacobian.
!
         ifreez = 0
!
!...       evaluate right hand side and its norm  and
!...       find the first newton correction
!
         CALL LSYSLV(msing,Xi,Xiold,Z,Dmz,Delz,Deldmz,G,W,V,Rhs,Dqdmz,  &
                   & Integs,Ipvtg,Ipvtw,rnold,1,FSUB,DFSUB,GSUB,DGSUB,  &
                   & GUESS)
!
         IF ( IPRint.LT.0 ) WRITE (IOUt,99001)
99001    FORMAT (/' FIXED JACOBIAN ITERATIONS,')
         IF ( IPRint.LT.0 ) WRITE (IOUt,99014) ITEr , rnold
         GOTO 400
      ELSE
!
!...       the linear case.
!...       set up and solve equations
!
         CALL LSYSLV(msing,Xi,Xiold,dummy,dummy,Z,Dmz,G,W,V,Rhs,dummy,  &
                   & Integs,Ipvtg,Ipvtw,rnorm,0,FSUB,DFSUB,GSUB,DGSUB,  &
                   & GUESS)
!
!...       check for a singular matrix
!
         IF ( msing.EQ.0 ) GOTO 1000
      ENDIF
 200  IF ( msing.LT.0 ) THEN
         IF ( IPRint.LT.1 ) WRITE (IOUt,99002)
!     ---------------------------------------------------------------
99002    FORMAT (//' THE GLOBAL BVP-MATRIX IS SINGULAR ')
         Iflag = 0
         RETURN
      ELSE
         IF ( IPRint.LT.1 ) WRITE (IOUt,99003)
99003    FORMAT (//' A LOCAL ELIMINATION MATRIX IS SINGULAR ')
         GOTO 1200
      ENDIF
!
!...       solve for the next iterate .
!...       the value of ifreez determines whether this is a full
!...       newton step (=0) or a fixed jacobian iteration (=1).
!
 300  IF ( IPRint.LT.0 ) WRITE (IOUt,99014) ITEr , rnorm
      rnold = rnorm
      CALL LSYSLV(msing,Xi,Xiold,Z,Dmz,Delz,Deldmz,G,W,V,Rhs,dummy,     &
                & Integs,Ipvtg,Ipvtw,rnorm,3+ifreez,FSUB,DFSUB,GSUB,    &
                & DGSUB,GUESS)
!
!...       check for a singular matrix
!
 400  IF ( msing.NE.0 ) GOTO 200
      IF ( ifreez.NE.1 ) THEN
!
!...       a full newton step
!
         ITEr = ITEr + 1
         ifrz = 0
      ENDIF
!
!...       update   z and dmz , compute new  rhs  and its norm
!
      DO i = 1 , NZ
         Z(i) = Z(i) + Delz(i)
      ENDDO
      DO i = 1 , NDMz
         Dmz(i) = Dmz(i) + Deldmz(i)
      ENDDO
      CALL LSYSLV(msing,Xi,Xiold,Z,Dmz,Delz,Deldmz,G,W,V,Rhs,dummy,     &
                & Integs,Ipvtg,Ipvtw,rnorm,2,FSUB,DFSUB,GSUB,DGSUB,     &
                & GUESS)
!
!...       check monotonicity. if the norm of  rhs  gets smaller,
!...       proceed with a fixed jacobian; else proceed cautiously,
!...       as if convergence has not been obtained before (iconv=0).
!
      IF ( rnorm.LT.PREcis ) GOTO 900
      IF ( rnorm.GT.rnold ) THEN
!
!...      convergence of fixed jacobian iteration failed.
!
         IF ( IPRint.LT.0 ) WRITE (IOUt,99014) ITEr , rnorm
         IF ( IPRint.LT.0 ) WRITE (IOUt,99004)
99004    FORMAT (/' SWITCH TO DAMPED NEWTON ITERATION,')
         iconv = 0
         relax = rstart
         DO i = 1 , NZ
            Z(i) = Z(i) - Delz(i)
         ENDDO
         DO i = 1 , NDMz
            Dmz(i) = Dmz(i) - Deldmz(i)
         ENDDO
!
!...       update old mesh
!
         np1 = N + 1
         DO i = 1 , np1
            Xiold(i) = Xi(i)
         ENDDO
         NOLd = N
!
         ITEr = 0
      ELSEIF ( ifreez.EQ.1 ) THEN
!
!...       verify that the linear convergence with fixed jacobian
!...       is fast enough.
!
         ifrz = ifrz + 1
         IF ( ifrz.GE.lmtfrz ) ifreez = 0
         IF ( rnold.LT.4.D0*rnorm ) ifreez = 0
!
!...       check convergence (iconv = 1).
!
         DO it = 1 , NTOl
            inz = LTOl(it)
            DO iz = inz , NZ , MSTar
               IF ( DABS(Delz(iz)).GT.TOLin(it)*(DABS(Z(iz))+1.D0) )    &
                  & GOTO 300
            ENDDO
         ENDDO
!
!...       convergence obtained
!
         IF ( IPRint.LT.1 ) WRITE (IOUt,99015) ITEr
         GOTO 1000
      ELSE
         ifreez = 1
         GOTO 300
      ENDIF
!
!...       no previous convergence has been obtained. proceed
!...       with the damped newton method.
!...       evaluate rhs and find the first newton correction.
!
 500  IF ( IPRint.LT.0 ) WRITE (IOUt,99005)
99005 FORMAT (/' FULL DAMPED NEWTON ITERATION,')
      CALL LSYSLV(msing,Xi,Xiold,Z,Dmz,Delz,Deldmz,G,W,V,Rhs,Dqdmz,     &
                & Integs,Ipvtg,Ipvtw,rnold,1,FSUB,DFSUB,GSUB,DGSUB,     &
                & GUESS)
!
!...       check for a singular matrix
!
      IF ( msing.NE.0 ) GOTO 200
!
!...       bookkeeping for first mesh
!
      IF ( IGUess.EQ.1 ) IGUess = 0
!
!...       find initial scaling
!
      CALL SKALE(N,MSTar,KD,Z,Xi,Scale,Dscale)
      GOTO 700
!
!...       main iteration loop
!
 600  rnold = rnorm
      IF ( ITEr.GE.LIMit ) THEN
!
!...       diagnostics for failure of nonlinear iteration.
!
         IF ( IPRint.LT.1 ) WRITE (IOUt,99006) ITEr
99006    FORMAT (/' NO CONVERGENCE AFTER ',I3,' ITERATIONS'/)
         GOTO 1100
      ELSE
!
!...       update scaling
!
         CALL SKALE(N,MSTar,KD,Z,Xi,Scale,Dscale)
!
!...       compute norm of newton correction with new scaling
!
         anscl = 0.D0
         DO i = 1 , NZ
            anscl = anscl + (Delz(i)*Scale(i))**2
         ENDDO
         DO i = 1 , NDMz
            anscl = anscl + (Deldmz(i)*Dscale(i))**2
         ENDDO
         anscl = DSQRT(anscl/DFLOAT(NZ+NDMz))
!
!...       find a newton direction
!
         CALL LSYSLV(msing,Xi,Xiold,Z,Dmz,Delz,Deldmz,G,W,V,Rhs,dummy,  &
                   & Integs,Ipvtg,Ipvtw,rnorm,3,FSUB,DFSUB,GSUB,DGSUB,  &
                   & GUESS)
!
!...       check for a singular matrix
!
         IF ( msing.NE.0 ) GOTO 200
!
!...       predict relaxation factor for newton step.
!
         andif = 0.D0
         DO i = 1 , NZ
            andif = andif + ((Dqz(i)-Delz(i))*Scale(i))**2
         ENDDO
         DO i = 1 , NDMz
            andif = andif + ((Dqdmz(i)-Deldmz(i))*Dscale(i))**2
         ENDDO
         andif = DSQRT(andif/DFLOAT(NZ+NDMz)+PREcis)
         relax = relax*anscl/andif
         IF ( relax.GT.1.D0 ) relax = 1.D0
      ENDIF
 700  rlxold = relax
      ipred = 1
      ITEr = ITEr + 1
!
!...       determine a new  z and dmz  and find new  rhs  and its norm
!
      DO i = 1 , NZ
         Z(i) = Z(i) + relax*Delz(i)
      ENDDO
      DO i = 1 , NDMz
         Dmz(i) = Dmz(i) + relax*Deldmz(i)
      ENDDO
 800  CALL LSYSLV(msing,Xi,Xiold,Z,Dmz,Dqz,Dqdmz,G,W,V,Rhs,dummy,Integs,&
                & Ipvtg,Ipvtw,rnorm,2,FSUB,DFSUB,GSUB,DGSUB,GUESS)
!
!...       compute a fixed jacobian iterate (used to control relax)
!
      CALL LSYSLV(msing,Xi,Xiold,Z,Dmz,Dqz,Dqdmz,G,W,V,Rhs,dummy,Integs,&
                & Ipvtg,Ipvtw,rnorm,4,FSUB,DFSUB,GSUB,DGSUB,GUESS)
!
!...       find scaled norms of various terms used to correct relax
!
      anorm = 0.D0
      anfix = 0.D0
      DO i = 1 , NZ
         anorm = anorm + (Delz(i)*Scale(i))**2
         anfix = anfix + (Dqz(i)*Scale(i))**2
      ENDDO
      DO i = 1 , NDMz
         anorm = anorm + (Deldmz(i)*Dscale(i))**2
         anfix = anfix + (Dqdmz(i)*Dscale(i))**2
      ENDDO
      anorm = DSQRT(anorm/DFLOAT(NZ+NDMz))
      anfix = DSQRT(anfix/DFLOAT(NZ+NDMz))
      IF ( icor.EQ.1 ) THEN
         IF ( IPRint.LT.0 ) WRITE (IOUt,99007) relax , anorm , anfix ,  &
                                 & rnold , rnorm
99007    FORMAT (' RELAXATION FACTOR CORRECTED TO RELAX = ',            &
                &D10.2/' NORM OF SCALED RHS CHANGES FROM ',D10.2,' TO', &
               & D10.2/' NORM   OF   RHS  CHANGES  FROM  ',D10.2,' TO', &
               & D10.2,D10.2)
      ELSE
         IF ( IPRint.LT.0 ) WRITE (IOUt,99008) ITEr , relax , anorm ,   &
                                 & anfix , rnold , rnorm
99008    FORMAT (' ITERATION = ',I3,'  RELAXATION FACTOR = ',           &
                &D10.2/' NORM OF SCALED RHS CHANGES FROM ',D10.2,' TO', &
               & D10.2/' NORM   OF   RHS  CHANGES  FROM  ',D10.2,' TO', &
               & D10.2,D10.2)
      ENDIF
      icor = 0
!
!...       check for monotonic decrease in  delz and deldmz.
!
      IF ( anfix.GE.PREcis .AND. rnorm.GE.PREcis ) THEN
         IF ( anfix.LE.anorm ) THEN
!
!...       we have a decrease.
!...       if  dqz  and dqdmz  small, check for convergence
!
            IF ( anfix.LE.check ) THEN
!
!...       check convergence (iconv = 0).
!
               DO it = 1 , NTOl
                  inz = LTOl(it)
                  DO iz = inz , NZ , MSTar
                     IF ( DABS(Dqz(iz)).GT.TOLin(it)*(DABS(Z(iz))+1.D0) &
                        & ) GOTO 600
                  ENDDO
               ENDDO
!
!...       convergence obtained
!
               IF ( IPRint.LT.1 ) WRITE (IOUt,99015) ITEr
!
!...       since convergence obtained, update  z and dmz  with term
!...       from the fixed jacobian iteration.
!
               DO i = 1 , NZ
                  Z(i) = Z(i) + Dqz(i)
               ENDDO
               DO i = 1 , NDMz
                  Dmz(i) = Dmz(i) + Dqdmz(i)
               ENDDO
               GOTO 900
!
!...       correct the predicted  relax  unless the corrected
!...       value is within 10 percent of the predicted one.
!
            ELSEIF ( ipred.NE.1 ) THEN
               GOTO 600
            ENDIF
         ENDIF
         IF ( ITEr.GE.LIMit ) THEN
            IF ( IPRint.LT.1 ) WRITE (IOUt,99006) ITEr
            GOTO 1100
         ELSE
!
!...       correct the relaxation factor.
!
            ipred = 0
            arg = (anfix/anorm-1.D0)/relax + 1.D0
            IF ( arg.LT.0.D0 ) GOTO 600
            IF ( arg.LE..25D0*relax+.125D0*relax**2 ) THEN
               IF ( relax.GE..9D0 ) GOTO 600
               relax = 1.D0
            ELSE
               factor = -1.D0 + DSQRT(1.D0+8.D0*arg)
               IF ( DABS(factor-1.D0).LT..1D0*factor ) GOTO 600
               IF ( factor.LT.0.5D0 ) factor = 0.5D0
               relax = relax/factor
            ENDIF
            icor = 1
            IF ( relax.LT.relmin ) THEN
               IF ( IPRint.LT.1 ) WRITE (IOUt,99009) relax , relmin
99009          FORMAT (/' NO CONVERGENCE.  RELAXATION FACTOR =',D10.3,  &
                      &' IS TOO SMALL (LESS THAN',D10.3,')'/)
               GOTO 1100
            ELSE
               fact = relax - rlxold
               DO i = 1 , NZ
                  Z(i) = Z(i) + fact*Delz(i)
               ENDDO
               DO i = 1 , NDMz
                  Dmz(i) = Dmz(i) + fact*Deldmz(i)
               ENDDO
               rlxold = relax
               GOTO 800
            ENDIF
         ENDIF
      ENDIF
 900  IF ( (anfix.LT.PREcis .OR. rnorm.LT.PREcis) .AND. IPRint.LT.1 )   &
         & WRITE (IOUt,99015) ITEr
      iconv = 1
      IF ( ICAre.EQ.(-1) ) ICAre = 0
!
!...       if full output has been requested, print values of the
!...       solution components   z  at the meshpoints.
!
 1000 IF ( IPRint.LT.0 ) THEN
         DO j = 1 , MSTar
            WRITE (IOUt,99010) j
99010       FORMAT (' MESH VALUES FOR Z(',I2,'),')
            WRITE (IOUt,99011) (Z(lj),lj=j,NZ,MSTar)
99011       FORMAT (' ',8D15.7)
         ENDDO
      ENDIF
!
!...       check for error tolerance satisfaction
!
      ifin = 1
      IF ( imesh.EQ.2 ) CALL ERRCHK(Xi,Z,Dmz,Valstr,ifin)
      IF ( imesh.EQ.1 .OR. ifin.EQ.0 .AND. ICAre.NE.2 ) GOTO 1200
      Iflag = 1
      RETURN
 1100 Iflag = -2
      noconv = noconv + 1
      IF ( ICAre.EQ.2 .AND. noconv.GT.1 ) RETURN
      IF ( ICAre.EQ.0 ) ICAre = -1
!
!...       update old mesh
!
 1200 np1 = N + 1
      DO i = 1 , np1
         Xiold(i) = Xi(i)
      ENDDO
      NOLd = N
!
!...       pick a new mesh
!...       check safeguards for mesh refinement
!
      imesh = 1
      IF ( iconv.EQ.0 .OR. MSHnum.GE.MSHlmt .OR. MSHalt.GE.MSHlmt )     &
         & imesh = 2
      IF ( MSHalt.GE.MSHlmt .AND. MSHnum.LT.MSHlmt ) MSHalt = 1
      CALL NEWMSH(imesh,Xi,Xiold,Z,Dmz,Valstr,Slope,Accum,Nfxpnt,Fixpnt)
!
!...       exit if expected n is too large (but may try n=nmax once)
!
      IF ( N.LE.NMAx ) THEN
         IF ( iconv.EQ.0 ) imesh = 1
         IF ( ICAre.EQ.1 ) iconv = 0
         GOTO 100
      ENDIF
      N = N/2
      Iflag = -1
      IF ( iconv.EQ.0 .AND. IPRint.LT.1 ) WRITE (IOUt,99012)
99012 FORMAT ('  (NO CONVERGENCE)')
      IF ( iconv.EQ.1 .AND. IPRint.LT.1 ) WRITE (IOUt,99013)
99013 FORMAT ('  (PROBABLY TOLERANCES TOO STRINGENT, OR NMAX TOO ',     &
             &'SMALL)')
      RETURN
99014 FORMAT (' ITERATION = ',I3,'  NORM (RHS) = ',D10.2)
99015 FORMAT (/' CONVERGENCE AFTER',I3,' ITERATIONS'/)
      END
