!*==NEWMSH.spg  processed by SPAG 6.72Dc at 06:28 on  8 Sep 2021
!----------------------------------------------------------------------
!                            p a r t  2
!          mesh selection, error estimation, (and related
!          constant assignment) routines -- see [3], [4]
!----------------------------------------------------------------------
!
      SUBROUTINE NEWMSH(Mode,Xi,Xiold,Z,Dmz,Valstr,Slope,Accum,Nfxpnt,  &
                      & Fixpnt)
!
!**********************************************************************
!
!c NEWMSH selects a mesh on which the solution is to be determined.
!
!   purpose
!            select a mesh on which a collocation solution is to be
!            determined
!
!                           there are 5 possible modes of action:
!            mode = 5,4,3 - deal mainly with definition of an initial
!                           mesh for the current boundary value problem
!                 = 2,1   - deal with definition of a new mesh, either
!                           by simple mesh halving or by mesh selection
!            more specifically, for
!            mode = 5  an initial (generally nonuniform) mesh is
!                      defined by the user and no mesh selection is to
!                      be performed
!                 = 4  an initial (generally nonuniform) mesh is
!                      defined by the user
!                 = 3  a simple uniform mesh (except possibly for some
!                      fixed points) is defined; n= no. of subintervals
!                 = 1  the automatic mesh selection procedure is used
!                      (see [3] for details)
!                 = 2  a simple mesh halving is performed
!
!**********************************************************************
!
!   variables
!
!            n      = number of mesh subintervals
!            nold   = number of subintervals for former mesh
!            xi     - mesh point array
!            xiold  - former mesh point array
!            mshlmt - maximum no. of mesh selections which are permitted
!                     for a given n before mesh halving
!            mshnum - no. of mesh selections which have actually been
!                     performed for the given n
!            mshalt - no. of consecutive times ( plus 1 ) the mesh
!                     selection has alternately halved and doubled n.
!                     if mshalt .ge. mshlmt then  contrl  requires
!                     that the current mesh be halved.
!            mshflg = 1  the mesh is a halving of its former mesh
!                       (so an error estimate has been calculated)
!                   = 0  otherwise
!            iguess - ipar(9) in subroutine colnew.  it is used
!                     here only for mode=5 and 4, where
!                   = 2 the subroutine sets xi=xiold.  this is
!                       used e.g. if continuation is being per-
!                       formed, and a mesh for the old differen-
!                       tial equation is being used
!                   = 3 same as for =2, except xi uses every other
!                       point of xiold (so mesh xiold is mesh xi
!                       halved)
!                   = 4 xi has been defined by the user, and an old
!                       mesh xiold is also available
!                       otherwise, xi has been defined by the user
!                       and we set xiold=xi in this subroutine
!            slope  - an approximate quantity to be equidistributed for
!                     mesh selection (see [3]), viz,
!                             .                        (k+mj)
!                     slope(i)=     max   (weight(l) *u      (xi(i)))
!                               1.le.l.le.ntol         j
!
!                     where j=jtol(l)
!            slphmx - maximum of slope(i)*(xiold(i+1)-xiold(i)) for
!                     i = 1 ,..., nold.
!            accum  - accum(i) is the integral of  slope  from  aleft
!                     to  xiold(i).
!            valstr - is assigned values needed in  errchk  for the
!                     error estimate.
!**********************************************************************
!
      IMPLICIT NONE
!*--NEWMSH84
!*** Start of declarations inserted by SPAG
      REAL*8 accl , accr , Accum , ACOl , ALEft , ARIght , ASAve ,      &
           & avrg , B , d1 , d2 , degequ , Dmz , dummy , dx , Fixpnt ,  &
           & hd6 , hiold , oneovh , PREcis
      REAL*8 ROOt , Slope , slphmx , temp , TOL , TOLin , tsum ,        &
           & Valstr , WGTerr , WGTmsh , x , Xi , Xiold , xleft ,        &
           & xright , Z , ZETa
      INTEGER i , ICAre , IDUm , iflip , IGUess , ileft , in , IOUt ,   &
            & IPRint , iright , ITEr , IZEta , j , jj , JTOl , jz , K , &
            & KD , kstore , l
      INTEGER lcarry , LIMit , lnew , lold , LTOl , M , MMAx , Mode ,   &
            & MSHalt , MSHflg , MSHlmt , MSHnum , MSTar , N , n2 ,      &
            & naccum , NCOmp , NDMz , nfxp1 , Nfxpnt
      INTEGER NMAx , nmax2 , nmin , nmx , NOLd , noldp1 , NONlin , np1 ,&
            & nregn , NTOl , NZ
!*** End of declarations inserted by SPAG
      DIMENSION d1(40) , d2(40) , Slope(1) , Accum(1) , Valstr(1)
      DIMENSION Xi(1) , Xiold(1) , Z(1) , Dmz(1) , Fixpnt(1) , dummy(1)
!
      COMMON /COLOUT/ PREcis , IOUt , IPRint
      COMMON /COLORD/ K , NCOmp , MSTar , KD , MMAx , M(20)
      COMMON /COLAPR/ N , NOLd , NMAx , NZ , NDMz
      COMMON /COLMSH/ MSHflg , MSHnum , MSHlmt , MSHalt
      COMMON /COLNLN/ NONlin , ITEr , LIMit , ICAre , IGUess
      COMMON /COLSID/ ZETa(40) , ALEft , ARIght , IZEta , IDUm
      COMMON /COLBAS/ B(28) , ACOl(28,7) , ASAve(28,4)
      COMMON /COLEST/ TOL(40) , WGTmsh(40) , WGTerr(40) , TOLin(40) ,   &
                    & ROOt(40) , JTOl(40) , LTOl(40) , NTOl
!
      nfxp1 = Nfxpnt + 1
      IF ( Mode.EQ.1 ) THEN
!
!...  mode=1  we do mesh selection if it is deemed worthwhile
!
         IF ( NOLd.EQ.1 ) GOTO 100
         IF ( NOLd.LE.2*Nfxpnt ) GOTO 100
!
!...  the first interval has to be treated separately from the
!...  other intervals (generally the solution on the (i-1)st and ith
!...  intervals will be used to approximate the needed derivative, but
!...  here the 1st and second intervals are used.)
!
         i = 1
         hiold = Xiold(2) - Xiold(1)
         CALL HORDER(1,d1,hiold,Dmz,NCOmp,K)
         hiold = Xiold(3) - Xiold(2)
         CALL HORDER(2,d2,hiold,Dmz,NCOmp,K)
         Accum(1) = 0.D0
         Slope(1) = 0.D0
         oneovh = 2.D0/(Xiold(3)-Xiold(1))
         DO j = 1 , NTOl
            jj = JTOl(j)
            jz = LTOl(j)
            Slope(1) = DMAX1(Slope(1),(DABS(d2(jj)-d1(jj))*WGTmsh(j)*   &
                     & oneovh/(1.D0+DABS(Z(jz))))**ROOt(j))
         ENDDO
         slphmx = Slope(1)*(Xiold(2)-Xiold(1))
         Accum(2) = slphmx
         iflip = 1
!
!...  go through the remaining intervals generating  slope
!...  and  accum .
!
         DO i = 2 , NOLd
            hiold = Xiold(i+1) - Xiold(i)
            IF ( iflip.EQ.-1 ) CALL HORDER(i,d1,hiold,Dmz,NCOmp,K)
            IF ( iflip.EQ.1 ) CALL HORDER(i,d2,hiold,Dmz,NCOmp,K)
            oneovh = 2.D0/(Xiold(i+1)-Xiold(i-1))
            Slope(i) = 0.D0
!
!...       evaluate function to be equidistributed
!
            DO j = 1 , NTOl
               jj = JTOl(j)
               jz = LTOl(j) + (i-1)*MSTar
               Slope(i) = DMAX1(Slope(i),(DABS(d2(jj)-d1(jj))*WGTmsh(j)*&
                        & oneovh/(1.D0+DABS(Z(jz))))**ROOt(j))
            ENDDO
!
!...       accumulate approximate integral of function to be
!...       equidistributed
!
            temp = Slope(i)*(Xiold(i+1)-Xiold(i))
            slphmx = DMAX1(slphmx,temp)
            Accum(i+1) = Accum(i) + temp
            iflip = -iflip
         ENDDO
         avrg = Accum(NOLd+1)/DFLOAT(NOLd)
         degequ = avrg/DMAX1(slphmx,PREcis)
!
!...  naccum=expected n to achieve .1x user requested tolerances
!
         naccum = Accum(NOLd+1) + 1.D0
         IF ( IPRint.LT.0 ) WRITE (IOUt,99001) degequ , naccum
99001    FORMAT (/' MESH SELECTION INFO,'/                              &
                &' DEGREE OF EQUIDISTRIBUTION = ',F8.5,                 &
                &' PREDICTION FOR REQUIRED N =',I8)
!
!...  decide if mesh selection is worthwhile (otherwise, halve)
!
         IF ( avrg.LT.PREcis ) GOTO 100
         IF ( degequ.GE..5D0 ) GOTO 100
!
!...  nmx assures mesh has at least half as many subintervals as the
!...  previous mesh
!
         nmx = MAX0(NOLd+1,naccum)/2
!
!...  this assures that halving will be possible later (for error est)
!
         nmax2 = NMAx/2
!
!...  the mesh is at most halved
!
         N = MIN0(nmax2,NOLd,nmx)
         GOTO 200
      ELSEIF ( Mode.EQ.2 ) THEN
         GOTO 100
      ELSEIF ( Mode.EQ.3 ) THEN
!
!...  mode=3   generate a (piecewise) uniform mesh. if there are
!...  fixed points then ensure that the n being used is large enough.
!
         IF ( N.LT.nfxp1 ) N = nfxp1
         np1 = N + 1
         Xi(1) = ALEft
         ileft = 1
         xleft = ALEft
!
!...  loop over the subregions between fixed points.
!
         DO j = 1 , nfxp1
            xright = ARIght
            iright = np1
            IF ( j.NE.nfxp1 ) THEN
               xright = Fixpnt(j)
!
!...       determine where the j-th fixed point should fall in the
!...       new mesh - this is xi(iright) and the (j-1)st fixed
!...       point is in xi(ileft)
!
               nmin = (xright-ALEft)/(ARIght-ALEft)*DFLOAT(N) + 1.5D0
               IF ( nmin.GT.N-Nfxpnt+j ) nmin = N - Nfxpnt + j
               iright = MAX0(ileft+1,nmin)
            ENDIF
            Xi(iright) = xright
!
!...       generate equally spaced points between the j-1st and the
!...       j-th fixed points.
!
            nregn = iright - ileft - 1
            IF ( nregn.NE.0 ) THEN
               dx = (xright-xleft)/DFLOAT(nregn+1)
               DO i = 1 , nregn
                  Xi(ileft+i) = xleft + DFLOAT(i)*dx
               ENDDO
            ENDIF
            ileft = iright
            xleft = xright
         ENDDO
         GOTO 300
      ELSEIF ( Mode.NE.4 ) THEN
!
!...  mode=5   set mshlmt=1 so that no mesh selection is performed
!
         MSHlmt = 1
      ENDIF
!
!...  mode=4   the user-specified initial mesh is already in place.
!
      IF ( IGUess.GE.2 ) THEN
!
!...  iguess=2, 3 or 4.
!
         noldp1 = NOLd + 1
         IF ( IPRint.LT.1 ) WRITE (IOUt,99002) NOLd ,                   &
                                 & (Xiold(i),i=1,noldp1)
99002    FORMAT (/' THE FORMER MESH (OF',I5,' SUBINTERVALS),',          &
               & 100(/8F12.6))
         IF ( IGUess.EQ.3 ) THEN
!
!...  if iread ( ipar(8) ) .ge. 1 and iguess ( ipar(9) ) .eq. 3
!...  then the first mesh is every second point of the
!...  mesh in  xiold .
!
            N = NOLd/2
            i = 0
            DO j = 1 , NOLd , 2
               i = i + 1
               Xi(i) = Xiold(j)
            ENDDO
         ENDIF
      ENDIF
      np1 = N + 1
      Xi(1) = ALEft
      Xi(np1) = ARIght
      GOTO 300
!
!...  mode=2  halve the current mesh (i.e. double its size)
!
 100  n2 = 2*N
!
!...  check that n does not exceed storage limitations
!
      IF ( n2.LE.NMAx ) THEN
!
!...  calculate the old approximate solution values at
!...  points to be used in  errchk  for error estimates.
!...  if  mshflg  =1 an error estimate was obtained for
!...  for the old approximation so half the needed values
!...  will already be in  valstr .
!
         IF ( MSHflg.EQ.0 ) THEN
!
!...  save in  valstr  the values of the old solution
!...  at the relative positions 1/6, 2/6, 4/6 and 5/6 in
!...  each subinterval.
!
            kstore = 1
            DO i = 1 , N
               x = Xi(i)
               hd6 = (Xi(i+1)-Xi(i))/6.D0
               DO j = 1 , 4
                  x = x + hd6
                  IF ( j.EQ.3 ) x = x + hd6
                  CALL APPROX(i,x,Valstr(kstore),ASAve(1,j),dummy,Xiold,&
                            & NOLd,Z,Dmz,K,NCOmp,MMAx,M,MSTar,4,dummy,0)
                  kstore = kstore + MSTar
               ENDDO
            ENDDO
         ELSE
!
!...  save in  valstr  the values of the old solution
!...  at the relative positions 1/6 and 5/6 in each subinterval.
!
            kstore = 1
            DO i = 1 , NOLd
               hd6 = (Xiold(i+1)-Xiold(i))/6.D0
               x = Xiold(i) + hd6
               CALL APPROX(i,x,Valstr(kstore),ASAve(1,1),dummy,Xiold,   &
                         & NOLd,Z,Dmz,K,NCOmp,MMAx,M,MSTar,4,dummy,0)
               x = x + 4.D0*hd6
               kstore = kstore + 3*MSTar
               CALL APPROX(i,x,Valstr(kstore),ASAve(1,4),dummy,Xiold,   &
                         & NOLd,Z,Dmz,K,NCOmp,MMAx,M,MSTar,4,dummy,0)
               kstore = kstore + MSTar
            ENDDO
         ENDIF
         MSHflg = 0
         MSHnum = 1
         Mode = 2
!
!...  generate the halved mesh.
!
         j = 2
         DO i = 1 , N
            Xi(j) = (Xiold(i)+Xiold(i+1))/2.D0
            Xi(j+1) = Xiold(i+1)
            j = j + 2
         ENDDO
         N = n2
         GOTO 300
!
!...  if possible, try with n=nmax. redistribute first.
!
      ELSEIF ( Mode.EQ.2 ) THEN
         IF ( IPRint.LT.1 ) WRITE (IOUt,99003)
99003    FORMAT (/'  EXPECTED N TOO LARGE ')
         N = n2
         RETURN
      ELSE
         N = NMAx/2
      ENDIF
 200  noldp1 = NOLd + 1
      IF ( N.LT.nfxp1 ) N = nfxp1
      MSHnum = MSHnum + 1
!
!...  if the new mesh is smaller than the old mesh set mshnum
!...  so that the next call to  newmsh  will produce a halved
!...  mesh. if n .eq. nold / 2 increment mshalt so there can not
!...  be an infinite loop alternating between n and n/2 points.
!
      IF ( N.LT.NOLd ) MSHnum = MSHlmt
      IF ( N.GT.NOLd/2 ) MSHalt = 1
      IF ( N.EQ.NOLd/2 ) MSHalt = MSHalt + 1
      MSHflg = 0
!
!...  having decided to generate a new mesh with n subintervals we now
!...  do so, taking into account that the nfxpnt points in the array
!...  fixpnt must be included in the new mesh.
!
      in = 1
      accl = 0.D0
      lold = 2
      Xi(1) = ALEft
      Xi(N+1) = ARIght
      DO i = 1 , nfxp1
         IF ( i.EQ.nfxp1 ) THEN
            accr = Accum(noldp1)
            lnew = noldp1
            nregn = N - in
         ELSE
            DO j = lold , noldp1
               lnew = j
               IF ( Fixpnt(i).LE.Xiold(j) ) GOTO 220
            ENDDO
 220        accr = Accum(lnew) + (Fixpnt(i)-Xiold(lnew))*Slope(lnew-1)
            nregn = (accr-accl)/Accum(noldp1)*DFLOAT(N) - .5D0
            nregn = MIN0(nregn,N-in-nfxp1+i)
            Xi(in+nregn+1) = Fixpnt(i)
         ENDIF
         IF ( nregn.NE.0 ) THEN
            temp = accl
            tsum = (accr-accl)/DFLOAT(nregn+1)
            DO j = 1 , nregn
               in = in + 1
               temp = temp + tsum
               DO l = lold , lnew
                  lcarry = l
                  IF ( temp.LE.Accum(l) ) GOTO 230
               ENDDO
 230           lold = lcarry
               Xi(in) = Xiold(lold-1) + (temp-Accum(lold-1))            &
                      & /Slope(lold-1)
            ENDDO
         ENDIF
         in = in + 1
         accl = accr
         lold = lnew
      ENDDO
      Mode = 1
 300  np1 = N + 1
      IF ( IPRint.LT.1 ) WRITE (IOUt,99004) N , (Xi(i),i=1,np1)
!----------------------------------------------------------------
99004 FORMAT (/' THE NEW MESH (OF',I5,' SUBINTERVALS), ',100(/8F12.6))
      NZ = MSTar*(N+1)
      NDMz = KD*N
      RETURN
      END
