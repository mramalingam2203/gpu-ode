c----------------------------------------------------------------------
c                            p a r t  2
c          mesh selection, error estimation, (and related
c          constant assignment) routines -- see [3], [4]
c----------------------------------------------------------------------
c
      subroutine newmsh (mode, xi, xiold, z, dmz, valstr,
     1                   slope, accum, nfxpnt, fixpnt)
c
c**********************************************************************
c
cc newmsh selects a mesh on which the solution is to be determined.
c
c   purpose
c            select a mesh on which a collocation solution is to be
c            determined
c
c                           there are 5 possible modes of action:
c            mode = 5,4,3 - deal mainly with definition of an initial
c                           mesh for the current boundary value problem
c                 = 2,1   - deal with definition of a new mesh, either
c                           by simple mesh halving or by mesh selection
c            more specifically, for
c            mode = 5  an initial (generally nonuniform) mesh is
c                      defined by the user and no mesh selection is to
c                      be performed
c                 = 4  an initial (generally nonuniform) mesh is
c                      defined by the user
c                 = 3  a simple uniform mesh (except possibly for some
c                      fixed points) is defined; n= no. of subintervals
c                 = 1  the automatic mesh selection procedure is used
c                      (see [3] for details)
c                 = 2  a simple mesh halving is performed
c
c**********************************************************************
c
c   variables
c
c            n      = number of mesh subintervals
c            nold   = number of subintervals for former mesh
c            xi     - mesh point array
c            xiold  - former mesh point array
c            mshlmt - maximum no. of mesh selections which are permitted
c                     for a given n before mesh halving
c            mshnum - no. of mesh selections which have actually been
c                     performed for the given n
c            mshalt - no. of consecutive times ( plus 1 ) the mesh
c                     selection has alternately halved and doubled n.
c                     if mshalt .ge. mshlmt then  contrl  requires
c                     that the current mesh be halved.
c            mshflg = 1  the mesh is a halving of its former mesh
c                       (so an error estimate has been calculated)
c                   = 0  otherwise
c            iguess - ipar(9) in subroutine colnew.  it is used
c                     here only for mode=5 and 4, where
c                   = 2 the subroutine sets xi=xiold.  this is
c                       used e.g. if continuation is being per-
c                       formed, and a mesh for the old differen-
c                       tial equation is being used
c                   = 3 same as for =2, except xi uses every other
c                       point of xiold (so mesh xiold is mesh xi
c                       halved)
c                   = 4 xi has been defined by the user, and an old
c                       mesh xiold is also available
c                       otherwise, xi has been defined by the user
c                       and we set xiold=xi in this subroutine
c            slope  - an approximate quantity to be equidistributed for
c                     mesh selection (see [3]), viz,
c                             .                        (k+mj)
c                     slope(i)=     max   (weight(l) *u      (xi(i)))
c                               1.le.l.le.ntol         j
c
c                     where j=jtol(l)
c            slphmx - maximum of slope(i)*(xiold(i+1)-xiold(i)) for
c                     i = 1 ,..., nold.
c            accum  - accum(i) is the integral of  slope  from  aleft
c                     to  xiold(i).
c            valstr - is assigned values needed in  errchk  for the
c                     error estimate.
c**********************************************************************
c
      implicit real*8 (a-h,o-z)
      dimension d1(40), d2(40), slope(1), accum(1), valstr(1)
      dimension xi(1), xiold(1), z(1), dmz(1), fixpnt(1), dummy(1)
c
      common /colout/ precis, iout, iprint
      common /colord/ k, ncomp, mstar, kd, mmax, m(20)
      common /colapr/ n, nold, nmax, nz, ndmz
      common /colmsh/ mshflg, mshnum, mshlmt, mshalt
      common /colnln/ nonlin, iter, limit, icare, iguess
      common /colsid/  zeta(40), aleft, aright, izeta, idum
      common /colbas/ b(28), acol(28,7), asave(28,4)
      common /colest/ tol(40), wgtmsh(40), wgterr(40), tolin(40),
     1                root(40), jtol(40), ltol(40), ntol
c
      nfxp1 = nfxpnt +1
      go to (180, 100, 50, 20, 10), mode
c
c...  mode=5   set mshlmt=1 so that no mesh selection is performed
c
   10 mshlmt = 1
c
c...  mode=4   the user-specified initial mesh is already in place.
c
   20 if ( iguess .lt. 2 )                          go to 40
c
c...  iguess=2, 3 or 4.
c
      noldp1 = nold + 1
      if (iprint .lt. 1)  write(iout,360) nold, (xiold(i), i=1,noldp1)
      if ( iguess .ne. 3 )                          go to 40
c
c...  if iread ( ipar(8) ) .ge. 1 and iguess ( ipar(9) ) .eq. 3
c...  then the first mesh is every second point of the
c...  mesh in  xiold .
c
      n = nold /2
      i = 0
      do 30 j = 1, nold, 2
           i = i + 1
   30 xi(i) = xiold(j)
   40 continue
      np1 = n + 1
      xi(1) = aleft
      xi(np1) = aright
      go to 320
c
c...  mode=3   generate a (piecewise) uniform mesh. if there are
c...  fixed points then ensure that the n being used is large enough.
c
   50 if ( n .lt. nfxp1 )  n = nfxp1
      np1 = n + 1
      xi(1) = aleft
      ileft = 1
      xleft = aleft
c
c...  loop over the subregions between fixed points.
c
      do 90 j = 1, nfxp1
           xright = aright
           iright = np1
           if ( j .eq. nfxp1 )                      go to 60
           xright = fixpnt(j)
c
c...       determine where the j-th fixed point should fall in the
c...       new mesh - this is xi(iright) and the (j-1)st fixed
c...       point is in xi(ileft)
c
           nmin = (xright-aleft) / (aright-aleft) * dfloat(n) + 1.5d0
           if (nmin .gt. n-nfxpnt+j)  nmin = n - nfxpnt + j
           iright = max0 (ileft+1, nmin)
   60      xi(iright) = xright
c
c...       generate equally spaced points between the j-1st and the
c...       j-th fixed points.
c
           nregn = iright - ileft - 1
           if ( nregn .eq. 0 )                      go to 80
           dx = (xright - xleft) / dfloat(nregn+1)
           do 70 i = 1, nregn
   70      xi(ileft+i) = xleft  +  dfloat(i) * dx
   80      ileft = iright
           xleft = xright
   90 continue
      go to 320
c
c...  mode=2  halve the current mesh (i.e. double its size)
c
  100 n2 = 2 * n
c
c...  check that n does not exceed storage limitations
c
      if ( n2 .le. nmax )                           go to 120
c
c...  if possible, try with n=nmax. redistribute first.
c
      if ( mode .eq. 2 )                            go to 110
      n = nmax / 2
      go to 220
  110 if ( iprint .lt. 1 )  write(iout,370)
      n = n2
      return
c
c...  calculate the old approximate solution values at
c...  points to be used in  errchk  for error estimates.
c...  if  mshflg  =1 an error estimate was obtained for
c...  for the old approximation so half the needed values
c...  will already be in  valstr .
c
  120 if ( mshflg .eq. 0 )                          go to 140
c
c...  save in  valstr  the values of the old solution
c...  at the relative positions 1/6 and 5/6 in each subinterval.
c
      kstore = 1
      do 130 i = 1, nold
          hd6 = (xiold(i+1) - xiold(i)) / 6.d0
          x = xiold(i) + hd6
          call approx (i, x, valstr(kstore), asave(1,1), dummy, xiold,
     1         nold, z, dmz, k, ncomp, mmax, m, mstar, 4, dummy, 0)
          x = x + 4.d0 * hd6
          kstore = kstore  +  3 * mstar
          call approx (i, x, valstr(kstore), asave(1,4), dummy, xiold,
     1         nold, z, dmz, k, ncomp, mmax, m, mstar, 4, dummy, 0)
          kstore = kstore  +  mstar
  130 continue
      go to 160
c
c...  save in  valstr  the values of the old solution
c...  at the relative positions 1/6, 2/6, 4/6 and 5/6 in
c...  each subinterval.
c
  140 kstore = 1
      do 150 i = 1, n
         x = xi(i)
         hd6 = (xi(i+1) - xi(i)) / 6.d0
         do 150 j = 1, 4
           x = x + hd6
           if ( j.eq.3 )  x = x + hd6
           call approx (i, x, valstr(kstore), asave(1,j), dummy, xiold,
     1          nold, z, dmz, k, ncomp, mmax, m, mstar, 4, dummy, 0)
           kstore = kstore  +  mstar
  150 continue
  160 mshflg = 0
      mshnum = 1
      mode = 2
c
c...  generate the halved mesh.
c
      j = 2
      do 170 i = 1, n
           xi(j) = (xiold(i) + xiold(i+1)) / 2.d0
           xi(j+1) = xiold(i+1)
  170 j = j + 2
      n = n2
      go to 320
c
c...  mode=1  we do mesh selection if it is deemed worthwhile
c
  180 if ( nold .eq. 1 )                            go to 100
      if ( nold .le. 2*nfxpnt )                     go to 100
c
c...  the first interval has to be treated separately from the
c...  other intervals (generally the solution on the (i-1)st and ith
c...  intervals will be used to approximate the needed derivative, but
c...  here the 1st and second intervals are used.)
c
      i = 1
      hiold = xiold(2) - xiold(1)
      call horder (1, d1, hiold, dmz, ncomp, k)
      hiold = xiold(3) - xiold(2)
      call horder (2, d2, hiold, dmz, ncomp, k)
      accum(1) = 0.d0
      slope(1) = 0.d0
      oneovh = 2.d0 / ( xiold(3) - xiold(1) )
      do 190 j = 1, ntol
           jj = jtol(j)
           jz = ltol(j)
  190 slope(1) = dmax1(slope(1),(dabs(d2(jj)-d1(jj))*wgtmsh(j)*
     1           oneovh / (1.d0 + dabs(z(jz)))) **root(j))
      slphmx = slope(1) * (xiold(2) - xiold(1))
      accum(2) = slphmx
      iflip = 1
c
c...  go through the remaining intervals generating  slope
c...  and  accum .
c
      do 210 i = 2, nold
           hiold = xiold(i+1) - xiold(i)
           if ( iflip .eq. -1 )
     1                   call horder ( i, d1, hiold, dmz, ncomp, k)
           if ( iflip .eq. 1 )
     1                   call horder ( i, d2, hiold, dmz, ncomp, k)
           oneovh = 2.d0 / ( xiold(i+1) - xiold(i-1) )
           slope(i) = 0.d0
c
c...       evaluate function to be equidistributed
c
           do 200 j = 1, ntol
             jj = jtol(j)
             jz = ltol(j)  +  (i-1) * mstar
  200      slope(i) = dmax1(slope(i),(dabs(d2(jj)-d1(jj))*wgtmsh(j)*
     1                oneovh / (1.d0 + dabs(z(jz)))) **root(j))
c
c...       accumulate approximate integral of function to be
c...       equidistributed
c
           temp = slope(i) * (xiold(i+1)-xiold(i))
           slphmx = dmax1(slphmx,temp)
           accum(i+1) = accum(i) + temp
  210 iflip = - iflip
      avrg = accum(nold+1) / dfloat(nold)
      degequ = avrg / dmax1(slphmx,precis)
c
c...  naccum=expected n to achieve .1x user requested tolerances
c
      naccum = accum(nold+1) + 1.d0
      if ( iprint .lt. 0 )  write(iout,350) degequ, naccum
c
c...  decide if mesh selection is worthwhile (otherwise, halve)
c
      if ( avrg .lt. precis )                       go to 100
      if ( degequ .ge. .5d0 )                       go to 100
c
c...  nmx assures mesh has at least half as many subintervals as the
c...  previous mesh
c
      nmx = max0 ( nold+1, naccum ) / 2
c
c...  this assures that halving will be possible later (for error est)
c
      nmax2 = nmax / 2
c
c...  the mesh is at most halved
c
      n = min0 ( nmax2, nold, nmx )
  220 noldp1 = nold + 1
      if ( n .lt. nfxp1 )  n = nfxp1
      mshnum = mshnum + 1
c
c...  if the new mesh is smaller than the old mesh set mshnum
c...  so that the next call to  newmsh  will produce a halved
c...  mesh. if n .eq. nold / 2 increment mshalt so there can not
c...  be an infinite loop alternating between n and n/2 points.
c
      if ( n .lt. nold )  mshnum = mshlmt
      if ( n .gt. nold/2 )  mshalt = 1
      if ( n .eq. nold/2 )  mshalt = mshalt + 1
      mshflg = 0
c
c...  having decided to generate a new mesh with n subintervals we now
c...  do so, taking into account that the nfxpnt points in the array
c...  fixpnt must be included in the new mesh.
c
      in = 1
      accl = 0.d0
      lold = 2
      xi(1) = aleft
      xi(n+1) = aright
      do 310 i = 1, nfxp1
           if ( i .eq. nfxp1 )                      go to 250
           do 230 j = lold, noldp1
             lnew = j
             if ( fixpnt(i) .le. xiold(j) )         go to 240
  230      continue
  240      continue
           accr = accum(lnew) + (fixpnt(i)-xiold(lnew))*slope(lnew-1)
           nregn = (accr-accl) / accum(noldp1) * dfloat(n) - .5d0
           nregn = min0(nregn, n - in - nfxp1 + i)
           xi(in + nregn + 1) = fixpnt(i)
           go to 260
  250      accr = accum(noldp1)
           lnew = noldp1
           nregn = n - in
  260      if ( nregn .eq. 0 )                      go to 300
           temp = accl
           tsum = (accr - accl) / dfloat(nregn+1)
           do 290 j = 1, nregn
             in = in + 1
             temp = temp + tsum
             do 270 l = lold, lnew
               lcarry = l
               if ( temp .le. accum(l) )            go to 280
  270        continue
  280        continue
             lold = lcarry
  290      xi(in) = xiold(lold-1) + (temp - accum(lold-1)) /
     1     slope(lold-1)
  300      in = in + 1
           accl = accr
           lold = lnew
  310 continue
      mode = 1
  320 continue
      np1 = n + 1
      if ( iprint .lt. 1 )  write(iout,340) n, (xi(i),i=1,np1)
      nz   = mstar * (n + 1)
      ndmz = kd * n
      return
c----------------------------------------------------------------
  340 format(/17h the new mesh (of,i5,16h subintervals), ,100(/8f12.6))
  350 format(/21h mesh selection info,/30h degree of equidistribution =
     1       , f8.5, 28h prediction for required n = , i8)
  360 format(/20h the former mesh (of,i5,15h subintervals),,
     1       100(/8f12.6))
  370 format (/23h  expected n too large  )
      end
      subroutine consts (k, rho, coef)
c
c**********************************************************************
c
cc consts assigns values to various array constants.
c
c   purpose
c            assign (once) values to various array constants.
c
c   arrays assigned during compilation:
c     cnsts1 - weights for extrapolation error estimate
c     cnsts2 - weights for mesh selection
c              (the above weights come from the theoretical form for
c              the collocation error -- see [3])
c
c   arrays assigned during execution:
c     wgterr - the particular values of cnsts1 used for current run
c              (depending on k, m)
c     wgtmsh - gotten from the values of cnsts2 which in turn are
c              the constants in the theoretical expression for the
c              errors. the quantities in wgtmsh are 10x the values
c              in cnsts2 so that the mesh selection algorithm
c              is aiming for errors .1x as large as the user
c              requested tolerances.
c     jtol   - components of differential system to which tolerances
c              refer (viz, if ltol(i) refers to a derivative of u(j),
c              then jtol(i)=j)
c     root   - reciprocals of expected rates of convergence of compo-
c              nents of z(j) for which tolerances are specified
c     rho    - the k collocation points on (0,1)
c     coef   -
c     acol  -  the runge-kutta coefficients values at collocation
c              points
c
c**********************************************************************
c
      implicit real*8 (a-h,o-z)
      dimension rho(7), coef(k,1), cnsts1(28), cnsts2(28), dummy(1)
c
      common /colord/ kdum, ncomp, mstar, kd, mmax, m(20)
      common /colbas/ b(28), acol(28,7), asave(28,4)
      common /colest/ tol(40), wgtmsh(40), wgterr(40), tolin(40),
     1                root(40), jtol(40), ltol(40), ntol
c
      data cnsts1 /    .25d0,     .625d-1,  7.2169d-2, 1.8342d-2,
     1     1.9065d-2, 5.8190d-2, 5.4658d-3, 5.3370d-3, 1.8890d-2,
     2     2.7792d-2, 1.6095d-3, 1.4964d-3, 7.5938d-3, 5.7573d-3,
     3     1.8342d-2, 4.673d-3,  4.150d-4,  1.919d-3,  1.468d-3,
     4     6.371d-3,  4.610d-3,  1.342d-4,  1.138d-4,  4.889d-4,
     5     4.177d-4,  1.374d-3,  1.654d-3,  2.863d-3  /
      data cnsts2 /   1.25d-1,   2.604d-3,  8.019d-3,  2.170d-5,
     1     7.453d-5,  5.208d-4,  9.689d-8,  3.689d-7,  3.100d-6,
     2     2.451d-5,  2.691d-10, 1.120d-9,  1.076d-8,  9.405d-8,
     3     1.033d-6,  5.097d-13, 2.290d-12, 2.446d-11, 2.331d-10,
     4     2.936d-9,  3.593d-8,  7.001d-16, 3.363d-15, 3.921d-14,
     5     4.028d-13, 5.646d-12, 7.531d-11, 1.129d-9  /
c
c...  assign weights for error estimate
c
      koff = k * ( k + 1 ) / 2
      iz = 1
      do 10 j = 1, ncomp
           mj = m(j)
           do 10 l = 1, mj
             wgterr(iz) = cnsts1(koff - mj + l)
             iz = iz + 1
   10 continue
c
c...  assign array values for mesh selection: wgtmsh, jtol, and root
c
      jcomp = 1
      mtot = m(1)
      do 40 i = 1, ntol
           ltoli = ltol(i)
   20      continue
           if ( ltoli .le. mtot )                   go to 30
           jcomp = jcomp + 1
           mtot = mtot + m(jcomp)
           go to 20
   30      continue
           jtol(i) = jcomp
           wgtmsh(i) = 1.d1 * cnsts2(koff+ltoli-mtot) / tolin(i)
           root(i) = 1.d0 / dfloat(k+mtot-ltoli+1)
   40 continue
c
c...  specify collocation points
c
      go to (50,60,70,80,90,100,110), k
   50 rho(1) = 0.d0
      go to 120
   60 rho(2) = .57735026918962576451d0
      rho(1) = - rho(2)
      go to 120
   70 rho(3) = .77459666924148337704d0
      rho(2) = .0d0
      rho(1) = - rho(3)
      go to 120
   80 rho(4) = .86113631159405257523d0
      rho(3) = .33998104358485626480d0
      rho(2) = - rho(3)
      rho(1) = - rho(4)
      go to 120
   90 rho(5) = .90617984593866399280d0
      rho(4) = .53846931010568309104d0
      rho(3) = .0d0
      rho(2) = - rho(4)
      rho(1) = - rho(5)
      go to 120
  100 rho(6) = .93246951420315202781d0
      rho(5) = .66120938646626451366d0
      rho(4) = .23861918608319690863d0
      rho(3) = -rho(4)
      rho(2) = -rho(5)
      rho(1) = -rho(6)
      go to 120
  110 rho(7) = .949107991234275852452d0
      rho(6) = .74153118559939443986d0
      rho(5) = .40584515137739716690d0
      rho(4) = 0.d0
      rho(3) = -rho(5)
      rho(2) = -rho(6)
      rho(1) = -rho(7)
  120 continue
c
c...  map (-1,1) to (0,1) by  t = .5 * (1. + x)
c
      do 130 j = 1, k
         rho(j) = .5d0 * (1.d0 + rho(j))
  130 continue
c
c...  now find runge-kutta coeffitients b, acol and asave
c...  the values of asave are to be used in  newmsh  and errchk .
c
      do 140 j = 1, k
         do 135 i = 1, k
  135      coef(i,j) = 0.d0
         coef(j,j) = 1.d0
         call vmonde (rho, coef(1,j), k)
  140 continue
      call rkbas ( 1.d0, coef, k, mmax, b, dummy, 0)
      do 150 i = 1, k
         call rkbas ( rho(i), coef, k, mmax, acol(1,i), dummy, 0)
  150 continue
      call rkbas ( 1.d0/6.d0, coef, k, mmax, asave(1,1), dummy, 0)
      call rkbas ( 1.d0/3.d0, coef, k, mmax, asave(1,2), dummy, 0)
      call rkbas ( 2.d0/3.d0, coef, k, mmax, asave(1,3), dummy, 0)
      call rkbas ( 5.d0/6.d0, coef, k, mmax, asave(1,4), dummy, 0)
      return
      end
      subroutine errchk (xi, z, dmz, valstr, ifin)
c
c**********************************************************************
c
cc errchk determines error estimates and tests error tolerances.
c
c      purpose
c               determine the error estimates and test to see if the
c               error tolerances are satisfied.
c
c      variables
c        xi     - current mesh points
c        valstr - values of the previous solution which are needed
c                 for the extrapolation- like error estimate.
c        wgterr - weights used in the extrapolation-like error
c                 estimate. the array values are assigned in
c                 subroutine  consts.
c        errest - storage for error estimates
c        err    - temporary storage used for error estimates
c        z      - approximate solution on mesh xi
c        ifin   - a 0-1 variable. on return it indicates whether
c                 the error tolerances were satisfied
c        mshflg - is set by errchk to indicate to newmsh whether
c                 any values of the current solution are stored in
c                 the array valstr. (0 for no, 1 for yes)
c
c**********************************************************************
c
      implicit real*8 (a-h,o-z)
      dimension err(40), errest(40), dummy(1)
      dimension xi(1), z(1), dmz(1), valstr(1)
c
      common /colout/ precis, iout, iprint
      common /colord/ k, ncomp, mstar, kd, mmax, m(20)
      common /colapr/ n, nold, nmax, nz, ndmz
      common /colmsh/ mshflg, mshnum, mshlmt, mshalt
      common /colbas/ b(28), acol(28,7), asave(28,4)
      common /colest/ tol(40), wgtmsh(40), wgterr(40), tolin(40),
     1                root(40), jtol(40), ltol(40), ntol
c
c...  error estimates are to be generated and tested
c...  to see if the tolerance requirements are satisfied.
c
      ifin = 1
      mshflg = 1
      do 10 j = 1, mstar
   10   errest(j) = 0.d0
      do 60 iback = 1, n
           i = n + 1 - iback
c
c...       the error estimates are obtained by combining values of
c...       the numerical solutions for two meshes.
c...       for each value of iback we will consider the two
c...       approximations at 2 points in each of
c...       the new subintervals.  we work backwards through
c...       the subinterval so that new values can be stored
c...       in valstr in case they prove to be needed later
c...       for an error estimate. the routine  newmsh
c...       filled in the needed values of the old solution
c...       in valstr.
c
           knew = ( 4 * (i-1) + 2 ) * mstar + 1
           kstore = ( 2 * (i-1) + 1 ) * mstar + 1
           x = xi(i) +  (xi(i+1)-xi(i)) * 2.d0 / 3.d0
           call approx (i, x, valstr(knew), asave(1,3), dummy, xi,
     1            n, z, dmz, k, ncomp, mmax, m, mstar, 4, dummy, 0)
           do 20 l = 1,mstar
             err(l) = wgterr(l) * dabs(valstr(knew) -
     1       valstr(kstore))
             knew = knew + 1
             kstore = kstore + 1
   20      continue
           knew = ( 4 * (i-1) + 1 ) * mstar + 1
           kstore = 2 * (i-1) * mstar + 1
           x = xi(i) +  (xi(i+1)-xi(i)) / 3.d0
           call approx (i, x, valstr(knew), asave(1,2), dummy, xi,
     1            n, z, dmz, k, ncomp, mmax, m, mstar, 4, dummy, 0)
           do 30 l = 1,mstar
             err(l) = err(l) + wgterr(l) * dabs(valstr(knew) -
     1       valstr(kstore))
             knew = knew + 1
             kstore = kstore + 1
   30      continue
c
c...       find component-wise maximum error estimate
c
           do 40 l = 1,mstar
             errest(l) = dmax1(errest(l),err(l))
   40      continue
c


c...       test whether the tolerance requirements are satisfied
c...       in the i-th interval.
c
           if ( ifin .eq. 0 )                       go to 60
           do 50 j = 1, ntol
             ltolj = ltol(j)
             ltjz = ltolj  +  (i-1) * mstar
           if ( err(ltolj) .gt.
     1          tolin(j) * (dabs(z(ltjz))+1.d0) )  ifin = 0
   50      continue
   60 continue
      if ( iprint .ge. 0 )                          return
      write(iout,130)
      lj = 1
      do 70 j = 1,ncomp
           mj = lj - 1 + m(j)
           write(iout,120) j, (errest(l), l= lj, mj)
           lj = mj + 1
   70 continue
      return
c--------------------------------------------------------------
  120 format (3h u(, i2, 3h) -,4d12.4)
  130 format (/26h the estimated errors are,)
      end
