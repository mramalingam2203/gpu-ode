c---------------------------------------------------------------------
c                            p a r t  3
c          collocation system setup routines
c---------------------------------------------------------------------
c
      subroutine lsyslv (msing, xi, xiold, z, dmz, delz, deldmz,
     1           g, w, v, rhs, dmzo, integs, ipvtg, ipvtw, rnorm,
     2           mode, fsub, dfsub, gsub, dgsub, guess )
c*********************************************************************
c
cc lsyslv controls the solution of linear systems.
c
c   purpose
c         this routine controls the set up and solution of a linear
c      system of collocation equations.
c         the matrix  g  is cast into an almost block diagonal
c      form by an appropriate ordering of the columns and solved
c      using the package of de boor-weiss [5]. the matrix is composed
c      of n blocks. the i-th block has the size
c                  integs(1,i) * integs(2,i).
c      it contains in its last rows the linearized collocation
c      equations, condensed as described in [2],
c      and the linearized side conditions corresponding to
c      the i-th subinterval.  integs(3,i)  steps of gaussian
c      elimination are applied to it to achieve a  partial plu
c      decomposition.  the right hand side vector is put into  rhs
c      and the solution vector is returned in  delz and deldmz.
c
c         lsyslv operates according to one of 5 modes:
c      mode = 0 - set up the collocation matrices  v , w , g
c                 and the right hand side  rhs ,  and solve.
c                 (for linear problems only.)
c      mode = 1 - set up the collocation matrices  v , w , g
c                 and the right hand sides  rhs  and  dmzo ,
c                 and solve. also set up  integs .
c                 (first iteration of nonlinear problems only).
c      mode = 2 - set up  rhs  only and compute its norm.
c      mode = 3 - set up  v, w, g  only and solve system.
c      mode = 4 - perform forward and backward substitution only
c                 (do not set up the matrices nor form the rhs).
c
c   variables
c
c      ig,izeta  - pointers to g,zeta respectively
c                       (necessary to keep track of blocks of g
c                       during matrix manipulations)
c      idmz,irhs,iv,iw - pointers to  rhs,v,w rspectively
c      df    - partial derivatives of f from dfsub
c      rnorm - euclidean norm of rhs
c      lside - number of side conditions in current and previous blocks
c      iguess = 1 when current soln is user specified via  guess
c             = 0 otherwise
c
c*********************************************************************
      implicit real*8 (a-h,o-z)
      dimension  z(1), dmz(1), delz(1), deldmz(1), xi(1), xiold(1)
      dimension  g(1), w(1), v(1),  rhs(1) , dmzo(1), dummy(1)
      dimension  integs(3,1), ipvtg(1), ipvtw(1)
      dimension  zval(40), f(40), dgz(40), dmval(20), df(800), at(28)
c
      common /colout/ precis, iout, iprint
      common /colloc/ rho(7), coef(49)
      common /colord/ k, ncomp, mstar, kd,  mmax, m(20)
      common /colsid/ zeta(40), aleft, aright, izeta, izsave
      common /colapr/ n, nold, nmax, nz, ndmz
      common /colnln/ nonlin, iter, limit, icare, iguess
      common /colbas/ b(28), acol(28,7), asave(28,4)
c
      external dfsub, dgsub
c
      m1 = mode + 1
      go to (10, 30, 30, 30, 310), m1
c
c...  linear problem initialization
c
   10 do 20 i=1,mstar
   20 zval(i) = 0.d0
c
c...  initialization
c
   30 idmz = 1
      idmzo = 1
      irhs = 1
      ig = 1
      iw = 1
      iv = 1
      izeta = 1
      lside = 0
      iold = 1
      ncol = 2 * mstar
      rnorm = 0.d0
      if ( mode .gt. 1 )                            go to 80
c
c...  build integs (describing block structure of matrix)
c
      do 70 i = 1,n
           integs(2,i) = ncol
           if (i .lt. n)                            go to 40
           integs(3,n) = ncol
           lside = mstar
           go to 60
   40      integs(3,i) = mstar
   50      if( lside .eq. mstar )                   go to 60
           if ( zeta(lside+1) .ge. xi(i)+precis )   go to 60
           lside = lside + 1
           go to 50
   60      nrow = mstar + lside
   70 integs(1,i) = nrow
   80 continue
      if ( mode .eq. 2 )                            go to 90
c
c...  zero the matrices to be computed
c
      lw = kd * kd * n
      do 84 l = 1, lw
   84   w(l) = 0.d0
c
c...  the do loop 290 sets up the linear system of equations.
c
  90  continue
      do 290 i=1, n
c
c...       construct a block of  a  and a corresponding piece of  rhs.
c
           xii = xi(i)
           h = xi(i+1) - xi(i)
           nrow = integs(1,i)
c
c...       go thru the ncomp collocation equations and side conditions
c...       in the i-th subinterval
c
  100      if ( izeta .gt. mstar )                  go to 140
           if ( zeta(izeta) .gt. xii + precis )      go to 140
c
c...       build equation for a side condition.
c
           if ( mode .eq. 0 )                       go to 110
           if ( iguess .ne. 1 )                     go to 102
c
c...       case where user provided current approximation
c
           call guess (xii, zval, dmval)
           go to 110
c
c...       other nonlinear case
c
  102      if ( mode .ne. 1 )                       go to 106
           call approx (iold, xii, zval, at, coef, xiold, nold,
     1          z, dmz, k, ncomp, mmax, m, mstar, 2, dummy, 0)
           go to 110
  106      call approx (i, xii, zval, at, dummy, xi, n, z, dmz,
     1                  k, ncomp, mmax, m, mstar, 1, dummy, 0)
  108      if ( mode .eq. 3 )                       go to 120
c
c...       find  rhs  boundary value.
c
  110      call gsub (izeta, zval, gval)
           rhs(ndmz+izeta) = -gval
           rnorm = rnorm + gval**2
           if ( mode .eq. 2 )                       go to 130
c
c...       build a row of  a  corresponding to a boundary point
c
  120      call gderiv (g(ig), nrow, izeta, zval, dgz, 1, dgsub)
  130      izeta = izeta + 1
           go to 100
c
c...       assemble collocation equations
c
  140      do 220 j = 1, k
             hrho = h * rho(j)
             xcol = xii + hrho
c
c...         this value corresponds to a collocation (interior)
c...         point. build the corresponding  ncomp  equations.
c
             if ( mode .eq. 0 )                     go to 200
             if ( iguess .ne. 1 )                   go to 160
c
c...         use initial approximation provided by the user.
c
             call guess (xcol, zval, dmzo(irhs) )
             go to 170
c
c...         find  rhs  values
c
  160        if ( mode .ne. 1 )                     go to 190
             call approx (iold, xcol, zval, at, coef, xiold, nold,
     1            z, dmz, k, ncomp, mmax, m, mstar, 2, dmzo(irhs), 1)
c
  170        call fsub (paramsl,paramsh,xcol, zval, f)
c        type (params) :: paramsl,paramsh

             do 180 jj = 1, ncomp
               value = dmzo(irhs) - f(jj)
               rhs(irhs) = - value
               rnorm = rnorm + value**2
               irhs = irhs + 1
  180        continue
             go to 210
c
c...         evaluate former collocation solution
c
  190        call approx (i, xcol, zval, acol(1,j), coef, xi, n,
     1            z, dmz, k, ncomp, mmax, m, mstar, 4, dummy, 0)
             if ( mode .eq. 3 )                     go to 210
c
c...         fill in  rhs  values (and accumulate its norm).
c
      call fsub (paramsl,paramsh,xcol, zval, f)
!        type (params) :: paramsl,paramsh

             do 195 jj = 1, ncomp
               value = dmz(irhs) - f(jj)
               rhs(irhs) = - value
               rnorm = rnorm + value**2
               irhs = irhs + 1
  195        continue
             go to 220
c
c...         the linear case
c
  200        call fsub (paramsl,paramsh,xcol, zval, rhs(irhs))
!        type (params) :: paramsl,paramsh

             irhs = irhs + ncomp
c
c...         fill in ncomp rows of  w and v
c
  210        call vwblok (xcol, hrho, j, w(iw), v(iv), ipvtw(idmz),
     1       kd, zval, df, acol(1,j), dmzo(idmzo), ncomp, dfsub, msing)
             if ( msing .ne. 0 )                    return
  220      continue
c
c...       build global bvp matrix  g
c
           if ( mode .ne. 2 )
     1      call gblock (h, g(ig), nrow, izeta, w(iw), v(iv), kd,
     2                  dummy, deldmz(idmz), ipvtw(idmz), 1 )
           if ( i .lt. n )                          go to 280
           izsave = izeta
  240      if ( izeta .gt. mstar )                  go to 290
c
c...       build equation for a side condition.
c
           if ( mode .eq. 0 )                       go to 250
           if ( iguess .ne. 1 )                     go to 245
c
c...       case where user provided current approximation
c
           call guess (aright, zval, dmval)
           go to 250
c
c...       other nonlinear case
c
  245      if ( mode .ne. 1 )                       go to 246
           call approx (nold+1, aright, zval, at, coef, xiold, nold,
     1          z, dmz, k, ncomp, mmax, m, mstar, 1, dummy, 0)
           go to 250
  246      call approx (n+1, aright, zval, at, coef, xi, n,
     1          z, dmz, k, ncomp, mmax, m, mstar, 1, dummy, 0)
  248      if ( mode .eq. 3 )                       go to 260
c
c...       find  rhs  boundary value.
c
  250      call gsub (izeta, zval, gval)
           rhs(ndmz+izeta) = - gval
           rnorm = rnorm + gval**2
           if ( mode .eq. 2 )                       go to 270
c
c...       build a row of  a  corresponding to a boundary point
c
  260      call gderiv (g(ig), nrow, izeta+mstar, zval, dgz, 2, dgsub)
  270      izeta = izeta + 1
           go to 240
c
c...       update counters -- i-th block completed
c
  280      ig = ig + nrow * ncol
           iv = iv + kd * mstar
           iw = iw + kd * kd
           idmz = idmz + kd
           if ( mode .eq. 1 )  idmzo = idmzo + kd
  290 continue
c
c...       assembly process completed
c
      if ( mode .eq. 0 .or. mode .eq. 3 )           go to 300
      rnorm = dsqrt(rnorm / dfloat(nz+ndmz))
      if ( mode .ne. 2 )                            go to 300
      return
c
c...  solve the linear system.
c
c...  matrix decomposition
c
  300 call fcblok (g, integs, n, ipvtg, df, msing)
c
c...  check for singular matrix
c
      msing = - msing
      if( msing .ne. 0 )                            return
c
c...  perform forward and backward substitution for mode=4 only.
c
  310 continue
      do 311 l = 1, ndmz
        deldmz(l) = rhs(l)
  311 continue
      iz = 1
      idmz = 1
      iw = 1
      izet = 1
      do 320 i=1, n
         nrow = integs(1,i)
         izeta = nrow + 1 - mstar
         if ( i .eq. n ) izeta = izsave
  322    if ( izet .eq. izeta )                     go to 324
           delz(iz-1+izet) = rhs(ndmz+izet)
           izet = izet + 1
         go to 322
  324    h = xi(i+1) - xi(i)
         call gblock (h, g(1), nrow, izeta, w(iw), v(1), kd,
     1                delz(iz), deldmz(idmz), ipvtw(idmz), 2 )
         iz = iz + mstar
         idmz = idmz + kd
         iw = iw + kd * kd
         if ( i .lt. n )                            go to 320
  326    if ( izet .gt. mstar )                     go to 320
           delz(iz-1+izet) = rhs(ndmz+izet)
           izet = izet + 1
         go to 326
  320 continue
c
c...  perform forward and backward substitution for mode=0,2, or 3.
c
      call sbblok (g, integs, n, ipvtg, delz)
c
c...  finaly find deldmz
c
      call dmzsol (kd, mstar, n, v, delz, deldmz)
c
      if ( mode .ne. 1 )                            return
      do 321 l = 1, ndmz
        dmz(l) = dmzo(l)
  321 continue
      iz = 1
      idmz = 1
      iw = 1
      izet = 1
      do 350 i=1, n
         nrow = integs(1,i)
         izeta = nrow + 1 - mstar
         if ( i .eq. n ) izeta = izsave
  330    if ( izet .eq. izeta )                     go to 340
           z(iz-1+izet) = dgz(izet)
           izet = izet + 1
         go to 330
  340    h = xi(i+1) - xi(i)
         call gblock (h, g(1), nrow, izeta, w(iw), df, kd,
     1                z(iz), dmz(idmz), ipvtw(idmz), 2 )
         iz = iz + mstar
         idmz = idmz + kd
         iw = iw + kd * kd
         if ( i .lt. n )                            go to 350
  342    if ( izet .gt. mstar )                     go to 350
            z(iz-1+izet) = dgz(izet)
            izet = izet + 1
         go to 342
  350 continue
      call sbblok (g, integs, n, ipvtg, z)
c
c...  finaly find dmz
c
      call dmzsol (kd, mstar, n, v, z, dmz)
c
      return
      end
      subroutine gderiv ( gi, nrow, irow, zval, dgz, mode, dgsub)
c
c**********************************************************************
c
cc gderiv constructs a row of the collocation matrix.
c
c   purpose:
c
c      construct a collocation matrix row according to mode:
c      mode = 1  -  a row corresponding to a initial condition
c                   (i.e. at the left end of the subinterval).
c      mode = 2  -  a row corresponding to a final condition.
c
c   variables:
c
c      gi     - the sub-block of the global bvp matrix in
c               which the equations are to be formed.
c      nrow   - no. of rows in gi.
c      irow   - the row in gi to be used for equations.
c      zval   - z(xi)
c      dg     - the derivatives of the side condition.
c
c**********************************************************************
      implicit real*8 (a-h,o-z)
      dimension gi(nrow,1), zval(1), dgz(1), dg(40)
c
      common /colord/ kdum, ndum, mstar, kd, mmax, m(20)
      common /colsid/ zeta(40), aleft, aright, izeta, idum
      common /colnln/ nonlin, iter, limit, icare, iguess
c
c...  zero jacobian dg
c
      do 10 j=1,mstar
   10   dg(j) = 0.d0
c
c...  evaluate jacobian dg
c
      call dgsub (izeta, zval, dg)
c
c...  evaluate  dgz = dg * zval  once for a new mesh
c
      if (nonlin .eq. 0 .or. iter .gt. 0)           go to 30
      dot = 0.d0
      do 20 j = 1, mstar
   20   dot = dot  +  dg(j) * zval(j)
      dgz(izeta) = dot
c
c...  branch according to  m o d e
c
   30 if ( mode .eq. 2 )                            go to 50
c
c...  provide coefficients of the j-th linearized side condition.
c...  specifically, at x=zeta(j) the j-th side condition reads
c...  dg(1)*z(1) + ... +dg(mstar)*z(mstar) + g = 0
c
c
c...  handle an initial condition
c
      do 40 j = 1, mstar
        gi(irow,j) =  dg(j)
   40 gi(irow,mstar+j) = 0.d0
      return
c
c...  handle a final condition
c
   50 do 60 j= 1, mstar
        gi(irow,j) = 0.d0
   60 gi(irow,mstar+j) = dg(j)
      return
      end
      subroutine vwblok (xcol, hrho, jj, wi, vi, ipvtw, kd, zval,
     1                   df, acol, dmzo, ncomp, dfsub, msing)
c
c**********************************************************************
c
cc vwblok constructs a group of ncomp rows of the wi and vi matrices.
c
c   purpose:
c
c      construct a group of  ncomp  rows of the matrices  wi  and  vi.
c      corresponding to an interior collocation point.
c
c
c   variables:
c
c      xcol   - the location of the collocation point.
c      jj     - xcol is the jj-th of k collocation points
c               in the i-th subinterval.
c      wi,vi  - the i-th block of the collocation matrix
c               before parameter condensation.
c      kd     - no. of rows in vi and wi .
c      zval   - z(xcol)
c      df     - the jacobian at xcol .
c      jcomp  - counter for the component being dealt with.
c
c**********************************************************************
      implicit real*8 (a-h,o-z)
      dimension wi(kd,1), vi(kd,1), zval(1), dmzo(1), df(ncomp,1)
      dimension ipvtw(1),  ha(7,4), acol(7,4), basm(5)
c
      common /colord/ k, ncdum, mstar, kdum, mmax, m(20)
      common /colnln/ nonlin, iter, limit, icare, iguess
c
c...  if jj = 1 initialize  wi .
c
      if ( jj .gt. 1 )                              go to 30
      do 10 id = 1, kd
        wi(id,id) = 1.d0
   10 continue
c
c...  calculate local basis
c
   30        fact = 1.d0
             do 150 l=1,mmax
                fact = fact * hrho / dfloat(l)
                basm(l) = fact
                do 150 j=1,k
                   ha(j,l) = fact * acol(j,l)
  150        continue
c
c... zero jacobian
c
      do 40 jcol = 1, mstar
        do 40 ir = 1, ncomp
   40 df(ir,jcol) = 0.d0
c
c...  build ncomp rows for interior collocation point x.
c...  the linear expressions to be constructed are:
c...   (m(id))
c...  u     -  df(id,1)*z(1) - ... - df(id,mstar)*z(mstar)
c...   id
c...  for id = 1 to ncomp.
c
      call dfsub (paramsl,paramsh,xcol, zval, df)
!      type (params) :: paramsl,paramsh

      i0 = (jj-1) * ncomp
      i1 = i0 + 1
      i2 = i0 + ncomp
c
c...  evaluate  dmzo = dmz - df * zval  once for a new mesh
c
      if (nonlin .eq. 0 .or. iter .gt. 0)          go to 60
      do 50 j = 1, mstar
        fact = - zval(j)
        do 50 id = 1, ncomp
          dmzo(i0+id) = dmzo(i0+id)  +  fact * df(id,j)
  50  continue
c
c...  loop over the  ncomp  expressions to be set up for the
c...  current collocation point.
c
   60 do 70 j = 1, mstar
        do 70 id = 1, ncomp
          vi(i0+id,j) = df(id,j)
   70 continue
      jn = 1
      do 140 jcomp = 1, ncomp
         mj = m(jcomp)
         jn = jn + mj
         do 130 l = 1, mj
            jv = jn - l
            jw = jcomp
            do 90 j = 1, k
              ajl = - ha(j,l)
              do 80 iw = i1, i2
                 wi(iw,jw) = wi(iw,jw)  +  ajl * vi(iw,jv)
   80         continue
   90       jw = jw + ncomp
            lp1 = l + 1
            if ( l .eq. mj )                        go to 130
            do 110 ll = lp1, mj
              jdf = jn - ll
              bl = basm(ll-l)
              do 100 iw = i1, i2
                vi(iw,jv) = vi(iw,jv)  +  bl * vi(iw,jdf)
  100         continue
  110       continue
  130    continue
  140 continue
      if ( jj .lt. k )                          return
c
c   ...decompose the wi block and solve for the mstar columns of vi
c
c
c...  do parameter condensation
c
      msing = 0
      call dgefa  (wi, kd, kd, ipvtw, msing)
c
c...   check for singularity
c
      if ( msing .ne. 0 )                         return
      do 250 j= 1,mstar
         call dgesl  (wi, kd, kd, ipvtw, vi(1,j), 0)
  250 continue
      return
      end
      subroutine gblock (h, gi, nrow, irow, wi, vi, kd,
     1                   rhsz, rhsdmz, ipvtw, mode)
c
c**********************************************************************
c
cc gblock constructs certain rows of the collocation matrix.
c
c   purpose:
c
c      construct collocation matrix rows according to mode:
c      mode = 1  -  a group of  mstar    rows corresponding
c                   an interior mesh interval.
c           = 2  -  corresponding right hand side
c
c   variables:
c
c      h      - the  local stepsize.
c      gi     - the sub-block of the collocation matrix in
c               which the equations are to be formed.
c      wi     - the sub-block of noncondensed collocation equations,
c               left-hand side part.
c      vi     - the sub-block of noncondensed collocation equations,
c               right-hand side part.
c      rhsdmz - the inhomogenous term of the uncondensed collocation
c               equations.
c      rhsz   - the inhomogenous term of the condensed collocation
c               equations.
c      nrow   - no. of rows in gi.
c      irow   - the first row in gi to be used for equations.
c
c**********************************************************************
      implicit real*8 (a-h,o-z)
      dimension hb(7,4), basm(5)
      dimension gi(nrow,1), wi(1), vi(kd,1)
      dimension rhsz(1), rhsdmz(1), ipvtw(1)
c
      common /colord/  k, ncomp, mstar, kdum, mmax, m(20)
      common /colbas/ b(7,4), acol(28,7), asave(28,4)
c
c...  compute local basis
c
      fact = 1.d0
      basm(1) = 1.d0
      do 30 l=1,mmax
         fact = fact * h / dfloat(l)
         basm(l+1) = fact
         do 20 j=1,k
   20       hb(j,l) = fact * b(j,l)
   30 continue
c
c...  branch according to  m o d e
c
      go to (40, 110), mode
c
c...  set right gi-block to identity
c
   40 continue
      do 60 j = 1, mstar
        do 50 ir = 1, mstar
          gi(irow-1+ir,j) = 0.d0
   50   gi(irow-1+ir,mstar+j) = 0.d0
   60 gi(irow-1+j,mstar+j) = 1.d0
c
c...  compute the block gi
c
      ir = irow
      do 100 icomp = 1, ncomp
         mj = m(icomp)
         ir = ir + mj
         do 90 l = 1, mj
            id = ir - l
            do 80 jcol = 1, mstar
               ind = icomp
               rsum = 0.d0
               do 70 j = 1, k
                  rsum = rsum  -  hb(j,l) * vi(ind,jcol)
   70          ind = ind + ncomp
               gi(id,jcol) = rsum
   80       continue
            jd = id - irow
            do 85 ll = 1, l
               gi(id,jd+ll) = gi(id,jd+ll) - basm(ll)
   85       continue
   90    continue
  100 continue
      return
c
c...  compute the appropriate piece of  rhsz
c
  110 continue
      call dgesl  (wi, kd, kd, ipvtw, rhsdmz, 0)
      ir = irow
      do 140 jcomp = 1, ncomp
         mj = m(jcomp)
         ir = ir + mj
         do 130 l = 1, mj
            ind = jcomp
            rsum = 0.d0
            do 120 j = 1, k
               rsum = rsum  +  hb(j,l) * rhsdmz(ind)
  120       ind = ind + ncomp
            rhsz(ir-l) = rsum
  130    continue
  140 continue
      return
      end
