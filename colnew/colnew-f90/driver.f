   subroutine contrl (xi, xiold, z, dmz, rhs, delz, deldmz,
     1           dqz, dqdmz, g, w, v, valstr, slope, scale, dscale,
     2           accum, ipvtg, integs, ipvtw, nfxpnt, fixpnt, iflag,
     3           fsub, dfsub, gsub, dgsub, guess )

c*********************************************************************72
c
cc contrl is the driver for colnew.
c
c   purpose
c     this subroutine is the actual driver.  the nonlinear iteration
c     strategy is controlled here ( see [4] ). upon convergence, errchk
c     is called to test for satisfaction of the requested tolerances.
c
c   variables
c
c     check  - maximum tolerance value, used as part of criteria for
c              checking for nonlinear iteration convergence
c     relax  - the relaxation factor for damped newton iteration
c     relmin - minimum allowable value for relax  (otherwise the
c              jacobian is considered singular).
c     rlxold - previous relax
c     rstart - initial value for relax when problem is sensitive
c     ifrz   - number of fixed jacobian iterations
c     lmtfrz - maximum value for ifrz before performing a reinversion
c     iter   - number of iterations (counted only when jacobian
c              reinversions are performed).
c     xi     - current mesh
c     xiold  - previous mesh
c     ipred  = 0  if relax is determined by a correction
c            = 1  if relax is determined by a prediction
c     ifreez = 0  if the jacobian is to be updated
c            = 1  if the jacobian is currently fixed (frozen)
c     iconv  = 0  if no previous convergence has been obtained
c            = 1  if convergence on a previous mesh has been obtained
c     icare  =-1  no convergence occurred (used for regular problems)
c            = 0  a regular problem
c            = 1  a sensitive problem
c            = 2  used for continuation (see description of ipar(10)
c                 in colnew).
c     rnorm  - norm of rhs (right hand side) for current iteration
c     rnold  - norm of rhs for previous iteration
c     anscl  - scaled norm of newton correction
c     anfix  - scaled norm of newton correction at next step
c     anorm  - scaled norm of a correction obtained with jacobian fixed
c     nz     - number of components of  z  (see subroutine approx)
c     ndmz   - number of components of  dmz  (see subroutine approx)
c     imesh  - a control variable for subroutines newmsh and errchk
c            = 1  the current mesh resulted from mesh selection
c                 or is the initial mesh.
c            = 2  the current mesh resulted from doubling the
c                 previous mesh
c
c**********************************************************************
c
      implicit real*8 (a-h,o-z)
      dimension xi(1), xiold(1), z(1), dmz(1), rhs(1)
      dimension g(1), w(1), v(1), valstr(1), slope(1), accum(1)
      dimension delz(1), deldmz(1), dqz(1), dqdmz(1) , fixpnt(1)
      dimension dummy(1), scale(1), dscale(1)
      dimension integs(1), ipvtg(1), ipvtw(1)
c
      common /colout/ precis, iout, iprint
      common /colord/ k, ncomp, mstar, kd, mmax, m(20)
      common /colapr/ n, nold, nmax, nz, ndmz
      common /colmsh/ mshflg, mshnum, mshlmt, mshalt
      common /colsid/ zeta(40), aleft, aright, izeta, idum
      common /colnln/ nonlin, iter, limit, icare, iguess
      common /colest/ tol(40), wgtmsh(40), wgterr(40), tolin(40),
     1                root(40), jtol(40), ltol(40), ntol
c
      external fsub, dfsub, gsub, dgsub, guess
c
c...  constants for control of nonlinear iteration
c
      relmin = 1.d-3
      rstart = 1.d-2
      lmtfrz = 4
c
c...  compute the maximum tolerance
c
      check = 0.d0
      do 10 i = 1, ntol
   10   check = dmax1 ( tolin(i), check )
      imesh = 1
      iconv = 0
      if ( nonlin .eq. 0 ) iconv = 1
      icor = 0
      noconv = 0
      msing = 0
c
c...  the main iteration begins here .
c...  loop 20 is executed until error tolerances are satisfied or
c...  the code fails (due to a singular matrix or storage limitations)
c
   20      continue
c
c...       initialization for a new mesh
c
           iter = 0
           if ( nonlin .gt. 0 )                     go to 50
c
c...       the linear case.
c...       set up and solve equations
c
           call lsyslv (msing, xi, xiold, dummy, dummy, z, dmz, g,
     1          w, v, rhs, dummy, integs, ipvtg, ipvtw, rnorm, 0,
     2          fsub, dfsub, gsub, dgsub, guess )
c
c...       check for a singular matrix
c
           if ( msing .eq. 0 )                      go to 400
   30      if ( msing .lt. 0 )                      go to 40
           if ( iprint .lt. 1 )  write (iout,495)
           go to 460
   40      if ( iprint .lt. 1 )  write (iout,490)
           iflag = 0
           return
c
c...       iteration loop for nonlinear case
c...       define the initial relaxation parameter (= relax)
c
   50      relax = 1.d0
c
c...       check for previous convergence and problem sensitivity
c
           if ( icare .eq. 1 .or. icare .eq. (-1) )  relax = rstart
           if ( iconv .eq. 0 )                      go to 160
c
c...       convergence on a previous mesh has been obtained.    thus
c...       we have a very good initial approximation for the newton
c...       process.    proceed with one full newton and then iterate
c...       with a fixed jacobian.
c
           ifreez = 0
c
c...       evaluate right hand side and its norm  and
c...       find the first newton correction
c
           call lsyslv (msing, xi, xiold, z, dmz, delz, deldmz, g,
     1          w, v, rhs, dqdmz, integs, ipvtg, ipvtw, rnold, 1,
     2          fsub, dfsub, gsub, dgsub, guess )
c
           if ( iprint .lt. 0 )  write(iout,530)
           if ( iprint .lt. 0 )  write (iout,510) iter, rnold
           go to 70
c
c...       solve for the next iterate .
c...       the value of ifreez determines whether this is a full
c...       newton step (=0) or a fixed jacobian iteration (=1).
c
   60      if ( iprint .lt. 0 )  write (iout,510) iter, rnorm
           rnold = rnorm
           call lsyslv (msing, xi, xiold, z, dmz, delz, deldmz, g,
     1          w, v, rhs, dummy, integs, ipvtg, ipvtw, rnorm,
     2          3+ifreez, fsub, dfsub, gsub, dgsub, guess )
c
c...       check for a singular matrix
c
   70      if ( msing .ne. 0 )                      go to 30
           if ( ifreez .eq. 1 )                     go to 80
c
c...       a full newton step
c
           iter = iter + 1
           ifrz = 0
   80      continue
c
c...       update   z and dmz , compute new  rhs  and its norm
c
           do 90 i = 1, nz
             z(i) = z(i) + delz(i)
   90      continue
           do 100 i = 1, ndmz
             dmz(i) = dmz(i) + deldmz(i)
  100      continue
           call lsyslv (msing, xi, xiold, z, dmz, delz, deldmz, g,
     1          w, v, rhs, dummy, integs, ipvtg, ipvtw, rnorm, 2,
     2          fsub, dfsub, gsub, dgsub, guess )
c
c...       check monotonicity. if the norm of  rhs  gets smaller,
c...       proceed with a fixed jacobian; else proceed cautiously,
c...       as if convergence has not been obtained before (iconv=0).
c
           if ( rnorm .lt. precis )                 go to 390
           if ( rnorm .gt. rnold )                  go to 130
           if ( ifreez .eq. 1 )                     go to 110
           ifreez = 1
           go to 60
c
c...       verify that the linear convergence with fixed jacobian
c...       is fast enough.
c
  110      ifrz = ifrz + 1
           if ( ifrz .ge. lmtfrz )       ifreez = 0
           if ( rnold .lt. 4.d0*rnorm )  ifreez = 0
c
c...       check convergence (iconv = 1).
c
           do 120 it = 1, ntol
             inz = ltol(it)
             do 120 iz = inz, nz, mstar
               if ( dabs(delz(iz)) .gt.
     1           tolin(it) * (dabs(z(iz)) + 1.d0))  go to 60
  120      continue
c
c...       convergence obtained
c
           if ( iprint .lt. 1 )  write (iout,560) iter
           go to 400
c
c...      convergence of fixed jacobian iteration failed.
c
  130      if ( iprint .lt. 0 )  write (iout,510) iter, rnorm
           if ( iprint .lt. 0 )  write (iout,540)
           iconv = 0
           relax = rstart
           do 140 i = 1, nz
             z(i) = z(i) - delz(i)
  140      continue
           do 150 i = 1, ndmz
             dmz(i) = dmz(i) - deldmz(i)
  150      continue
c
c...       update old mesh
c
           np1 = n + 1
           do 155 i = 1, np1
  155        xiold(i) = xi(i)
           nold = n
c
           iter = 0
c
c...       no previous convergence has been obtained. proceed
c...       with the damped newton method.
c...       evaluate rhs and find the first newton correction.
c
  160      if(iprint .lt. 0)  write (iout,500)
           call lsyslv (msing, xi, xiold, z, dmz, delz, deldmz, g,
     1          w, v, rhs, dqdmz, integs, ipvtg, ipvtw, rnold, 1,
     2          fsub, dfsub, gsub, dgsub, guess )
c
c...       check for a singular matrix
c
           if ( msing .ne. 0 )                       go to 30
c
c...       bookkeeping for first mesh
c
           if ( iguess .eq. 1 )  iguess = 0
c
c...       find initial scaling
c
           call skale (n, mstar, kd, z, xi, scale, dscale)
           go to 220
c
c...       main iteration loop
c
  170      rnold = rnorm
           if ( iter .ge. limit )                   go to 430
c
c...       update scaling
c
           call skale (n, mstar, kd, z, xi, scale, dscale)
c
c...       compute norm of newton correction with new scaling
c
           anscl = 0.d0
           do 180 i = 1, nz
             anscl = anscl + (delz(i) * scale(i))**2
  180      continue
           do 190 i = 1, ndmz
             anscl = anscl + (deldmz(i) * dscale(i))**2
  190      continue
           anscl = dsqrt(anscl / dfloat(nz+ndmz))
c
c...       find a newton direction
c
           call lsyslv (msing, xi, xiold, z, dmz, delz, deldmz, g,
     1          w, v, rhs, dummy, integs, ipvtg, ipvtw, rnorm, 3,
     2          fsub, dfsub, gsub, dgsub, guess )
c
c...       check for a singular matrix
c
           if ( msing .ne. 0 )                      go to 30
c
c...       predict relaxation factor for newton step.
c
           andif = 0.d0
           do 200 i = 1, nz
             andif = andif + ((dqz(i) - delz(i)) * scale(i))**2
  200      continue
           do 210 i = 1, ndmz
             andif = andif + ((dqdmz(i) - deldmz(i)) * dscale(i))**2
  210      continue
           andif = dsqrt(andif/dfloat(nz+ndmz) + precis)
           relax = relax * anscl / andif
           if ( relax .gt. 1.d0 )  relax = 1.d0
  220      rlxold = relax
           ipred = 1
           iter = iter + 1
c
c...       determine a new  z and dmz  and find new  rhs  and its norm
c
           do 230 i = 1, nz
             z(i) = z(i)  +  relax * delz(i)
  230      continue
           do 240 i = 1, ndmz
             dmz(i) = dmz(i)  +  relax * deldmz(i)
  240      continue
  250      call lsyslv (msing, xi, xiold, z, dmz, dqz, dqdmz, g,
     1          w, v, rhs, dummy, integs, ipvtg, ipvtw, rnorm, 2,
     2          fsub, dfsub, gsub, dgsub, guess )
c
c...       compute a fixed jacobian iterate (used to control relax)
c
           call lsyslv (msing, xi, xiold, z, dmz, dqz, dqdmz, g,
     1          w, v, rhs, dummy, integs, ipvtg, ipvtw, rnorm, 4,
     2          fsub, dfsub, gsub, dgsub, guess )
c
c...       find scaled norms of various terms used to correct relax
c
           anorm = 0.d0
           anfix = 0.d0
           do 260 i = 1, nz
             anorm = anorm  +  (delz(i) * scale(i))**2
             anfix = anfix  +  (dqz(i) * scale(i))**2
  260      continue
           do 270 i = 1, ndmz
             anorm = anorm  +  (deldmz(i) * dscale(i))**2
             anfix = anfix  +  (dqdmz(i) * dscale(i))**2
  270      continue
           anorm = dsqrt(anorm / dfloat(nz+ndmz))
           anfix = dsqrt(anfix / dfloat(nz+ndmz))
           if ( icor .eq. 1 )                         go to 280
           if (iprint .lt. 0)  write (iout,520) iter, relax, anorm,
     1           anfix, rnold, rnorm
           go to 290
  280      if (iprint .lt. 0) write (iout,550) relax, anorm, anfix,
     1           rnold, rnorm
  290      icor = 0
c
c...       check for monotonic decrease in  delz and deldmz.
c
           if (anfix.lt.precis .or. rnorm.lt.precis)  go to 390
           if ( anfix .gt. anorm )                    go to 300
c
c...       we have a decrease.
c...       if  dqz  and dqdmz  small, check for convergence
c
           if ( anfix .le. check )                    go to 350
c
c...       correct the predicted  relax  unless the corrected
c...       value is within 10 percent of the predicted one.
c
           if ( ipred .ne. 1 )                        go to 170
  300      if ( iter .ge. limit )                     go to 430
c
c...       correct the relaxation factor.
c
           ipred = 0
           arg = (anfix/anorm - 1.d0) / relax + 1.d0
           if ( arg .lt. 0.d0 )                       go to 170
           if (arg .le. .25d0*relax+.125d0*relax**2)  go to 310
           factor = -1.d0 + dsqrt (1.d0+8.d0 * arg)
           if ( dabs(factor-1.d0) .lt. .1d0*factor )  go to 170
           if ( factor .lt. 0.5d0 )  factor = 0.5d0
           relax = relax / factor
           go to 320
  310      if ( relax .ge. .9d0 )                     go to 170
           relax = 1.d0
  320      icor = 1
           if ( relax .lt. relmin )                   go to 440
           fact = relax - rlxold
           do 330 i = 1, nz
            z(i) = z(i)  +  fact * delz(i)
  330      continue
           do 340 i = 1, ndmz
             dmz(i) = dmz(i)  +  fact * deldmz(i)
  340      continue
           rlxold = relax
           go to 250
c
c...       check convergence (iconv = 0).
c
  350      continue
           do 360 it = 1, ntol
             inz = ltol(it)
             do 360 iz = inz, nz, mstar
               if ( dabs(dqz(iz)) .gt.
     1         tolin(it) * (dabs(z(iz)) + 1.d0) )   go to 170
  360      continue
c
c...       convergence obtained
c
           if ( iprint .lt. 1 )  write (iout,560) iter
c
c...       since convergence obtained, update  z and dmz  with term
c...       from the fixed jacobian iteration.
c
           do 370 i = 1, nz
             z(i) = z(i)  +  dqz(i)
  370      continue
           do 380 i = 1, ndmz
             dmz(i) = dmz(i)  +  dqdmz(i)
  380      continue
  390      if ( (anfix .lt. precis .or. rnorm .lt. precis)
     1          .and. iprint .lt. 1 )  write (iout,560) iter
           iconv = 1
           if ( icare .eq. (-1) )  icare = 0
c
c...       if full output has been requested, print values of the
c...       solution components   z  at the meshpoints.
c
  400      if ( iprint .ge. 0 )                     go to 420
           do 410 j = 1, mstar
             write(iout,610) j
  410      write(iout,620) (z(lj), lj = j, nz, mstar)
c
c...       check for error tolerance satisfaction
c
  420      ifin = 1
           if (imesh .eq. 2) call errchk (xi, z, dmz, valstr, ifin)
           if ( imesh .eq. 1 .or.
     1          ifin .eq. 0 .and. icare .ne. 2)     go to 460
           iflag = 1
           return
c
c...       diagnostics for failure of nonlinear iteration.
c
  430      if ( iprint .lt. 1 )  write (iout,570) iter
           go to 450
  440      if( iprint .lt. 1 )  write(iout,580) relax, relmin
  450      iflag = -2
           noconv = noconv + 1
           if ( icare .eq. 2 .and. noconv .gt. 1 )  return
           if ( icare .eq. 0 )  icare = -1
c
c...       update old mesh
c
  460      np1 = n + 1
           do 470 i = 1, np1
  470        xiold(i) = xi(i)
           nold = n
c
c...       pick a new mesh
c...       check safeguards for mesh refinement
c
           imesh = 1
           if ( iconv .eq. 0 .or. mshnum .ge. mshlmt
     1                       .or. mshalt .ge. mshlmt )  imesh = 2
           if ( mshalt .ge. mshlmt .and.
     1          mshnum .lt. mshlmt )  mshalt = 1
           call newmsh (imesh, xi, xiold, z, dmz, valstr,
     1                  slope, accum, nfxpnt, fixpnt)
c
c...       exit if expected n is too large (but may try n=nmax once)
c
           if ( n .le. nmax )                       go to 480
           n = n / 2
           iflag = -1
           if ( iconv .eq. 0 .and. iprint .lt. 1 )  write (iout,590)
           if ( iconv .eq. 1 .and. iprint .lt. 1 )  write (iout,600)
           return
  480      if ( iconv .eq. 0 )  imesh = 1
           if ( icare .eq. 1 )  iconv = 0
      go to 20
c     ---------------------------------------------------------------
  490 format(//35h the global bvp-matrix is singular )
  495 format(//40h a local elimination matrix is singular )
  500 format(/30h full damped newton iteration,)
  510 format(13h iteration = , i3, 15h  norm (rhs) = , d10.2)
  520 format(13h iteration = ,i3,22h  relaxation factor = ,d10.2
     1       /33h norm of scaled rhs changes from ,d10.2,3h to,d10.2
     2       /33h norm   of   rhs  changes  from  ,d10.2,3h to,d10.2,
     2       d10.2)
  530 format(/27h fixed jacobian iterations,)
  540 format(/35h switch to damped newton iteration,)
  550 format(40h relaxation factor corrected to relax = , d10.2
     1       /33h norm of scaled rhs changes from ,d10.2,3h to,d10.2
     2       /33h norm   of   rhs  changes  from  ,d10.2,3h to,d10.2
     2       ,d10.2)
  560 format(/18h convergence after , i3,11h iterations /)
  570 format(/22h no convergence after , i3, 11h iterations/)
  580 format(/37h no convergence.  relaxation factor =,d10.3
     1       ,24h is too small (less than, d10.3, 1h)/)
  590 format(18h  (no convergence) )
  600 format(50h  (probably tolerances too stringent, or nmax too
     1       ,6hsmall) )
  610 format( 19h mesh values for z(, i2, 2h), )
  620 format(1h , 8d15.7)
      end
      subroutine skale (n, mstar, kd, z, xi, scale, dscale)
c
c**********************************************************************
c
cc skale provides a scaling for the state variables.
c
c   purpose
c            provide a proper scaling of the state variables, used
c            to control the damping factor for a newton iteration [2].
c
c   variables
c
c            n      = number of mesh subintervals
c            mstar  = number of unknomns in z(u(x))
c            kd     = number of unknowns in dmz
c            z      = the global unknown vector
c            xi     = the current mesh
c            scale  = scaling vector for z
c            dscale = scaling vector for dmz
c
c**********************************************************************
c
      implicit real*8 (a-h,o-z)
      dimension z(mstar,1), scale(mstar,1), dscale(kd,1)
      dimension xi(1), basm(5)
c
      common /colord/ k, ncomp, id1, id2, mmax, m(20)
c
      basm(1) = 1.d0
      do 50 j=1,n
        iz = 1
        h = xi(j+1) - xi(j)
        do 10 l = 1, mmax
          basm(l+1) = basm(l) * h / dfloat(l)
  10    continue
        do 40 icomp = 1, ncomp
          scal = (dabs(z(iz,j)) + dabs(z(iz,j+1))) * .5d0 + 1.d0
          mj = m(icomp)
          do 20 l = 1, mj
            scale(iz,j) = basm(l) / scal
            iz = iz + 1
  20      continue
          scal = basm(mj+1) / scal
          do 30 idmz = icomp, kd, ncomp
            dscale(idmz,j) = scal
  30      continue
  40    continue
  50  continue
      np1 = n + 1
      do 60 iz = 1, mstar
        scale(iz,np1) = scale(iz,n)
  60  continue
      return
      end
