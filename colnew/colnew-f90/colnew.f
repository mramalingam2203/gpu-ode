c  This package solves boundary value problems for
c  ordinary differential equations, as described below.
c
c  colnew is a modification of the package colsys by ascher,
c  christiansen and russell [1]. it incorporates a new basis
c  representation replacing b-splines, and improvements for
c  the linear and nonlinear algebraic equation solvers.
c  the package can be referenced as either colnew or colsys.
c
c----------------------------------------------------------------------
c                            p a r t  1
c        main storage allocation and program control subroutines
c----------------------------------------------------------------------
c
      subroutine colnew (ncomp, m, aleft, aright, zeta, ipar, ltol,
     1                   tol, fixpnt, ispace, fspace, iflag,
     2                   fsub, dfsub, gsub, dgsub, guess)

c*********************************************************************72
c
cc colnew solves a multipoint boundary value problem using collocation.
c
c     written by
c                  u. ascher,
c                            department of computer science,
c                            university of british columbia,
c                            vancouver, b. c., canada   v6t 1w5
c                  g. bader,
c                            institut f. angewandte mathematik
c                            university of heidelberg
c                            im neuenheimer feld 294
c                            d-6900 heidelberg 1
c
c**********************************************************************
c
c     purpose
c
c     this package solves a multi-point boundary value
c     problem for a mixed order system of ode-s given by
c
c          (m(i))
c         u       =  f  ( x; z(u(x)) )      i = 1, ... ,ncomp
c          i          i
c
c                                          aleft .lt. x .lt. aright,
c
c
c         g  ( zeta(j); z(u(zeta(j))) ) = 0   j = 1, ... ,mstar
c          j
c                                    mstar = m(1)+m(2)+...+m(ncomp),
c
c
c         where                          t
c               u = (u , u , ... ,u     )  is the exact solution vector
c                     1   2        ncomp
c
c                (mi)
c               u     is the mi=m(i) th  derivative of u
c                i                                      i
c
c                                  (1)        (m1-1)       (mncomp-1)
c               z(u(x)) = ( u (x),u  (x),...,u    (x),...,u      (x) )
c                            1     1          1            ncomp
c
c                f (x,z(u))   is a (generally) nonlinear function of
c                 i
c                             z(u)=z(u(x)).
c
c                g (zeta(j);z(u))  is a (generally) nonlinear function
c                 j
c                               used to represent a boundary condition.
c
c         the boundary points satisfy
c               aleft .le. zeta(1) .le. .. .le. zeta(mstar) .le. aright
c
c         the orders mi of the differential equations satisfy
c                            1 .le. m(i) .le. 4.
c
c
c**********************************************************************
c
c     method
c
c        the method used to approximate the solution u is
c     collocation at gaussian points, requiring m(i)-1 continuous
c     derivatives in the i-th component, i = 1, ..., ncomp.
c     here, k is the number of collocation points (stages) per
c     subinterval and is chosen such that k .ge. max m(i).
c     a runge-kutta-monomial solution representation is utilized.
c
c     references
c
c     [1] u. ascher, j. christiansen and r.d. russell,
c         collocation software for boundary-value odes,
c         acm trans. math software 7 (1981), 209-222.
c         this paper contains examples where use of the code
c         is demonstrated.
c
c     [2] g. bader and u. ascher,
c         a new basis implementation for a mixed order
c         boundary value ode solver,
c         siam j. scient. stat. comput. (1987).
c
c     [3] u. ascher, j. christiansen and r.d. russell,
c         a collocation solver for mixed order
c         systems of boundary value problems,
c         math. comp. 33 (1979), 659-679.
c
c     [4] u. ascher, j. christiansen and r.d. russell,
c         colsys - a collocation code for boundary
c         value problems,
c         lecture notes comp.sc. 76, springer verlag,
c         b. childs et. al. (eds.) (1979), 164-185.
c
c     [5] c. deboor and r. weiss,
c         solveblok: a package for solving almost block diagonal
c         linear systems,
c         acm trans. math. software 6 (1980), 80-87.
c
c**********************************************************************
c
c     ***************     input to colnew     ***************
c
c     variables
c
c     ncomp - no. of differential equations   (ncomp .le. 20)
c
c     m(j) - order of the j-th differential equation
c            ( mstar = m(1) + ... + m(ncomp) .le. 40 )
c
c     aleft - left end of interval
c
c     aright - right end of interval
c
c     zeta(j) - j-th side condition point (boundary point). must
c               have  zeta(j) .le. zeta(j+1). all side condition
c               points must be mesh points in all meshes used,
c               see description of ipar(11) and fixpnt below.
c
c     ipar - an integer array dimensioned at least 11.
c            a list of the parameters in ipar and their meaning follows
c            some parameters are renamed in colnew; these new names are
c            given in parentheses.
c
c     ipar(1)     ( = nonlin )
c             = 0 if the problem is linear
c             = 1 if the problem is nonlinear
c
c     ipar(2) = no. of collocation points per subinterval  (= k )
c               where max m(i) .le.  k .le. 7 . if ipar(2)=0 then
c               colnew sets  k = max ( max m(i)+1, 5-max m(i) )
c
c     ipar(3) = no. of subintervals in the initial mesh  ( = n ).
c               if ipar(3) = 0 then colnew arbitrarily sets n = 5.
c
c     ipar(4) = no. of solution and derivative tolerances.  ( = ntol )
c               we require  0 .lt. ntol .le. mstar.
c
c     ipar(5) = dimension of fspace.     ( = ndimf )
c
c     ipar(6) = dimension of ispace.     ( = ndimi )
c
c     ipar(7) -  output control ( = iprint )
c              = -1 for full diagnostic printout
c              = 0 for selected printout
c              = 1 for no printout
c
c     ipar(8)     ( = iread )
c             = 0 causes colnew to generate a uniform initial mesh.
c             = 1 if the initial mesh is provided by the user.  it
c                 is defined in fspace as follows:  the mesh
c                 aleft=x(1).lt.x(2).lt. ... .lt.x(n).lt.x(n+1)=aright
c                 will occupy  fspace(1), ..., fspace(n+1). the
c                 user needs to supply only the interior mesh
c                 points  fspace(j) = x(j), j = 2, ..., n.
c             = 2 if the initial mesh is supplied by the user
c                 as with ipar(8)=1, and in addition no adaptive
c                 mesh selection is to be done.
c
c     ipar(9)     ( = iguess )
c             = 0 if no initial guess for the solution is
c                 provided.
c             = 1 if an initial guess is provided by the user
c                 in subroutine  guess.
c             = 2 if an initial mesh and approximate solution
c                 coefficients are provided by the user in  fspace.
c                 (the former and new mesh are the same).
c             = 3 if a former mesh and approximate solution
c                 coefficients are provided by the user in fspace,
c                 and the new mesh is to be taken twice as coarse;
c                 i.e.,every second point from the former mesh.
c             = 4 if in addition to a former initial mesh and
c                 approximate solution coefficients, a new mesh
c                 is provided in fspace as well.
c                 (see description of output for further details
c                 on iguess = 2, 3, and 4.)
c
c     ipar(10)= 0 if the problem is regular
c             = 1 if the first relax factor is =rstart, and the
c                 nonlinear iteration does not rely on past covergence
c                 (use for an extra sensitive nonlinear problem only).
c             = 2 if we are to return immediately upon  (a) two
c                 successive nonconvergences, or  (b) after obtaining
c                 error estimate for the first time.
c
c     ipar(11)= no. of fixed points in the mesh other than aleft
c               and aright. ( = nfxpnt , the dimension of fixpnt)
c               the code requires that all side condition points
c               other than aleft and aright (see description of
c               zeta ) be included as fixed points in fixpnt.
c
c     ltol  -  an array of dimension ipar(4). ltol(j) = l  specifies
c              that the j-th tolerance in  tol  controls the error
c              in the l-th component of z(u).   also require that
c              1.le.ltol(1).lt.ltol(2).lt. ... .lt.ltol(ntol).le.mstar
c
c     tol    - an array of dimension ipar(4). tol(j) is the
c              error tolerance on the ltol(j) -th component
c              of z(u). thus, the code attempts to satisfy
c              for j=1,...,ntol  on each subinterval
c              abs(z(v)-z(u))       .le. tol(j)*abs(z(u))       +tol(j)
c                            ltol(j)                     ltol(j)
c
c              if v(x) is the approximate solution vector.
c
c     fixpnt - an array of dimension ipar(11).   it contains
c              the points, other than aleft and aright, which
c              are to be included in every mesh.
c
c     ispace - an integer work array of dimension ipar(6).
c              its size provides a constraint on nmax,
c              the maximum number of subintervals. choose
c              ipar(6) according to the formula
c                      ipar(6)  .ge.  nmax*nsizei
c                where
c                      nsizei = 3 + kdm
c                with
c                      kdm = kd + mstar  ;  kd = k * ncomp ;
c                      nrec = no. of right end boundary conditions.
c
c
c     fspace - a real work array of dimension ipar(5).
c              its size provides a constraint on nmax.
c              choose ipar(5) according to the formula
c                      ipar(5)  .ge.  nmax*nsizef
c                where
c                      nsizef = 4 + 3 * mstar + (5+kd) * kdm +
c                              (2*mstar-nrec) * 2*mstar.
c
c
c     iflag - the mode of return from colnew.
c           = 1 for normal return
c           = 0 if the collocation matrix is singular.
c           =-1 if the expected no. of subintervals exceeds storage
c               specifications.
c           =-2 if the nonlinear iteration has not converged.
c           =-3 if there is an input data error.
c
c
c**********************************************************************
c
c     *************    user supplied subroutines   *************
c
c
c     the following subroutines must be declared external in the
c     main program which calls colnew.
c
c
c     fsub  - name of subroutine for evaluating f(x,z(u(x))) =
c                            t
c             (f ,...,f     )  at a point x in (aleft,aright).  it
c               1      ncomp
c             should have the heading
c
c                       subroutine fsub (x , z , f)
c
c             where f is the vector containing the value of fi(x,z(u))
c             in the i-th component and                            t
c                                       z(u(x))=(z(1),...,z(mstar))
c             is defined as above under  purpose .
c
c
c     dfsub - name of subroutine for evaluating the jacobian of
c             f(x,z(u)) at a point x.  it should have the heading
c
c                       subroutine dfsub (x , z , df)
c
c             where z(u(x)) is defined as for fsub and the (ncomp) by
c             (mstar) array df should be filled by the partial deriv-
c             atives of f, viz, for a particular call one calculates
c                                df(i,j) = dfi / dzj, i=1,...,ncomp
c                                                     j=1,...,mstar.
c
c
c     gsub  - name of subroutine for evaluating the i-th component of
c             g(x,z(u(x))) = g (zeta(i),z(u(zeta(i)))) at a point x =
c                             i
c             zeta(i) where 1.le.i.le.mstar. it should have the heading
c
c                       subroutine gsub (i , z , g)
c
c             where z(u) is as for fsub, and i and g=g  are as above.
c                                                     i
c             note that in contrast to f in  fsub , here
c             only one value per call is returned in g.
c
c
c     dgsub - name of subroutine for evaluating the i-th row of
c             the jacobian of g(x,u(x)).  it should have the heading
c
c                       subroutine dgsub (i , z , dg)
c
c             where z(u) is as for fsub, i as for gsub and the mstar-
c             vector dg should be filled with the partial derivatives
c             of g, viz, for a particular call one calculates
c                   dg(i,j) = dgi / dzj      j=1,...,mstar.
c
c
c     guess - name of subroutine to evaluate the initial
c             approximation for  z(u(x)) and for dmval(u(x))= vector
c             of the mj-th derivatives of u(x). it should have the
c             heading
c
c                       subroutine guess (x , z , dmval)
c
c             note that this subroutine is needed only if using
c             ipar(9) = 1, and then all  mstar  components of z
c             and  ncomp  components of  dmval  should be specified
c             for any x,  aleft .le. x .le. aright .
c
c
c**********************************************************************
c
c     ************   use of output from colnew   ************
c
c                 ***   solution evaluation   ***
c
c     on return from colnew, the arrays fspace and ispace
c     contain information specifying the approximate solution.
c     the user can produce the solution vector  z( u(x) )  at
c     any point x, aleft .le. x .le. aright, by the statement,
c
c           call appsln (x, z, fspace, ispace)
c
c     when saving the coefficients for later reference, only
c     ispace(1),...,ispace(7+ncomp)    and
c     fspace(1),...,fspace(ispace(7))    need to be saved as
c     these are the quantities used by appsln.
c
c
c                 ***   simple continuation   ***
c
c
c     a formerly obtained solution can easily be used as the
c     first approximation for the nonlinear iteration for a
c     new problem by setting   (iguess =) ipar(9) = 2, 3 or 4.
c
c     if the former solution has just been obtained then the
c     values needed to define the first approximation are
c     already in ispace and fspace.
c     alternatively, if the former solution was obtained in a
c     previous run and its coefficients were saved then those
c     coefficients must be put back into
c     ispace(1),..., ispace(7+ncomp)    and
c     fspace(1),..., fspace(ispace(7)).
c
c     for ipar(9) = 2 or 3 set ipar(3) = ispace(1) ( = the
c     size of the previous mesh ).
c
c     for ipar(9) = 4 the user specifies a new mesh of n subintervals
c     as follows.
c     the values in  fspace(1),...,fspace(ispace(7))  have to be
c     shifted by n+1 locations to  fspace(n+2),..,fspace(ispace(7)+n+1)
c     and the new mesh is then specified in fspace(1),..., fspace(n+1).
c     also set ipar(3) = n.
c
c
c**********************************************************************
c
c     ***************      package subroutines      ***************
c
c     the following description gives a brief overview of how the
c     procedure is broken down into the subroutines which make up
c     the package called  colnew . for further details the
c     user should refer to documentation in the various subroutines
c     and to the references cited above.
c
c     the subroutines fall into four groups:
c
c part 1 - the main storage allocation and program control subr
c
c     colnew - tests input values, does initialization and breaks up
c              the work areas, fspace and ispace, into the arrays
c              used by the program.
c     colsys - another name for colnew
c
c     contrl - is the actual driver of the package. this routine
c              contains the strategy for nonlinear equation solving.
c
c     skale  - provides scaling for the control
c              of convergence in the nonlinear iteration.
c
c
c part 2 - mesh selection and error estimation subroutines
c
c     consts - is called once by  colnew  to initialize constants
c              which are used for error estimation and mesh selection.
c
c     newmsh - generates meshes. it contains the test to decide
c              whether or not to redistribute a mesh.
c
c     errchk - produces error estimates and checks against the
c              tolerances at each subinterval
c
c
c part 3 - collocation system set-up subroutines
c
c     lsyslv - controls the set-up and solution of the linear
c              algebraic systems of collocation equations which
c              arise at each newton iteration.
c
c     gderiv - is used by lsyslv to set up the equation associated
c              with a side condition point.
c
c     vwblok - is used by lsyslv to set up the equation(s) associated
c              with a collocation point.
c
c     gblock - is used by lsyslv to construct a block of the global
c              collocation matrix or the corresponding right hand
c              side.
c
c
c part 4 - service subroutines
c
c     appsln - sets up a standard call to  approx .
c
c     approx - evaluates a piecewise polynomial solution.
c
c     rkbas  - evaluates the mesh independent runge-kutta basis
c
c     vmonde - solves a vandermonde system for given right hand
c              side
c
c     horder - evaluates the highest order derivatives of the
c              current collocation solution used for mesh refinement.
c
c
c part 5 - linear algebra  subroutines
c
c     to solve the global linear systems of collocation equations
c     constructed in part 3,  colnew  uses a column oriented version
c     of the package  solveblok originally due to de boor and weiss.
c
c     to solve the linear systems for static parameter condensation
c     in each block of the collocation equations, the linpack
c     routines  dgefa and  dgesl  are included. but these
c     may be replaced when solving problems on vector processors
c     or when solving large scale sparse jacobian problems.
c
c

	  use collocation
      use colnew
      use driver
      use lingebra
      use mesh
      use numethods
      use vectors
      implicit real*8 (a-h,o-z)

      

      type params
		real::r, delta, v, sigmap, sigmaa, theta, lambda,l, gamma,phi, ksi, etaa, etap,chia, chip, rho, alpha, mu
	end type

	type (params) :: paramsl,paramsh

      dimension m(1), zeta(1), ipar(1), ltol(1), tol(1), dummy(1),
     1          fixpnt(1), ispace(1), fspace(1)
c
      common /colout/ precis, iout, iprint
      common /colloc/ rho(7), coef(49)
      common /colord/ k, nc, mstar, kd, mmax, mt(20)
      common /colapr/ n, nold, nmax, nz, ndmz
      common /colmsh/ mshflg, mshnum, mshlmt, mshalt
      common /colsid/ tzeta(40), tleft, tright, izeta, idum
      common /colnln/ nonlin, iter, limit, icare, iguess
      common /colest/ ttl(40), wgtmsh(40), wgterr(40), tolin(40),
     1                root(40), jtol(40), lttol(40), ntol
c
      external fsub, dfsub, gsub, dgsub, guess
c
c     this subroutine can be called either colnew or colsys
c
      entry      colsys (ncomp, m, aleft, aright, zeta, ipar, ltol,
     1                   tol, fixpnt, ispace, fspace, iflag,
     2                   fsub, dfsub, gsub, dgsub, guess)
c
c*********************************************************************
c
c     the actual subroutine colnew serves as an interface with
c     the package of subroutines referred to collectively as
c     colnew. the subroutine serves to test some of the input
c     parameters, rename some of the parameters (to make under-
c     standing of the coding easier), to do some initialization,
c     and to break the work areas fspace and ispace up into the
c     arrays needed by the program.
c
c**********************************************************************
c
c...  specify machine dependent output unit  iout  and compute machine
c...  dependent constant  precis = 100 * machine unit roundoff
c
      if ( ipar(7) .le. 0 )  write(6,99)
  99  format(//,33h version *colnew* of colsys .    ,//)
c
      iout = 6
      precis = 1.d0
   10 precis = precis / 2.d0
      precp1 = precis + 1.d0
      if ( precp1 .gt. 1.d0 )                       go to 10
      precis = precis * 100.d0
c
c...  in case incorrect input data is detected, the program returns
c...  immediately with iflag=-3.
c
      iflag = -3
      if ( ncomp .lt. 1 .or. ncomp .gt. 20 )        return
      do 20 i=1,ncomp
         if ( m(i) .lt. 1 .or. m(i) .gt. 4 )        return
   20 continue
c
c...  rename some of the parameters and set default values.
c
      nonlin = ipar(1)
      k = ipar(2)
      n = ipar(3)
      if ( n .eq. 0 )  n = 5
      iread = ipar(8)
      iguess = ipar(9)
      if ( nonlin .eq. 0 .and. iguess .eq. 1 )  iguess = 0
      if ( iguess .ge. 2 .and. iread .eq. 0 )   iread = 1
      icare = ipar(10)
      ntol = ipar(4)
      ndimf = ipar(5)
      ndimi = ipar(6)
      nfxpnt = ipar(11)
      iprint = ipar(7)
      mstar = 0
      mmax = 0
      do  30 i = 1, ncomp
         mmax = max0 ( mmax, m(i) )
         mstar = mstar + m(i)
         mt(i) = m(i)
   30 continue
      if ( k .eq. 0 )   k = max0( mmax + 1 , 5 - mmax )
      do 40 i = 1, mstar
   40 tzeta(i) = zeta(i)
      do 50 i = 1, ntol
         lttol(i) = ltol(i)
   50 tolin(i) = tol(i)
      tleft = aleft
      tright = aright
      nc = ncomp
      kd = k * ncomp
c
c...  print the input data for checking.
c
      if ( iprint .gt. -1 )                         go to 80
      if ( nonlin .gt. 0 )                          go to 60
      write (iout,260) ncomp, (m(ip), ip=1,ncomp)
      go to 70
   60 write (iout,270) ncomp, (m(ip), ip=1,ncomp)
   70 write (iout,280) (zeta(ip), ip=1,mstar)
      if ( nfxpnt .gt. 0 )
     1   write (iout,340) nfxpnt, (fixpnt(ip), ip=1,nfxpnt)
      write (iout,290) k
      write (iout,300) (ltol(ip), ip=1,ntol)
      write (iout,310) (tol(ip), ip=1,ntol)
      if (iguess .ge. 2) write (iout,320)
      if (iread .eq. 2) write (iout,330)
   80 continue
c
c...  check for correctness of data
c
      if ( k .lt. 0 .or. k .gt. 7 )                 return
      if ( n .lt. 0 )                               return
      if ( iread .lt. 0 .or. iread .gt. 2 )         return
      if ( iguess .lt. 0 .or. iguess .gt. 4 )       return
      if ( icare .lt. 0 .or. icare .gt. 2 )         return
      if ( ntol .lt. 0 .or. ntol .gt. mstar )       return
      if ( nfxpnt .lt. 0 )                          return
      if ( iprint .lt. (-1) .or. iprint .gt. 1 )    return
      if ( mstar .lt. 0 .or. mstar .gt. 40 )        return
      ip = 1
      do 100 i = 1, mstar
      if ( dabs(zeta(i) - aleft) .lt. precis .or.
     1     dabs(zeta(i) - aright) .lt. precis )     go to 100
   90 if ( ip .gt. nfxpnt )                         return
        if ( zeta(i) - precis .lt. fixpnt(ip) )     go to 95
        ip = ip + 1
      go to 90
   95 if ( zeta(i) + precis .lt. fixpnt(ip) )       return
  100 continue
c
c...  set limits on iterations and initialize counters.
c...  limit = maximum number of newton iterations per mesh.
c...  see subroutine  newmsh  for the roles of  mshlmt , mshflg ,
c...  mshnum , and  mshalt .
c
      mshlmt = 3
      mshflg = 0
      mshnum = 1
      mshalt = 1
      limit = 40
c
c...  compute the maxium possible n for the given sizes of
c...  ispace  and  fspace.
c
      nrec = 0
      do 110 i = 1, mstar
           ib = mstar + 1 - i
           if ( zeta(ib) .ge. aright )  nrec = i
  110 continue
      nfixi = mstar
      nsizei = 3 + kd + mstar
      nfixf = nrec * (2*mstar) + 5 * mstar + 3
      nsizef = 4 + 3 * mstar + (kd+5) * (kd+mstar) +
     1(2*mstar-nrec) * 2*mstar
      nmaxf = (ndimf - nfixf) / nsizef
      nmaxi = (ndimi - nfixi) / nsizei
      if ( iprint .lt. 1 )  write(iout,350) nmaxf, nmaxi
      nmax = min0( nmaxf, nmaxi )
      if ( nmax .lt. n )                            return
      if ( nmax .lt. nfxpnt+1 )                     return
      if (nmax .lt. 2*nfxpnt+2 .and. iprint .lt. 1)  write(iout,360)
c
c...  generate pointers to break up  fspace  and  ispace .
c
      lxi = 1
      lg = lxi + nmax + 1
      lxiold = lg + 2*mstar * (nmax * (2*mstar-nrec) + nrec)
      lw     = lxiold + nmax + 1
      lv     = lw + kd**2 * nmax
      lz     = lv + mstar * kd * nmax
      ldmz   = lz + mstar * (nmax + 1)
      ldelz  = ldmz + kd * nmax
      ldeldz = ldelz + mstar * (nmax + 1)
      ldqz   = ldeldz + kd * nmax
      ldqdmz = ldqz + mstar * (nmax + 1)
      lrhs   = ldqdmz + kd * nmax
      lvalst = lrhs   + kd * nmax + mstar
      lslope = lvalst + 4 * mstar * nmax
      laccum = lslope + nmax
      lscl   = laccum + nmax + 1
      ldscl  = lscl + mstar * (nmax + 1)
      lpvtg = 1
      lpvtw = lpvtg + mstar * (nmax + 1)
      linteg = lpvtw + kd * nmax
c
c...  if  iguess .ge. 2, move  xiold, z, and  dmz  to their proper
c...  locations in  fspace.
c
      if ( iguess .lt. 2 )                          go to 160
      nold = n
      if (iguess .eq. 4)  nold = ispace(1)
      nz = mstar * (nold + 1)
      ndmz = kd * nold
      np1 = n + 1
      if ( iguess .eq. 4 )  np1 = np1 + nold + 1
      do 120 i=1,nz
  120 fspace( lz+i-1 )  =  fspace( np1+i )
      idmz = np1 + nz
      do 125 i=1,ndmz
  125 fspace( ldmz+i-1 )  =  fspace( idmz+i )
      np1 = nold + 1
      if ( iguess .eq. 4 )                          go to 140
      do 130 i=1,np1
  130 fspace( lxiold+i-1 )  =  fspace( lxi+i-1 )
      go to 160
  140 do 150 i=1,np1
  150 fspace( lxiold+i-1 )  =  fspace( n+1+i )
  160 continue
c
c...  initialize collocation points, constants, mesh.
c
      call consts ( k, rho, coef )
      call newmsh (3+iread, fspace(lxi), fspace(lxiold), dummy,
     1             dummy, dummy, dummy, dummy, nfxpnt, fixpnt)
c
c...  determine first approximation, if the problem is nonlinear.
c
      if (iguess .ge. 2)                            go to 230
      np1 = n + 1
      do 210 i = 1, np1
  210 fspace( i + lxiold - 1 ) = fspace( i + lxi - 1 )
      nold = n
      if ( nonlin .eq. 0  .or. iguess .eq. 1 )      go to 230
c
c...  system provides first approximation of the solution.
c...  choose z(j) = 0  for j=1,...,mstar.
c
      do 220 i=1, nz
  220 fspace( lz-1+i ) = 0.d0
      do 225 i=1, ndmz
  225 fspace( ldmz-1+i ) = 0.d0
  230 continue
      if (iguess .ge. 2)  iguess = 0
      call contrl (fspace(lxi),fspace(lxiold),fspace(lz),fspace(ldmz),
     1     fspace(lrhs),fspace(ldelz),fspace(ldeldz),fspace(ldqz),
     2     fspace(ldqdmz),fspace(lg),fspace(lw),fspace(lv),
     3     fspace(lvalst),fspace(lslope),fspace(lscl),fspace(ldscl),
     4     fspace(laccum),ispace(lpvtg),ispace(linteg),ispace(lpvtw),
     5     nfxpnt,fixpnt,iflag,fsub,dfsub,gsub,dgsub,guess )
c
c...  prepare output
c
      ispace(1) = n
      ispace(2) = k
      ispace(3) = ncomp
      ispace(4) = mstar
      ispace(5) = mmax
      ispace(6) = nz + ndmz + n + 2
      k2 = k * k
      ispace(7) = ispace(6) + k2 - 1
      do 240 i = 1, ncomp
  240 ispace(7+i) = m(i)
      do 250 i = 1, nz
  250 fspace( n+1+i ) = fspace( lz-1+i )
      idmz = n + 1 + nz
      do 255 i = 1, ndmz
  255 fspace( idmz+i ) = fspace( ldmz-1+i )
      ic = idmz + ndmz
      do 258 i = 1, k2
  258 fspace( ic+i ) = coef(i)
      return
c
  260 format(/// 37h the number of (linear) diff eqns is , i3/ 1x,
     1       16htheir orders are, 20i3)
  270 format(/// 40h the number of (nonlinear) diff eqns is , i3/ 1x,
     1       16htheir orders are, 20i3)
  280 format(27h side condition points zeta, 8f10.6, 4( / 27x, 8f10.6))
  290 format(37h number of colloc pts per interval is, i3)
  300 format(39h components of z requiring tolerances -,8(7x,i2,1x),
     1       4(/38x,8i10))
  310 format(33h corresponding error tolerances -,6x,8d10.2,
     1       4(/39x,8d10.2))
  320 format(44h initial mesh(es) and z,dmz provided by user)
  330 format(27h no adaptive mesh selection)
  340 format(10h there are ,i5,27h fixed points in the mesh - ,
     1       10(6f10.6/))
  350 format(44h the maximum number of subintervals is min (, i4,
     1       23h (allowed from fspace),,i4, 24h (allowed from ispace) ))
  360 format(/53h insufficient space to double mesh for error estimate)
      end
   