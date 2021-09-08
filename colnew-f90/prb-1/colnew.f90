module colnewy
	contains

!*==COLNEW.spg  processed by SPAG 6.72Dc at 03:32 on  8 Sep 2021
!  this package solves boundary value problems for
!  ordinary differential equations, as described below.
!
!  COLNEW is a modification of the package COLSYS by ascher,
!  christiansen and russell [1]. It incorporates a new basis
!  representation replacing b-splines, and improvements for
!  the linear and nonlinear algebraic equation solvers.
!  the package can be referenced as either COLNEW or COLSYS.
!
!----------------------------------------------------------------------
!                            p a r t  1
!        main storage allocation and program control subroutines
!----------------------------------------------------------------------
!
      SUBROUTINE COLNEW(Ncomp,M,Aleft,Aright,Zeta,Ipar,Ltol,Tol,Fixpnt, &
     &                  Ispace,Fspace,Iflag,FSUB,DFSUB,GSUB,DGSUB,GUESS)
 
!*********************************************************************72
!
!c COLNEW solves a multipoint boundary value problem using collocation.
!
!     written by
!                  u. ascher,
!                            department of computer science,
!                            university of british columbia,
!                            vancouver, b. c., canada   v6t 1w5
!                  g. bader,
!                            institut f. angewandte mathematik
!                            university of heidelberg
!                            im neuenheimer feld 294
!                            d-6900 heidelberg 1
!
!**********************************************************************
!
!     purpose
!
!     this package solves a multi-point boundary value
!     problem for a mixed order system of ode-s given by
!
!          (m(i))
!         u       =  f  ( x; z(u(x)) )      i = 1, ... ,ncomp
!          i          i
!
!                                          aleft .lt. x .lt. aright,
!
!
!         g  ( zeta(j); z(u(zeta(j))) ) = 0   j = 1, ... ,mstar
!          j
!                                    mstar = m(1)+m(2)+...+m(ncomp),
!
!
!         where                          t
!               u = (u , u , ... ,u     )  is the exact solution vector
!                     1   2        ncomp
!
!                (mi)
!               u     is the mi=m(i) th  derivative of u
!                i                                      i
!
!                                  (1)        (m1-1)       (mncomp-1)
!               z(u(x)) = ( u (x),u  (x),...,u    (x),...,u      (x) )
!                            1     1          1            ncomp
!
!                f (x,z(u))   is a (generally) nonlinear function of
!                 i
!                             z(u)=z(u(x)).
!
!                g (zeta(j);z(u))  is a (generally) nonlinear function
!                 j
!                               used to represent a boundary condition.
!
!         the boundary points satisfy
!               aleft .le. zeta(1) .le. .. .le. zeta(mstar) .le. aright
!
!         the orders mi of the differential equations satisfy
!                            1 .le. m(i) .le. 4.
!
!
!**********************************************************************
!
!     method
!
!        the method used to approximate the solution u is
!     collocation at gaussian points, requiring m(i)-1 continuous
!     derivatives in the i-th component, i = 1, ..., ncomp.
!     here, k is the number of collocation points (stages) per
!     subinterval and is chosen such that k .ge. max m(i).
!     a runge-kutta-monomial solution representation is utilized.
!
!     references
!
!     [1] u. ascher, j. christiansen and r.d. russell,
!         collocation software for boundary-value odes,
!         acm trans. math software 7 (1981), 209-222.
!         this paper contains EXAMPLES where use of the code
!         is demonstrated.
!
!     [2] g. bader and u. ascher,
!         a new basis implementation for a mixed order
!         boundary value ode solver,
!         siam j. scient. stat. comput. (1987).
!
!     [3] u. ascher, j. christiansen and r.d. russell,
!         a collocation solver for mixed order
!         systems of boundary value problems,
!         math. comp. 33 (1979), 659-679.
!
!     [4] u. ascher, j. christiansen and r.d. russell,
!         colsys - a collocation code for boundary
!         value problems,
!         lecture notes comp.sc. 76, springer verlag,
!         b. childs et. al. (eds.) (1979), 164-185.
!
!     [5] c. deboor and r. weiss,
!         solveblok: a package for solving almost block diagonal
!         linear systems,
!         acm trans. math. software 6 (1980), 80-87.
!
!**********************************************************************
!
!     ***************     input to colnew     ***************
!
!     variables
!
!     ncomp - no. of differential equations   (ncomp .le. 20)
!
!     m(j) - order of the j-th differential equation
!            ( mstar = m(1) + ... + m(ncomp) .le. 40 )
!
!     aleft - left end of interval
!
!     aright - right end of interval
!
!     zeta(j) - j-th side condition point (boundary point). must
!               have  zeta(j) .le. zeta(j+1). all side condition
!               points must be mesh points in all meshes used,
!               see description of ipar(11) and fixpnt below.
!
!     ipar - an integer array dimensioned at least 11.
!            a list of the parameters in ipar and their meaning follows
!            some parameters are renamed in colnew; these new names are
!            given in parentheses.
!
!     ipar(1)     ( = nonlin )
!             = 0 if the problem is linear
!             = 1 if the problem is nonlinear
!
!     ipar(2) = no. of collocation points per subinterval  (= k )
!               where max m(i) .le.  k .le. 7 . if ipar(2)=0 then
!               colnew sets  k = max ( max m(i)+1, 5-max m(i) )
!
!     ipar(3) = no. of subintervals in the initial mesh  ( = n ).
!               if ipar(3) = 0 then colnew arbitrarily sets n = 5.
!
!     ipar(4) = no. of solution and derivative tolerances.  ( = ntol )
!               we require  0 .lt. ntol .le. mstar.
!
!     ipar(5) = dimension of fspace.     ( = ndimf )
!
!     ipar(6) = dimension of ispace.     ( = ndimi )
!
!     ipar(7) -  output control ( = iprint )
!              = -1 for full diagnostic printout
!              = 0 for selected printout
!              = 1 for no printout
!
!     ipar(8)     ( = iread )
!             = 0 causes colnew to generate a uniform initial mesh.
!             = 1 if the initial mesh is provided by the user.  it
!                 is defined in fspace as follows:  the mesh
!                 aleft=x(1).lt.x(2).lt. ... .lt.x(n).lt.x(n+1)=aright
!                 will occupy  fspace(1), ..., fspace(n+1). the
!                 user needs to supply only the interior mesh
!                 points  fspace(j) = x(j), j = 2, ..., n.
!             = 2 if the initial mesh is supplied by the user
!                 as with ipar(8)=1, and in addition no adaptive
!                 mesh selection is to be done.
!
!     ipar(9)     ( = iguess )
!             = 0 if no initial guess for the solution is
!                 provided.
!             = 1 if an initial guess is provided by the user
!                 in subroutine  guess.
!             = 2 if an initial mesh and approximate solution
!                 coefficients are provided by the user in  fspace.
!                 (the former and new mesh are the same).
!             = 3 if a former mesh and approximate solution
!                 coefficients are provided by the user in fspace,
!                 and the new mesh is to be taken twice as coarse;
!                 i.e.,every second point from the former mesh.
!             = 4 if in addition to a former initial mesh and
!                 approximate solution coefficients, a new mesh
!                 is provided in fspace as well.
!                 (see description of output for further details
!                 on iguess = 2, 3, and 4.)
!
!     ipar(10)= 0 if the problem is regular
!             = 1 if the first relax factor is =rstart, and the
!                 nonlinear iteration does not rely on past covergence
!                 (use for an extra sensitive nonlinear problem only).
!             = 2 if we are to return immediately upon  (a) two
!                 successive nonconvergences, or  (b) after obtaining
!                 error estimate for the first time.
!
!     ipar(11)= no. of fixed points in the mesh other than aleft
!               and aright. ( = nfxpnt , the dimension of fixpnt)
!               the code requires that all side condition points
!               other than aleft and aright (see description of
!               zeta ) be included as fixed points in fixpnt.
!
!     ltol  -  an array of dimension ipar(4). ltol(j) = l  specifies
!              that the j-th tolerance in  tol  controls the error
!              in the l-th component of z(u).   also require that
!              1.le.ltol(1).lt.ltol(2).lt. ... .lt.ltol(ntol).le.mstar
!
!     tol    - an array of dimension ipar(4). tol(j) is the
!              error tolerance on the ltol(j) -th component
!              of z(u). thus, the code attempts to satisfy
!              for j=1,...,ntol  on each subinterval
!              abs(z(v)-z(u))       .le. tol(j)*abs(z(u))       +tol(j)
!                            ltol(j)                     ltol(j)
!
!              if v(x) is the approximate solution vector.
!
!     fixpnt - an array of dimension ipar(11).   it contains
!              the points, other than aleft and aright, which
!              are to be included in every mesh.
!
!     ispace - an integer work array of dimension ipar(6).
!              its size provides a constraint on nmax,
!              the maximum number of subintervals. choose
!              ipar(6) according to the formula
!                      ipar(6)  .ge.  nmax*nsizei
!                where
!                      nsizei = 3 + kdm
!                with
!                      kdm = kd + mstar  ;  kd = k * ncomp ;
!                      nrec = no. of right end boundary conditions.
!
!
!     fspace - a real work array of dimension ipar(5).
!              its size provides a constraint on nmax.
!              choose ipar(5) according to the formula
!                      ipar(5)  .ge.  nmax*nsizef
!                where
!                      nsizef = 4 + 3 * mstar + (5+kd) * kdm +
!                              (2*mstar-nrec) * 2*mstar.
!
!
!     iflag - the mode of return from colnew.
!           = 1 for normal return
!           = 0 if the collocation matrix is singular.
!           =-1 if the expected no. of subintervals exceeds storage
!               specifications.
!           =-2 if the nonlinear iteration has not converged.
!           =-3 if there is an input data error.
!
!
!**********************************************************************
!
!     *************    user supplied subroutines   *************
!
!
!     the following subroutines must be declared external in the
!     main program which calls colnew.
!
!
!     fsub  - name of subroutine for evaluating f(x,z(u(x))) =
!                            t
!             (f ,...,f     )  at a point x in (aleft,aright).  it
!               1      ncomp
!             should have the heading
!
!                       subroutine fsub (x , z , f)
!
!             where f is the vector containing the value of fi(x,z(u))
!             in the i-th component and                            t
!                                       z(u(x))=(z(1),...,z(mstar))
!             is defined as above under  purpose .
!
!
!     dfsub - name of subroutine for evaluating the jacobian of
!             f(x,z(u)) at a point x.  it should have the heading
!
!                       subroutine dfsub (x , z , df)
!
!             where z(u(x)) is defined as for fsub and the (ncomp) by
!             (mstar) array df should be filled by the partial deriv-
!             atives of f, viz, for a particular call one calculates
!                                df(i,j) = dfi / dzj, i=1,...,ncomp
!                                                     j=1,...,mstar.
!
!
!     gsub  - name of subroutine for evaluating the i-th component of
!             g(x,z(u(x))) = g (zeta(i),z(u(zeta(i)))) at a point x =
!                             i
!             zeta(i) where 1.le.i.le.mstar. it should have the heading
!
!                       subroutine gsub (i , z , g)
!
!             where z(u) is as for fsub, and i and g=g  are as above.
!                                                     i
!             note that in contrast to f in  fsub , here
!             only one value per call is returned in g.
!
!
!     dgsub - name of subroutine for evaluating the i-th row of
!             the jacobian of g(x,u(x)).  it should have the heading
!
!                       subroutine dgsub (i , z , dg)
!
!             where z(u) is as for fsub, i as for gsub and the mstar-
!             vector dg should be filled with the partial derivatives
!             of g, viz, for a particular call one calculates
!                   dg(i,j) = dgi / dzj      j=1,...,mstar.
!
!
!     guess - name of subroutine to evaluate the initial
!             approximation for  z(u(x)) and for dmval(u(x))= vector
!             of the mj-th derivatives of u(x). it should have the
!             heading
!
!                       subroutine guess (x , z , dmval)
!
!             note that this subroutine is needed only if using
!             ipar(9) = 1, and then all  mstar  components of z
!             and  ncomp  components of  dmval  should be specified
!             for any x,  aleft .le. x .le. aright .
!
!
!**********************************************************************
!
!     ************   use of output from colnew   ************
!
!                 ***   solution evaluation   ***
!
!     on return from colnew, the arrays fspace and ispace
!     contain information specifying the approximate solution.
!     the user can produce the solution vector  z( u(x) )  at
!     any point x, aleft .le. x .le. aright, by the statement,
!
!           call appsln (x, z, fspace, ispace)
!
!     when saving the coefficients for later reference, only
!     ispace(1),...,ispace(7+ncomp)    and
!     fspace(1),...,fspace(ispace(7))    need to be saved as
!     these are the quantities used by appsln.
!
!
!                 ***   simple continuation   ***
!
!
!     a formerly obtained solution can easily be used as the
!     first approximation for the nonlinear iteration for a
!     new problem by setting   (iguess =) ipar(9) = 2, 3 or 4.
!
!     if the former solution has just been obtained then the
!     values needed to define the first approximation are
!     already in ispace and fspace.
!     alternatively, if the former solution was obtained in a
!     previous run and its coefficients were saved then those
!     coefficients must be put back into
!     ispace(1),..., ispace(7+ncomp)    and
!     fspace(1),..., fspace(ispace(7)).
!
!     for ipar(9) = 2 or 3 set ipar(3) = ispace(1) ( = the
!     size of the previous mesh ).
!
!     for ipar(9) = 4 the user specifies a new mesh of n subintervals
!     as follows.
!     the values in  fspace(1),...,fspace(ispace(7))  have to be
!     shifted by n+1 locations to  fspace(n+2),..,fspace(ispace(7)+n+1)
!     and the new mesh is then specified in fspace(1),..., fspace(n+1).
!     also set ipar(3) = n.
!
!
!**********************************************************************
!
!     ***************      package subroutines      ***************
!
!     the following description gives a brief overview of how the
!     procedure is broken down into the subroutines which make up
!     the package called  colnew . for further details the
!     user should refer to documentation in the various subroutines
!     and to the references cited above.
!
!     the subroutines fall into four groups:
!
! part 1 - the main storage allocation and program control subr
!
!     colnew - tests input values, does initialization and breaks up
!              the work areas, fspace and ispace, into the arrays
!              used by the program.
!     colsys - another name for colnew
!
!     contrl - is the actual driver of the package. this routine
!              contains the strategy for nonlinear equation solving.
!
!     skale  - provides scaling for the control
!              of convergence in the nonlinear iteration.
!
!
! part 2 - mesh selection and error estimation subroutines
!
!     consts - is called once by  colnew  to initialize constants
!              which are used for error estimation and mesh selection.
!
!     newmsh - generates meshes. it contains the test to decide
!              whether or not to redistribute a mesh.
!
!     errchk - produces error estimates and checks against the
!              tolerances at each subinterval
!
!
! part 3 - collocation system set-up subroutines
!
!     lsyslv - controls the set-up and solution of the linear
!              algebraic systems of collocation equations which
!              arise at each newton iteration.
!
!     gderiv - is used by lsyslv to set up the equation associated
!              with a side condition point.
!
!     vwblok - is used by lsyslv to set up the equation(s) associated
!              with a collocation point.
!
!     gblock - is used by lsyslv to construct a block of the global
!              collocation matrix or the corresponding right hand
!              side.
!
!
! part 4 - service subroutines
!
!     appsln - sets up a standard call to  approx .
!
!     approx - evaluates a piecewise polynomial solution.
!
!     rkbas  - evaluates the mesh independent runge-kutta basis
!
!     vmonde - solves a vandermonde system for given right hand
!              side
!
!     horder - evaluates the highest order derivatives of the
!              current collocation solution used for mesh refinement.
!
!
! part 5 - linear algebra  subroutines
!
!     to solve the global linear systems of collocation equations
!     constructed in part 3,  colnew  uses a column oriented version
!     of the package  solveblok originally due to de boor and weiss.
!
!     to solve the linear systems for static parameter condensation
!     in each block of the collocation equations, the linpack
!     routines  dgefa and  dgesl  are included. but these
!     may be replaced when solving problems on vector processors
!     or when solving large scale sparse jacobian problems.
!
!
      IMPLICIT NONE
!*--COLNEW463
!*** Start of declarations inserted by SPAG
      REAL*8 Aleft , Aright , COEf , DFSUB , DGSUB , dummy , Fixpnt ,   &
     &       Fspace , FSUB , GSUB , GUESS , PREcis , precp1 , RHO ,     &
     &       ROOt , TLEft , Tol , TOLin , TRIght , TTL
      REAL*8 TZEta , WGTerr , WGTmsh , Zeta
      INTEGER i , ib , ic , ICAre , idmz , IDUm , Iflag , IGUess ,      &
     &        IOUt , ip , Ipar , IPRint , iread , Ispace , ITEr ,       &
     &        IZEta , JTOl , K , k2 , KD
      INTEGER laccum , ldeldz , ldelz , ldmz , ldqdmz , ldqz , ldscl ,  &
     &        lg , LIMit , linteg , lpvtg , lpvtw , lrhs , lscl ,       &
     &        lslope , Ltol , LTTol , lv , lvalst , lw
      INTEGER lxi , lxiold , lz , M , MMAx , MSHalt , MSHflg , MSHlmt , &
     &        MSHnum , MSTar , MT , N , NC , Ncomp , ndimf , ndimi ,    &
     &        NDMz , nfixf , nfixi , nfxpnt
      INTEGER NMAx , nmaxf , nmaxi , NOLd , NONlin , np1 , nrec ,       &
     &        nsizef , nsizei , NTOl , NZ
!*** End of declarations inserted by SPAG
 
!      STRUCTURE /PARAMS/
!         REAL :: R , DELTA , V , SIGMAP , SIGMAA , THETA , LAMBDA , L , &
!     &           GAMMA , PHI
!         REAL :: KSI , ETAA , ETAP , CHIA , CHIP , RHO , ALPHA , MU
!      END STRUCTURE
 
!      RECORD /PARAMS/ ::Paramsl , Paramsh
 
      DIMENSION M(1) , Zeta(1) , Ipar(1) , Ltol(1) , Tol(1) , dummy(1) ,&
     &          Fixpnt(1) , Ispace(1) , Fspace(1)
!
      COMMON /COLOUT/ PREcis , IOUt , IPRint
      COMMON /COLLOC/ RHO(7) , COEf(49)
      COMMON /COLORD/ K , NC , MSTar , KD , MMAx , MT(20)
      COMMON /COLAPR/ N , NOLd , NMAx , NZ , NDMz
      COMMON /COLMSH/ MSHflg , MSHnum , MSHlmt , MSHalt
      COMMON /COLSID/ TZEta(40) , TLEft , TRIght , IZEta , IDUm
      COMMON /COLNLN/ NONlin , ITEr , LIMit , ICAre , IGUess
      COMMON /COLEST/ TTL(40) , WGTmsh(40) , WGTerr(40) , TOLin(40) ,   &
     &                ROOt(40) , JTOl(40) , LTTol(40) , NTOl
!
      EXTERNAL FSUB , DFSUB , GSUB , DGSUB , GUESS
!
!     this subroutine can be called either COLNEW or COLSYS
!
      ENTRY COLSYS(Ncomp,M,Aleft,Aright,Zeta,Ipar,Ltol,Tol,Fixpnt,      &
     &             Ispace,Fspace,Iflag,FSUB,DFSUB,GSUB,DGSUB,GUESS)
!
!*********************************************************************
!
!     the actual subroutine colnew serves as an interface with
!     the package of subroutines referred to collectively as
!     colnew. the subroutine serves to test some of the input
!     parameters, rename some of the parameters (to make under-
!     standing of the coding easier), to do some initialization,
!     and to break the work areas fspace and ispace up into the
!     arrays needed by the program.
!
!**********************************************************************
!
!...  specify machine dependent output unit  iout  and compute machine
!...  dependent constant  precis = 100 * machine unit roundoff
!
      IF ( Ipar(7).LE.0 ) WRITE (6,99001)
99001 FORMAT (//,' VERSION *COLNEW* OF COLSYS .    ',//)
!
      IOUt = 6
      PREcis = 1.D0
 100  PREcis = PREcis/2.D0
      precp1 = PREcis + 1.D0
      IF ( precp1.GT.1.D0 ) GOTO 100
      PREcis = PREcis*100.D0
!
!...  in case incorrect input data is detected, the program returns
!...  immediately with iflag=-3.
!
      Iflag = -3
      IF ( Ncomp.LT.1 .OR. Ncomp.GT.20 ) RETURN
      DO i = 1 , Ncomp
         IF ( M(i).LT.1 .OR. M(i).GT.4 ) RETURN
      ENDDO
!
!...  rename some of the parameters and set default values.
!
      NONlin = Ipar(1)
      K = Ipar(2)
      N = Ipar(3)
      IF ( N.EQ.0 ) N = 5
      iread = Ipar(8)
      IGUess = Ipar(9)
      IF ( NONlin.EQ.0 .AND. IGUess.EQ.1 ) IGUess = 0
      IF ( IGUess.GE.2 .AND. iread.EQ.0 ) iread = 1
      ICAre = Ipar(10)
      NTOl = Ipar(4)
      ndimf = Ipar(5)
      ndimi = Ipar(6)
      nfxpnt = Ipar(11)
      IPRint = Ipar(7)
      MSTar = 0
      MMAx = 0
      DO i = 1 , Ncomp
         MMAx = MAX0(MMAx,M(i))
         MSTar = MSTar + M(i)
         MT(i) = M(i)
      ENDDO
      IF ( K.EQ.0 ) K = MAX0(MMAx+1,5-MMAx)
      DO i = 1 , MSTar
         TZEta(i) = Zeta(i)
      ENDDO
      DO i = 1 , NTOl
         LTTol(i) = Ltol(i)
         TOLin(i) = Tol(i)
      ENDDO
      TLEft = Aleft
      TRIght = Aright
      NC = Ncomp
      KD = K*Ncomp
!
!...  print the input data for checking.
!
      IF ( IPRint.LE.-1 ) THEN
         IF ( NONlin.GT.0 ) THEN
            WRITE (IOUt,99002) Ncomp , (M(ip),ip=1,Ncomp)
99002       FORMAT (///' THE NUMBER OF (NONLINEAR) DIFF EQNS IS ',I3/1X,&
     &              'THEIR ORDERS ARE',20I3)
         ELSE
            WRITE (IOUt,99003) Ncomp , (M(ip),ip=1,Ncomp)
!
99003       FORMAT (///' THE NUMBER OF (LINEAR) DIFF EQNS IS ',I3/1X,   &
     &              'THEIR ORDERS ARE',20I3)
         ENDIF
         WRITE (IOUt,99004) (Zeta(ip),ip=1,MSTar)
99004    FORMAT (' SIDE CONDITION POINTS ZETA',8F10.6,4(/27X,8F10.6))
         IF ( nfxpnt.GT.0 ) WRITE (IOUt,99005) nfxpnt ,                 &
     &                             (Fixpnt(ip),ip=1,nfxpnt)
99005    FORMAT (' THERE ARE',I5,' FIXED POINTS IN THE MESH -',         &
     &           10(6F10.6/))
         WRITE (IOUt,99006) K
99006    FORMAT (' NUMBER OF COLLOC PTS PER INTERVAL IS',I3)
         WRITE (IOUt,99007) (Ltol(ip),ip=1,NTOl)
99007    FORMAT (' COMPONENTS OF Z REQUIRING TOLERANCES -',8(7X,I2,1X), &
     &           4(/38X,8I10))
         WRITE (IOUt,99008) (Tol(ip),ip=1,NTOl)
99008    FORMAT (' CORRESPONDING ERROR TOLERANCES -',6X,8D10.2,         &
     &           4(/39X,8D10.2))
         IF ( IGUess.GE.2 ) WRITE (IOUt,99009)
99009    FORMAT (' INITIAL MESH(ES) AND Z,DMZ PROVIDED BY USER')
         IF ( iread.EQ.2 ) WRITE (IOUt,99010)
99010    FORMAT (' NO ADAPTIVE MESH SELECTION')
      ENDIF
!
!...  check for correctness of data
!
      IF ( K.LT.0 .OR. K.GT.7 ) RETURN
      IF ( N.LT.0 ) RETURN
      IF ( iread.LT.0 .OR. iread.GT.2 ) RETURN
      IF ( IGUess.LT.0 .OR. IGUess.GT.4 ) RETURN
      IF ( ICAre.LT.0 .OR. ICAre.GT.2 ) RETURN
      IF ( NTOl.LT.0 .OR. NTOl.GT.MSTar ) RETURN
      IF ( nfxpnt.LT.0 ) RETURN
      IF ( IPRint.LT.(-1) .OR. IPRint.GT.1 ) RETURN
      IF ( MSTar.LT.0 .OR. MSTar.GT.40 ) RETURN
      ip = 1
      DO i = 1 , MSTar
         IF ( DABS(Zeta(i)-Aleft).GE.PREcis .AND. DABS(Zeta(i)-Aright)  &
     &        .GE.PREcis ) THEN
 120        IF ( ip.GT.nfxpnt ) RETURN
            IF ( Zeta(i)-PREcis.GE.Fixpnt(ip) ) THEN
               ip = ip + 1
               GOTO 120
            ELSEIF ( Zeta(i)+PREcis.LT.Fixpnt(ip) ) THEN
               RETURN
            ENDIF
         ENDIF
      ENDDO
!
!...  set limits on iterations and initialize counters.
!...  limit = maximum number of newton iterations per mesh.
!...  see subroutine  newmsh  for the roles of  mshlmt , mshflg ,
!...  mshnum , and  mshalt .
!
      MSHlmt = 3
      MSHflg = 0
      MSHnum = 1
      MSHalt = 1
      LIMit = 40
!
!...  compute the maxium possible n for the given sizes of
!...  ispace  and  fspace.
!
      nrec = 0
      DO i = 1 , MSTar
         ib = MSTar + 1 - i
         IF ( Zeta(ib).GE.Aright ) nrec = i
      ENDDO
      nfixi = MSTar
      nsizei = 3 + KD + MSTar
      nfixf = nrec*(2*MSTar) + 5*MSTar + 3
      nsizef = 4 + 3*MSTar + (KD+5)*(KD+MSTar) + (2*MSTar-nrec)*2*MSTar
      nmaxf = (ndimf-nfixf)/nsizef
      nmaxi = (ndimi-nfixi)/nsizei
      IF ( IPRint.LT.1 ) WRITE (IOUt,99011) nmaxf , nmaxi
99011 FORMAT (' THE MAXIMUM NUMBER OF SUBINTERVALS IS MIN (',I4,        &
     &        ' (ALLOWED FROM FSPACE),',I4,' (ALLOWED FROM ISPACE) )')
      NMAx = MIN0(nmaxf,nmaxi)
      IF ( NMAx.LT.N ) RETURN
      IF ( NMAx.LT.nfxpnt+1 ) RETURN
      IF ( NMAx.LT.2*nfxpnt+2 .AND. IPRint.LT.1 ) WRITE (IOUt,99012)
99012 FORMAT (/' INSUFFICIENT SPACE TO DOUBLE MESH FOR ERROR ESTIMATE')
!
!...  generate pointers to break up  fspace  and  ispace .
!
      lxi = 1
      lg = lxi + NMAx + 1
      lxiold = lg + 2*MSTar*(NMAx*(2*MSTar-nrec)+nrec)
      lw = lxiold + NMAx + 1
      lv = lw + KD**2*NMAx
      lz = lv + MSTar*KD*NMAx
      ldmz = lz + MSTar*(NMAx+1)
      ldelz = ldmz + KD*NMAx
      ldeldz = ldelz + MSTar*(NMAx+1)
      ldqz = ldeldz + KD*NMAx
      ldqdmz = ldqz + MSTar*(NMAx+1)
      lrhs = ldqdmz + KD*NMAx
      lvalst = lrhs + KD*NMAx + MSTar
      lslope = lvalst + 4*MSTar*NMAx
      laccum = lslope + NMAx
      lscl = laccum + NMAx + 1
      ldscl = lscl + MSTar*(NMAx+1)
      lpvtg = 1
      lpvtw = lpvtg + MSTar*(NMAx+1)
      linteg = lpvtw + KD*NMAx
!
!...  if  iguess .ge. 2, move  xiold, z, and  dmz  to their proper
!...  locations in  fspace.
!
      IF ( IGUess.GE.2 ) THEN
         NOLd = N
         IF ( IGUess.EQ.4 ) NOLd = Ispace(1)
         NZ = MSTar*(NOLd+1)
         NDMz = KD*NOLd
         np1 = N + 1
         IF ( IGUess.EQ.4 ) np1 = np1 + NOLd + 1
         DO i = 1 , NZ
            Fspace(lz+i-1) = Fspace(np1+i)
         ENDDO
         idmz = np1 + NZ
         DO i = 1 , NDMz
            Fspace(ldmz+i-1) = Fspace(idmz+i)
         ENDDO
         np1 = NOLd + 1
         IF ( IGUess.EQ.4 ) THEN
            DO i = 1 , np1
               Fspace(lxiold+i-1) = Fspace(N+1+i)
            ENDDO
         ELSE
            DO i = 1 , np1
               Fspace(lxiold+i-1) = Fspace(lxi+i-1)
            ENDDO
         ENDIF
      ENDIF
!
!...  initialize collocation points, constants, mesh.
!
      CALL CONSTS(K,RHO,COEf)
      CALL NEWMSH(3+iread,Fspace(lxi),Fspace(lxiold),dummy,dummy,dummy, &
     &            dummy,dummy,nfxpnt,Fixpnt)
!
!...  determine first approximation, if the problem is nonlinear.
!
      IF ( IGUess.LT.2 ) THEN
         np1 = N + 1
         DO i = 1 , np1
            Fspace(i+lxiold-1) = Fspace(i+lxi-1)
         ENDDO
         NOLd = N
         IF ( NONlin.NE.0 .AND. IGUess.NE.1 ) THEN
!
!...  system provides first approximation of the solution.
!...  choose z(j) = 0  for j=1,...,mstar.
!
            DO i = 1 , NZ
               Fspace(lz-1+i) = 0.D0
            ENDDO
            DO i = 1 , NDMz
               Fspace(ldmz-1+i) = 0.D0
            ENDDO
         ENDIF
      ENDIF
      IF ( IGUess.GE.2 ) IGUess = 0
      CALL CONTRL(Fspace(lxi),Fspace(lxiold),Fspace(lz),Fspace(ldmz),   &
     &            Fspace(lrhs),Fspace(ldelz),Fspace(ldeldz),Fspace(ldqz)&
     &            ,Fspace(ldqdmz),Fspace(lg),Fspace(lw),Fspace(lv),     &
     &            Fspace(lvalst),Fspace(lslope),Fspace(lscl),           &
     &            Fspace(ldscl),Fspace(laccum),Ispace(lpvtg),           &
     &            Ispace(linteg),Ispace(lpvtw),nfxpnt,Fixpnt,Iflag,FSUB,&
     &            DFSUB,GSUB,DGSUB,GUESS)
!
!...  prepare output
!
      Ispace(1) = N
      Ispace(2) = K
      Ispace(3) = Ncomp
      Ispace(4) = MSTar
      Ispace(5) = MMAx
      Ispace(6) = NZ + NDMz + N + 2
      k2 = K*K
      Ispace(7) = Ispace(6) + k2 - 1
      DO i = 1 , Ncomp
         Ispace(7+i) = M(i)
      ENDDO
      DO i = 1 , NZ
         Fspace(N+1+i) = Fspace(lz-1+i)
      ENDDO
      idmz = N + 1 + NZ
      DO i = 1 , NDMz
         Fspace(idmz+i) = Fspace(ldmz-1+i)
      ENDDO
      ic = idmz + NDMz
      DO i = 1 , k2
         Fspace(ic+i) = COEf(i)
      ENDDO
      RETURN
      END
      
end module colnewy