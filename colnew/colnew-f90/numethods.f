c
c----------------------------------------------------------------------
c                             p a r t  4
c               polynomial and service routines
c----------------------------------------------------------------------
c
      subroutine appsln (x, z, fspace, ispace)
c
c*****************************************************************
c
cc appsln sets up a call to approximate the solution.
c
c     purpose
c
c           set up a standard call to  approx  to evaluate the
c           approximate solution  z = z( u(x) )  at a point x
c           (it has been computed by a call to  colnew ).
c           the parameters needed for  approx  are retrieved
c           from the work arrays  ispace  and  fspace .
c
c*****************************************************************
c
      implicit real*8 (a-h,o-z)
      dimension z(1), fspace(1), ispace(1), a(28), dummy(1)
      is6 = ispace(6)
      is5 = ispace(1) + 2
      is4 = is5 + ispace(4) * (ispace(1) + 1)
      i = 1
      call approx (i, x, z, a, fspace(is6), fspace(1), ispace(1),
     1             fspace(is5), fspace(is4), ispace(2), ispace(3),
     2             ispace(5), ispace(8), ispace(4), 2, dummy, 0)
      return
      end
      subroutine approx (i, x, zval, a, coef, xi, n, z, dmz, k,
     1                   ncomp, mmax, m, mstar, mode, dmval, modm )
c
c**********************************************************************
c
cc approx evaluates the computes solution at a point.
c
c   purpose
c                                    (1)       (m1-1)     (mncomp-1)
c           evaluate z(u(x))=(u (x),u (x),...,u  (x),...,u  (x)      )
c                              1     1         1          mncomp
c           at one point x.
c
c   variables
c     a      - array of mesh independent rk-basis coefficients
c     basm   - array of mesh dependent monomial coefficients
c     xi     - the current mesh (having n subintervals)
c     z      - the current solution vector
c     dmz    - the array of mj-th derivatives of the current solution
c     mode   - determines the amount of initialization needed
c            = 4  forms z(u(x)) using z, dmz and ha
c            = 3  as in =4, but computes local rk-basis
c            = 2  as in =3, but determines i such that
c                       xi(i) .le. x .lt. xi(i+1) (unless x=xi(n+1))
c            = 1  retrieve  z=z(u(x(i)))  directly
c
c**********************************************************************
c
      implicit real*8 (a-h,o-z)
      dimension zval(1), dmval(1), xi(1), m(1), a(7,1), dm(7)
      dimension z(1), dmz(1), bm(4), coef(1)
c
      common /colout/ precis, iout, iprint
c
      go to (10, 30, 80, 90), mode
c
c...  mode = 1 , retrieve  z( u(x) )  directly for x = xi(i).
c
   10 x  = xi(i)
      iz = (i-1) * mstar
      do 20 j = 1, mstar
        iz = iz + 1
        zval(j) = z(iz)
   20 continue
      return
c
c...  mode = 2 ,  locate i so  xi(i) .le. x .lt. xi(i+1)
c
   30 continue
      if ( x .ge. xi(1)-precis .and. x .le. xi(n+1)+precis )
     1                                              go to 40
      if (iprint .lt. 1)  write(iout,900) x, xi(1), xi(n+1)
      if ( x .lt. xi(1) )  x = xi(1)
      if ( x .gt. xi(n+1) )  x = xi(n+1)
   40 if ( i .gt. n .or. i .lt. 1 )  i = (n+1) / 2
      ileft = i
      if ( x .lt. xi(ileft) )                       go to 60
      do 50 l = ileft, n
           i = l
           if ( x .lt. xi(l+1) )                    go to 80
   50 continue
      go to 80
   60 iright = ileft - 1
      do 70 l = 1, iright
           i = iright + 1 - l
           if ( x .ge. xi(i) )                      go to 80
   70 continue
c
c...  mode = 2 or 3 , compute mesh independent rk-basis.
c
   80 continue
      s = (x - xi(i)) / (xi(i+1) - xi(i))
      call rkbas ( s, coef, k, mmax, a, dm, modm )
c
c...  mode = 2, 3, or 4 , compute mesh dependent rk-basis.
c
   90 continue
      bm(1) = x - xi(i)
      do 95 l = 2, mmax
         bm(l) = bm(1) / dfloat(l)
   95 continue
c
c...  evaluate  z( u(x) ).
c
  100 ir = 1
      iz = (i-1) * mstar + 1
      idmz = (i-1) * k * ncomp
      do 140 jcomp = 1, ncomp
          mj = m(jcomp)
          ir = ir + mj
          iz = iz + mj
          do 130 l = 1, mj
             ind = idmz + jcomp
             zsum = 0.d0
             do 110 j = 1, k
               zsum = zsum  +  a(j,l) * dmz(ind)
  110        ind = ind + ncomp
             do 120 ll = 1, l
               lb = l + 1 - ll
  120          zsum = zsum * bm(lb)  +  z(iz-ll)
  130     zval(ir-l) = zsum
  140 continue
      if ( modm .eq. 0 )                            return
c
c...  for modm = 1 evaluate  dmval(j) = mj-th derivative of uj.
c
      do 150 jcomp = 1, ncomp
  150 dmval(jcomp) = 0.d0
      idmz = idmz + 1
      do 170 j = 1, k
         fact = dm(j)
         do 160 jcomp = 1, ncomp
            dmval(jcomp) = dmval(jcomp)  +  fact * dmz(idmz)
            idmz = idmz + 1
  160    continue
  170 continue
      return
c--------------------------------------------------------------------
  900 format(37h ****** domain error in approx ******
     1       /4h x =,d20.10, 10h   aleft =,d20.10,
     2       11h   aright =,d20.10)
      end
      subroutine rkbas (s, coef, k, m, rkb, dm, mode)
c
c**********************************************************************
c
cc rkbas evaluates a mesh-independent runge-kutta basis.
c
c   purpose
c           evaluate mesh independent runge-kutta basis for given s
c
c   variables
c     s      - argument, i.e. the relative position for which
c              the basis is to be evaluated ( 0. .le. s .le. 1. ).
c     coef   - precomputed derivatives of the basis
c     k      - number of collocatin points per subinterval
c     m      - maximal order of the differential equation
c     rkb    - the runge-kutta basis (0-th to (m-1)-th derivatives )
c     dm     - basis elements for m-th derivative
c
c**********************************************************************
c
      implicit real*8 (a-h,o-z)
      dimension coef(k,1), rkb(7,1), dm(1), t(10)
c
      if ( k .eq. 1 )                            go to 70
      kpm1 = k + m - 1
      do 10 i = 1, kpm1
   10   t(i) = s / dfloat(i)
      do 40 l = 1, m
         lb = k + l + 1
         do 30 i = 1, k
           p = coef(1,i)
           do 20 j = 2, k
             p = p * t(lb-j)  + coef(j,i)
   20      continue
           rkb(i,l) = p
   30    continue
   40 continue
      if ( mode .eq. 0 )                         return
      do 60 i = 1, k
         p = coef(1,i)
         do 50 j = 2, k
   50       p = p * t(k+1-j) + coef(j,i)
         dm(i) = p
   60 continue
      return
   70 rkb(1,1) = 1.0d0
      dm(1) = 1.0d0
      return
      end
