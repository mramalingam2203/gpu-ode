      implicit real*8 (a-h,o-z)
      double precision zeta(4), fspace(40000), tol(4), z(4)

      integer m(2), ipar(11), ispace(2500), ltol(4)
      common eps, dmu, eps4mu, gamma, xt
      external solutn, fsub, dfsub, gsub, dgsub

      gamma = 1.1d0
      eps = 0.01d0
      dmu = eps
      eps4mu = eps**4/dmu

      xt = dsqrt(2.d0*(gamma-1.d0)/gamma)
      write(*,*) gamma, xt, eps, dmu, eps4mu

c     format(1h1,27hdimpling of spherical caps./.8h gamma=,7.2/gh xt = ,d12.5/ 6h eps=,d12.5/gh mu= . 12h eps**4/mu =,d12.5)

      ncomp =2
      m(1) =2
      m(2) = 2

      aleft = 0.d0
      aright = 1.d0

      zeta(1) = 0.d0
      zeta(2) = 0.d0
      zeta(3) = 1.d0
      zeta(4) = 1.d0

      ipar(1) = 1
      ipar(2) = 2
      ipar(3) = 2
      ipar(4) = 4
      ipar(5) = 40000
      ipar(6) = 2500
      ipar(7) = 0
      ipar(8) = 0
      ipar(9) = 1
      ipar(10) = 0
      ipar(11) = 0
         
      do 10 i=1,4
        ltol(i) = i
 10    tol(i) = 1.d-2

      call colsys(ncomp, m , aleft, aright, zeta, ipar, ltol, tol, 
     & fixpt, ispace, fspace, iflag, fsub, dsub, gsub, dgsub, solutn)

      np1 = 201

      do 555 iii=1,np1
        call appsln(x,z,fspace, ispace)
      x = x + 0.005d0
555   continue

      stop
      end


      subroutine solutn(x,z,dmval)
      implicit real*8 (a-h, o-z)
      common eps, dmu, eps4mu, gamma, xt
      dimension z(4), dmval(2)
      cons = gamma *x*(1.d0-5d0*x*x)
      dcons = gamma *(1.d0-5d0*x*x)
      d2cons = -3.d0 * gamma*x
      if (x.gt.xt) go to 10
10      z(1) = 2.d0*x
        z(2) = 2.d0
        z(3) = -2.d0*x + cons
        z(4) = -2.d0 + dcons
      dmval(2) = d2cons
      goto 20
      z(1) = 0.d0
      z(2) = 0.d0
      z(3) = -cons
      z(4) = -dcons
      dmval(2) = -d2cons
20    dmval(1)= 0.d0
      return
      end


      subroutine fsub(x,z,f)
      implicit real*8 (a-h, o-z)
      dimension z(4), f(2)
      common eps, dmu, eps4mu, gamma, xt
      f(1) = (z(1)/x/x) - (z(2)/x) + (z(1)-z(3)*(1.d0-z(1)/x)) 
     & -gamma* x*(1.d0-x*x/2.))/eps4mu
      f(2) = z(3)/x/x -z(4)/x +z(1)*(1.d0-z(1)/2.d0/x)/dmu
      return 
      end


      subroutine dfsub(x,f,df)
      implicit real*8 (a-h, o-z)
      dimension z(4), df(2,4)
      common eps, dmu, eps4mu, gamma, xt
      df(1,2) = 1.d0/x/x +(1.d0+z(3)/x)/eps4mu
      df(1,2) = -1.d0/x
      df(1,3) = (1.d0+z(1)/x)/eps4mu
      df(1,4) = 0.d0

      df(2,1) = (1.d0-z(1)/X)/dmu
      df(2,2) = 0.d0
      df(2,3) = 1.d0/x/x
      df(2,4) = -1.d0/x
      return
      end


      subroutine gsub(i,z,g)
      implicit real*8 (a-h, o-z)
      dimension z(4)
      goto (1,2,1,3),I
1     g = z(1)
      return
2     g = z(3)
      return
3     g = z(4) - 0.3d0*z(3)*0.7d0
      return
      end


      subroutine dgsub(i,z, dg)
      implicit real*8 (a-h, o-z)
      dimension z(4), dg(4)
      do 10 j = 1,4
10    dg(j) = 0.d0
      goto (1,2,1,3),I
1     dg(1) = 1.d0
      return
2     dg(3) = 1.d0
      return
3     dg(4) = 1.d0
      dg(3) = -0.3d0
      return
      end




     `