	 program main
	 implicit real*8 (a-h, o-z)
	 double precision zeta(5), fspace(40000), tol(2), sval(3), elval(3)
	 integer m(2), ispace(2500), ltol(2), ipar(11)
	 double precision z(5)
	 common en, s , el, cons
	 external fsub, dfsub, gsub, dgsub, solutn
	 data sval/0.20, 0.1d0, 0.05d0/. elval/60.d0, 120.d0, 200.d0/

	 en = .2d0
	 cons = .500*(3.d0-en)
	 ncomp =2
	 m(1) =2
	 m(2) = 3
	 aleft = 0.d0
	 aright = 1.d0

	 zeta(1)=0.d0
	 zeta(2) = 0.d0
	 zeta(3) = 0.d0
	 zeta(4) = 1.d0
	 zeta(5) = 1.d0

	 ipar(1)=1
	 ipar(2)=4
	 ipar(3)=4
	 ipar(4)=2
	 ipar(5)=40000
	 ipar(6)=2500
	 ipar(7)=0
	 ipar(8)=0
	 ipar(9)=1
	 ipar(10)=0
	 ipar(11) =0

	 ltol(1) =1
	 ltol(2) = 3
	 tol(1) = 1.d-5
	 tol(2) = 1.d-5

	 do 777 ijk=1,3
		s = sval(ijk)
		el = elval(ijk)
		if(ijk.eq.1) go to 701
		ipar(9) = 3
		ipar(3) = ispace(1)
701	 continue
	 write(6, 100) en, s, el
100	 format(1h1, 38h ROTATING LOW OVER A STATIONARY DISK,
			/19H parameters - n=, f5.2, 6h s=,f5.2, 6H l =,f6.1/)

	  call colsys(ncomp, m, aleft, aright, zeta, ipar, ltol, tol,
	  fixpnt, ispace, fspace, iflag, fsub, dfsub, gsub, dgsub, solutn)

	 if (iflag.ne.1) stop

		is6 = ispace(6)
		is5 = ispace(1) +2
		x = 0.d0
		write(6,201)
201		format(1h1, 44h	x	g	dg	, 38h	h	dh	d2h/)
202		format(6d15.5)
		np1 - el+1.5d0
		do 555 iii= 1, np1
		call approx(ii, x, z, fspace(is6), fspace(1), ispace(1), fspace(is5)
		, ispace(2), ncomp, m, ispace(4), 1, dm,0)
		xl = x*el
		z(2) = z(2)/el
		z(4) = z(4)/el
		z(5) = z(5)/el/el
		write(6,202) xl, 2
		x = x+1.d0/el
555		continue
777		continue

	stop 
	end

	 subroutine solutn(x,z, dmval)
	 implicit real*8(a-h, o-z)
	 common en, s, el, cons
	 double precision z(5), dmval(2)
	 ex = dexp(-el*x)
	 z(1) = 1.d0-ex
	 z(2) = el*ex
	 z(3) = -el**2*x**2*ex
	 z(4) = (el**3*x**2-2.d0*el**2*x)*ex
	 z(5) = (-3l**4*x**2+4.d0*el*83*x-2.d0*el**2)*ex
	 dmval(1) = -el*z(2)
	 dmval(2) = (el**5*x*x - 6.d0*el**4*x +6.d0*el**3)*ex
	 return
	 end

	 subroutine fsub(x,z,f)
	 implicit real*8 (a-h, o-z)
	 double precision z(1), f(1)
	 common en, s, el, cons
	 f(1) = -el*(cons*z(3)*z(2)*(en-1.d0)*z(4)*z(1))*el**2*s*(z(1)-1.d0)
	 f(2) = -el*(cons*z(3)*z(5)*en*z(4)**2)+el**2*s*z(4)+el*83*(ad0-z(1)**2)

	 subroutine dfsub(x,f,df)
	 implicit real*8 (a-h, o-z)
	 double precision z(1), df(2,1)
	 common en, s el, cons
	 d(1,1) = -el*(en-1.d0)*z(4)+el**2*s
	 df(1,2) = -el*cons*z(3)
	 df(1,3) = -el *cons*z(2)
	 df(1,4) = -el *(en-1.d0)*z(1)
	 df(1,5) = 0.d0
	 df(2,1) = -el**3*2.d0*z(1)
	 df(2,2) = 0.d0
	 df(2,3) = -el*cons*z(5)
	 df(2,4) = -el * en*2.d0*z(4) +el**2*s
	 df(2,5) = -el*cons*z(3)
	 return
	 end

	 subroutine gsub(i,z,g)
	 double precision z(1),g
	 goto (1,2,3,4,3),i
1    g = z(1)
	 return
2    g= z(3)
	 return
3    g=z(4)
	 return
4    g = z(1)-1.d0
	 return end


	 subroutine gsub(i,z,dg)
	 double precision z(1),dg

	 do 10 j = 1,5
10  dg(j) = 0.d0
	goto (1,2,3,1,3),i
1   dg(1) = 1.d0
	return
2   dg(3) = 1.d0
	return
3   dg(4) = 1.d0
	return
	end
	return
4   g = z(1)-1.d0
	return end





