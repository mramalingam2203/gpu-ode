
program main
use auxiliary
!use adds

implicit none
integer a,b,c
real k(8), ksi(2), k_g(7), Kappa_g(6)

real r,delta,lambda,v,etaA,etaP,kappaH,kappaL,& 
	ksiH,ksiL,mu,alpha,sigmaP,sigmaA,rho,chiA,chiP,& 
	theta,gamma,phiH,phiL,lH,lL


type (params)::paraL,paraH



double precision aleft
double precision aright
external dfsub
external dgsub
double precision fixpnt
double precision fspace(52000)
!external fsub
external gsub
external guess
integer i
integer iflag
integer ipar(11)
integer ispace(3000)
integer j
integer ltol(3)
integer m(2)
integer ncomp
integer n_out
double precision tol(3)
double precision x
double precision z(6)
double precision zeta(5)


ipar(1) = 1
ipar(2) = 2
ipar(3) = 2
ipar(4) = 4
ipar(5) = 40000
ipar(6) = 2500
ipar(7) = 0
ipar(8) = 1
ipar(9) = 1
ipar(10) = 0
ipar(11) = 0


ltol(1) = 1
ltol(2) = 2
ltol(3) = 4

tol(1) = 1.0D-07
tol(2) = 1.0D-07
tol(3) = 1.0D-07



ncomp = 2
m(1) = 1
m(2) = 4
aleft = 0.0D+00
aright = 1.0D+00


r = 0.05
delta = 0.15
lambda = 0.015
v = 0.15
etaA = -0.10
etaP = 0.10
kappaH = log(3.0)
kappaL = -log(3.0)

ksiH = 0.1
ksiL = 0.9

mu = 0.0729
alpha = 0.10
sigmaP = 0.7559
sigmaA = 0.1280
rho =  -0.0370
chiA =  -0.0077
chiP = 0.2017

theta = 1.80
gamma = 0.06
phiH = 0.001
phiL = 0.50
lH = 0.8
lL = 0.4




 call make_params(k, ksi, k_g, Kappa_g,r,delta,lambda,v,etaA,etaP,&
kappaH,kappaL,ksiH,ksiL,mu,alpha,sigmaP,sigmaA,rho,chiA,chiP,& 
 theta,gamma,phiH,phiL,lH,lL)

 call parms_type(paraL, paraH,k, ksi, k_g, Kappa_g)

 call fsub(paraL,paraH)
 end







subroutine dfsub()

end

subroutine gsub


end

subroutine dgsub 
end

subroutine guess
end
