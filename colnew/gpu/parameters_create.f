!module parameters_create
!contains

!subroutine pars(params,r,delta,lambda,v,etaA,etaP,kappaH,kappaL,& 
!ksiH,ksiL,mu,alpha,sigmaP,sigmaA,rho,chiA,chiP,& 
!theta,gamma, gamma_H,phiH,phiL,lH,lL)


program main

integer i

real r,delta,lambda,v,etaA,etaP,kappaH,kappaL,& 
ksiH,ksiL,mu,alpha,sigmaP,sigmaA,rho,chiA,chiP,& 
theta,gamma, gamma_H,phiH,phiL,lH,lL

real k(8), ksi(2), k_g(7), Kappa_g(6)

type ragged_array
real,allocatable::vec(:)
end type ragged_array

type(ragged_array),allocatable::params(:)
allocate(params(1)%vec(8))
allocate(params(2)%vec(2))
allocate(params(3)%vec(7))
allocate(params(4)%vec(6))

r = 1.0
delta = 1.0
lambda = 1.0
v = 1.0
etaA = 1.0
etaP = 1.0
kappaH = 1.0
kappaL = 1.0 
ksiH = 1.0
ksiL = 1.0
mu = 1.0
alpha = 1.0
sigmaP = 1.0 
sigmaA = 1.0
rho =1.0
chiA = 1.0
chiP = 1.0
theta = 1.0
gamma = 1.0
gamma_H =1.0
phiH = 1.0
phiL = 1.0
lH = 1.0
lL = 1.0


	k(1) = r
	k(2) = delta 
	k(3) = lambda
	k(4) = v
	k(5) = etaA
	k(6) = etaP 
	k(7) = kappaH
	k(8) = kappaL

	ksi(1) = ksiH
	ksi(2) = ksiL
	
	k_g(1) = mu
	k_g(2) = alpha
	k_g(3) = sigmaP
	k_g(4) = sigmaA
	k_g(5) = rho
	k_g(6) = chiA
	k_g(7) = chiP

	Kappa_g(1) = theta
	Kappa_g(2) = gamma
	Kappa_g(3) = phiH
	Kappa_g(4) = phiL
	Kappa_g(5) = lH
	Kappa_g(6) = lL

do i = 1, 8
	params(i)%vec(i) = k(i)
enddo


do i = 1, 2
	params(i)%vec(i) = ksi(i)
enddo

do i = 1, 7
	params(i)%vec(i) = k_g(i)
enddo

do i = 1, 7
	params(i)%vec(i) = Kappa_g(i)
enddo


deallocate(params)

end program main


!end subroutine pars
!end module parameters_create



