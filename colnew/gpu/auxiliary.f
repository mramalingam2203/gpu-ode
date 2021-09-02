module auxiliary

implicit none

type params
		real::r, delta, v, sigmaP, sigmaA, theta, lambda,l, gamma,phi, ksi, etaA, etaP,chiA, chiP, rho, alpha, mu
	end type


! declarations
contains

subroutine make_params(k, ksi, k_g, Kappa_g,r,delta,lambda,v,etaA,etaP,kappaH,kappaL,& 
	ksiH,ksiL,mu,alpha,sigmaP,sigmaA,rho,chiA,chiP,& 
	theta,gamma,phiH,phiL,lH,lL)


real r,delta,lambda,v,etaA,etaP,kappaH,kappaL,& 
	ksiH,ksiL,mu,alpha,sigmaP,sigmaA,rho,chiA,chiP,& 
	theta,gamma,phiH,phiL,lH,lL
	real k(8), ksi(2), k_g(7), Kappa_g(6)


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

	end subroutine make_params



	! Create the high and low range parameters
	! returns a 
	
	subroutine parms_type(paramsL, paramsH,k, ksi, k_g, Kappa_g)
	real k(8), ksi(2), k_g(7), Kappa_g(6)

	type (params) :: paramsL,paramsH

	paramsH%r = k(1)
	paramsH%delta = k(2)
	paramsH%v = k(4)
	paramsH%sigmaP = k_g(3)
	paramsH%sigmaA = k_g(4)
	paramsH%theta = Kappa_g(1)
	paramsH%lambda = k(3)
	paramsH%l = Kappa_g(5)
	paramsH%gamma = Kappa_g(2)
 	paramsH%phi = Kappa_g(3)
	paramsH%ksi = exp(k(7)) * ksi(1)
	paramsH%etaA = k(5)
	paramsH%etaP = k(6)
	paramsH%chiA = k_g(6)
	paramsH%chiP = k_g(7)
	paramsH%rho = k_g(5)
	paramsH%alpha = k_g(2) - paramsH%sigmaA * paramsH%etaA * paramsH%chiA
	paramsH%mu = k_g(1) - paramsH%sigmaP * paramsH%etaP * paramsH%chiP

	paramsL%r = k(1)
	paramsL%delta = k(2)
	paramsL%v = k(4)
	paramsL%sigmaP = k_g(3)
	paramsL%sigmaA = k_g(4)
	paramsL%theta = Kappa_g(1)
	paramsL%lambda = k(3)
	paramsL%l = Kappa_g(6)
	paramsL%gamma = Kappa_g(2)
	paramsL%phi = Kappa_g(4)
	paramsL%ksi = exp(k(8)) * ksi(2)
	paramsL%etaA = k(5)
	paramsL%etaP = k(6)
	paramsL%chiA = k_g(6)
	paramsL%chiP = k_g(7)
	paramsL%rho = k_g(5)
	paramsL%alpha = k_g(2) - paramsH%sigmaA * paramsH%etaA * paramsH%chiA
	paramsL%mu = k_g(1) - paramsH%sigmaP * paramsH%etaP * paramsH%chiP

	end


subroutine  fsub(parL,parH)
	type (params)::parL,parH
	write(*,*)paraL
end





end module auxiliary