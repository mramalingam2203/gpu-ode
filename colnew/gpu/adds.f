
	! Create the high and low range parameters
	! returns a 
	
	subroutine parameters(paramsL, params,k, ksi, k_g, Kappa_g)
	real k(8), ksi(2), k_g(7), Kappa_g(6)

	type params
		real::r, delta, v, sigmaP, sigmaA, theta, lambda,l, gamma,phi, ksi, etaA, etaP,chiA, chiP, rho, alpha, mu
	end type


	type (params) :: paramsL,paramsH

	paramsH%r = k(1);
	paramsH%delta = k(2);
	paramsH%v = k(4);
	paramsH%sigmaP = k_g(3);
	paramsH%sigmaA = k_g(4);
	paramsH%theta = Kappa_g(1);
	paramsH%lambda = k(3);
	paramsH%l = Kappa_g(5);
	paramsH%gamma = Kappa_g(2);
 	paramsH%phi = Kappa_g(3);
	paramsH%ksi = exp(k(7)) * ksi(1);
	paramsH%etaA = k(5);
	paramsH%etaP = k(6);
	paramsH%chiA = k_g(6);
	paramsH%chiP = k_g(7);
	paramsH%rho = k_g(5);
	paramsH%alpha = k_g(2) - paramsH%sigmaA * paramsH%etaA * paramsH%chiA;
	paramsH%mu = k_g(1) - paramsH%sigmaP * paramsH%etaP * paramsH%chiP;

	paramsL%r = k(1);
	paramsL%delta = k(2);
	paramsL%v = k(4);
	paramsL%sigmaP = k_g(3);
	paramsL%sigmaA = k_g(4); 
	paramsL%theta = Kappa_g(1);
	paramsL%lambda = k(3);
	paramsL%l = Kappa_g(6);
	paramsL%gamma = Kappa_g(2);
	paramsL%phi = Kappa_g(4);
	paramsL%ksi = exp(k(8)) * ksi(2);
	paramsL%etaA = k(5);
	paramsL%etaP = k(6);
	paramsL%chiA = k_g(6);
	paramsL%chiP = k_g(7);
	paramsL%rho = k_g(5);
	paramsL%alpha = k_g(2) - paramsH%sigmaA * paramsH%etaA * paramsH%chiA;
	paramsL%mu = k_g(1) - paramsH%sigmaP * paramsH%etaP * paramsH%chiP;

		

	return
	end