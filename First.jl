include("Funcs.jl")                

# Parameters(Float64) = 0.0 ;
#=Parameters common to all units (known)
r::Float64
delta::Float64
lambda::Float64
v::Float64
etaA::Float64
etaP::Float64
kappaH::Float64
kappaL::Float64
ksiH::Float64
ksiL::Float64
mu::Float64
alpha::Float64
sigmaP::Float64
sigmaA::Float64
rho::Float64
chiA::Float64
chiP::Float64
theta::Float64
gamma::Float64
phiH::Float64
phiL::Float64
lH::Float64
lL::Float64
=#

#Parameters common to all units (known)
r = Float32(0.05);
delta = Float32(0.15); 
lambda = Float32(0.015);
v = Float32(0.15);
etaA = Float32(-0.1);
etaP = Float32(0.10);
kappaH = log(3);
kappaL = -log(3);


k = [r, delta, lambda, v, etaA, etaP, kappaH, kappaL];

#Parameters common to all units (unknown)
ksiH = Float32(0.1);
ksiL = Float32(0.9);

ksi = [ksiH, ksiL];

# Parameters specific to each unit (known)
mu = Float32(0.0729);
alpha = Float32(0.10);
sigmaP = Float32(0.7559);
sigmaA = Float32(0.1280);
rho =  Float32(-0.0370);
chiA =  Float32(-0.0077);
chiP = Float32(0.2017);


k_g = [mu, alpha, sigmaP, sigmaA, rho, chiA,chiP]

# Parameters specific to each unit (unknown)
theta = Float32(1.80);
gamma = Float32(0.06);
phiH = Float32(0.001);
phiL = Float32(0.50);
lH = Float32(0.8);
lL = Float32(0.4);

Kappa_g = [theta, gamma, phiH, phiL, lH, lL];


# Initial guess for the unknown boundaries 
# x0 = [x_lH m_H m_L x_uH x_uL]

x0 = [0.02 0.15 0.20 0.36 0.41];

b = boundarySol(k, ksi, k_g, Kappa_g, x0);


# Find boundary solutions 

