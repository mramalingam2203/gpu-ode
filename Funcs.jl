mutable struct params
    r::Float64;delta::Float64;v::Float64;sigmaP::Float64;sigmaA::Float64;theta::Float64;lambda::Float64;
    l::Float64;gamma::Float64;phi::Float64;ksi::Float64;etaA::Float64;etaP::Float64;chiA::Float64;chiP::Float64;
    rho::Float64;alpha::Float64;mu::Float64
end

function parameters(k, ksi, k_g, Kappa_g)
    #Create parameter structures for each state s = H, L
    paramsH =  params(k[1],k[2],k[4],k_g[3],k_g[4],Kappa_g[1],k[3],Kappa_g[5],Kappa_g[2],Kappa_g[3],
        exp(k[7])*ksi[1],k[5],k[6],k_g[6],k_g[7],k_g[5],k_g[2] - sigmaA * etaA * chiA,k_g[1] - sigmaP *etaP * chiP)
     
    paramsL = params(k[1],k[2],k[4],k_g[3],k_g[4],Kappa_g[1],k[3],Kappa_g[6],Kappa_g[2],Kappa_g[4],
        exp(k[8]) * ksi[2],k[5],k[6],k_g[6],k_g[7],k_g[5],k_g[2] - sigmaA * etaA * chiA,k_g[1] - sigmaP * etaP * chiP)
    
    return [paramsL paramsH]
end


function boundarySol(k, ksi, k_g, Kappa_g, x0)
    # This function finds the unknown ODE boundaries using a direct serach 
    # method. The boundary solutions minimize the sum of squared distances
    # from the boundary conditions.
      
    parms =parameters(k, ksi, k_g, Kappa_g);
    print(parms)
   
    A = [-1 0 0 0 0; 1 -1 0 0 0; 0 0 -1 0 0; 0 1 -1 0 0;
        0 0 1 -1 0; 0 0 0 1 -1];
    b = [0.01 0.01 0.01 0.01 0.01 0.01];
    
    #fun = sysODE(parms[1], parms[2]);

    #options = optimoptions('patternsearch','Cache','on',...
    #    'UseParallel', true, 'UseCompletePoll', true, ...
    #    'UseVectorized', false, 'PlotFcn', {@psplotbestf, @psplotfuncount},...
    #    'TolFun', 1e-6, 'CompleteSearch', 'on', 'SearchMethod', @searchlhs, ...
    #    'OutputFcn', @stopfn);s
    
    #x = patternsearch(fun, x0, A, -b, [], [], [], [], [], options);
    #print(x)
end



function sysODE(paramsH, paramsL)
    # Computes the cost function based on how well the solution for given
    # boundaries sastifies the boundary conditions
    
    d = f_solver(paramsH, paramsL, x)
    #d =[1,2, 3, 4, 5]
    F = [0.0,0.0,0.0,0.0,0.0]
    
    F[1] = d[1] - (1 + paramsH.gamma)
    F[2] = d[2] - (1 + paramsH.gamma)
    F[3] = d[3] - (1 + paramsL.gamma)
    F[4] = d[4] - 1
    F[5] = d[5] - 1
    #Cost = sum(F.^2)^(1 / 2)
    return nothing
end

#boundarySol -> sysODE -> f_solver -> two_ode         


function d = f_solver(paramsH, paramsL, p)
    # This function provides the values of the first derivatives at
    # boundaries in each state s = H, L, for given boundaries in p
    
    ode = @(c, f, region) twoode(c, f, region, paramsH, paramsL)
    
    f_init = [1 1 1 1]
    x_init = [0 p[1] p[1] p[2] p[2] p[3] p[3] p[4] p[4] p[5]
    solinit = bvpinit(x_init, f_init)
        
    bc = @(yl, yr) twobc(yl, yr, p, paramsH, paramsL)
    sol = bvp4c(ode, bc, solinit)
    x = [0 p[1] p[2] p[3] p[4] p[5]]
    y = deval(sol, x)
    
    d[1] = y[2 2] # x_lH (lower boundary for state H)
    d[2] = y[2 3] # m_H (return point for state H)
    d[3] = y[4 4] # m_L (return point for state L)
    d[4] = y[2 5] # x_uH (upper boundary for state H)
    d[5] = y[4 6] # x_uL (upper boundary for state L)

end

    
function dfdc = twoode(c, f, region, paramsH, paramsL)
    # ODE equations
    a13 = 0.5 * (paramsH.sigmaP^2 * c^2 - 2 * paramsH.sigmaP *paramsH.sigmaA * paramsH.rho * c + paramsH.sigmaA^2)

    a23 = 0.5 * (paramsL.sigmaP^2 * c^2 - 2 * paramsL.sigmaP *paramsL.sigmaA * paramsL.rho * c + paramsL.sigmaA^2)
    
    i1 = (1 / paramsH.theta) * ((f[1] / f[2]) - c - 1) + paramsH.v
    g_i1 = 0.5 * paramsH.theta * (i1 - paramsH.v)^2
    
    
    i2 = (1 / paramsL.theta) * ((f[3] / f[4]) - c - 1) + paramsL.v
    g_i2 = 0.5 * paramsL.theta * (i2 - paramsL.v)^2
    
    
    a2 = (paramsH.r * f[1] - ((paramsH.r - paramsH.lambda) * c +
        paramsH.alpha - i1 - g_i1) * f[2] - (i1 - paramsH.delta
        + paramsH.mu) * (f(1) - c * f(2)) - paramsH.ksi *
        (f[3] - f[1])) / a13
            
    a4 = (paramsL.r * f[3] - ((paramsL.r - paramsL.lambda) * c + 
        paramsL.alpha - i2 - g_i2) * f[4] - (i2 - paramsL.delta
        + paramsL.mu) * (f[3] - c * f[4]) - paramsL.ksi *
        (f[1] - f[3])) / a23
    
    dfdc = [f[2] a2 f[4] a4]
    
    return dfdc
end
            
#=
function twobc(yl, yr, p, paramsH, paramsL)
    % Boundary conditions 

    i1 = (1 / paramsH.theta) * (yr[1 4] - p[4] - 1) + paramsH.v;
    g_i1 = 0.5 * paramsH.theta * (i1 - paramsH.v)^2;
    
    
    i2 = (1 / paramsL.theta) * (yr[3 5] - p[5] - 1) + paramsL.v;
    g_i2 = 0.5 * paramsL.theta * (i2 - paramsL.v)^2;
    
    res = [ 
            yl[3 1] - yr[3 3] + paramsL.phi +
            (1 + paramsL.gamma) * p[3] 
            yl[1 2] - yr[1 2] + paramsH.phi +
            (1 + paramsH.gamma) * (p[2] - p[1])  
            
            yl[1 2] - yr[1 1] 
            yl[2 2] - yr[2 1] 
            yl[3 2] - yr[3 1] 
            yl[4 2] - yr[4 1] 
            
            # Conditions at m_H
            yl[1 3] - yr[1 2] 
            yl[2 3] - yr[2 2] 
            yl[3 3] - yr[3 2] 
            yl[4 3] - yr[4 2] 
            
            # Conditions at m_L
            yl[1 4] - yr[1 3] 
            yl[2 4] - yr[2 3] 
            yl[3 4] - yr[3 3] 
            yl[4 4] - yr[4 3]  
            
            # Conditions at upper boundary 1
            yl[2 5] - yr[2 4] 
            yl[3 5] - yr[3 4]  
            yl[4 5] - yr[4 4] 
            
            # Conditions at upper boundary 2
            yr[1 5] - yl[1 5] - p[5] + p[4];
            
            paramsH.r * yr[1 4] - (i1 - paramsH.delta + paramsH.mu) *(yr(1, 4) - p(4)) - ((paramsH.r - paramsH.lambda) * p(4) +
            paramsH.alpha - i1 - g_i1) - paramsH.ksi *(yr[3 4] - yr[1 4])
            
            
            paramsL.r * yr[3 5] - (i2 - paramsL.delta + paramsL.mu) *(yr[3 5] - p[5])- ((paramsL.r - paramsL.lambda) * p[5] +
            paramsL.alpha - i2 - g_i2) - paramsL.ksi * (yr[1 5] - yr[3 5])
        ]            
          
end
=#