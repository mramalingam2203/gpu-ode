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