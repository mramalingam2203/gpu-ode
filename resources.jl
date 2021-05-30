#=
function juliaSetPixel(z0, c)
    z = z0
    niter = 255
    for i in 1:niter
        abs2(z)> 4.0 && return (i - 1)%UInt8
        z = z*z + c
    end
    return niter%UInt8
end


function calcColumn!(pic, c, n, j)
    x = -2.0 + (j-1)*4.0/(n-1)
    for i in 1:n
        y = -2.0 + (i-1)*4.0/(n-1)
        @inbounds pic[i,j] = juliaSetPixel(x+im*y, c)
    end
    nothing
end

function juliaSetCalc!(pic, c, n)
    for j in 1:n
        calcColumn!(pic, c, n, j)
    end
    nothing
end

function juliaSet(x, y, n=1000, method = juliaSetCalc!, extra...)
    c = x + y*im
    pic = Array{UInt8,2}(undef,n,n)
    method(pic, c, n, extra...)
    return pic
end

frac = juliaSet(-0.79,0.15)
using Plots
plot(heatmap(1:size(frac,1),1:size(frac,2), frac, color=:Spectral))


using BenchmarkTools
@btime juliaSet(-0.79,0.15)


import Base.Threads.@threads
function juliaSetCalcThread!(pic, c, n)
    @threads for j in 1:n
        calcColumn!(pic, c, n, j)
    end
    nothing
end


fracP1 = juliaSet(-0.79,0.15,100,juliaSetCalcThread!)
fracP1 == frac


function juliaSetCalcRecSpawn!(pic, c, n, lo=1, hi=n, ntasks=16)
    if hi - lo > n/ntasks-1
        mid = (lo+hi)>>>1
        finish = Threads.@spawn juliaSetCalcRecSpawn!(pic, c, n, lo, mid, ntasks)
        juliaSetCalcRecSpawn!(pic, c, n, mid+1, hi, ntasks)
        wait(finish)
        return
    end
    for j in lo:hi
        calcColumn!(pic, c, n, j)
    end
    nothing
end

fracP2 = juliaSet(-0.79,0.15,1000,juliaSetCalcRecSpawn!)
fracP2 == frac

@btime juliaSet(-0.79,0.15,1000,juliaSetCalcRecSpawn!)

@btime juliaSet(-0.79,0.15,1000,juliaSetCalcThread!)

=#


using BoundaryValueDiffEq
const g = 9.81
L = 1.0
tspan = (0.0,pi/2)
function simplependulum(du,u,p,t)
    θ  = u[1]
    dθ = u[2]
    du[1] = dθ
    du[2] = -(g/L)*sin(θ)
end
function bc1(residual, u, p, t)
    residual[1] = u[end÷2][1] + pi/2 # the solution at the middle of the time span should be -pi/2
    residual[2] = u[end][1] - pi/2 # the solution at the end of the time span should be pi/2
end
bvp1 = BVProblem(simplependulum, bc1, [pi/2,pi/2], tspan)
sol1 = solve(bvp1, GeneralMIRK4(), dt=0.05)

