using DifferentialEquations
using LinearAlgebra
using ForwardDiff
using Plots
include("stability.jl")
include("rhsFunctions.jl")

ΔE = 0.2;
ΔI = 0.2;
η_0E = 1.5;
η_0I = 1.5;
τE = 12;
τI = 15;

σE,σI,τxE,τxI,
κSEE,κSIE,κSEI,κSII,
αEE,αIE,αEI,αII,
κVEE,κVIE,κVEI,κVII,
VsynEE,VsynIE,VsynEI,VsynII = GenRandParams()


p = κSEE,κSIE,κSEI,κSII,αEE,αIE,αEI,αII,
    κVEE,κVIE,κVEI,κVII,
    VsynEE,VsynIE,VsynEI,VsynII,ΔE,ΔI,η_0E,η_0I,τE,τI


EIGS = get_EigenValues(p)

if maximum(real(EIGS)) >0
    println("should be unstable")
else
    println("should be stable")
end
u0 = init_conds_SS(p) .+0.001
tspan = (0.0,5000.0)

prob2 = ODEProblem(fdeterministic,u0,tspan,p)
sol2 = solve(prob2,saveat = 0.1)

p1 = plot(sol2[1,:])
p2 = plot(sol2[2,:])

plot(p1,p2,layout=2)