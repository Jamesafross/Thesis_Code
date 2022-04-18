using DifferentialEquations
using SparseArrays
using LinearAlgebra
include("Banded_Matrices.jl")


function rhs_func(du,u,p,t)
    du[:] = -k * Dxx2 * u[:]
end

function w(x,σ)
    return (exp(-abs(x)/σ))/(2*σ)
end
k = 0.1
α = 1

dx = 0.1
max_X = 10
X = Int(max_X / dx)
Dxx2 = D2xTest(1, X, X, 1)



σ = 1
x = collect(-20:0.01:20)

alfa = 1

W = w.(x,σ)



u0 = zeros(X)
u0[45:55] .= 1


tspan = (0.0, 20.0)
params  = [k,α]
prob = ODEProblem(rhs_func, u0,tspan,params)
print("Solving...")
sol = solve(prob,saveat = 0.02,progress = true)
print("\n Done!")
