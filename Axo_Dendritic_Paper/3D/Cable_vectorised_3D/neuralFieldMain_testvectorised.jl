using ProgressMeter
using Plots
using SparseArrays
using LinearAlgebra
using DifferentialEquations
using NPZ


#Load functions
include("firingrate.jl")
include("banded_matrices.jl")
include("rhsFuncs.jl")
include("Guassian_Delta.jl")
include("params.jl")
#GRID DISCRETISATION
#R (r) - soma space
#X (x) - cable space
dtdt = dt * dt
dr1dr1 = dr1 * dr1
dr2dr2 = dr2 *dr2
dxdx = dx * dx
dtdx = dt * dx




R1 = Int(round(max_R1 / dr1))
R2 = Int(round(max_R2 / dr2))
R = R1 * R2
T = Int(round(max_T / dt))
X = Int(round((max_X - min_X) / dx)) + 1
RX = R * X
linspace_X = collect(min_X:dx:max_X) #linear grid for X
xzero = findfirst(linspace_X .== 0)


linspace_X = collect(min_X:dx:max_X) #linear grid for X
xzero = Int(findfirst(linspace_X .== 0)) #this is the grid number for x=0
d = Int(round.(xzero + pos / dx))#this is the actual grid number for the position

#deltaF = zeros(RX)
#deltaF[(d - 1) * R + 1:(d) * R] .= 1/dx
deltaF = repeat(smooth_Delta(linspace_X .- pos,dx),inner = R)

Dr1r1 = (1 / dr1dr1) .* D2r1(R1, R, X, RX)
Dr2r2 = (1 / dr2dr2) .* D2r2(R1, R2, R, X, RX)
Dxx = (1 / dxdx) .* D2x(R, X, RX)
#Dx = (1 / 2*dx) .* D1x(R, X, RX)
Dx = (1 / dx) .* D1x(R, X, RX)

g = zeros(R1,R2, X)
psi = zeros(R1,R2, X)

st =Int(R1/2) - 4
en = Int(R1/2) + 4
g[st:en,st:en,xzero] .= 1
psi[st:en,st:en,xzero] .= 1


g = reshape(g, RX)
psi = reshape(g, RX)
V = psi
u0 = zeros(5 * RX)
u0[RX + 1:2 * RX] .= psi
u0[4*RX+1:5*RX] .= V
u0[3 * RX + 1:4 * RX] .= g

tspan = (0.0, max_T)
params = [type,β, α, tau, k, v,W, h, beta, deltaF, Dxx, Dx, Dr1r1,Dr2r2]
condition(u,t,integrator) = abs(maximum(u[3*RX+1:4*RX]) - minimum(u[3*RX+1:4*RX])) <0.1
affect!(integrator) = terminate!(integrator)
cb = DiscreteCallback(condition,affect!)
prob = ODEProblem(rhsFun, u0, tspan, params)
print("Solving...")
sol = solve(prob,saveat = 0.5,progress = true,callback=cb)
print("\n Done!")

g_sol = reshape(sol[3 * RX + 1:4 * RX,:], R1,R2, X, size(sol, 2)) #conductance
V_sol = reshape(sol[4 * RX + 1:5 * RX,:], R1,R2, X, size(sol, 2)) #voltage
 
