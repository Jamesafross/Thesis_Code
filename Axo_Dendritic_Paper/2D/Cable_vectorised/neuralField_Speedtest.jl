


using ProgressMeter
using Plots
using DelimitedFiles
using SparseArrays
using SuiteSparse
using LinearAlgebra
using DifferentialEquations
using NPZ
using Sundials

#Load functions
include("firingrate.jl")
include("banded_matrices.jl")
include("rhsFuncs.jl")
include("speedtracker.jl")
include("SpeedEquation.jl")
include("Guassian_Delta.jl")
#GRID DISCRETISATION
#R (r) - soma space
#X (x) - cable space
min_T = 0;
max_T = 20
min_R = 0;
max_R = 40
min_X = -3
max_X = 3
dt = 0.002
dx = 0.25;
dr = 0.2;
dtdt = dt * dt
drdr = dr * dr
dxdx = dx * dx
dtdx = dt * dx
R = Int(round(max_R / dr))
T = Int(round(max_T / dt))
X = Int(round((max_X - min_X) / dx)) + 1
RX = R * X
linspace_X = collect(min_X:dx:max_X) #linear grid for X
xzero = findfirst(linspace_X .== 0) #finds the grid number for origin on the cable
#PARAMETERS
v = 15 #connectivity speed
sigma = 1.0; #
h = 0.01; #threshold
α = 1 / sigma;
β = 1 / v
τ = 1.0   #time constant for cable equation
D = 0.2#diffusion coefficient cable equation
#diffusion coefficient for x spatial interaction for wave PDE
alfa = 1.0
W = 1.0
type = 0 #type of firing rate function (0 for heaviside, 1 for sigmoid)
beta = 20 #sigmoid steepness (if using...)
pos = 1 #position on cable (multiples of dx)

V_plus = 1.0
#d = Int(round.(xzero + pos / dx))#this is the actual grid number for the position

#deltaF = zeros(RX)
#deltaF[(d - 1) * R + 1:(d) * R] .= 1/dx


Drr = (1 / drdr) .* D2z(R, X, RX)
Dxx = (1 / dxdx) .* D2xCent(R, X, RX,xzero)
Dx = (1 / 2*dx) .* D1xF(R, X, RX) + D1xB(R,X,RX)
Dx = (1 / dx) .* D1xB(R, X, RX)
g = zeros(R, X)
psi = zeros(R, X)
st = 98
en = 102
g[st:en,xzero] .= 1
psi[st:en,xzero] .= 1
W = 1
tau = 1


g = reshape(g, RX)
psi = reshape(g, RX)
V = psi
u0 = zeros(5 * RX)
u0[RX + 1:2 * RX] .= psi
u0[4*RX+1:5*RX] .= V
u0[3 * RX + 1:4 * RX] .= g
h_vec = collect(0.01:0.01:0.2)
k_vec = collect(0:0.01:0.2)
c_numeric_vector_k = zeros(size(k_vec,1))
pos_vec = [-1.0, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1]
d_vec = Int.(round.(xzero .+ pos_vec ./ dx))
T_vec = [60.0, 30.0, 30.0, 20.0, 20.0, 20.0, 30.0, 30.0, 60.0]



for j = 1:size(pos_vec,1)

    for i = 1:size(k_vec,1)
    global deltaF = repeat(smooth_Delta(linspace_X .- pos_vec[j],dx),inner = R)
    global  k = k_vec[i]
    global d_vec
    global T_vec
    tspan = (0.0, T_vec[j])
    params = [type,β, α, tau, k, v,W, h, beta, deltaF, Dxx, Dx, Drr]
    prob = ODEProblem(rhsFun, u0, tspan, params)
    print("\n Solving... (",i,"/",size(k_vec,1),")"," at position ",pos_vec[j])
    sol = solve(prob,saveat = 0.001,progress = true)
    print("\n Done!")

    #psi_sol = reshape(sol[1 * RX + 1:2 * RX,:], R, X, size(sol, 2))
    g_sol = reshape(sol[3 * RX + 1:4 * RX,:], R, X, size(sol, 2))
    V_sol = reshape(sol[4 * RX + 1:5 * RX,:], R, X, size(sol, 2))
    c_numeric_vector_k[i] = trackSpeed(g_sol[:,d_vec[j],:],Int(25/dr),Int(35/dr),0.01,0.1,0.5,size(g_sol,3))

    end

    p = pos_vec[j]
    npzwrite("WaveSpeedk_$p.npy", c_numeric_vector_k)
end
#c_analytic = get_speed(sigma,alfa,tau,D,v,h,pos,7)

#print("\n Analytic Speed = ",c_analytic,"\n Numerical Speed = ",c_numeric)
