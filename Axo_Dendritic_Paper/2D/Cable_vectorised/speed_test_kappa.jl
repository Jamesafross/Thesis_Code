using ProgressMeter
using Plots
using SparseArrays
using LinearAlgebra
using DifferentialEquations
using NPZ
#Load functions
include("firingrate.jl")
include("Banded_Matrices.jl")
include("rhsFuncs2.jl")
include("speedtracker.jl")
include("SpeedEquation.jl")
include("Guassian_Delta.jl")
const max_T = 200
#mesh paramaters
const min_R = 0;
const max_R = 80
const min_X = -0.3
const max_X = 0.3
const dt = 0.004
const dx = 0.02;
const dr = 0.2;
#sqrt(D*tau) \approx 0.1 - 1 mm (D*tau)
const W = 5.0
const tau = 1.0  #milliseconds #x10
#finds the grid number for origin on the cable
#PARAMETERS
const v = 18 #mm/s  #connectivity speed
const sigma = 1.0;   #x10^-1 mm #
#const h = 0.15; #threshold
const α = 1 / sigma;
const β = 1 / v
const τ = 1   #time constant for cable equation
const D = (0.1)^2 #diffusion coefficient cable equation#d #milliseconds #x10^-3 seconds
const type = 0#type of firing rate function (0 for heaviside, 1 for sigmoid)
const beta = 45 #sigmoid steepness (if using...)
#const pos = 0.0 #position on cable (multiples of dx)
const alfa = 1
const h = 0.01
const V_plus = 1
#GRID DISCRETISATION
#R (r) - soma space
#X (x) - cable space
const drdr = dr * dr
const dxdx = dx * dx
const R = Int(round(max_R / dr))

const X = Int(round((max_X - min_X) / dx)) + 1
const RX = R * X
const linspace_X = collect(min_X:dx:max_X) #linear grid for X #linear grid for X
const xzero = Int(findfirst(linspace_X .== 0)) #this is the grid number for x=0
#const d = Int(round.(xzero + pos / dx))#this is the actual grid number for the position
#deltaF = zeros(RX)
#deltaF[(d - 1) * R + 1:(d) * R] .= 1/dx

const Drr = (1 / drdr) .* D2z2(R, X, RX)
const Dxx = (1 / dxdx) .* D2xCent(R, X, RX,xzero)
#Dxx22 = (1 / dxdx) .*D2xF(R,X,RX)
#Dx1 = (1 / dx) .* D1xF(R, X, RX)
const Dx2 = (1 / dx) .* D1xB(R, X, RX)
const Dx = Dx2
#Dxx2[1:R*(d-5),:] =  Dxx[1:R*(d-5),:]
g = zeros(R, X)
psi = zeros(R, X)
st =Int(R/2) - 4
en = Int(R/2) + 4
g[1:10,xzero] .= 1
g = reshape(g, RX)
u0 = zeros(5 * RX)

#u0[RX + 1:2*RX] .= g
u0[2 * RX + 1:3 * RX] .= g
u0[3 * RX + 1:4 * RX] .= g

u0[4 * RX + 1:5*RX] .= g./maximum(g)
kappa_vec = collect(0:0.005:0.15)

c_numeric_vector_k = zeros(size(kappa_vec,1))
pos_vec = [0.0,0.02,0.04,0.06]
global counter = 0
for jj = 1:4
    global counter = counter + 1
    global pos = pos_vec[jj]
    global deltaF = repeat(smooth_Delta(linspace_X .- (pos),dx),inner = R)
    global d = Int(round.(xzero + pos / dx))
for i = 1:size(kappa_vec,1)
    if i == 1
        print(pos)
    end
    global u0
    global kappa_vec
    global k = kappa_vec[i]

    print("κ = ",k, "\n")
tspan = (0.0, max_T)
params = [type,β, α, tau, k, v,W, h, beta, deltaF, Dxx, Dx, Drr]
prob = ODEProblem(rhsFun,u0, tspan, params)
print("Solving...")
condition(u,t,integrator) = reshape(u[3*RX+1:4*RX],R,X)[Int(70/dr),d]>0.1
affect!(integrator) = terminate!(integrator)
cb = DiscreteCallback(condition,affect!)
sol = solve(prob,RK4(),saveat = dt,progress = true,callback = cb)
print("\n Done! \n")
g_sol =  reshape(sol[3 * RX + 1:4 * RX,:], R, X, size(sol, 2))
c_numeric_vector_k[i] = trackSpeed(g_sol[:,d,:],Int(55/dr),Int(60/dr),dt,dr,h,size(g_sol,3),0.1)
print("c = ",c_numeric_vector_k[i],"\n")
end

npzwrite("/tmp/c_numeric_kappa_$counter.npy",c_numeric_vector_k)
end
