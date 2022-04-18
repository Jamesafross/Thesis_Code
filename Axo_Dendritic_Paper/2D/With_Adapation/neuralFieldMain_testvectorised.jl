using ProgressMeter
using Plots
using SparseArrays
using LinearAlgebra
using DifferentialEquations
using NPZ
#Load functions
include("firingrate.jl")
include("Banded_Matrices.jl")
include("rhsFuncs2shunts.jl")
include("speedtracker.jl")
include("SpeedEquation.jl")
include("Guassian_Delta.jl")
const max_T = 150
#mesh paramaters
const min_R = 0;
const max_R =80
const min_X = -0.8
const max_X = 0.8
const dt = 1
const dx = 0.01;
const dr = 0.1;
#sqrt(D*tau) \approx 0.1 - 1 mm (D*tau)
const W = 1.0
const tau = 1  #milliseconds #x10
#finds the grid number for origin on the cable
#PARAMETERS
const v = 8 #mm/s  #connectivity speed
const sigma = 1.0;   #x10^-1 mm #
const h = 0.15; #threshold
const α = 1 / sigma;
const β = 1 / v
const τ = 1   #time constant for cable equation
const D = (0.02)^2 #diffusion coefficient cable equation
const k = 0.1
const alfa = 1.0
 #milliseconds #x10^-3 seconds
const type = 1#type of firing rate function (0 for heaviside, 1 for sigmoid)
const beta = 100#sigmoid steepness (if using...)
const pos = 0.02 #position on cable (multiples of dx)

const V_plus = 70
#GRID DISCRETISATION
#R (r) - soma space|
#X (x) - cable space
const drdr = dr * dr
const dxdx = dx * dx
const R = Int(round(max_R / dr))

const X = Int(round((max_X - min_X) / dx)) + 1
const RX = R * X
const linspace_X = collect(min_X:dx:max_X) #linear grid for X
const xzero = findfirst(linspace_X .== 0) #linear grid for X
const xzero = Int(findfirst(linspace_X .== 0)) #this is the grid number for x=0
const d = Int(round.(xzero + pos / dx))#this is the actual grid number for the position
#deltaF = zeros(RX)
#deltaF[(d - 1) * R + 1:(d) * R] .= 1/dx
deltaF = repeat(smooth_Delta(linspace_X .- (pos),dx),inner = R)
const Drr = (1 / drdr) .* D2z(R, X, RX)
const Dxx = (1 / dxdx) .* D2xCent(R, X, RX,xzero)
#Dxx22 = (1 / dxdx) .*D2xF(R,X,RX)
const Dx1 = (1 / dx) .* D1xF(R, X, RX)
const Dx2 = (1 / dx) .* D1xB(R, X, RX)
const Dx = (Dx2)
#Dxx2[1:R*(d-5),:] =  Dxx[1:R*(d-5),:]
g = zeros(R, X)
psi = zeros(R, X)
st =Int(R/2) - 4
en = Int(R/2) + 4
g[st:en,xzero] .= 1
g = reshape(g, RX)
u0 = zeros(6 * RX)

#adaption parameters
tau_a = 150
beta_a = 0
α_a =1

u0[RX + 1:2*RX] .= g
u0[2 * RX + 1:3 * RX] .= g
u0[3 * RX + 1:4 * RX] .= g
u0[4 * RX + 1:5*RX] .= g./maximum(g)
u0[5 * RX + 1:6*RX] .= 0
h_vec = collect(0.01:0.01:0.2)
k_vec = collect(0:0.002:0.2);
c_numeric_vector_k = zeros(size(k_vec,1))
counter = 1
tspan = (0.0, max_T)
params = [type,β, α, tau, k, v,W, h, beta, deltaF, Dxx, Dx, Drr]
prob = ODEProblem(rhsFun,u0, tspan, params)
print("Solving...")
sol = solve(prob,saveat = dt,progress = true)
print("\n Done!")

#z_sol = reshape(sol[1:RX,:], R, X, size(sol, 2))
#psi_sol = reshape(sol[RX + 1:2 * RX,:], R, X, size(sol, 2))
g_sol = reshape(sol[3 * RX + 1:4 * RX,:], R, X, size(sol, 2))
#z2_sol =  reshape(sol[2 * RX + 1:3 * RX,:], R, X, size(sol, 2))
V_sol = reshape(sol[4 * RX + 1:5 * RX,:], R, X, size(sol, 2)) #voltage
#c_numeric_vector_k = trackSpeed(g_sol[:,d,:],Int(80/dr),Int(90/dr),0.1,dr,0.5,size(g_sol,3),0.1)
#npzwrite("/home/james/c_numeric_alfa.npy",c_numeric_vector_k)

npzwrite("/home/james/PhD_Work/Python_Code/Axo_Dendritic_Paper/PythonGraphingScripts/DataGraphing/Conductance_2D.npy",g_sol)
npzwrite("/home/james/PhD_Work/Python_Code/Axo_Dendritic_Paper/PythonGraphingScripts/DataGraphing/Voltage_2D.npy",V_sol)
@gif for i = 1:1:size(g_sol,3)
    heatmap(g_sol[:,:,i]',levels=LinRange(minimum(g_sol),maximum(g_sol),100))
end
