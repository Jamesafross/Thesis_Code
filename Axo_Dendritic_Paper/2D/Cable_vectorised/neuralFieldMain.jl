using ProgressMeter
using Plots
using SparseArrays
using LinearAlgebra
using DifferentialEquations
using NPZ
#Load functions

include("Banded_Matrices.jl")
include("rhsFuncs2shunts.jl")
  max_T = 500
#mesh paramaters
  min_R = -2pi;
  max_R =2pi
  min_X = -0.5
  max_X = 0.5
  dt = 0.1
  dx = 0.1;
  dr = pi/(2^6);
#sqrt(D*tau) \approx 0.1 - 1 mm (D*tau)
  W = 2.0
  tau = 1  #milliseconds #x10
#finds the grid number for origin on the cable
#PARAMETERS
  v = 8 #mm/s  #connectivity speed
  σ = 1.0;   #x10^-1 mm #
  h = 0; #threshold
  β = 1 / v
  τ = 1   #time constant for cable equation
  D = (0.0)^2 #diffusion coefficient cable equation
  k = 0.0
  α = 1.0
  V_plus = 70
  g_c = 0
V_ss = (g_c*V_plus)/((1/τ) + g_c)

 #milliseconds #x10^-3 seconds

  type = 2 #type of firing rate function (0 for heaviside, 1 for sigmoid 2 for tanh)
  beta = 100#sigmoid steepness (if using...)
  pos = 0.2 #position on cable (multiples of dx)


#GRID DISCRETISATION
#R (r) - soma space|
#X (x) - cable space
  drdr = dr * dr
  dxdx = dx * dx
  R = Int(round(max_R / dr))

  X = Int(round((max_X - min_X) / dx)) + 1
  RX = R * X
  linspace_X = collect(min_X:dx:max_X) #linear grid for X
  xzero = findfirst(linspace_X .== 0) #linear grid for X
  xzero = Int(findfirst(linspace_X .== 0)) #this is the grid number for x=0
  d = Int(round.(xzero + pos / dx))#this is the actual grid number for the position
#deltaF = zeros(RX)
#deltaF[(d - 1) * R + 1:(d) * R] .= 1/dx
deltaF = repeat(smooth_Delta(linspace_X .- (pos),dx),inner = R)
  Drr = (1 / drdr) .* D2z(R, X, RX)
  Dxx = (1 / dxdx) .* D2xCent(R, X, RX,xzero)
#Dxx22 = (1 / dxdx) .*D2xF(R,X,RX)
  Dx1 = (1 / dx) .* D1xF(R, X, RX)
  Dx2 = (1 / dx) .* D1xB(R, X, RX)
  Dx = (Dx2)
#Dxx2[1:R*(d-5),:] =  Dxx[1:R*(d-5),:]

perturb =zeros(R, X) .+0.00*cos.(LinRange(min_R,max_R,R))
perturb = reshape(perturb,RX)

psi = zeros(R, X)
g = zeros(R,X)
g = reshape(g, RX) .+ perturb
u0 = zeros(5 * RX)

g_ext1 =g_ext.(pos,linspace_X,σ,k,W,dx,type,beta,V_ss)*0
g_extin = repeat(g_ext1,inner = R)

psi0 = -g_extin .- g_c


u0[1:RX] .= α*psi0 - α*W*deltaF.*firingrate.(zeros(RX),type,beta)
u0[RX + 1:2*RX] .= psi0
u0[2 * RX + 1:3 * RX] .= g .+ perturb
u0[3 * RX + 1:4 * RX] .= g .+ perturb
u0[4 * RX + 1:5*RX] .= 0 .+ perturb
h_vec = collect(0.01:0.01:0.2)
k_vec = collect(0:0.002:0.2);
c_numeric_vector_k = zeros(size(k_vec,1))
counter = 1
tspan = (0.0, max_T)
params = [g_c, g_extin, type, β, σ, α, tau, k, v, W, h, beta, deltaF, Dxx, Dx, Drr]
prob = ODEProblem(rhsFun,u0, tspan, params)
print("Solving...")
sol = solve(prob,saveat = dt,progress = true)
print("\n Done!")

z_sol = reshape(sol[1:RX,:], R, X, size(sol, 2))
psi_sol = reshape(sol[RX + 1:2 * RX,:], R, X, size(sol, 2))
g_sol = reshape(sol[3 * RX + 1:4 * RX,:], R, X, size(sol, 2))
z2_sol =  reshape(sol[2 * RX + 1:3 * RX,:], R, X, size(sol, 2))
V_sol = reshape(sol[4 * RX + 1:5 * RX,:], R, X, size(sol, 2)) #voltage
#c_numeric_vector_k = trackSpeed(g_sol[:,d,:],Int(80/dr),Int(90/dr),0.1,dr,0.5,size(g_sol,3),0.1)
#npzwrite("/home/james/c_numeric_alfa.npy",c_numeric_vector_k)
@gif for i = 1:1:size(g_sol,3)
    heatmap(g_sol[:,:,i]')
end
#heatmap(V_sol[:,:,end])
#npzwrite("/tmp/Conductance_2D_kappa10.npy",g_sol)
#npzwrite("/tmp/Voltage_2D_kappa10.npy",V_sol)
