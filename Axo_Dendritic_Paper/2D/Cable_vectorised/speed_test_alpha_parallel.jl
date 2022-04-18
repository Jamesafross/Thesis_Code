using ProgressMeter
using Plots
using SparseArrays
using LinearAlgebra
using DifferentialEquations
using NPZ
using
#Load functions
include("firingrate.jl")
include("Banded_Matrices.jl")
include("rhsFuncs2.jl")
include("speedtracker.jl")
include("SpeedEquation.jl")
include("Guassian_Delta.jl")
const max_T = 150
#mesh paramaters
const min_R = 0;
const max_R = 800
const min_X = -0.5
const max_X = 0.5
const dt = 0.1
const dx = 0.1;
const dr = 0.2;
#sqrt(D*tau) \approx 0.1 - 1 mm (D*tau)
const W = 3.0
const tau = 1.0  #milliseconds #x10
#finds the grid number for origin on the cable
#PARAMETERS
const v = 100 #mm/s  #connectivity speed
const sigma = 10.0;   #x10^-1 mm #
#const h = 0.15; #threshold
const α = 1 / sigma;
const β = 1 / v
const τ = 1   #time constant for cable equation
const D = (0.1)^2 #diffusion coefficient cable equation
const k = 0.005#d #milliseconds #x10^-3 seconds
const type = 0#type of firing rate function (0 for heaviside, 1 for sigmoid)
const beta = 45 #sigmoid steepness (if using...)
const pos = 0.2 #position on cable (multiples of dx)
const h = 0.1
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
const d = Int(round.(xzero + pos / dx))#this is the actual grid number for the position
#deltaF = zeros(RX)
#deltaF[(d - 1) * R + 1:(d) * R] .= 1/dx
deltaF = repeat(smooth_Delta(linspace_X .- (pos),dx),inner = R)
const Drr = (1 / drdr) .* D2z(R, X, RX)
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
g[st:en,xzero] .= 1
g = reshape(g, RX)
u0 = zeros(5 * RX)

#u0[RX + 1:2*RX] .= g
u0[2 * RX + 1:3 * RX] .= g
u0[3 * RX + 1:4 * RX] .= g
u0[4 * RX + 1:5*RX] .= g./maximum(g)
alfa_vec = collect(1:0.2:10)

c_numeric_vector_k = zeros(size(alfa_vec,1))

for i = 1:size(alfa_vec,1)
    global u0
    global alfa_vec
    global alfa = alfa_vec[i]

    print("α = ",alfa, "\n")
tspan = (0.0, max_T)
params = [type,β, α, tau, k, v,W, h, beta, deltaF, Dxx, Dx, Drr]
prob = ODEProblem(rhsFun,u0, tspan, params)
print("Solving...")
sol = solve(prob,RK4(),saveat = dt,progress = true)
print("\n Done! \n")
g_sol =  reshape(sol[3 * RX + 1:4 * RX,:], R, X, size(sol, 2))
c_numeric_vector_k[i] = trackSpeed(g_sol[:,d,:],Int(740/dr),Int(750/dr),dt,dr,h,size(g_sol,3),0.1)
print("c = ",c_numeric_vector_k[i],"\n")
end

npzwrite("/tmp/c_numeric_alfa.npy",c_numeric_vector_k)
