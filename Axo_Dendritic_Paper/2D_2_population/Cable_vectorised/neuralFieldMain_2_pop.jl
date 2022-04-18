using ProgressMeter
using Plots
using SparseArrays
using LinearAlgebra
using DifferentialEquations
using NPZ
#Load functions
include("firingrate.jl")
include("Banded_Matrices.jl")
include("rhsFuncs2_noshunts.jl")
include("Guassian_Delta.jl")
kc=2.49

#options
save_data = 1
max_T = 1500
#mesh paramaters
  min_R = -4*pi;
  max_R = 4*pi
  min_X = -1.6
  max_X = 1.6
  nX = 20
  dt = 0.5
  dx = 0.05
  dr = pi/2^5;


#sqrt(D*tau) \approx 0.1 - 1 mm (D*tau)
#PARAMETERS
α_EE = 0.5
α_EI = 0.5
α_IE = 0.7
α_II = 0.2

sigma_I = 1
sigma_E = sigma_I
σ_EE = sigma_E*1.5  #x10^-1 mm #
σ_EI = sigma_I*0.2   #x10^-1 mm #
σ_IE = sigma_E*1.5   #x10^-1 mm #
σ_II = sigma_I*0.2 #x10^-1 mm #

τ_E = 2

τ_I =2

κ_EE = .01
κ_EI = .01
κ_IE = .01
κ_II = .01

d = 0.1
d_EE = d
d_EI = d
d_IE = d
d_II = d

Wc = 1
W0_EE =Wc*1.0
W0_EI = Wc*2.2
W0_IE = Wc*2.0
W0_II = Wc*2.7

D_E = 0.05
D_I = 0.05


h_E= 0.0
h_I= 0.0
γ = 0.4 #good bif param!
γ_E = γ
γ_I = γ
v = .1
g_EXT = 50

#k_vec = LinRange(0,5,100000)

V_plus = 1
V_neg = -V_plus
  type = "tanh"
  beta = γ
#GRID DISCRETISATION
#R (r) - soma space
#X (x) - cable space
v_EE = v
v_EI = v
v_IE = v
v_II = v
  drdr = dr * dr
  dxdx = dx * dx
  R = Int(round((max_R-min_R) / dr))

  X = Int(round((max_X - min_X) / dx)) + 1
  RX = R * X
  linspace_X = collect(min_X:dx:max_X) #linear grid for X
  xzero = Int(findfirst(linspace_X .== 0)) #this is the grid number for x=0
  pos_EE = Int(round.(xzero + d_EE / dx))#this is the actual grid number for the position
  pos_EI = Int(round.(xzero + d_EI / dx))#this is the actual grid number for the position
  pos_IE = Int(round.(xzero + d_IE / dx))#this is the actual grid number for the position
  pos_II = Int(round.(xzero + d_II / dx))#this is the actual grid number for the position
#deltaF = zeros(RX)
#deltaF[(d - 1) * R + 1:(d) * R] .= 1/dx
deltaF_EE = repeat(smooth_Delta(linspace_X .- (d_EE),dx),inner = R)
deltaF_EI = repeat(smooth_Delta(linspace_X .- (d_EI),dx),inner = R)
deltaF_IE = repeat(smooth_Delta(linspace_X .- (d_IE),dx),inner = R)
deltaF_II = repeat(smooth_Delta(linspace_X .- (d_II),dx),inner = R)
deltaF_zero = repeat(smooth_Delta(linspace_X,dx),inner = R)

function V_ss_func(tau,g1,g2,V_p,V_n)
  return (V_p*g1 + V_n*g2)/((1/tau) + g1 +g2)
end
V_ssE = V_ss_func(τ_E,g_EXT,g_EXT,V_plus,V_neg)
V_ssI = V_ss_func(τ_I,g_EXT,g_EXT,V_plus,V_neg)



hE = V_ssE;
hI = V_ssI;


  Drr = (1 / drdr) .* D2z(R, X, RX)
  Dxx = (1 / dxdx) .* D2xCent(R, X, RX,xzero)
#Dxx22 = (1 / dmakegif(gEE_sol)xdx) .*D2xF(R,X,RX)
  Dx1 = (1 / dx) .* D1xF(R, X, RX)
  Dx2 = (1 / dx) .* D1xB(R, X, RX)
  Dx = (Dx1)
#Dxx2[1:R*(d-5),:] =  Dxx[1:R*(d-5),:]


#perturbh = reshape(zeros(R, X)  .+ 1*exp.(-(LinRange(min_R,max_R,R)).^2),RX)
perturbh = reshape(zeros(R, X)  .+ .1*cos.(kc*LinRange(min_R,max_R,R)),RX)
perturb_g_EE = perturbh
perturb_g_EI = perturbh
perturb_g_IE = perturbh
perturb_g_II = perturbh
g = zeros(R, X) .+ g_EXT
gEE = reshape(g, RX) .+ perturb_g_EE
gEI = reshape(g, RX) .+ perturb_g_EI
gIE = reshape(g, RX) .+ perturb_g_IE
gII = reshape(g, RX) .+ perturb_g_II
hE = V_ssE .+ firingrate.(deltaF_zero,"heaviside",maximum(deltaF_zero)-0.1).*perturbh;
hI = V_ssI .+ firingrate.(deltaF_zero,"heaviside",maximum(deltaF_zero)-0.1).*perturbh;
u0 = zeros(18 * RX)



psiEE0 = 0
psiEI0 = 0
psiIE0 = 0
psiII0 = 0

u0[1:RX] .= (psiEE0)
u0[RX + 1:2 *RX] .=  (psiEI0)
u0[2 * RX + 1:3 * RX] .=  (psiIE0)
u0[3 * RX + 1:4 * RX] .= (psiII0)
u0[4 * RX + 1:5 * RX] .= psiEE0
u0[5 * RX + 1:6 * RX] .= psiEI0
u0[6 * RX + 1:7 * RX] .= psiIE0
u0[7 * RX + 1:8 * RX] .= psiII0
u0[8 * RX + 1:9 * RX] .= gEE
u0[9 * RX + 1:10 * RX] .= gEI
u0[10 * RX + 1:11 * RX] .= gIE
u0[11 * RX + 1:12 * RX] .= gII
u0[12 * RX + 1:13 * RX] .= gEE
u0[13 * RX + 1:14 * RX] .= gEI
u0[14 * RX + 1:15 * RX] .= gIE
u0[15 * RX + 1:16 * RX] .= gII
u0[16 * RX + 1:17 * RX] .= hE
u0[17 * RX + 1:18 * RX] .= hI
h_vec = collect(0.01:0.01:0.2)
k_vec = collect(0:0.002:0.2);
c_numeric_vector_k = zeros(size(k_vec,1))
counter = 1
tspan = (0.0, max_T)
params = [g_EXT,type,α_EE,α_EI,α_IE,α_II, τ_E,τ_I,κ_EE,κ_EI,κ_IE,κ_II, v_EE,v_EI, v_IE, v_II,W0_EE,W0_EI,W0_IE,W0_II, hE,hI,beta,d_EE,d_EI,d_IE,d_II, deltaF_EE,deltaF_EI,deltaF_IE,deltaF_II, Dxx, Dx, Drr,linspace_X]
prob = ODEProblem(rhsFun_2pop,u0, tspan, params)
print("\n Solving... \n")
sol = solve(prob,RK4(),saveat=dt,progress = true,reltol=1e-8,abstol=1e-8)
print("\n Done! \n")

#z_sol = reshape(sol[1:RX,:], R, X, size(sol, 2))
#psi_sol = reshape(sol[RX + 1:2 * RX,:], R, X, size(sol, 2))

AEE_sol = reshape(sol[1:1 * RX,:], R, X, size(sol, 2))
AEI_sol = reshape(sol[1 * RX + 1:2 * RX,:], R, X, size(sol, 2))
AIE_sol = reshape(sol[2 * RX + 1:3 * RX,:], R, X, size(sol, 2))
AII_sol = reshape(sol[3 * RX + 1:4 * RX,:], R, X, size(sol, 2))

psiEE_sol = reshape(sol[4 * RX + 1:5 * RX,:], R, X, size(sol, 2))
psiEI_sol = reshape(sol[5 * RX + 1:6 * RX,:], R, X, size(sol, 2))
psiIE_sol = reshape(sol[6 * RX + 1:7 * RX,:], R, X, size(sol, 2))
psiII_sol = reshape(sol[7 * RX + 1:8 * RX,:], R, X, size(sol, 2))
YEE_sol = reshape(sol[8 * RX + 1:9 * RX,:], R, X, size(sol, 2))
YEI_sol = reshape(sol[9 * RX + 1:10 * RX,:], R, X, size(sol, 2))
YIE_sol = reshape(sol[10 * RX + 1:11 * RX,:], R, X, size(sol, 2))
YII_sol = reshape(sol[11 * RX + 1:12 * RX,:], R, X, size(sol, 2))
gEE_sol = reshape(sol[12 * RX + 1:13 * RX,:], R, X, size(sol, 2))
gEI_sol = reshape(sol[13 * RX + 1:14 * RX,:], R, X, size(sol, 2))
gIE_sol = reshape(sol[14 * RX + 1:15 * RX,:], R, X, size(sol, 2))
gII_sol = reshape(sol[15 * RX + 1:16 * RX,:], R, X, size(sol, 2))

VE_sol = reshape(sol[16 * RX + 1:17 * RX,:], R, X, size(sol, 2))
VI_sol = reshape(sol[17 * RX + 1:18 * RX,:], R, X, size(sol, 2)) #voltage

function makegif(sol)
  @gif for i = 1:1:size(sol,3)
    heatmap(sol[:,:,i]')
  end
end
p1 = heatmap(VE_sol[:,pos_EE-1,1:end])
p2 = heatmap(VI_sol[:,pos_EE,1:end])
timespan = LinRange(0,max_T,size(VE_sol,3))

#plot([gII_sol[20,pos_EE,1:end],gEI_sol[20,pos_EE,1:end],gIE_sol[20,pos_EE,1:end],gEE_sol[20,pos_EE,1:end]])
plot(p1,p2)
if save_data == 1
    npzwrite("/home/james/PhD_Work/Python_Code/Axo_Dendritic_Paper/2_Pop_plots/data/turing_hopf_sim.npy", VE_sol[:,:,end-500:end])
end
makegif(VE_sol[:,:,1:end])
