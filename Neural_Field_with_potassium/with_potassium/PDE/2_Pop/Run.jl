using Plots
using SparseArrays
using DifferentialEquations
using NPZ
using NLsolve
using JLD

include("functions/banded_matrices.jl")
include("functions/RHS_funcs.jl")
include("functions/miscfunctions.jl")

#moreK : A1 = 2.0
#midK : A1 = 1.0
#lessK : A1 = 0.1

#lowdecay = 0.09
#highdecay = 0.8

#lowdiffusion = 0.01
#highdiffusion = 1.5

DIRSAVE = "/1D_potassium"
file_name = "potassium_1D_highdiffusion"
function savedata(save_data)
    if save_data == 1
        print("SAVING")
        npzwrite("/home/james/PhD_Work/Python_Code/ExtraCellularPotassium/data$DIRSAVE/ZE_$file_name.npy", ZE)
        npzwrite("/home/james/PhD_Work/Python_Code/ExtraCellularPotassium/data$DIRSAVE/ZI_$file_name.npy", ZI)
        npzwrite("/home/james/PhD_Work/Python_Code/ExtraCellularPotassium/data$DIRSAVE/K_$file_name.npy", K)

    end
end

function bump(X,xmid,size)
    bump = zeros(X)
end


use_rand_init_conds = 1 # custom initial conditions if not perturbing steady state
attached_steady_state = 0
restart_solve = 0
load_data = 0
save_data = 1

# space nd tim parameters
X_max = 2.0pi;
dx = pi/2^6
T_max =6000.0
dxdx = dx * dx;
saveat = collect(T_max-100:0.1:T_max)


#options
dimension = 1  # number of spatial dimensions (1 or 2)

Add_Potassium = 1 #option to add or remove potassium dynamics (0 or 1)
gaps_option = "I-I"  #("ALL-ALL" or "I-I" or "OFF")
options = [Add_Potassium, gaps_option, dimension]


kc= 1.1#wave length of perturbations
amp_perturb = 0.001
run_sim = 1

#perturb = .0*exp.(-(X_space.^2)/10)
if dimension == 1
    X = Int(2*X_max / dx)
    X_space = LinRange(-X_max,X_max,X)
    ∇ = (1.0 / dxdx) * D2x(X)
    ∇ψ = ∇
    perturb = amp_perturb*cos.(kc*X_space)
    
elseif dimension == 2
    X1 = Int(2*X_max / dx)
    X = X1*X1
    D2x1 = (1 / dxdx) * D2r1(X1,X)
    D2x2 = (1 / dxdx) * D2r2(X1,X)
    ∇ = D2x1 + D2x2
    ∇ψ = (3.0/2.0)*(D2x1 .+ D2x2)
    X1_space = LinRange(-X_max,X_max,X1)
    kc = kc/(sqrt(2.0))
    X2_space = LinRange(-X_max,X_max,X1)'
    #perturb = reshape(amp_perturb*cos.(kc*X1_space)*cos.(kc*X2_space),X)
    perturb = reshape(real(amp_perturb*init_conds(3,kc,[0,2pi/3,4pi/3],[1,1,1],X1_space,X2_space)),X)
    
end


function init_conditions(u0,X,Add_Potassium)
    u0[1:20*X] .+= 0.1randn(20X)
    u0[1:X] .+= abs.(u0[1:X])
    u0[2X+1:3X] .+= abs.(u0[2X+1:3X])
    if Add_Potassium == 1
        #reshape(twoDimGuassian(2pi,2pi,.1,X1_space,5),X)
        u0[20*X + 1:21*X] .+= 0.1randn(X)
        u0[21*X + 1:22*X] .+= 0.1randn(X)
        u0[22*X + 1:23*X] .+= 0.1randn(X)
        u0[23*X + 1:24*X] .+= 0.1rand(X)
    end
return 5u0
end

# # # # MODEL PARAMETERS # # # #

# # # # MODEL PARAMETERS # # # #
ΔA = 0.20
ΔE = ΔA
ΔI = ΔA
τE = 9.0
τI = 11.0
VsynEE = 12.0
VsynEI = VsynEE*(-10.0/16.0)
VsynIE = VsynEE*(10.0/16.0)
VsynII = VsynEE*(-11.0/16.0)

#g and ψ params
ss = 1
σEE = ss * 0.2
σEI = ss * 1.5
σIE = 0.2
σII = 1.5

#conduction velocity
v =.08

#synaptic time constant
αEE = 0.7
αEI = 0.5
αIE = 0.7
αII = 0.8

#synaptic connection strengths
κSc = 1.0
κSEE = κSc*3.5
κSEI = κSc*κSEE*(4.9/5.0)
κSIE = κSc*κSEE*(4.9/5.0)
κSII =  κSc*κSEE*(3.0/5.0)

#Drive params
ηE = 1
ηI = 1
# gaps Strengths
κV = .0
κVEE = 0.0
κVEI = 0.0
κVIE = κVEI
κVII = 0.0
τ0EE = 0
τ0EI = 0
τ0IE = 0
τ0II = 0

# Potassium params # # #
# sigmoid params
A1 = 1.0
A2 = 2.0 #steepness
A3 = 0.0
βK = 1.0 #time constant
# other potassiumn params
δ = 0.2 #decay rate
A4 = 1.5 #diffusion rate
fK(x,A1,A2,A3) = x*f("sigmoid",x,A1,A2,A3)
dfK(x,A1,A2,A3) = f("sigmoid",x,A1,A2,A3) + x*dfdx("sigmoid",x,A1,A2,A3)
# # # # # # # # # # # #
# Drive params # # # #
#Excitatory
B1 =15.0
B2 = 8.0# steepness
B3 = 1.0 #threshold
B4 = 0.0
βηE = 25 #time constant
fηE(x,B1,B2,B3) = f("sigmoid",x,B1,B2,B3)
dfηE(x,B1,B2,B3) = dfdx("sigmoid",x,B1,B2,B3)
#Inhibitory
C1 = 15.0
C2 = 8.0 # steepness
C3 = 1.0 # threshold
C4 = 0.0
βηI = 25 #time constant

fηI(x,C1,C2,C3) = f("sigmoid",x,C1,C2,C3)
dfηI(x,C1,C2,C3) = dfdx("sigmoid",x,C1,C2,C3)
# # # # # # # # # # # #
# gaps Strength params
D1 = 0.5
D2 = -2.0 # steepness
D3 = 1.0
βκV = 4.0 #time constant
fκV(x,D1,D2,D3) = f("sigmoid",x,D1,D2,D3)
dfκV(x,D1,D2,D3) = dfdx("sigmoid",x,D1,D2,D3)
params = [ΔE, ΔI, κVEE,κVEI, κVIE,κVII,ηE,ηI,τE,τI,
          σEE, σEI,σIE, σII, v, κSEE,κSEI,κSIE, κSII,
          αEE, αEI, αIE, αII, δ,
          A1, A2, A3, B1, B2, B3, C1, C2, C3,
          D1, D2, D3, βK, βηE, βηI,βκV,∇,∇ψ,options]


u0,ssMat = setup(params,use_rand_init_conds,restart_solve,X)


if attached_steady_state == 1 && size(ssMat,1) > 1
    if dimension == 1
        attach_ss!(ssMat,u0,2,3,Add_Potassium,X,X_space,dimension)
    elseif dimension == 2
        attach_ss!(ssMat,u0,2,3,Add_Potassium,X,X1_space,dimension)
    end
end



# # # SOLVE # # # #
tspan = (0.0, T_max) # TIME SPAN
# Set up ODE problem
println("setting up Problem")
prob = ODEProblem(rhsFun, u0, tspan, params)

println(" Solving...")
#SOLVER
sol = solve(prob,BS3(),maxiters=10e20,saveat = saveat, progress = true)
print("\n Done!")

    RE = sol[1:X,:]
    VE = sol[X + 1:2*X,:]
    RI = sol[2*X+1:3*X,:]
    VI = sol[3*X + 1:4*X,:]
    ψEE = sol[4*X+1:5*X,:]
    ψEI = sol[6*X+1:7*X,:]
    ψIE = sol[8*X+1:9*X,:]
    ψII = sol[10*X+1:11*X,:]
    gEE  = sol[12*X+1:13*X,:]
    gEI  = sol[14*X+1:15*X,:]
    gIE  = sol[16*X+1:17*X,:]
    gII  = sol[18*X+1:19*X,:]
    if Add_Potassium == 1
        K = sol[20*X + 1:21*X,:]
        ηEsol = sol[21*X + 1:22*X,:]
        ηIsol = sol[22*X + 1:23*X,:]
        κV = sol[23*X + 1:24*X,:]
    end



WE = τE*pi*RE + im*VE
ZE = (1.0 .-conj.(WE)) ./(1.0 .+conj.(WE))
syncE = abs.(ZE)

WI = τI*pi*RI + im*VI
ZI = (1.0 .-conj.(WI)) ./(1.0 .+conj.(WI))
syncI = abs.(ZI)

t = LinRange(0,T_max,size(sol,2))

M = size(sol,2)



if dimension == 2
    RE = reshape(RE,X1,X1,size(sol,2))
    RI = reshape(RI,X1,X1,size(sol,2))
    ZE = reshape(ZE,X1,X1,size(sol,2))
    ZI = reshape(ZI,X1,X1,size(sol,2))
    K = reshape(K,X1,X1,size(sol,2))
end


savedata(save_data)
if dimension == 1
    if Add_Potassium == 1
        plot(heatmap(syncE),heatmap(syncI),heatmap(ηEsol),heatmap(ηIsol),heatmap(K),heatmap(κV) )
    else
        plot( heatmap(syncE),heatmap(syncI) )
    end
elseif dimension == 2
    N = size(sol,2)
    levels = LinRange(minimum(syncE[:,1:N]),maximum(syncE[:,1:N]),200)
    @gif for i = 1:4:N; contourf(reshape(syncE[:,i],X1,X1),line_smoothing=0.85,linewidth=0);end
end


#@gif for i = 1:2000; plot(plot(syncE[:,i],ylims=[minimum(syncE[:,i]),maximum(syncE[:,i])]),plot(syncI[:,i],ylims=[minimum(syncI[:,i]),maximum(syncI[:,i])]));end


 