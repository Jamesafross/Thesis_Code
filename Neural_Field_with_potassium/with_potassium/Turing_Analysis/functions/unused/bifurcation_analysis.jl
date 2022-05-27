using Plots
using SparseArrays
using ProgressMeter
using DifferentialEquations
using NPZ
using LinearAlgebra
using NLsolve
include("TuringHopf_IFT_functions.jl")
include("hopf_IFT_functions.jl")

#options
delay = 1
Add_Potassium = 0
k_vec = LinRange(-5, 5, 20000)
# #PARAMETERS # #
ΔA = 0.5
ΔE = ΔA
ΔI = ΔA
τE = 1
τI = 20
VsynEE = 10
VsynEI = -5
VsynIE = 6
VsynII = -8

#g and ψ params
ss = 1
σEE = ss * 1
σEI = ss * 1.5
σIE = σEE
σII = ss * σEI

#conduction velocity
v =10

#synaptic time constant
αEE = 0.5
αEI = 0.5
αIE = 0.7
αII = 0.8

#synaptic connection strengths
κSc = 1
κSEE = κSc*1
κSEI = κSc*3.8
κSIE = κSc*1.8
κSII = 6.5

#Drive params
ηE = -5.3
ηI = 9

# gaps Strengths
κVc = .1
κVEE = κVc
κVEI = κVc
κVIE = κVc
τ0EE = 0
κVII = κVc
τ0EI = 0
τ0IE = 0
τ0II = 0

# # # # # # # # # # # # # # # # # # # # # # # # # # # #
#Potassium Things
f1(x1,x2,A1,A2,A3) = f("sigmoid",x1,x2,A1,A2,A3)
df1(x1,x2,A1,A2,A3) = dfdx("sigmoid",x1,x2,A1,A2,A3)
A1 = 1
A2 = 100
A3 = 0
# other params
δ = 1
A4 = 0.1
β1 = 1
# # # # # # # # # # # #
# Drive params # # # #
#Excitatory
f2(x2,B1,B2,B3) = f("tanh",1,x2,B1,B2,B3)
df2(x2,B1,B2,B3) = dfdx("tanh",1,x2,B1,B2,B3)
B1 = 10
B2 = 1
B3 = 0.1
β2 = 1
#Inhibitory
f3(x2,C1,C2,C3) = f("tanh",1,x2,C1,C2,C3)
df3(x2,C1,C2,C3) = dfdx("tanh",1,x2,C1,C2,C3)
C1 = 20
C2 = 1
C3 = 0.1
β3 = 1
# # # # # # # # # # # #
# gaps Strength params (D_II)
f4(x2,C1,C2,C3) = f("sigmoid",1,x2,D1,D2,D3)
df4(x2,D1,D2,D3) = dfdx("tanh",1,x2,D1,D2,D3)
D1 = 0.2
D2 = -10
D3= 0.1
β4 =  10
#make matrix

σ = [σEE σEI σIE σII]
κV = [κVEE κVEI κVIE κVII]
κS = [κSEE κSEI κSIE κSII]
η0 = [ηE, ηI]
τ0 = [τ0EE τ0EI τ0IE τ0II]
α = [αEE αEI αIE αII]
Vsyn = [VsynEE VsynEI VsynIE VsynII]
τ = [τE τI]
β = [β1 β2 β3 β4]

τ0EE,τ0EI,τ0IE,τ0II = [0;0;0;0]

#bif param
M = 1000
bifparam1 = LinRange(0.1,0.2,M)
b1_hopf = zeros(M)
b2_hopf = zeros(M)

init_hopf = [0.1,0.5]
for i =1:100
    v = bifparam1[i]
    b1_hopf[i],b2_hopf[i] = SteadyState_g_hopf_IFT(init_hopf)
    global init_hopf = [b1_hopf[i],b2_hopf[i]]

end

init_turinghopf = [0.34, -27.696384819240954, 2.7063642874415934]
b1_turinghopf = zeros(M)
b2_turinghopf = zeros(M)
b3_turinghopf = zeros(M)
for i =1:M
    v = bifparam1[i]
    print("\n v = ",v,"\n")
    if i == 1; iters = 1000;else iters = 1000; end
    b1_turinghopf[i],b2_turinghopf[i],b3_turinghopf[i] = SteadyState_turinghopf_IFT(init_turinghopf,iters)
    print("\n bifparam = ",b1_turinghopf[i],"\n")
    global init = [b1_turinghopf[i],b2_turinghopf[i],b3_turinghopf[i]]

end

plot(v,b1_hopf)
