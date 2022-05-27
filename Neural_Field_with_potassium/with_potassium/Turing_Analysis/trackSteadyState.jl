using Plots
using SparseArrays
using DifferentialEquations
using NPZ
using LinearAlgebra
using NLsolve
using Parameters
include("functions/functions.jl")

@with_kw struct SteadyStates
    S 
    A
end

@with_kw struct SteadyStates2
    S 
    A
    Stability
    eigs
end




#options
file_save = "1D_potassium/spectrum_1D_turinghopf"
save_data = "N" #(Y/N)
Add_Potassium = 1
dimensions = 1
gaps_option = "ALL-ALL" # ( I-I or ALL-ALL or OFF)
k_vec = LinRange(-10,10, 5000)
shunts = 1

options = [Add_Potassium,gaps_option]

# # # # MODEL PARAMETERS # # # #
ΔA = 0.20
ΔE = ΔA
ΔI = ΔA
τE = 9.0
τI = 11.0
VsynEE = 11.0
VsynEI = VsynEE*(-10.0/16.0)
VsynIE = VsynEE*(10.0/16.0)
VsynII = VsynEE*(-11.0/16.0)

#g and ψ params
ss = 1
σEE = ss * 0.2
σEI = ss * 1.6
σIE = 0.2
σII = 1.6

#conduction velocity
v =.08

#synaptic time constant
αEE = 0.7
αEI = 0.5
αIE = 0.7
αII = 0.8

#synaptic connection strengths
κSc = 1.0
κSEE = κSc*2.0
κSEI = κSc*κSEE*(4.6/5.0)
κSIE = κSc*κSEE*(4.6/5.0)
κSII =  κSc*κSEE*(3.1/5.0)

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
A1 = 2.0
A2 = 2.0 #steepness
A3 = 0.0
βK = 1.0 #time constant
# other potassiumn params
δ = 0.2 #decay rate
A4 = 0.3 #diffusion rate
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

#make matrix
Nsteps = 2000
rEss = zeros(Nsteps,3)
vEss = zeros(Nsteps,3)
rIss = zeros(Nsteps,3)
vIss = zeros(Nsteps,3)
A1_vec = LinRange(3.2,.1,Nsteps)
for j = 1:3
    A1 = 3.2
    p = ΔE, ΔI, κVEE,κVEI,κVIE,κVII,ηE,ηI,τE,τI,
    σEE,σEI,σIE, σII, v, κSEE,κSEI,κSIE, κSII,
    αEE, αEI, αIE, αII,τ0EE,τ0EI,τ0IE,τ0II, δ,
    A1, A2, A3, A4,B1, B2, B3, C1, C2, C3,
    D1, D2, D3,βK,βηE,βηI,βκV
    ssMat1 = find_steady_state(p,Add_Potassium,gaps_option)
    while size(ssMat1,1) < 3
        ssMat1 = find_steady_state(p,Add_Potassium,gaps_option)
    end
    ssMat1 = ssMat1[sortperm(ssMat1[:,1]),:]
    ssINIT = ssMat1[j,:]
    τ0EE,τ0EI,τ0IE,τ0II = [0. 0. 0. 0.]
   
    
    ssMat = []
    for i = 1:Nsteps
        A1 = A1_vec[i]
        #println(A1)
        p = ΔE, ΔI, κVEE,κVEI,κVIE,κVII,ηE,ηI,τE,τI,
        σEE,σEI,σIE, σII, v, κSEE,κSEI,κSIE, κSII,
        αEE, αEI, αIE, αII,τ0EE,τ0EI,τ0IE,τ0II, δ,
        A1, A2, A3, A4,B1, B2, B3, C1, C2, C3,
        D1, D2, D3,βK,βηE,βηI,βκV
        #println(i)
            if i == 1 
                ssMat = cat(ssMat,SteadyStates(S = find_steady_state2(p,Add_Potassium,gaps_option,ssINIT), A = A1),dims=1)
            else
                ssMat = cat(ssMat,SteadyStates(S = find_steady_state2(p,Add_Potassium,gaps_option,ssMat[i-1].S), A = A1),dims=1)
            end
    end

    
    for i = 1:Nsteps
        rEss[i,j] = ssMat[i].S[1]
        vEss[i,j] = ssMat[i].S[2]
        rIss[i,j] = ssMat[i].S[3]
        vIss[i,j] = ssMat[i].S[4]
    end
end

rEssN = zeros(Nsteps,3)
rEssN .= NaN
vEssN = zeros(Nsteps,3)
vEssN .= NaN
rIssN = zeros(Nsteps,3)
rIssN .= NaN
vIssN = zeros(Nsteps,3)
vIssN .= NaN
for j = 1:3
    #print(j)
    for i = 1:Nsteps
        bc = false
        if i == 1
            rEssN[i,j] = rEss[i,j]
            vEssN[i,j] = vEss[i,j]
            rIssN[i,j] = rIss[i,j]
            vIssN[i,j] = vIss[i,j]

        else
            if abs.(rEss[i,j] - rEss[i-1,j]) > 0.01
                bc = true
                break
            else
                rEssN[i,j] = rEss[i,j]
                vEssN[i,j] = vEss[i,j]
                rIssN[i,j] = rIss[i,j]
                vIssN[i,j] = vIss[i,j]
    
            end
        end
        if bc == true
            break
        end
    end
end

SSMAT = Array{SteadyStates2}(undef,Nsteps,3)
for j = 1:3
    println(j)
    for i = 1:Nsteps
        A1 = A1_vec[i]
        p = ΔE, ΔI, κVEE,κVEI,κVIE,κVII,ηE,ηI,τE,τI,
        σEE,σEI,σIE, σII, v, κSEE,κSEI,κSIE, κSII,
        αEE, αEI, αIE, αII,τ0EE,τ0EI,τ0IE,τ0II, δ,
        A1, A2, A3, A4,B1, B2, B3, C1, C2, C3,
        D1, D2, D3,βK,βηE,βηI,βκV
        
    RbarE, VbarE, RbarI, VbarI = rEssN[i,j],vEssN[i,j],rIssN[i,j],vIssN[i,j]
     V = [VbarE, VbarI]
    R = [RbarE, RbarI]
    if Add_Potassium == 1
        Kbar = (1/δ)*fK(RbarE+RbarI,A1,A2,A3)
        κbarV = fκV(Kbar,D1,D2,D3)
        ηbarE = fηE(Kbar,B1,B2,B3)
         ηbarI = fηI(Kbar,C1,C2,C3)
        d_fK_dRE = dfK(RbarE+RbarI,A1,A2,A3)
        d_fK_dRI = dfK(RbarE+RbarI,A1,A2,A3)
        d_fK_dVE = 0.
        d_fK_dVI = 0.
        d_fηE_dK = dfηE(Kbar,B1,B2,B3)
        d_fηI_dK = dfηI(Kbar,C1,C2,C3)
        d_fκV_dK = dfκV(Kbar,D1,D2,D3)
         derivatives  = d_fK_dRE,d_fK_dRI,d_fK_dVE,d_fK_dVI,d_fηE_dK,d_fηI_dK,
        d_fκV_dK 
    else
        derivatives = 0

    end
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    k_vec = LinRange(0,5,500)
    reL_vec = zeros(length(k_vec))
    imL_vec = zeros(length(k_vec))

    eigs = find_all_eigs(V,R,κbarV,k_vec,p,options,derivatives)

        if maximum(real(eigs)) > 0.001
            stability = "unstable"
        else
            stability = "stable"
        end
        SSMAT[i,j] = SteadyStates2(S = [RbarE, VbarE, RbarI, VbarI],A = A1_vec[i], Stability = stability,eigs = eigs)

    end

   

end

(1/δ)*fK(RbarE+RbarI,A1,A2,A3)

firstBranch = zeros(Nsteps,3)
secondBranch = zeros(Nsteps,3)
thirdBranch = zeros(Nsteps,3)
for i = 1:Nsteps
    
firstBranch[i,1] = (1/δ)*fK(SSMAT[i,1].S[1]+SSMAT[i,1].S[3],SSMAT[i,1].A,A2,A3)
firstBranch[i,2] = SSMAT[i,1].A
secondBranch[i,1] = (1/δ)*fK(SSMAT[i,2].S[1]+SSMAT[i,2].S[3],SSMAT[i,2].A,A2,A3)
secondBranch[i,2] = SSMAT[i,2].A
thirdBranch[i,1] = (1/δ)*fK(SSMAT[i,3].S[1]+SSMAT[i,3].S[3],SSMAT[i,3].A,A2,A3)
thirdBranch[i,2] = SSMAT[i,3].A
    if SSMAT[i,1].Stability == "unstable"
        firstBranch[i,3] = 1
    else
        firstBranch[i,3] = 0
    end
    if SSMAT[i,2].Stability == "unstable"
        secondBranch[i,3] = 1
    else
        secondBranch[i,3] = 0
    end
    if SSMAT[i,3].Stability == "unstable"
        thirdBranch[i,3] = 1
    else
        thirdBranch[i,3] = 0
    end
end