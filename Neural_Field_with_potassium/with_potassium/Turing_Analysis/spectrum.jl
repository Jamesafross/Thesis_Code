using Plots
using SparseArrays
using DifferentialEquations
using NPZ
using LinearAlgebra
using NLsolve
include("functions/functions.jl")

#options
file_save = "1D_potassium/spectrum_1D_turinghopf"
save_data = "N" #(Y/N)
Add_Potassium = 1
dimensions = 1
gaps_option = "I-I" # ( I-I or ALL-ALL or OFF)
k_vec = LinRange(-10,10, 10000)

options = [Add_Potassium,gaps_option]

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
κSEE = κSc*5.0
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
δ = 0.09 #decay rate
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
τ0EE,τ0EI,τ0IE,τ0II = [0. 0. 0. 0.]

    p = ΔE, ΔI, κVEE,κVEI,κVIE,κVII,ηE,ηI,τE,τI,
    σEE,σEI,σIE, σII, v, κSEE,κSEI,κSIE, κSII,
    αEE, αEI, αIE, αII,τ0EE,τ0EI,τ0IE,τ0II, δ,
    A1, A2, A3, A4,B1, B2, B3, C1, C2, C3,
    D1, D2, D3,βK,βηE,βηI,βκV

# # # # # # # # # STEADY STATE(S) # # # # # # # #


ssMat = find_steady_state(p,Add_Potassium,gaps_option)

ssMat = ssMat[sortperm(ssMat[:,1]),:] #sort ssMat from smallest RE to largest RE
numSS = size(ssMat,1)
RE_ss = ssMat[:,1]
VE_ss = ssMat[:,2]
RI_ss = ssMat[:,3]
VI_ss = ssMat[:,4]

WE = pi*τE*RE_ss .+ im*VE_ss
WI = pi*τI*RI_ss .+ im*VI_ss
ZE = (1 .-conj.(WE)) ./(1 .+conj.(WE))
ZI = (1 .-conj.(WI)) ./(1 .+conj.(WI))
syncE = abs.(ZE)
syncI = abs.(ZI)


print("",numSS," STEADY STATE(S) FOUND!!!")
for i = 1:numSS
    if Add_Potassium == 0
        print("\n [",i,"]: RE_ss = ",RE_ss[i]," VE_ss = ",VE_ss[i],
         " RI_ss = ",RI_ss[i]," VI_ss = ",VI_ss[i]," \n    |ZE_ss| = ",syncE[i], " |ZI_ss| = ",syncI[i],"\n" )

    elseif Add_Potassium == 1
        K_ss = (1/δ)*fK.(RE_ss+RI_ss,A1,A2,A3)
        ηE_ss = fηE.(K_ss,B1,B2,B3)
        ηI_ss = fηI.(K_ss,C1,C2,C3)
        κVss = fκV.(K_ss,D1,D2,D3)
        print("\n [",i,"]: ",RE_ss[i]," ",VE_ss[i]," ",RI_ss[i]," ",VI_ss[i],
        " ",K_ss[i]," ",ηE_ss[i]," ",ηI_ss[i]," ",κVss[i])

    end
end


if numSS > 1
print("\n Choose Steady State to use: 1 - ", numSS,"\n")
chooseSS = Int(parse(Float64,readline()))
Steady_State = ssMat[chooseSS,:]
else
    Steady_State =ssMat[1,:]
end

RbarE, VbarE, RbarI, VbarI = Steady_State
V = [VbarE, VbarI]
R = [RbarE, RbarI]
if Add_Potassium == 1
    Kbar = (1/δ)*fK(RbarE+RbarI,A1,A2,A3)
    κbarV = fκV(Kbar,D1,D2,D3)
    ηbarE = fηE(Kbar,B1,B2,B3)
    ηbarI = fηI(Kbar,C1,C2,C3)
    d_fK_dRE = dfK(RbarE+RbarI,A1,A2,A3)
    d_fK_dRI = dfK(RbarE+RbarI,A1,A2,A3)
    d_fK_dVE = 0
    d_fK_dVI = 0
    d_fηE_dK = dfηE(Kbar,B1,B2,B3)
    d_fηI_dK = dfηI(Kbar,C1,C2,C3)
    d_fκV_dK = dfκV(Kbar,D1,D2,D3)
    derivatives  = d_fK_dRE,d_fK_dRI,d_fK_dVE,d_fK_dVI,d_fηE_dK,d_fηI_dK,
    d_fκV_dK 
else
    derivatives = 0

end
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


reL_vec = zeros(length(k_vec))
imL_vec = zeros(length(k_vec))

reL_vec, imL_vec = makeSpectrum(reL_vec,imL_vec,k_vec,[0.5,0.5])

k_max = k_vec[findmax(reL_vec)[2]]
ω_max = imL_vec[findmax(reL_vec)[2]]

#print("\n dif = ",findmax(reL_vec)[1] - findmin(reL_vec)[1],"\n")
if maximum(reL_vec) > 0.0
    stability = "unstable"
else
    stability = "stable"
end
print("\n maximum μ at |k| = ", abs(k_max)," (",stability,")\n")

k_vec_save = cat(k_vec,k_vec,dims=1)
reL_vec_save = cat(reL_vec,reL_vec,dims=1)
imL_vec_save = cat(imL_vec,-imL_vec,dims=1)
spectrum_data = [k_vec_save reL_vec_save imL_vec_save]


if save_data == "Y" 
npzwrite("/home/james/PhD_Work/Python_Code/ExtraCellularPotassium/data/$file_save", spectrum_data)
end

plot(plot(k_vec, [zeros(length(reL_vec)), reL_vec]), plot(imL_vec))






#plot(scatter((reL+imL)[1:50:Int(length(reL)/2)],color = :bluesreds),plot(reL_vec))
