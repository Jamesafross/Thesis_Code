using Plots
using LinearAlgebra
using NLsolve
include("functions.jl")

#options
k_vec = LinRange(0,10,20000)
save_data = 0
# params
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
γ = 0.31 #good bif param!
γ_E = γ
γ_I = γ
v = 1
g_EXT = 50


V_plus = 1
V_neg = -V_plus
kc = 1



p = [α_EE,γ_E,σ_EE,κ_EE,τ_E,v,D_E,d_EE,W0_EE,
     α_EI,γ_I,σ_EI,κ_EI,v,D_E,d_EI,W0_EI,
     α_IE,γ_E,σ_IE,κ_IE,v,D_I,d_IE,W0_IE,
     α_II,γ_I,σ_II,κ_II,τ_I,v,D_I,d_II,W0_II ]



reλprevious =-0.01
imλprevious = 0.27


for j = 1:length(k_vec)
    #print("\n",j,"\n")
        global k = k_vec[j]
        z = nlsolve(f!, [ reλprevious; imλprevious])
        while z.f_converged == false
            z = nlsolve(f!, [ reλprevious + randn(); imλprevious+randn()])
        end
        reλ = z.zero[1]
        imλ = z.zero[2]
        global reλprevious = reλ
        global imλprevious = imλ
        #print("\n k = ",k,"zeros =",zeros,"\n")
        if j == 1
            global s1 = [reλ imλ  k]
        else
            global s1 = vcat(s1,[reλ imλ  k])
        end

end

reλprevious =-0.01
imλprevious = -0.27


for j = 1:length(k_vec)
    #print("\n",j,"\n")
        global k = -k_vec[j]
        z = nlsolve(f!, [ reλprevious; imλprevious])
        while z.f_converged == false
            z = nlsolve(f!, [ reλprevious - 0.01*randn(); imλprevious-0.01*randn()])
        end
        reλ = z.zero[1]
        imλ = z.zero[2]
        global reλprevious = reλ
        global imλprevious = imλ
        #print("\n k = ",k,"zeros =",zeros,"\n")
        if j == 1
            global s2 = [reλ imλ  k]
        else
            global s2 = vcat(s2,[reλ imλ  k])
        end

end





p1 = plot([real(s1[1:end,3]) real(s2[1:end,3])],[[zeros(size(s2,1)),zeros(size(s2,1))], [real(s1[1:end,1]) real(s2[1:end,1])]],lab="Re(λ_2)")
p2 = plot(real(s1[1:end,3]),[real(s1[1:end,2]),real(s2[1:end,2])],lab="im(λ_2)")
k_save = cat(reverse(s2[1:end,3]),s1[1:end,3],dims=1)
reL = cat(reverse(s2[1:end,1]),s1[1:end,1],dims=1)
imL = cat(reverse(-s1[1:end,2]),s1[1:end,2],dims=1)
data = [k_save reL imL]

k_max = k_save[findmax(reL)[2]]
ω_max = imL[findmax(reL)[2]]




if save_data == 1
  #  npzwrite("/home/james/PhD_Work/Python_Code/Axo_Dendritic_Paper/2_Pop_plots/data/turing_hopf_spectrum.npy", data)
end
plot(p1,p2)

#scatter(reL + im*imL)
