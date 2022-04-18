using Plots
using LinearAlgebra
using NLsolve
include("functions.jl")
α_EE = 0.5
α_EI = 0.5
α_IE = 0.7
α_II = 0.2

sigma_I = 1
sigma_E = 1
σ_EE = sigma_E*2.0;   #x10^-1 mm #
σ_EI = sigma_I*1.0;   #x10^-1 mm #
σ_IE = sigma_E*2.0;   #x10^-1 mm #
σ_II = sigma_I*1.0;   #x10^-1 mm #

τ_EE = 1
τ_EI = 1
τ_IE = 1
τ_II = 1

κ_EE = .0
κ_EI = .0
κ_IE = .0
κ_II = .0

d = 0.1
d_EE = d
d_EI = d
d_IE = d
d_II = d

Wc = 1
W0_EE =Wc*1
W0_EI = Wc*1
W0_IE = Wc*1
W0_II = Wc*1

D_E = 0.8^2
D_I = 0.8^2


h_E= 0.0
h_I= 0.0
d_vec = LinRange(.05,.4,55)
γ = 20
γ_E = γ
γ_I = γ
v = 1

k_vec = LinRange(0,10,10000)

V_plus = 1
V_neg = -V_plus
kc = 1

p1 =[]
p2 =[]
p3 =[]
p4 =[]


p = [α_EE,γ_E,σ_EE,κ_EE,τ_EE,v,D_E,d_EE,W0_EE,
     α_EI,γ_I,σ_EI,κ_EI,τ_EI,v,D_E,d_EI,W0_EI,
     α_IE,γ_E,σ_IE,κ_IE,τ_IE,v,D_I,d_IE,W0_IE,
     α_II,γ_I,σ_II,κ_II,τ_II,v,D_I,d_II,W0_II ]


for i = 1:length(d_vec)
print("\n",i,"\n")
global d = d_vec[i]

d_EE = d
d_EI = d
d_IE = d
d_II = d
reλprevious = 2
imλprevious = -1.0


s1 = [1 1 1]
for j = 1:length(k_vec)
    #print("\n",j,"\n")
        global k = k_vec[j]
        zeros = nlsolve(f!, [ reλprevious; imλprevious]).zero
        reλ = zeros[1]
        imλ = zeros[2]
         reλprevious = reλ
         imλprevious = imλ
        #print("\n k = ",k,"zeros =",zeros,"\n")
        s1 = vcat(s1,[reλ imλ  k])

end

reλprevious = 0.1
imλprevious = -1.0

s2 = [1 1 1]
for j = 1:length(k_vec)
    #print("\n",j,"\n")
        global k = k_vec[j]
        zeros = nlsolve(f!, [ reλprevious; imλprevious]).zero
        reλ = zeros[1]
        imλ = zeros[2]
         reλprevious = reλ
         imλprevious = imλ
        #print("\n k = ",k,"zeros =",zeros,"\n")
        s2 = vcat(s2,[reλ imλ  k])
end




p1 = vcat(p1,plot(real(s1[2:end,3]),real(s2[2:end,1]),lab = "re(λ_2)"))
p2 = vcat(p2,plot(real(s1[2:end,3]),real(s2[2:end,2]),lab="im(λ_2)"))
p3 = vcat(p3,plot(real(s1[2:end,3]),real(s1[2:end,1]),lab = "re(λ_1)"))
p4 = vcat(p4,plot(real(s1[2:end,3]),real(s1[2:end,2]),lab="im(λ_1)"))
end

@gif for i = 1:length(d_vec)
plot(p1[i],p2[i],p3[i],p4[i])
end
