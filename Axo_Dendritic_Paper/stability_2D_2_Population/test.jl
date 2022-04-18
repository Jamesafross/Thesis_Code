using Plots
using LinearAlgebra
using NLsolve

function f!(F, x)
    #gEE
    F[1] = gEXT + tanh(γ*x[5])*C - x[1]
    #gEI
    F[2] = gEXT + tanh(γ*x[6])*C - x[2]
    #gIE
    F[3] = gEXT + tanh(γ*x[5])*C - x[3]
    #gII
    F[4] = gEXT + tanh(γ*x[6])*C - x[4]
    #VE
    F[5] = (Vplus*(x[1] - x[2]))/((1/tau) + x[1] + x[2]) - x[5]
    #VI
    F[6] = (Vplus*(x[3] - x[4]))/((1/tau) + x[3] + x[4]) - x[6]
end


x0 = [0.5;1;1;1;2;2]

γ = 18
gEXT=100
C = 10
tau=1
Vplus = 1
tau = 1
_zeros_ = nlsolve(f!,x0).zero
print(_zeros_)
