using NLsolve
using LinearAlgebra
using Plots
using DifferentialEquations
using FFTW
include("GenerateInitialPopulation.jl")
include("../massModelFunctions/stability.jl")
include("../massModelFunctions/rhsFunctions.jl")


ΔE = 0.2;
ΔI = 0.2;
η_0E = 2.0;
η_0I = 1.8;
τE = 12;
τI = 15;

constParams = ΔE,ΔI,η_0E,η_0I,τE,τI
pop = GenerateGoodParameters(constParams)
print("found params")
σE,σI,τxE,τxI,
κSEE,κSIE,κSEI,κSII,
αEE,αIE,αEI,αII,
κVEE,κVIE,κVEI,κVII,
VsynEE,VsynIE,VsynEI,VsynII  = pop

params = [κSEE,κSIE,κSEI,κSII,αEE,αIE,αEI,αII,
κVEE,κVIE,κVEI,κVII,
VsynEE,VsynIE,VsynEI,VsynII,ΔE,ΔI,η_0E,η_0I,τE,τI]

params2 = [κSEE,κSIE,κSEI,κSII,αEE,αIE,αEI,αII,
κVEE,κVIE,κVEI,κVII,
VsynEE,VsynIE,VsynEI,VsynII,ΔE,ΔI,η_0E+randn(),η_0I+randn(),τE,τI]


tspan = (0.0,4000.0)


SS1 = init_conds_SS(params) .+ 0.000
SS2 = init_conds_SS(params2) .+ 0.000

prob1 = ODEProblem(fdeterministic,u01,tspan,params)
sol1 = solve(prob1,Tsit5(),saveat = 0.0:0.1:4000)



prob2 = ODEProblem(fdeterministic,u02,tspan,params2)
sol2 = solve(prob2,Tsit5(),saveat = 1000:0.1:4000)

p1 = plot(sol1[1,:])
p2 = plot(sol2[1,:])

t0 = 0
tmax = 3
# time coordinate
t = t0:0.1/1000:tmax

# signal
signal = sol2[1,:] # sin (2π f t)

# Fourier Transform of it
F = fft(signal) |> fftshift
freqs = fftfreq(length(t), 1.0/(0.1/1000)) |> fftshift

# plots
time_domain = plot(t, signal, title = "Signal")
freq_domain = plot(freqs, abs.(F), title = "Spectrum", xlim=(5, +60))
plot(p1,time_domain, freq_domain, layout = 3)

#plot(p1,p2)




#scatter(EIGS)
