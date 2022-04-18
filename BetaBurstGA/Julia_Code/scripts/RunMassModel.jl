using DifferentialEquations
using NLsolve
using LinearAlgebra
using Parameters
using MAT
using JLD
using ProgressBars
include("../massModelFunctions/FuncRunMassModel.jl")
include("../GAFunctions/GenerateInitialPopulation.jl")
include("../GAFunctions/structures.jl")
include("../massModelFunctions/stability.jl")
include("../massModelFunctions/rhsFunctions.jl")
include("config.jl")



@time Pop = load("$JulDIR/Generation.jld", "GenPop")

PopCurrent = zeros(length(tRange),MCtrials,sizePop)

for i in ProgressBar(1:sizePop)
    #println("Running MC trials for phenotype: ", i)
    PhenotypeParams = Pop[i].P
    σE = PhenotypeParams[1]
    σI = PhenotypeParams[2]
    τxE = PhenotypeParams[3]
    τxI = PhenotypeParams[4]
    κSEE = PhenotypeParams[5]
    κSIE = PhenotypeParams[6]
    κSEI = PhenotypeParams[7]
    κSII = PhenotypeParams[8]
    αEE = PhenotypeParams[9]
    αIE = PhenotypeParams[10]
    αEI = PhenotypeParams[11]
    αII = PhenotypeParams[12]
    κVEE= PhenotypeParams[13]
    κVIE= PhenotypeParams[14]
    κVEI = PhenotypeParams[15]
    κVII = PhenotypeParams[16]
    VsynEE = PhenotypeParams[17]
    VsynIE = PhenotypeParams[18]
    VsynEI = PhenotypeParams[19]
    VsynII = PhenotypeParams[20]

    params = σE,σI,τxE,τxI,
    κSEE,κSIE,κSEI,κSII,
    αEE,αIE,αEI,αII,
    κVEE,κVIE,κVEI,κVII,
    VsynEE,VsynIE,VsynEI,VsynII,ΔE,ΔI,η_0E,η_0I,τE,τI

    u0 =  init_conds_SS(params[5:end])
    u0 = cat(u0,[randn(); randn()],dims=1)

    for j = 1:MCtrials
        #println("Running MC trial ",j)
        PopCurrent[:,j,i] = runMassModel(u0, dt, saveat, params,tspan)
    end
end

#save("$JulDIR/Generation.jld", "GenPop", Pop)

file = matopen("$MatDIR/PopCurrent.mat", "w") # save current for use in HMM-MAR
write(file, "popcurrent", 10*PopCurrent[Int(bufferPeriod/saveat):end,:,:])
close(file)  