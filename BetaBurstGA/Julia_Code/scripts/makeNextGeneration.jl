using DifferentialEquations
using NLsolve
using LinearAlgebra
using Parameters
using MAT
using JLD
using Random

include("../massModelFunctions/FuncRunMassModel.jl")
include("../GAFunctions/GenerateInitialPopulation.jl")
include("../GAFunctions/tournament.jl")
include("../GAFunctions/structures.jl")
include("../massModelFunctions/stability.jl")
include("../massModelFunctions/rhsFunctions.jl")
include("../GAFunctions/Genetic_Operators.jl")
include("config.jl")


GenPhenotype = load("$JulDIR/Generation.jld", "GenPop")
file = matopen("$JulDIR/stats.mat")
phenotype_stats = read(file, "stats")
close(file)

file = matopen("$JulDIR/realDataStats.mat")
realdata_stats = read(file, "realdatastats")
close(file)



assignStats(GenPhenotype,phenotype_stats)
assignFitness(GenPhenotype,realdata_stats)

winnersTour = tournament(GenPhenotype,tourRounds,minParents)
bestInGen,bestFitness,bestindx = findBest(winnersTour)
bestFitnessParams = bestInGen.P
println("best phenotype in this generation is: \n $bestFitnessParams. \n \n with total fitness: $bestFitness \n \n at index ",bestindx)

open("$JulDIR/BestFitnessLog.txt", "w") do f #outputs best fitness phenotype to a .txt file
    write(f, "best phenotype in this generation is: $bestInGen. \n with total fitness: $bestFitness \n")
    end

nextGen = makeNextGen(winnersTour,bestInGen,sizePop,numChild,constParams)
save("$JulDIR/LastGeneration.jld", "GenPop", GenPhenotype)
save("$JulDIR/Generation.jld", "GenPop", nextGen)

# get the fitness of each phenotype


