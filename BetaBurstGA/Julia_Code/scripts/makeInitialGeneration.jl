using DifferentialEquations
using NLsolve
using LinearAlgebra
using Parameters
using MAT
using JLD
include("../massModelFunctions/FuncRunMassModel.jl")
include("../GAFunctions/GenerateInitialPopulation.jl")
include("../GAFunctions/structures.jl")
include("../massModelFunctions/stability.jl")
include("../massModelFunctions/rhsFunctions.jl")
include("config.jl")


constParams = ΔE,ΔI,η_0E,η_0I,τE,τI

println("Generating Initial Population of ", sizePop, " Phenotypes ... ")
InitPop = GenPop(sizePop,constParams)
println("Done.")

save("$JulDIR/Generation.jld", "GenPop", InitPop)