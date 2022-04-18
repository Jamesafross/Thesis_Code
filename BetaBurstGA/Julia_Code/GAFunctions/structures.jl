using NLsolve
using LinearAlgebra
using Parameters
#mutable struct GenPhenotypes
#
#end



@with_kw mutable struct phenotype
    P
    Stats = 0
    Fitness = 0.0
    Rank = 0
end
