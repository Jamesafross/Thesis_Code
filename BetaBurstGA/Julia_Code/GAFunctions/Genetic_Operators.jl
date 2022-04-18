using Random
#genetic operations
function crossover(numParams,parent1,parent2,max_crossoverpoints)
    #function performs a crossover genetic operator with 2 or 3 parents
    #output is one child
    child1 = zeros(numParams); child2 = zeros(numParams)
    child1 .= parent1
    child2 .= parent2
    child = zeros(numParams,2)
    numCrossoverPoints = Int.(round(max_crossoverpoints*rand()))
    CrossoverPoints = Int.(randperm(numParams))[1:numCrossoverPoints]
    crsov1 = parent1[CrossoverPoints]
    crsov2 = parent2[CrossoverPoints]

    child1[CrossoverPoints] .= crsov2
    child2[CrossoverPoints] .= crsov1

    child[:,1] = child1
    child[:,2] = child2

    return child #returns two children

end


function mutate(numParams,parent)
    #mutate genes in the parameter set
    numMutations = Int(round(numParams*rand()))
    mutant = parent
    mInd = Int.(randperm(numParams)[1:Int(numMutations)]) # mutation indicies
    for i = 1:numMutations
            mutant[mInd[i]] = mutant[mInd[i]] .+ (-0.5 + 1*rand())*mutant[mInd[i]] # changes a gene
    end
 return mutant
end

function crossover(Parent1,Parent2)
    numCrossover = rand(1:4)
    crossoverPoints = randperm(4)[1:numCrossover]

    crossoverParamsParent1 = Parent1.P[1:4]
    crossoverParamsParent2 = Parent2.P[1:4]
    
    for i = 1:numCrossover
        crossoverParamsParent1[i] = Parent2.P[i]
        crossoverParamsParent2[i] = Parent1.P[i]
    end
    child1 = phenotype(P = [crossoverParamsParent1;Parent1.P[5:end]])
    child2 = phenotype(P = [crossoverParamsParent2;Parent2.P[5:end]])
    #println(crossoverPoints)
    return child1,child2

end

function doCrossover(parents,n)
    nextGen = []
    i = 1
    sizeP = size(parents,1)
    while size(nextGen,1) < n
        r1,r2 = randperm(sizeP)[1:2]
        child1, child2 = crossover(parents[r1],parents[r2]);
        if i == 1
            nextGen = child1
            nextGen = cat(nextGen,child2,dims=1)
        else
            nextGen = cat(nextGen,child1,dims=1)
            nextGen = cat(nextGen,child2,dims=1)
        end

        i += 1

        #println(size(nextGen,1))
    
    end

    return nextGen
end

function makeNextGen(phenotypes,bestPhenotype,sizePop,n,constParams)
    
    children = doCrossover(phenotypes,n)
    numNew = sizePop - size(children,1) - 1
    return cat(bestPhenotype,children,GenPop(numNew,constParams),dims=1)
end

function normaliseFitness(phenotypes)
    fitnessSum = 0
    for i in size(phenotyes,1)
        fitnessSum += phenotypes[i].Fitness
    end
end



    

