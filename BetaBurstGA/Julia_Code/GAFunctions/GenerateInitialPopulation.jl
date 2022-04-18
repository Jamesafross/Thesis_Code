#functions for generating the initial population of parameters

function randNum(a,b,y)
    # creates a uniform random number between a and b
    return a .+ (b-a)*rand(y)
end

function GenerateGoodParameters(constParams)
   
    ΔE,ΔI,η_0E,η_0I,τE,τI =  constParams
    p = GenRandParams()
    σE,σI,τxE,τxI,
    κSEE,κSIE,κSEI,κSII,
    αEE,αIE,αEI,αII,
    κVEE,κVIE,κVEI,κVII,
    VsynEE,VsynIE,VsynEI,VsynII  = p

    params = κSEE,κSIE,κSEI,κSII,αEE,αIE,αEI,αII,
    κVEE,κVIE,κVEI,κVII,
    VsynEE,VsynIE,VsynEI,VsynII,ΔE,ΔI,η_0E,η_0I,τE,τI

    params2 = κSEE,κSIE,κSEI,κSII,
    αEE,αIE,αEI,αII,
    κVEE,κVIE,κVEI,κVII,
    VsynEE,VsynIE,VsynEI,VsynII,ΔE,ΔI,η_0E+randn(),η_0I+randn(),τE,τI

    EIGS = get_EigenValues(params)
    EIGS2 = get_EigenValues(params2)

    #println(maximum(real.(EIGS)))

    breakcond = false
    while  maximum(real.(EIGS)) > 0.0 || maximum(real.(EIGS2)) < 0.0
        #println(maximum(real.(EIGS)))
  

           p = GenRandParams()

            σE,σI,τxE,τxI,
            κSEE,κSIE,κSEI,κSII,
            αEE,αIE,αEI,αII,
            κVEE,κVIE,κVEI,κVII,
            VsynEE,VsynIE,VsynEI,VsynII = p

            params = κSEE,κSIE,κSEI,κSII,αEE,αIE,αEI,αII,
            κVEE,κVIE,κVEI,κVII,
            VsynEE,VsynIE,VsynEI,VsynII,ΔE,ΔI,η_0E,η_0I,τE,τI

            params2 = κSEE,κSIE,κSEI,κSII,
            αEE,αIE,αEI,αII,
            κVEE,κVIE,κVEI,κVII,
            VsynEE,VsynIE,VsynEI,VsynII,ΔE,ΔI,η_0E+randn(),η_0I+randn(),τE,τI

            EIGS = get_EigenValues(params)
            EIGS2 = get_EigenValues(params2)
    
    end
    #println( maximum(real.(EIGS)))
    #println( maximum(real.(get_EigenValues(params))))

    return p

end


function GenRandParams()
    InitPop = zeros(20)
    #row 1 : sigE
    #row 2 : sigI
    #row 3 : tauxE
    #row 4 : tauxI
    #row 5-8 : kappaS
    #row 9-12 : alpha
    #row 13-16 : kappaV
    #row 17-20 : VsynAB

    InitPop[1:2,:] .= randNum(0.001,0.1,2)
    InitPop[3:4,:] .= randNum(40,100,2)
    InitPop[5:8,:] .= randNum(0.1,5,4)
    InitPop[9:12,:] .= randNum(0.1,1.0,4)
    InitPop[13:16,:] .= randNum(0.0,0.2,4)
    InitPop[14,:] = InitPop[15,:]
    InitPop[17:18,:] .= randNum(1,30,2)
    InitPop[19:20,:] .= randNum(-1,-30,2)
    return InitPop
end


function GenPop(SizePop,constParams)
    phenotypes = []
    for i = 1:SizePop
        #println(i)
        if i == 1
            phenotypes = phenotype(P = GenerateGoodParameters(constParams))
        else
            phenotypes = cat(phenotypes,phenotype(P = GenerateGoodParameters(constParams)),dims=1)
        end
    end
    return phenotypes
end
