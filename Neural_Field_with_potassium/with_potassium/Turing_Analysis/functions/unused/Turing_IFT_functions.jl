function turing_IFT!(F, x)
    global κSII=x[1]
    σ = [σEE σEI σIE σII]
    global κV = [κVEE κVEI κVIE κVII]
    κS = [κSEE κSEI κSIE κSII]
    η0 = [ηE, ηI]
    τ0 = [τ0EE τ0EI τ0IE τ0II]
    α = [αEE αEI αIE αII]
    Vsyn = [VsynEE VsynEI VsynIE VsynII]
    τ = [τE τI]
    β = [β1 β2 β3 β4]
    rE,vE,rI,vI = SteadyStateNoPotassium()
    V = [vE, vI]
    R = [rE, rI]

    gEE = κSEE * rE
    gEI = κSEI * rI
    gIE = κSIE * rE
    gII = κSII * rI
    print("\n" ,SteadyStateNoPotassium(),"\n")



    ℰ = det(ℰfunc(V,R,α,τ0,η0,v,τ,Vsyn,σ,κV,κS,x[2],0,0))
    dℰdω = dℰdωfunc(V,R,α,τ0,η0,v,τ,Vsyn,σ,κV,κS,x[2],0,0)
    dℰdk = dℰdkfunc(V,R,α,τ0,η0,v,τ,Vsyn,σ,κV,κS,x[2],0,0)

    F[1] = real(dℰdk)*imag(dℰdω) - real(dℰdω)*imag(dℰdk)
    F[2] = ℰ

end

function SteadyState_turing_IFT(init)
    X = init
    counter = 0
    SS = nlsolve(turing_IFT!, [ X[1],X[2] ])
    while SS.f_converged == false
        counter += 1
        X = init.+0.01*randn(2)
        SS = nlsolve(turing_IFT!, [ X[1],X[2]])
    end
    if counter == 10000
        print("Could not find zeros")
        return 0
    else
        return SS.zero
    end
end
