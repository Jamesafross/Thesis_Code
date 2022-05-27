function hopf_IFT!(F, x)
    global κSII = x[1]
    #print(κSII)
    σ = [σEE σEI σIE σII]
    κV = [κVEE κVEI κVIE κVII]
    κS = [κSEE κSEI κSIE κSII]
    η0 = [ηE, ηI]
    τ0 = [τ0EE τ0EI τ0IE τ0II]
    α = [αEE αEI αIE αII]
    Vsyn = [VsynEE VsynEI VsynIE VsynII]
    τ = [τE τI]
    rE,vE,rI,vI = SteadyStateNoPotassium()
    global V = [vE, vI]
    global R = [rE, rI]
    global gEE = κSEE * rE
    global gEI = κSEI * rI
    global gIE = κSIE * rE
    global gII = κSII * rI

    ℰ = det(ℰfunc(V,R,α,τ0,η0,v,τ,Vsyn,σ,κV,κS,0,0,x[2]))

    F[1] = real(ℰ)
    F[2] = imag(ℰ)
end

function SteadyState_g_hopf_IFT(init)
    X = init
    stopcond = 5000
    counter = 0
    SS = nlsolve(hopf_IFT!, [ X[1];X[2] ],iterations=60000)
    while SS.f_converged == false
            counter += 1
            X = init+0.1.*randn(2)
            SS = nlsolve(hopf_IFT!, [ X[1];X[2] ],iterations=60000)

    end
    if SS.f_converged == true
        print("\n converged! \n")
    end
        return SS.zero
end
