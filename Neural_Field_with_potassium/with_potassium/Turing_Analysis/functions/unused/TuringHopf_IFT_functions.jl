function turinghopf_IFT!(F, x)
    global κSII=x[1]
    σ = [σEE σEI σIE σII]
    κV = [κVEE κVEI κVIE κVII]
    global κS = [κSEE κSEI κSIE κSII]
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
    #print("\n" ,SteadyStateNoPotassium(),"\n")

    ℰ = det(ℰfunc(V,R,α,τ0,η0,v,τ,Vsyn,σ,κV,κS,x[2],0,x[3]))
    dℰdω = dℰdωfunc(V,R,α,τ0,η0,v,τ,Vsyn,σ,κV,κS,x[2],0,x[3])
    dℰdk = dℰdkfunc(V,R,α,τ0,η0,v,τ,Vsyn,σ,κV,κS,x[2],0,x[3])

    F[1] = real(dℰdk)*imag(dℰdω) - real(dℰdω)*imag(dℰdk)
    F[2] = real(ℰ)
    F[3] = imag(ℰ)
end

function SteadyState_turinghopf_IFT(init,iters)
    X = init
    counter = 0
    stop_cond = 100
    SS = nlsolve(turinghopf_IFT!, [ X[1];X[2];X[3]],iterations=iters,ftol=10e-5)
    while (SS.f_converged == false || SS.zero[1] < 0 || abs(SS.zero[2]) < 0.01 ) && counter <100
        if abs(SS.zero[2]) < 0.01
            print("\n found hopf \n")
        end
        counter += 1
        #print("\n",counter,"\n")
        #X = init.+0.001*randn(3)
        SS = nlsolve(turinghopf_IFT!, [ X[1];X[2];X[3] ],iterations=iters,ftol=10e-5)
    end
    if counter == stop_cond
        print("Could not find zeros")
        return [NaN,NaN,NaN]
    else
        return SS.zero
    end
end
