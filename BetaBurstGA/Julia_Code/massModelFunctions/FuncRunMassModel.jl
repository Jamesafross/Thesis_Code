include("rhsFunctions.jl")
function runMassModel(u0, dt,saveat, p,tspan)
    # This function runs the mass model using
    # the SingleNodeEM solver.
    # The output is the sum of excitatory and
    # inhibitory current
    σE,σI,τxE,τxI,
    κSEE,κSIE,κSEI,κSII,
    αEE,αIE,αEI,αII,
    κVEE,κVIE,κVEI,κVII,
    VsynEE,VsynIE,VsynEI,VsynII,ΔE,ΔI,η_0E,η_0I,τE,τI = p

    prob = SDEProblem(f,g,u0,tspan,p)

    sol = solve(prob,SOSRA(),maxiters=10^20,saveat = saveat)

    rE = sol[1,:]
    rI = sol[2,:]
    vE = sol[3,:]
    vI = sol[4,:]
    gEE = sol[6,:]
    gIE = sol[8,:]
    gEI = sol[10,:]
    gII = sol[12,:]
    noiseE = sol[13,:]
    noiseI = sol[14,:]
    WE = (1 .-conj.(sol[1,:]))./(1 .+conj.(sol[1,:]));
    WI = (1 .-conj.(sol[2,:]))./(1 .+conj.(sol[2,:]));
    currentE =  real(gEE.*(VsynEE .- vE) .+ gEI.*(VsynEI .- vE));
    currentI =  real(gIE.*(VsynIE .- vI) .+ gII.*(VsynII .- vI));
    current = currentE .+ currentI

    return current
end
