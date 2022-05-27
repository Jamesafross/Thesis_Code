@with_kw mutable struct parameters
    ΔE = 0.20
    ΔI = 0.20
    τE = 9.0
    τI = 11.0
    VsynEE = 11.0
    VsynEI = VsynEE*(-10.0/16.0)
    VsynIE = VsynEE*(10.0/16.0)
    VsynII = VsynEE*(-11.0/16.0)

    #g and ψ params
    σEE = 0.2
    σEI = 1.6
    σIE = 0.2
    σII = 1.6

    #conduction velocity
    v =.04

    #synaptic time constant
    αEE = 0.7
    αEI = 0.5
    αIE = 0.7
    αII = 0.8

    #synaptic connection strengths

    κSEE = 2.28
    κSEI = κSEE*(4.6/5.0)
    κSIE = κSEE*(4.6/5.0)
    κSII =  κSEE*(3.1/5.0)

    #Drive params
    ηE = 1
    ηI = 1
    # gaps Strengths
    κV = .0
    κVEE = 0.0
    κVEI = 0.0
    κVIE = κVEI
    κVII = 0.0
    τ0EE = 0
    τ0EI = 0
    τ0IE = 0
    τ0II = 0

    # Potassium params # # #
    # sigmoid params
    A1 = 1.0
    A2 = 2.0 #steepness
    A3 = 0.0
    βK = 1.0 #time constant
    # other potassiumn params
    δ = 0.2 #decay rate
    A4 = 0.3 #diffusion rate

    # # # # # # # # # # # #
    # Drive params # # # #
    #Excitatory
    B1 =15.0
    B2 = 8.0# steepness
    B3 = 1.0 #threshold
    B4 = 0.0
    βηE = 25 #time constant

    #Inhibitory
    C1 = 15.0
    C2 = 8.0 # steepness
    C3 = 1.0 # threshold
    C4 = 0.0
    βηI = 25 #time constant


    # # # # # # # # # # # #
    # gaps Strength params
    D1 = 0.5
    D2 = -2.0 # steepness
    D3 = 1.0
    βκV = 4.0 #time constant

end