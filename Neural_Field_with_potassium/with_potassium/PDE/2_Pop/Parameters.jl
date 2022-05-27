using Parameters


struct ModelSetup
    dimension::Int64
    gaps_option::String
    Add_Potassium::Int64
    ∇
    ∇ψ 
end

@with_kw struct ModelParams{R}

ΔE::R = 0.2
ΔI::R = 0.2
τE::R = 9.0
τI::R = 11.0
VsynEE::R = 12.0
VsynEI::R = VsynEE*(-10.0/16.0)
VsynIE::R = VsynEE*(10.0/16.0)
VsynII::R = VsynEE*(-11.0/16.0)

#g and ψ params
σEE::R = 0.2
σEI::R =  1.5
σIE::R = 0.2
σII::R = 1.5

#conduction velocity
v::R =.15 #(.15 2D)

#synaptic time constant
αEE::R = 0.7
αEI::R = 0.5
αIE::R = 0.7
αII::R = 0.8

#synaptic connection strengths
κSEE::R = 4.5 #4.5(2D)
κSEI::R = κSEE*(4.9/5.0)
κSIE::R = κSEE*(4.9/5.0)
κSII::R =  κSEE*(3.0/5.0)

#Drive params
ηE::R = 1.
ηI::R = 1.
# gaps Strengths (for model w/o potassium)
κV::R = .0
κVEE::R = 0.0
κVEI::R = 0.0
κVIE::R = κVEI
κVII::R = 0.0
τ0EE::R = 0.
τ0EI::R = 0.
τ0IE::R = 0.
τ0II::R = 0.

# Potassium params # # #
# sigmoid params
A1::R = 2.0
A2::R = 2.0 #steepness
A3::R = 0.0
βK::R = 1.0 #time constant
# other potassiumn params
δ::R = 0.2 #decay rate
A4::R = 0.3 #diffusion rate

# Drive params # # # #
#Excitatory
B1::R =15.0
B2::R = 8.0# steepness
B3::R = 1.0 #threshold
B4::R = 0.0
βηE::R = 25. #time constant

#Inhibitory
C1::R = 15.0
C2::R = 8.0 # steepness
C3::R = 1.0 # threshold
C4::R = 0.0
βηI::R = 25. #time constant

#gaps Strength params
D1::R = 0.5
D2::R = -2.0 # steepness
D3::R = 1.0
βκV::R = 4.0 #time constant

end


@with_kw struct ModelParams2{R}

    ΔE::R = 0.5
    ΔI::R = 0.5
    τE::R = 10.0
    τI::R = 11.0
    VsynEE::R = 16.0
    VsynEI::R = -10.
    VsynIE::R = 10.
    VsynII::R = -15.5
    
    #g and ψ params
    σEE::R = 0.4
    σEI::R =  0.8
    σIE::R = 0.4
    σII::R = 0.8
    
    #conduction velocity
    v::R =.1 #(.15 2D)
    
    #synaptic time constant
    αEE::R = 0.7
    αEI::R = 0.5
    αIE::R = 0.7
    αII::R = 0.8
    
    #synaptic connection strengths
    κSEE::R = 4.8#4.5(2D)
    κSEI::R = 5.5
    κSIE::R = 5.5
    κSII::R =  4.1
    
    #Drive params
    ηE::R = 1.
    ηI::R = 1.
    # gaps Strengths (for model w/o potassium)
    κV::R = .0
    κVEE::R = 0.44
    κVEI::R = 0.44
    κVIE::R = κVEI
    κVII::R = 0.44
    τ0EE::R = 0.
    τ0EI::R = 0.
    τ0IE::R = 0.
    τ0II::R = 0.
    
    # Potassium params # # #
    # sigmoid params
    A1::R = 2.0
    A2::R = 2.0 #steepness
    A3::R = 0.0
    βK::R = 1.0 #time constant
    # other potassiumn params
    δ::R = 0.2 #decay rate
    A4::R = 0.3 #diffusion rate
    
    # Drive params # # # #
    #Excitatory
    B1::R =15.0
    B2::R = 8.0# steepness
    B3::R = 1.0 #threshold
    B4::R = 0.0
    βηE::R = 25. #time constant
    
    #Inhibitory
    C1::R = 15.0
    C2::R = 8.0 # steepness
    C3::R = 1.0 # threshold
    C4::R = 0.0
    βηI::R = 25. #time constant
    
    #gaps Strength params
    D1::R = 0.5
    D2::R = -2.0 # steepness
    D3::R = 1.0
    βκV::R = 4.0 #time constant
    
    end
thesis_potassium_sims = ModelParams()
thesis_nopotassium_sims = ModelParams2()