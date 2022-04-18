function rhsFun(du, u, p, t)

    println(t)

    ΔE, ΔI, κVEE,κVEI,κVIE,κVII,ηE,ηI,τE,τI,
    σEE,σEI,σIE, σII, v, κSEE,κSEI,κSIE, κSII,
    αEE, αEI, αIE, αII, δ,
    A1, A2, A3, B1, B2, B3, C1, C2, C3,
    D1, D2, D3,βK,βηE,βηI,βκV,∇,∇ψ,options= p

    Add_Potassium,gaps_option,dimension = options

    @views RE   = u[1:X]
    @views VE   = u[X+1:2*X]
    @views RI   = u[2*X+1:3*X]
    @views VI   = u[3*X+1:4*X]

    @views ψEE = u[4*X+1:5*X]
    @views AuxEE = u[5*X+1:6*X]
    @views ψEI = u[6*X+1:7*X]
    @views AuxEI = u[7*X+1:8*X]
    @views ψIE = u[8*X+1:9*X]
    @views AuxIE = u[9*X+1:10*X]
    @views ψII = u[10*X+1:11*X]
    @views AuxII = u[11*X+1:12*X]

    @views gEE  = u[12*X+1:13*X]
    @views pEE  = u[13*X+1:14*X]
    @views gEI  = u[14*X+1:15*X]
    @views pEI  = u[15*X+1:16*X]
    @views gIE  = u[16*X+1:17*X]
    @views pIE  = u[17*X+1:18*X]
    @views gII  = u[18*X+1:19*X]
    @views pII  = u[19*X+1:20*X]

    if Add_Potassium == 1
        @views K  = u[20*X+1:21*X]
        @views ηE = u[21*X+1:22*X]
        @views ηI = u[22*X+1:23*X]
        @views κV = u[23*X+1:24*X]
        if gaps_option == "ALL-ALL"
            gaps_term_RE = -κV.*RE .- κV .* RE
            gaps_term_RI = .-κV.*RI .- κV .* RI
            gaps_term_VE = κV.*(VI.-VE)
            gaps_term_VI =  κV.*(VE.-VI)
        elseif gaps_option == "I-I"
            gaps_term_RE = 0.0
            gaps_term_RI = -κV.*RI
            gaps_term_VE = 0.0
            gaps_term_VI = 0.0
        elseif gaps_option == "OFF"
            gaps_term_RE = 0.0
            gaps_term_RI = 0.0
            gaps_term_VE = 0.0
            gaps_term_VI = 0.0
        end
    elseif Add_Potassium == 0
        gaps_term_RE = -κVEE.*RE .- κVEI .* RE
        gaps_term_RI = -κVII.*RI .- κVIE .* RI
        gaps_term_VE =  κVEI.*(VI-VE)
        gaps_term_VI =  κVIE.*(VE-VI)
    end



    dREdt =(1.0/τE)*(gaps_term_RE .-gEE.*RE .- gEI.*RE .+ 2.0.* RE .* VE .+ (ΔE / (τE*pi)))
    dVEdt =(1.0/τE).*(gaps_term_VE .+ gEE.*(VsynEE .- VE) .+ gEI.*(VsynEI .- VE) .- (τE^2)*(pi^2.0) * (RE.^2.0) .+  VE.^2.0 .+ ηE )
    dRIdt =(1.0/τI).*(gaps_term_RI .-gII.*RI .-  gIE.*RI .+ 2.0 .* RI .* VI .+ (ΔI / (τI*pi)))
    dVIdt =(1.0/τI).*(gaps_term_VI .+ gIE.*(VsynIE .-VI) .+ gII.*(VsynII .- VI) .- (pi^2).*(τI^2.0) .* (RI.^2.0) .+ VI.^2.0 .+ ηI )

    if dimension == 1
     ψ_auxEE = v.*(-(1.0/σEE).*AuxEE .+ ∇ψ*ψEE .+ (1.0/σEE^2.0)*RE .+ (1.0/(σEE*v)).*dREdt)
     ψ_auxEI = v.*(-(1.0/σEI).*AuxEI .+ ∇ψ*ψEI .+ (1.0/σEI^2.0)*RI .+ (1.0/(σEI*v)).*dRIdt)
     ψ_auxIE = v.*(-(1.0/σIE).*AuxIE .+ ∇ψ*ψIE .+ (1.0/σIE^2.0)*RE .+ (1.0/(σIE*v)).*dREdt)
     ψ_auxII = v.*(-(1.0/σII).*AuxII .+ ∇ψ*ψII.+ (1.0/σII^2.0)*RI .+ (1.0/(σII*v)).*dRIdt)
    elseif dimension == 2
     ψ_auxEE = v.*(-(1.0/σEE).*AuxEE .+ ∇ψ*ψEE .+ (1.0/σEE^2.0)*RE )
     ψ_auxEI = v.*(-(1.0/σEI).*AuxEI .+ ∇ψ*ψEI .+ (1.0/σEI^2.0)*RI )
     ψ_auxIE = v.*(-(1.0/σIE).*AuxIE .+ ∇ψ*ψIE .+ (1.0/σIE^2.0)*RE )
     ψ_auxII = v.*(-(1.0/σII).*AuxII .+ ∇ψ*ψII .+ (1.0/σII^2.0)*RI )
    end



# R & V
    du[1:X] = dREdt
    du[X+1:2*X] =dVEdt

    du[2*X+1:3*X] =dRIdt
    du[3*X+1:4*X] =dVIdt


#  ψ and Aux variables
    du[4*X+1:5*X]   = v.*((-1.0/σEE).* ψEE .+ AuxEE) # ψEE
    du[5*X+1:6*X]   =  ψ_auxEE                     #aux  ψEE variable
    du[6*X+1:7*X]   = v.*((-1.0/σEI).* ψEI .+ AuxEI) # ψEI
    du[7*X+1:8*X]   =  ψ_auxEI                     #aux  ψEI variable
    du[8*X+1:9*X]   = v.*((-1.0/σIE).* ψIE .+ AuxIE) # ψIE
    du[9*X+1:10*X]   =  ψ_auxIE                     #aux  ψIE variable
    du[10*X+1:11*X]   = v.*((-1.0/σII).* ψII .+ AuxII) # ψI
    du[11*X+1:12*X]   =  ψ_auxII                     #aux  ψI variable

    # g and Aux variables
    du[12*X+1:13*X]   = αEE.*(-gEE .+ pEE) #gEE
    du[13*X+1:14*X]   = αEE.*(-pEE .+ κSEE* ψEE)#pEE

    du[14*X+1:15*X]  = αEI.*(-gEI .+ pEI)#gEI
    du[15*X+1:16*X] = αEI.*(-pEI .+ κSEI* ψEI)#pEI

    du[16*X+1:17*X] = αIE.*(-gIE .+ pIE)#gIE
    du[17*X+1:18*X] = αIE.*(-pIE .+ κSIE* ψIE)#pIE

    du[18*X+1:19*X] = αII.*(-gII .+ pII)#gII
    du[19*X+1:20*X] = αII.*(-pII .+ κSII* ψII)#pII


    # potassium,drive and gaps strength variables
    if Add_Potassium == 1
    du[20*X+1:21*X] = (1.0/βK)*(.-δ*K .+  fK.(RE.+RI,A1,A2,A3) .+ A4.*∇ *K)
    du[21*X+1:22*X] = (1.0/βηE)*(.-ηE .+ fηE.(K,B1,B2,B3))
    du[22*X+1:23*X] = (1.0/βηI)*(.-ηI .+ fηI.(K,C1,C2,C3))
    #gap junctions
    du[23*X+1:24*X] = (1.0/βκV)*(.-κV .+ fκV.(K,D1,D2,D3))
    end
end

function noiseRHS(du,u,p,t)
    ΔE, ΔI, κVEE,κVEI,κVIE,κVII,ηE,ηI,τE,τI,
    σEE,σEI,σIE, σII, v, κSEE,κSEI,κSIE, κSII,
    αEE, αEI, αIE, αII, δ,
    A1, A2, A3, B1, B2, B3, C1, C2, C3,
    D1, D2, D3,βK,βηE,βηI,βκV,∇,∇ψ,options= p

end
