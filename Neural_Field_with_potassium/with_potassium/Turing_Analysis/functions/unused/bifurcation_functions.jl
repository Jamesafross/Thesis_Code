dηABdω(λ,α_ab) = (-2*im*α_ab^2)/((1+λ)^3)

W_AB(v,λ,k,σAB,κSAB) = (κSAB*v*(λ*σAB + v))/((v + σAB*(λ - im*k*v))*(v + σAB*(λ + im*k*v)))

ηAB(λ,α_ab) = (1+((λ)/α_ab))^-2

function dW_ABdk(v,λ,k,σAB,κSAB)
    tFrac = (-2*k*κSAB*σAB^2*v^3)*(v+σAB*λ)
    bFrac = ((v+σAB*(λ-im*k*v))^2)*((v+σAB*(λ+im*k*v))^2)
    return tFrac/bFrac
end

function dW_ABdω(v,λ,k,σAB,κSAB)
    tFrac = -0.5*im*κSAB*σAB*v
    bFrac1 = (v+σAB*(λ-im*k*v)^2)*((v+σAB*(λ+im*k*v))^2)
    bFrac2 = (v+σAB*(λ-im*k*v)^2)*((v+σAB*(λ-im*k*v))^2)
    return tFrac/bFrac1 + tFrac/bFrac2

end

function ℰfunc(V,R,α,τ0,η0,v,τ,Vsyn,σ,κV,κS,k,μ,ω)
    λ= μ+im*ω
    VbarE,VbarI,RbarE,RbarI = V[1],V[2],R[1],R[2]
    τ0EE,τ0EI,τ0IE,τ0II = τ0
    τE,τI = τ
    VsynEE,VsynEI,VsynIE,VsynII = Vsyn
    σEE,σEI,σIE,σII = σ
    κVEE,κVEI,κVIE,κVII = κV
    κSEE,κSEI,κSIE,κSII = κS
    αEE,αEI,αIE,αII = α
    gbarEE,gbarEI,gbarIE,gbarII = κSEE*RbarE,κSEI*RbarI,κSIE*RbarE,κSII*RbarI


    W_EE = W_AB(v,λ,k,σEE,κSEE)
    W_EI = W_AB(v,λ,k,σEI,κSEI)
    W_IE = W_AB(v,λ,k,σIE,κSIE)
    W_II = W_AB(v,λ,k,σII,κSII)

    ηEE = ηAB(λ,αEE)
    ηEI = ηAB(λ,αEI)
    ηIE = ηAB(λ,αIE)
    ηII = ηAB(λ,αII)

    AEE = 2*VbarE - (gbarEE + κVEE) - (gbarEI + κVEI) - RbarE*ηEE*W_EE
    AEI = -RbarE*ηEI*W_EI
    AIE = -RbarI*ηIE*W_IE
    AII = 2*VbarI - (gbarIE + κVIE) - (gbarII + κVII) - RbarI*ηII*W_II

    BEE = 2*RbarE
    BEI = 0
    BIE = 0
    BII = 2*RbarI

    CEE = -2*(τE^2)*(pi^2)*RbarE + (VsynEE - VbarE)*ηEE*W_EE
    CEI = (VsynEI - VbarE)*ηEI*W_EI
    CIE = (VsynIE - VbarI)*ηIE*W_IE
    CII = -2*(τI^2)*(pi^2)*RbarI + (VsynII - VbarI)*ηII*W_II

    DEE = 2*VbarE - (κVEE + gbarEE) - (κVEI + gbarEI) + κVEE
    DEI = κVEI
    DIE = κVIE
    DII = 2*VbarI - (κVIE + gbarIE) - (κVII + gbarII) + κVII

    return [AEE - τE*λ           AEI                  BEE                  BEI;
            AIE                  AII - τI*λ           BIE                  BII;
            CEE                  CEI                  DEE - τE*λ           DEI;
            CIE                  CII                  DIE                  DII - τI*λ]
end

function dℰdωfunc(V,R,α,τ0,η0,v,τ,Vsyn,σ,κV,κS,k,μ,ω)

    #ω = im(λ)
    λ = μ + im*ω

    VbarE,VbarI,RbarE,RbarI = V[1],V[2],R[1],R[2]
    τ0EE,τ0EI,τ0IE,τ0II = τ0
    τE,τI = τ
    VsynEE,VsynEI,VsynIE,VsynII = Vsyn
    σEE,σEI,σIE,σII = σ
    κVEE,κVEI,κVIE,κVII = κV
    κSEE,κSEI,κSIE,κSII = κS
    αEE,αEI,αIE,αII = α
    gbarEE,gbarEI,gbarIE,gbarII = κSEE*RbarE,κSEI*RbarI,κSIE*RbarE,κSII*RbarI


    W_EE = W_AB(v,λ,k,σEE,κSEE)
    W_EI = W_AB(v,λ,k,σEI,κSEI)
    W_IE = W_AB(v,λ,k,σIE,κSIE)
    W_II = W_AB(v,λ,k,σII,κSII)


    dW_EEdω = dW_ABdω(v,λ,k,σEE,κSEE)
    dW_EIdω = dW_ABdω(v,λ,k,σEI,κSEI)
    dW_IEdω = dW_ABdω(v,λ,k,σIE,κSIE)
    dW_IIdω = dW_ABdω(v,λ,k,σII,κSII)


    ηEE = ηAB(λ,αEE)
    ηEI = ηAB(λ,αEI)
    ηIE = ηAB(λ,αIE)
    ηII = ηAB(λ,αII)


    dηEEdω = dηABdω(λ,αEE)
    dηEIdω = dηABdω(λ,αEI)
    dηIEdω = dηABdω(λ,αIE)
    dηIIdω = dηABdω(λ,αII)

    dAEEdω = -RbarE*(dηEEdω*W_EE+ηEE*dW_EEdω)
    dAEIdω = -RbarE*(dηEIdω*W_EI+ηEI*dW_EIdω)
    dAIEdω = -RbarI*(dηIEdω*W_IE+ηIE*dW_IEdω)
    dAIIdω = -RbarI*(dηIIdω*W_II+ηII*dW_IIdω)

    dBEEdω = 0
    dBEIdω = 0
    dBIEdω = 0
    dBIIdω = 0

    dCEEdω = (VsynEE - VbarE)*(dηEEdω*W_EE+ηEE*dW_EEdω)
    dCEIdω = (VsynEI - VbarE)*(dηEIdω*W_EI+ηEI*dW_EIdω)
    dCIEdω = (VsynIE - VbarI)*(dηIEdω*W_IE+ηIE*dW_IEdω)
    dCIIdω = (VsynII - VbarI)*(dηIIdω*W_II+ηII*dW_IIdω)

    dDEEdω = 0
    dDEIdω = 0
    dDIEdω = 0
    dDIIdω = 0

    MAT_1  = ℰfunc(V,R,α,τ0,η0,v,τ,Vsyn,σ,κV,κS,k,μ,ω)

    MAT_2 = [dAEEdω - τE*(im)      dAEIdω                  dBEEdω                  dBEIdω;
              dAIEdω                  dAIIdω - τI*(im)      dBIEdω                  dBIIdω;
              dCEEdω                  dCEIdω                  dDEEdω - τE*(im)      dDEIdω;
              dCIEdω                  dCIIdω                  dDIEdω                  dDIIdω - τI*(im) ]

    dℰdω = tr((det(MAT_1))*(inv(MAT_1))*MAT_2)

    return dℰdω

end


function dℰdkfunc(V,R,α,τ0,η0,v,τ,Vsyn,σ,κV,κS,k,μ,ω)
    #ω = im(λ)
    λ = μ + im*ω
    VbarE,VbarI,RbarE,RbarI = V[1],V[2],R[1],R[2]
    τ0EE,τ0EI,τ0IE,τ0II = τ0
    τE,τI = τ
    VsynEE,VsynEI,VsynIE,VsynII = Vsyn
    σEE,σEI,σIE,σII = σ
    κVEE,κVEI,κVIE,κVII = κV
    κSEE,κSEI,κSIE,κSII = κS
    αEE,αEI,αIE,αII = α
    gbarEE,gbarEI,gbarIE,gbarII = κSEE*RbarE,κSEI*RbarI,κSIE*RbarE,κSII*RbarI

    W_EE = W_AB(v,λ,k,σEE,κSEE)
    W_EI = W_AB(v,λ,k,σEI,κSEI)
    W_IE = W_AB(v,λ,k,σIE,κSIE)
    W_II = W_AB(v,λ,k,σII,κSII)

    ηAB(λ,α_ab) = (1+((λ)/α_ab))^-2
    ηEE = ηAB(λ,αEE)
    ηEI = ηAB(λ,αEI)
    ηIE = ηAB(λ,αIE)
    ηII = ηAB(λ,αII)

    dW_EEdk = dW_ABdk(v,λ,k,σEE,κSEE)
    dW_EIdk = dW_ABdk(v,λ,k,σEI,κSEI)
    dW_IEdk = dW_ABdk(v,λ,k,σIE,κSIE)
    dW_IIdk = dW_ABdk(v,λ,k,σII,κSII)

    MAT_1 = ℰfunc(V,R,α,τ0,η0,v,τ,Vsyn,σ,κV,κS,k,μ,ω)

    dAEEdk = -RbarE*(ηEE*dW_EEdk)
    dAEIdk = -RbarE*(ηEI*dW_EIdk)
    dAIEdk = -RbarI*(ηIE*dW_IEdk)
    dAIIdk = -RbarI*(ηII*dW_IIdk)

    dBEEdk = 0
    dBEIdk = 0
    dBIEdk = 0
    dBIIdk = 0

    dCEEdk = (VsynEE - VbarE)*(ηEE*dW_EEdk)
    dCEIdk = (VsynEI - VbarE)*(ηEI*dW_EIdk)
    dCIEdk = (VsynIE - VbarI)*(ηIE*dW_IEdk)
    dCIIdk = (VsynII - VbarI)*(ηII*dW_IIdk)

    dDEEdk = 0
    dDEIdk = 0
    dDIEdk = 0
    dDIIdk = 0

    MAT_2 = [dAEEdk      dAEIdk                  dBEEdk                  dBEIdk;
              dAIEdk                  dAIIdk       dBIEdk                  dBIIdk;
              dCEEdk                  dCEIdk                  dDEEdk      dDEIdk;
              dCIEdk                  dCIIdk                  dDIEdk                  dDIIdk  ]

    dℰdk = tr((det(MAT_1))*(inv(MAT_1))*MAT_2)

    return dℰdk

end
