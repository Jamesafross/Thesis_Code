include("functions.jl")

function  dηABdω(λ,αAB) 
    return (-2*im*αAB^2)/((1+λ)^3)
end

function W_AB(v,λ,k,σAB,κSAB) 
    return (κSAB*v*(λ*σAB + v))/((v + σAB*(λ - im*k*v))*(v + σAB*(λ + im*k*v)))
end

function ηAB(λ,αAB) 
    return (1+((λ)/αAB))^-2
end

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

function dFABdk(λ,k,σAB,κSAB,v)
    tFrac = (-2*k*κSAB*σAB^2*v^3)*(v+σAB*λ)
    bFrac = ((v+σAB*(λ-im*k*v))^2)*((v+σAB*(λ+im*k*v))^2)
    return tFrac/bFrac
end

function dFABdω(λ,k,αAB,σAB,κSAB,v)
    return ηAB(λ,αAB)*dW_ABdω(v,λ,k,σAB,κSAB) + dηABdω(λ,αAB)*W_AB(v,λ,k,σAB,κSAB)
end 

function dAdω(λ,k,RbarE,RbarI,p)
    dAEEdω = -RbarE*dFABdω(λ,k,αEE,σEE,κSEE,v)
    dAEIdω = -RbarE*dFABdω(λ,k,αEI,σEI,κSEI,v)
    dAIEdω = -RbarI*dFABdω(λ,k,αIE,σIE,κSIE,v)
    dAIIdω = -RbarI*dFABdω(λ,k,αII,σII,κSII,v)
    return [dAEEdω  dAEIdω;dAIEdω dAIIdω   ]
end

function dBdω()
    return zeros(2,2)
end

function dCdω(λ,k,VbarE,VbarI,p)
    dCEEdω = (VsynEE - VbarE)*dFABdω(λ,k,αEE,σEE,κSEE,v)
    dCEIdω = (VsynEI - VbarE)*dFABdω(λ,k,αEI,σEI,κSEI,v)
    dCIEdω = (VsynIE - VbarI)*dFABdω(λ,k,αIE,σIE,κSIE,v)
    dCIIdω = (VsynII - VbarI)*dFABdω(λ,k,αII,σII,κSII,v)
    return [dCEEdω  dCEIdω;dCIEdω dCIIdω   ]
end

function dDdω()
    return zeros(2,2)
end


function dAdk(λ,k,RbarE,RbarI,p)
    dAEEdk = -RbarE*dFABdk(λ,k,σEE,κSEE,v)
    dAEIdk = -RbarE*dFABdk(λ,k,σEI,κSEI,v)
    dAIEdk = -RbarI*dFABdk(λ,k,σIE,κSIE,v)
    dAIIdk = -RbarI*dFABdk(λ,k,σII,κSII,v)
    return [dAEEdk  dAEIdk;dAIEdk dAIIdk   ]
end

function dBdk()
    return zeros(2,2)
end

function dCdk(λ,k,VbarE,VbarI,p)
    dCEEdk = (VsynEE - VbarE)*dFABdk(λ,k,σEE,κSEE,v)
    dCEIdk = (VsynEI - VbarE)*dFABdk(λ,k,σEI,κSEI,v)
    dCIEdk = (VsynIE - VbarI)*dFABdk(λ,k,σIE,κSIE,v)
    dCIIdk = (VsynII - VbarI)*dFABdk(λ,k,σII,κSII,v)
    return [dCEEdk  dCEIdk;dCIEdk dCIIdk   ]
end

function dDdk()
    return zeros(2,2)
end


function dEdk(gaps_option)
    if gaps_option == "ALL-ALL" || gaps_option == "I-I"
        return zeros(2,4)
    elseif gaps_option == "OFF"
        return zeros(2,3)
    else
        @error "gaps_option must be I-I, ALL-ALL or OFF"
    end
end

dEdω() =  dEdk(gaps_option)


function dFdk(gaps_option)
    if gaps_option == "ALL-ALL" || gaps_option == "I-I" 
        return zeros(2,4)
    elseif gaps_option == "OFF"
        return  zeros(2,3)
    else
        @error "gaps_option must be I-I, ALL-ALL or OFF"
    end
end

dFdω() =  dFdk(gaps_option)


function dGdω(gaps_option)
    if gaps_option == "ALL-ALL" || gaps_option == "I-I"
        return zeros(4,8)
    elseif gaps_option == "OFF"
        return zeros(3,7)
    else
        @error "gaps_option must be I-I, ALL-ALL or OFF"
    end
end

function dGdk(k,p,gaps_option)

    ΔE, ΔI, κVEE,κVEI,κVIE,κVII,ηE,ηI,τE,τI,
    σEE,σEI,σIE, σII, v, κSEE,κSEI,κSIE, κSII,
    αEE, αEI, αIE, αII,τ0EE,τ0EI,τ0IE,τ0II, δ,
    A1, A2, A3, A4,B1, B2, B3, C1, C2, C3,
    D1, D2, D3,βK,βηE,βηI,βκV = p

    if gaps_option == "ALL-ALL"
        return [0. 0. 0. 0. -2*A4*k 0. 0. 0.;
        0. 0. 0. 0. 0. 0. 0. 0. ;
        0. 0. 0. 0. 0.  0. 0. 0. ;
        0. 0. 0. 0. 0. 0. 0. 0. ]
    elseif gaps_option == "I-I"
        return [0. 0. 0. 0. -2*A4*k 0. 0. 0.;
        0. 0. 0. 0. 0. 0. 0. 0. ;
        0. 0. 0. 0. 0.  0. 0. 0. ;
        0. 0. 0. 0. 0. 0. 0. 0. ]
    elseif gaps_option == "OFF"
        return [0. 0. 0. 0. -2*A4*k 0. 0.;
        0. 0. 0. 0. 0. 0. 0. ;
        0. 0. 0. 0. 0.  0. 0. ]
    else
        @error "gaps_option must be I-I, ALL-ALL or OFF"
    end
end



function dεdωfunc(V,R,κbarV,k,μ,ω,p,options,derivatives)

    ΔE, ΔI, κVEE,κVEI,κVIE,κVII,ηE,ηI,τE,τI,
    σEE,σEI,σIE, σII, v, κSEE,κSEI,κSIE, κSII,
    αEE, αEI, αIE, αII,τ0EE,τ0EI,τ0IE,τ0II, δ,
    A1, A2, A3, A4,B1, B2, B3, C1, C2, C3,
    D1, D2, D3,βK,βηE,βηI,βκV = p
    Add_Potassium, gaps_option = options
    if Add_Potassium == 1
    d_fK_dRE,d_fK_dRI,d_fK_dVE,d_fK_dVI,d_fηE_dK,d_fηI_dK,
    d_fκV_dK = derivatives 
    end
    λ = μ + im*ω

    M1 = JacobianMatrix(V,R,κbarV,k,μ,ω,p,options,derivatives)
    τMat = [τE 0. ; 0. τI]
    
    M2 = [dAdω(λ,k,RbarE,RbarI,p).- τMat.*im dBdω() dEdω();
    dCdω(λ,k,VbarE,VbarI,p) dDdω().-τMat*im dFdω() ;
     dGdω(gaps_option) .- β_mat(p,gaps_option).*im]

    tr((det(M1))*(inv(M1))*M2)
end

function dεdkfunc(V,R,κbarV,k,μ,ω,p,options,derivatives)

    ΔE, ΔI, κVEE,κVEI,κVIE,κVII,ηE,ηI,τE,τI,
    σEE,σEI,σIE, σII, v, κSEE,κSEI,κSIE, κSII,
    αEE, αEI, αIE, αII,τ0EE,τ0EI,τ0IE,τ0II, δ,
    A1, A2, A3, A4,B1, B2, B3, C1, C2, C3,
    D1, D2, D3,βK,βηE,βηI,βκV = p
    Add_Potassium, gaps_option = options
    if Add_Potassium == 1
    d_fK_dRE,d_fK_dRI,d_fK_dVE,d_fK_dVI,d_fηE_dK,d_fηI_dK,
    d_fκV_dK = derivatives 
    end

    M1 = JacobianMatrix(V,R,κbarV,k,μ,ω,p,options,derivatives)
    λ = μ + im*ω

    
    M2 = [dAdk(λ,k,RbarE,RbarI,p)  dBdk() dEdk(gaps_option);
    dCdk(λ,k,VbarE,VbarI,p) dDdk()  dFdk(gaps_option);
     dGdk(k,p,gaps_option) ]

    tr((det(M1))*(inv(M1))*M2)
end


function IFT(V,R,κbarV,k,μ,ω,p,options,derivatives)
    dϵdω =  dεdωfunc(V,R,κbarV,k,μ,ω,p,options,derivatives)
    dϵdk =  dεdkfunc(V,R,κbarV,k,μ,ω,p,options,derivatives)

    return imag(dϵdω)*real(dϵdk)  - real(dϵdω)*imag(dϵdk)
end

function bif!(F,x)
    p = ΔE, ΔI, κVEE,κVEI,κVIE,κVII,ηE,ηI,τE,τI,
    σEE,σEI,σIE, σII, v, x[2],x[2]*(4.6/5.0),x[2]*(4.6/5.0), x[2]*(3.1/5.0),
    αEE, αEI, αIE, αII,τ0EE,τ0EI,τ0IE,τ0II, δ,
    A1, A2, A3, A4,B1, B2, B3, C1, C2, C3,
    D1, D2, D3,βK,βηE,βηI,βκV

    rE,vE,rI,vI = find_steady_state(p,Add_Potassium,gaps_option)
    K = (1/δ)*fK(rE+rI,A1,A2,A3)
    κV = fκV(K,D1,D2,D3)
    F[1] = IFT([vE,vI],[rE,rI],κV,x[3],0.0,x[1],p,options,derivatives)
    detM = det(JacobianMatrix([vE,vI],[rE,rI],κV,x[3],0.0,x[1],p,options,derivatives))
    F[2] = real(detM)
    F[3] = imag(detM)

    
end






dεdkfunc(V,R,κbarV,0,0,0,p,options,derivatives)
dεdωfunc(V,R,κbarV,2,2,-1,p,options,derivatives)


X1,X2,X3 = randn(3)
while nlsolve(bif!, [X1;X2;X3]).f_converged == false 
    X1,X2,X3 = X1+.1randn(),X2+.1randn(),X2+.1randn()
end
nlsolve(bif!, [X1;X2;X3])