
function  W_AB(v,λ,k,σAB,κSAB,dimensions)
    if dimensions == 1
        return (κSAB*v*(λ*σAB + v))/((v + σAB*(λ - im*k*v))*(v + σAB*(λ + im*k*v)))
    elseif dimensions == 2
        return (κSAB*v*(λ*σAB + v))/(pi*σAB*((v + σAB*(λ - im*k*v))*(v + σAB*(λ + im*k*v))))
    end
end

function ηAB(λ,α_ab) 
    return (1+((λ)/α_ab))^-2
end

function A_mat(SteadyState,FAB,p,Add_Potassium,gaps_option)

    ΔE, ΔI, κVEE,κVEI,κVIE,κVII,ηE,ηI,τE,τI,
    σEE,σEI,σIE, σII, v, κSEE,κSEI,κSIE, κSII,
    αEE, αEI, αIE, αII,τ0EE,τ0EI,τ0IE,τ0I, δ,
    A1, A2, A3, A4,B1, B2, B3, C1, C2, C3,
    D1, D2, D3,βK,βηE,βηI,βκV= p

    VbarE,VbarI,RbarE,RbarI,κbarV = SteadyState
    gbarEE,gbarEI,gbarIE,gbarII = κSEE*RbarE,κSEI*RbarI,κSIE*RbarE,κSII*RbarI

    FEE,FEI,FIE,FII = FAB

    if Add_Potassium == 0 
        AEE = 2*VbarE - (gbarEE + κVEE) - (gbarEI + κVEI) - RbarE*FEE
        AEI = -RbarE*FEI
        AIE = -RbarI*FIE
        AII = 2*VbarI - (gbarIE + κVIE) - (gbarII + κVII) - RbarI*FII
    elseif Add_Potassium == 1
        if gaps_option == "ALL-ALL"
            AEE = 2*VbarE-(gbarEE+κbarV)-(gbarEI+κbarV)-RbarE*FEE
            AEI = -RbarE*FEI
            AIE = -RbarI*FIE
            AII = 2*VbarI-(gbarIE+κbarV)-(gbarII-κbarV)-RbarI*FII
        elseif gaps_options == "I-I"
            AEE = 2*VbarE-(gbarEE)-(gbarEI)-RbarE*FEE
            AEI = -RbarE*FEI
            AIE = -RbarI*FIE
            AII = 2*VbarI-(gbarIE)-(gbarII-κbarV)-RbarI*FII
        elseif gaps_option == "OFF"
            AEE = 2*VbarE-(gbarEE)-(gbarEI)-RbarE*FEE
            AEI = -RbarE*FEI
            AIE = -RbarI*FIE
            AII = 2*VbarI-(gbarIE)-(gbarII)-RbarI*FII
        else
            @error "gaps_option must be I-I, ALL-ALL or OFF"
        end
    else
        @error "Add_Pottasium must be 0 or 1"
    end

    return [AEE AEI; AIE AII]
end

function B_mat(SteadyState)
    VbarE,VbarI,RbarE,RbarI,κbarV = SteadyState
    BEE = 2*RbarE
    BEI = 0
    BIE = 0
    BII = 2*RbarI
    return [BEE BEI; BIE BII]
end


function C_mat(SteadyState,FAB,p)
    ΔE, ΔI, κVEE,κVEI,κVIE,κVII,ηE,ηI,τE,τI,
    σEE,σEI,σIE, σII, v, κSEE,κSEI,κSIE, κSII,
    αEE, αEI, αIE, αII,τ0EE,τ0EI,τ0IE,τ0I, δ,
    A1, A2, A3, A4,B1, B2, B3, C1, C2, C3,
    D1, D2, D3,βK,βηE,βηI,βκV= p

    VbarE,VbarI,RbarE,RbarI,κbarV = SteadyState
    gbarEE,gbarEI,gbarIE,gbarII = κSEE*RbarE,κSEI*RbarI,κSIE*RbarE,κSII*RbarI

    FEE,FEI,FIE,FII = FAB
   

    CEE = -2.0*(τE^2.0)*(pi^2.0)*RbarE + (VsynEE - VbarE)*FEE
    CEI = (VsynEI - VbarE)*FEI
    CIE = (VsynIE - VbarI)*FIE
    CII = -2.0*(τI^2.0)*(pi^2.0)*RbarI + (VsynII - VbarI)*FII

    return [CEE CEI; CIE CII]


end

function D_mat(SteadyState,p,Add_Potassium,gaps_option)
    VbarE,VbarI,RbarE,RbarI,κbarV = SteadyState
    gbarEE,gbarEI,gbarIE,gbarII = κSEE*RbarE,κSEI*RbarI,κSIE*RbarE,κSII*RbarI

    ΔE, ΔI, κVEE,κVEI,κVIE,κVII,ηE,ηI,τE,τI,
    σEE,σEI,σIE, σII, v, κSEE,κSEI,κSIE, κSII,
    αEE, αEI, αIE, αII,τ0EE,τ0EI,τ0IE,τ0I, δ,
    A1, A2, A3, A4,B1, B2, B3, C1, C2, C3,
    D1, D2, D3,βK,βηE,βηI,βκV= p

    if Add_Potassium == 0 
        DEE = 2.0*VbarE - (κVEE + gbarEE) - (κVEI + gbarEI) + κVEE
        DEI = κVEI
        DIE = κVIE
        DII = 2.0*VbarI - (κVIE + gbarIE) - (κVII + gbarII) + κVII
    elseif Add_Potassium == 1
        if gaps_option == "ALL-ALL"
            DEE = 2.0*VbarE-κbarV-gbarEE-gbarEI ;
            DEI = κbarV
            DIE = κbarV
            DII = 2.0*VbarI-κbarV-gbarIE-gbarII
        elseif gaps_option == "I-I"
            DEE = 2.0*VbarE-gbarEE-gbarEI ;
            DEI = 0.
            DIE = 0.
            DII = 2.0*VbarI-gbarIE-gbarII
        elseif gaps_options == "OFF"
            DEE = 2.0*VbarE-gbarEE-gbarEI ;
            DEI = 0.
            DIE = 0.
            DII = 2.0*VbarI-gbarIE-gbarII
        else
            @error "gaps_option must be I-I, ALL-ALL or OFF"
        end
    else
        @error "Add_Pottasium must be 0 or 1"
    end
    return [DEE DEI; DIE DII]
end

function E_mat(SteadyState,gaps_options)
    VbarE,VbarI,RbarE,RbarI,κbarV = SteadyState
    if gaps_option == "ALL-ALL"
        return [0. 0. 0. -2.0*RbarE ; 0. 0. 0. -2.0*RbarI ]
    elseif gaps_option == "I-I"
        return [0. 0. 0. 0. ; 0. 0. 0. -RbarI ]
    elseif gaps_options == "OFF"
        return [0. 0. 0. ; 0. 0. 0.]
    else
        @error "gaps_option must be I-I, ALL-ALL or OFF"
    end
end

function F_mat(SteadyState,gaps_option)
    VbarE,VbarI,RbarE,RbarI,κbarV = SteadyState
    if gaps_option == "ALL-ALL"
        return [0. 1. 0. VbarI-VbarE; 0. 0. 1. VbarE-VbarI ]
    elseif gaps_option == "I-I"
        return  [0. 1. 0. 0. ; 0. 0. 1. 0. ]
    elseif gaps_options == "OFF"
        return  [0. 1. 0.; 0. 0. 1.]
    else
        @error "gaps_option must be I-I, ALL-ALL or OFF"
    end
end

function G_mat(derivatives,k,p)
    d_fK_dRE,d_fK_dRI,d_fK_dVE,d_fK_dVI,d_fηE_dK,d_fηI_dK,
    d_κV_dK = derivatives 

    ΔE, ΔI, κVEE,κVEI,κVIE,κVII,ηE,ηI,τE,τI,
    σEE,σEI,σIE, σII, v, κSEE,κSEI,κSIE, κSII,
    αEE, αEI, αIE, αII,τ0EE,τ0EI,τ0IE,τ0I, δ,
    A1, A2, A3, A4,B1, B2, B3, C1, C2, C3,
    D1, D2, D3,βK,βηE,βηI,βκV= p

    if gaps_option == "ALL-ALL"
        return [d_fK_dRE d_fK_dRI d_fK_dVE d_fK_dVI -δ-A4*k^2. 0. 0. 0.;
        0. 0. 0. 0. d_fηE_dK -1. 0. 0. ;
        0. 0. 0. 0. d_fηI_dK  0. -1. 0. ;
        0. 0. 0. 0. d_κV_dK  0. 0. -1. ]
    elseif gaps_option == "I-I"
        return [d_fK_dRE d_fK_dRI d_fK_dVE d_fK_dVI -δ-A4*k^2. 0. 0. 0.;
        0. 0. 0. 0. d_fηE_dK -1. 0. 0. ;
        0. 0. 0. 0. d_fηI_dK  0. -1. 0. ;
        0. 0. 0. 0. d_κV_dK  0. 0. -1. ]
    elseif gaps_options == "OFF"
        return [d_fK_dRE d_fK_dRI d_fK_dVE d_fK_dVI -δ-A4*k^2. 0. 0. ;
                        0. 0. 0. 0. d_fηE_dK -1. 0.  ;
                        0. 0. 0. 0. d_fηI_dK  0. -1. ] 
    else
        @error "gaps_option must be I-I, ALL-ALL or OFF"
    end
end

function β_mat(p,gaps_option)
    ΔE, ΔI, κVEE,κVEI,κVIE,κVII,ηE,ηI,τE,τI,
    σEE,σEI,σIE, σII, v, κSEE,κSEI,κSIE, κSII,
    αEE, αEI, αIE, αII,τ0EE,τ0EI,τ0IE,τ0I, δ,
    A1, A2, A3, A4,B1, B2, B3, C1, C2, C3,
    D1, D2, D3,βK,βηE,βηI,βκV= p 
    if gaps_option == "ALL-ALL"
        return [0. 0. 0. 0. βK 0. 0. 0.;
                0. 0. 0. 0. 0. βηE 0. 0. ;
                0. 0. 0. 0. 0. 0. βηI 0. ;
                0. 0. 0. 0. 0. 0. 0. βκV ]
    elseif gaps_option == "I-I"
        return     [0. 0. 0. 0. βK 0. 0. 0.;
                    0. 0. 0. 0. 0. βηE 0. 0. ;
                    0. 0. 0. 0. 0. 0. βηI 0. ;
                    0. 0. 0. 0. 0. 0. 0. βκV ]
        return 
    elseif gaps_options == "OFF"
        return     [0. 0. 0. 0. βK 0. 0.;
                    0. 0. 0. 0. 0. βηE 0. ;
                    0. 0. 0. 0. 0. 0. βηI ]
    else
        @error "gaps_option must be I-I, ALL-ALL or OFF"
    end
end



function JacobianMatrix(V,R,κbarV,k,μ,ω,p,options,derivatives)

    ΔE, ΔI, κVEE,κVEI,κVIE,κVII,ηE,ηI,τE,τI,
    σEE,σEI,σIE, σII, v, κSEE,κSEI,κSIE, κSII,
    αEE, αEI, αIE, αII,τ0EE,τ0EI,τ0IE,τ0I = p
    d_fK_dRE,d_fK_dRI,d_fK_dVE,d_fK_dVI,d_fηE_dK,d_fηI_dK,
    d_κV_dK = derivatives 
    options = Add_Potassium, gaps_option
    
    λ= μ+im*ω
    VbarE,VbarI,RbarE,RbarI = V[1],V[2],R[1],R[2]
    SteadyState = VbarE,VbarI,RbarE,RbarI,κbarV
    #println(dimensions)
    FEE = W_AB(v,λ,k,σEE,κSEE,dimensions)*ηAB(λ,αEE)
    FEI = W_AB(v,λ,k,σEI,κSEI,dimensions)*ηAB(λ,αEI)
    FIE = W_AB(v,λ,k,σIE,κSIE,dimensions)*ηAB(λ,αIE)
    FII = W_AB(v,λ,k,σII,κSII,dimensions)*ηAB(λ,αII)

    FAB = FEE,FEI,FIE,FII
    τMAT = [τE 0; 0 τI]

    A = A_mat(SteadyState,FAB,p,Add_Potassium,gaps_option)
    B = B_mat(SteadyState)
    C = C_mat(SteadyState,FAB,p)
    D = D_mat(SteadyState,p,Add_Potassium,gaps_option)
    if Add_Potassium == 0 
        return [A .- τMat.*λ B;  C D .- τMat.*λ]
    elseif Add_Potassium == 1
       E =  E_mat(SteadyState,gaps_options)
       F = F_mat(SteadyState,gaps_option)
       G = G_mat(derivatives,k,p)
       βMat =  β_mat(p,gaps_option)

       return [A.- τMat.*λ B E;C D.-λ.*τMat F; G .- βMat.*λ]
    else
        @error "Add_Pottasium must be 0 or 1"
    end
end

function JacobianMatrix_Potassium(μ,ω,k,Rbar,Vbar,κbarV,
    d_fK_dRE,d_fK_dRI,d_fK_dVE,d_fK_dVI,d_fηE_dK,d_fηI_dK,
    d_κV_dK,τ0,p)

    ΔE, ΔI, κVEE,κVEI,κVIE,κVII,ηE,ηI,τE,τI,
    σEE,σEI,σIE, σII, v, κSEE,κSEI,κSIE, κSII,
    αEE, αEI, αIE, αII,τ0EE,τ0EI,τ0IE,τ0II, δ,
    A1, A2, A3, A4,B1, B2, B3, C1, C2, C3,
    D1, D2, D3,βK,βηE,βηI,βκV = p
        ##MAKE THE MATRIX FOR SPECTRUM##
        λ = μ+im*ω
        τE,τI = τ
      
        RbarE,RbarI = Rbar
        VbarE,VbarI = Vbar
    
        τ0EE,τ0EI,τ0IE,τ0II = τ0
        gbarEE = κSEE*RbarE
        gbarEI = κSEI*RbarI
        gbarIE = κSIE*RbarE
        gbarII = κSII*RbarI
        FEE = W_AB(v,λ,k,σEE,κSEE,dimensions)*ηAB(λ,αEE)
        FEI = W_AB(v,λ,k,σEI,κSEI,dimensions)*ηAB(λ,αEI)
        FIE = W_AB(v,λ,k,σIE,κSIE,dimensions)*ηAB(λ,αIE)
        FII = W_AB(v,λ,k,σII,κSII,dimensions)*ηAB(λ,αII)
    
    
        τ_mat = [τE 0; 0 τI]

        B_mat = [2*RbarE 0;0 2*RbarI]

        C_mat = [-2*(pi^2)*(τE^2)*RbarE+FEE*(VsynEE-VbarE) FEI*(VsynEI-VbarI);
        FIE*(VsynIE - VbarE) -2*(pi^2)*(τI^2)*RbarI+FII*(VsynII-VbarI) ]

        if gaps_option == "ALL-ALL"
            A_mat = [2*VbarE-(gbarEE+κbarV)-(gbarEI+κbarV)-RbarE*FEE -RbarE*FEI;
                    -RbarI*FIE 2*VbarI-(gbarIE+κbarV)-(gbarII-κbarV)-RbarI*FII]
            E_mat = [0 0 0 -2*RbarE ; 0 0 0 -2*RbarI ]
            D_mat = [2*VbarE-κbarV-gbarEE-gbarEI κbarV;
                        κbarV 2*VbarI-κbarV-gbarIE-gbarII]
            F_mat = [0 1 0 VbarI-VbarE ;
                     0 0 1 VbarE-VbarI ]
            G_mat = [d_fK_dRE d_fK_dRI d_fK_dVE d_fK_dVI -δ-A4*k^2 0 0 0;
                    0 0 0 0 d_fηE_dK -1 0 0 ;
                    0 0 0 0 d_fηI_dK  0 -1 0 ;
                    0 0 0 0 d_κV_dK  0 0 -1 ]
            β_mat = [0 0 0 0 βK 0 0 0;
                    0 0 0 0 0 βηE 0 0 ;
                    0 0 0 0 0 0 βηI 0 ;
                    0 0 0 0 0 0 0 βκV ]
        elseif gaps_option == "I-I"
            A_mat = [2*VbarE-(gbarEE)-(gbarEI)-RbarE*FEE -RbarE*FEI;
                    -RbarI*FIE 2*VbarI-(gbarIE)-(gbarII+κbarV)-RbarI*FII]
            E_mat = [0 0 0 0 ; 0 0 0 -RbarI ]
            D_mat = [2*VbarE-gbarEE-gbarEI  0 ;
                     0     2*VbarI-gbarIE-gbarII]
            F_mat = [0 1 0 0 ;
                    0 0 1 0 ]
            G_mat = [d_fK_dRE d_fK_dRI d_fK_dVE d_fK_dVI -δ-A4*k^2 0 0 0;
                     0 0 0 0 d_fηE_dK -1 0 0 ;
                     0 0 0 0 d_fηI_dK  0 -1 0 ;
                     0 0 0 0 d_κV_dK  0 0 -1 ]
            β_mat = [0 0 0 0 βK 0 0 0;
                    0 0 0 0 0 βηE 0 0 ;
                    0 0 0 0 0 0 βηI 0 ;
                    0 0 0 0 0 0 0 βκV ]
        elseif gaps_option == "OFF"
            A_mat = [2*VbarE-(gbarEE)-(gbarEI)-RbarE*FEE -RbarE*FEI;
                    -RbarI*FIE 2*VbarI-(gbarIE)-(gbarII)-RbarI*FII]
            E_mat = [0 0 0  ; 0 0 0]
            D_mat = [2*VbarE-gbarEE-gbarEI  0 ;
                     0     2*VbarI-gbarIE-gbarII]
            F_mat = [0 1 0  ;
                    0 0 1  ]
            G_mat = [d_fK_dRE d_fK_dRI d_fK_dVE d_fK_dVI -δ-A4*k^2 0 0 ;
                    0 0 0 0 d_fηE_dK -1 0  ;
                    0 0 0 0 d_fηI_dK  0 -1  ]
            β_mat = [0 0 0 0 βK 0 0 ;
                    0 0 0 0 0 βηE 0  ;
                    0 0 0 0 0 0 βηI ]
        end

      

        return [A_mat.-λ.*τ_mat B_mat E_mat;
                C_mat D_mat.-λ.*τ_mat F_mat;
                G_mat .- β_mat.*λ]
end






function g(F, x, V,R,κbarV,k,μ,ω,p,options,derivatives)
    detM = det(JacobianMatrix(V,R,κbarV,k,μ,ω,p,options,derivatives))
    F[1] = real(detM)
    F[2] = imag(detM)
end




function g_potassium!(F, x)
    detM = det(JacobianMatrix_Potassium(x[1],x[2],k,R,V,κbarV,
    d_fK_dRE,d_fK_dRI,d_fK_dVE,d_fK_dVI,d_fηE_dK,d_fηI_dK,
    d_fκV_dK,τ0,p))
    F[1] = real(detM)
    F[2] = imag(detM)
end



function f_nopotassium!(F, x)
    rE  = x[1]
    vE  = x[2]
    rI  = x[3]
    vI  = x[4]
    gEE = κSEE * rE
    gEI = κSEI * rI
    gIE = κSIE * rE
    gII = κSII * rI


    F[1]=(-gEE*rE - gEI*rE - κVEE*rE - κVEI*rE + 2*rE*vE + (ΔE/(τE*pi)))

    F[2]=(κVEI*(vI-vE) + gEE*(VsynEE-vE)+gEI*(VsynEI-vE)-(τE^2)*(pi^2)*(rE^2)+(vE^2)+ηE)

    F[3]=(-gII*rI - gIE*rI - κVIE*rI - κVII*rI + 2*rI*vI + (ΔI/(τI*pi)))

    F[4]=(κVIE*(vE-vI) + gIE*(VsynIE-vI)+gII*(VsynII-vI)-(τI^2)*(pi^2)*(rI^2)+(vI^2)+ηI)
end

function f_potassium!(F, x)
    gEE = κSEE*x[1]
    gEI = κSEI*x[3]
    gIE = κSIE*x[1]
    gII = κSII*x[3]
    K = (1/δ)*fK(x[1]+x[3],A1,A2,A3)
    ηE = fηE(K,B1,B2,B3)
    ηI = fηI(K,C1,C2,C3)
    κV = fκV(K,D1,D2,D3)
    if gaps_option == "ALL-ALL"
        gaps_term_RE = -κV.*x[1] .- κV .* x[1]
        gaps_term_RI = .-κV.*x[3] .- κV .* x[3]
        gaps_term_VE = κV.*(x[4]-x[2])
        gaps_term_VI =  κV.*(x[2]-x[4])
    elseif gaps_option == "I-I"
        gaps_term_RE = 0
        gaps_term_RI = -κV.*x[3]
        gaps_term_VE = 0
        gaps_term_VI = 0
    elseif gaps_option == "OFF"
        gaps_term_RE = 0
        gaps_term_RI = 0
        gaps_term_VE = 0
        gaps_term_VI = 0
    end
        F[1] = (gaps_term_RE -gEE*x[1] -gEI*x[1] +2 * x[1] * x[2] + (ΔE / (τE*pi)))

        F[2] = (gaps_term_VE + gEE*(VsynEE - x[2]) + gEI*(VsynEI - x[2]) - (τE^2)*(pi^2) * (x[1]^2) +  x[2]^2 + ηE )

        F[3] =(gaps_term_RI -gII*x[3] - gIE*x[3] + 2 * x[3] * x[4] + (ΔI / (τI*pi)))

        F[4] =(gaps_term_VI + gIE*(VsynIE -x[4]) + gII*(VsynII - x[4]) - (pi^2)*(τI^2) * (x[3]^2) + x[4]^2 + ηI )

end



function SteadyStateNoPotassium()
    σEE, σEI, σIE, σII = σ
    κVEE, κVEI, κVIE, κVII = κV
    κSEE, κSEI, κSIE, κSII = κS
    ηE, ηI = η0
    αEE, αEI, αIE, αII = α
    VsynEE, VsynEI, VsynIE, VsynII = Vsyn
    τE,τI = τ

    X = rand(4)
    SS = nlsolve(f_nopotassium!, [ X[1];X[2];X[3];X[4]])
    while SS.f_converged == false || ((SS.zero[1] < 0) == true)
        X = rand(4)
        SS = nlsolve(f_nopotassium!, [ X[1];X[2];X[3];X[4]])
    end
    return SS.zero
end

function SteadyStatePotassium()
    X = rand(4)
    SS = nlsolve(f_potassium!, [ X[1];X[2];X[3];X[4]])
    while SS.f_converged == false
        X = rand(4)
        SS = nlsolve(f_potassium!, [ X[1];X[2];X[3];X[4]])
    end
    return SS.zero
end




function F_Func(k,λ,τAB0,σAB,v,αAB,κSAB)
        ηFunc = (1+λ/αAB)^-2
        w_hat = (v*(λ*σAB + v))/(λ^2*σAB^2 + 2*λ*σAB*v + (1+k^2*σAB^2)*v^2)
        return exp(-λ*τAB0)*ηFunc*w_hat
end


function f(type,x,C1,C2,C3)
        if type == "sigmoid"
                return (C1)/(1+exp(-C2*(x - C3)))
        elseif type == "tanh"
                return C1*tanh(C2*(x-C3))
        elseif type == "constant"
                return C1*(x-C3)
        else
                error("type must be 'sigmoid', 'tanh' or 'constant'")
        end
end

function dfdx(type,x,C1,C2,C3)
        if type == "sigmoid"
                return ((C1*C2*(x)*exp(-C2*(C3+(x))))/((1+exp(-C2*(-C3+(x))))^2))
        elseif type == "tanh"
                return C1*C2*sech(C2*(-C3 + x))^2
        elseif type == "constant"
                return C1*x
        else
                error("type must be 'sigmoid', 'tanh' or 'constant'")
        end
end


function find_steady_state()
ssMat = []
for i = 1:1000
    global flagSS = false
    if Add_Potassium == 0
        steadyStateTest =  SteadyStateNoPotassium()'
        #print("\n",steadyStateTest,"\n")
        while steadyStateTest[1] < 0 || steadyStateTest[3] < 0
            steadyStateTest =  SteadyStateNoPotassium()'
        end
    elseif Add_Potassium == 1
        steadyStateTest = SteadyStatePotassium()'
        while steadyStateTest[1] < 0 || steadyStateTest[3] < 0
            steadyStateTest =  SteadyStatePotassium()'
        end
    end
    if i == 1
            ssMat = cat(ssMat,steadyStateTest,dims=1)
    else
        for j = 1:size(ssMat,1)
            if sum(abs.((steadyStateTest' .- ssMat[j,:]))) < 10^-1
                global flagSS = true
            end
        end
        if flagSS  == false
            if steadyStateTest[1] > 0
                if  steadyStateTest[3] > 0
                    ssMat =  cat(ssMat,steadyStateTest,dims=1)
                end
            end
        end
    end
end
return ssMat
end
