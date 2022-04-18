
function sigmoid(x, a1, a2, a3)
    return a1 ./ (1.0 .+ exp.(a2 .* (x .- a3)))
end

function drivesTanh(x, a1, a2, a3,a4)
    return a1*tanh(a2*(x-a3)) + a4
end

function findSteadyState(X1,X2,X3,X4,X5,X6,X7,X8)
        return nlsolve(f!, [ X1;X2;X3;X4;X5;X6;X7;X8]).zero
end

function f(type,x,C1,C2,C3)
        if type == "sigmoid"
                return (C1)/(1.0+exp(-C2*(x - C3)))
        elseif type == "tanh"
                return C1*tanh(C2*(x-C3))
        elseif type == "constant"
                return C1*(x-C3)
        else
                error("type must be 'sigmoid', 'tanh' or 'constant'")
        end
end

function init_conds(nterms,kc,phi_vec,c_vec,X,Y)
    p = 0
    for i = 1:nterms
        k1,k2 = kc*[cos(phi_vec[i]) sin(phi_vec[i]) ]
        p = p .+ c_vec[i]*exp.(im*k1*X .+ im*k2*Y) + conj(c_vec[i])*exp.(-im*k1*X .- im*k2*Y)
    end
    return p
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

        F[2] = (gaps_term_VE + gEE*(VsynEE - x[2]) + gEI*(VsynEI - x[2]) - (τE^2)*(pi^2) * (x[1]^2) +  x[2]^2 + ηE)

        F[3] =(gaps_term_RI -gII*x[3] - gIE*x[3] + 2 * x[3] * x[4] + (ΔI / (τI*pi)))

        F[4] =(gaps_term_VI + gIE*(VsynIE -x[4]) + gII*(VsynII - x[4]) - (pi^2)*(τI^2) * (x[3]^2) + x[4]^2 + ηI)

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

    F[1]=-gEE*rE - gEI*rE - κVEE*rE - κVEI*rE + 2*rE*vE + (ΔE/(τE*pi))

    F[2]=κVEI*(vI-vE) + gEE*(VsynEE-vE)+gEI*(VsynEI-vE)-(τE^2)*(pi^2)*(rE^2)+(vE^2)+ηE

    F[3]=-gII*rI - gIE*rI - κVIE*rI - κVII*rI + 2*rI*vI + (ΔI/(τI*pi))

    F[4]=κVIE*(vE-vI) + gIE*(VsynIE-vI)+gII*(VsynII-vI)-(τI^2)*(pi^2)*(rI^2)+(vI^2)+ηI

end

function SteadyStateNoPotassium()
    X = rand(4)
    SS = nlsolve(f_nopotassium!, [ X[1];X[2];X[3];X[4]],ftol=1e-12,iterations=50000)
    while SS.f_converged == false
        X = rand(4)
        SS = nlsolve(f_nopotassium!, [ X[1];X[2];X[3];X[4]],ftol=1e-12,iterations=50000)
    end
    return SS.zero
end

function SteadyStatePotassium(gaps_option)
    X = rand(4)
    SS = nlsolve(f_potassium!, [ X[1];X[2];X[3];X[4]])
    while SS.f_converged == false
        X = rand(4)
        SS = nlsolve(f_potassium!, [ X[1];X[2];X[3];X[4]])
    end
    return SS.zero
end

#spatial and time parameters and set up
function attach_ss!(ssMat,u0,a,b,Add_Potassium,X,X_space,dimension)
    wp = 1
    RE1 = ssMat[a,1]
    VE1 = ssMat[a,2]
    RI1 = ssMat[a,3]
    VI1 = ssMat[a,4]
    RE2 = ssMat[b,1]
    VE2 = ssMat[b,2]
    RI2 = ssMat[b,3]
    VI2 = ssMat[b,4]
    K1 = (1/δ)*fK(RE1+RI1,A1,A2,A3)
    K2 = (1/δ)*fK(RE2+RI2,A1,A2,A3)
    difference = (ssMat[b,:] .- ssMat[a,:])

    differenceP = K2-K1
    differenceηE =  fηE(K2,B1,B2,B3) - fηE(K1,B1,B2,B3)
    differenceηI =  fηI(K2,C1,C2,C3) - fηI(K1,C1,C2,C3)
    differenceDII =  fκV(K2,D1,D2,D3) - fκV(K1,D1,D2,D3)
    if dimension == 1
        ff  = exp.(-(X_space.^2)/wp)
    elseif dimension == 2
        ff = reshape(exp.(-(X_space.^2/wp .+ X_space'.^2/wp)),X)
    end
    u0[1:X] .= RE1 .+ difference[1].*ff
    u0[X+1:2X] .= VE1 .+  difference[2].*ff
    u0[2X+1:3X] .=  RI1 .+ difference[3].*ff
    u0[3X+1:4X] .=  VI1 .+ difference[4].*ff
    u0[4X+1:5X]  .=  RE1 .+ difference[1].*ff
    u0[5X+1:6X]   .=  RE1/σEE .+ difference[1].*ff/σEE
    u0[6X+1:7X]   .= RI1 .+ difference[3].*ff
    u0[7X+1:8X]  .= RI1/σEI .+  difference[3].*ff/σEI
    u0[8X+1:9X]  .= RE1 .+  difference[1].*ff
    u0[9X+1:10X]  .= RE1/σIE .+ difference[1].*ff/σIE
    u0[10X+1:11X]  .= RI1 .+  difference[3].*ff
    u0[11X+1:12X]  .= RI1/σII .+  difference[3].*ff/σII
    u0[12X+1:13X]  .=κSEE*RE1 .+  κSEE*difference[1].*ff
    u0[13X+1:14X] .=  κSEE*RE1 .+ κSEE*difference[1].*ff
    u0[14X+1:15X] .= κSEI*RI1 .+  κSEI*difference[3].*ff
    u0[15X+1:16X] .= κSEI*RI1 .+  κSEI*difference[3].*ff
    u0[16X+1:17X] .=  κSIE*RE1 .+ κSIE*difference[1].*ff
    u0[17X+1:18X] .=κSIE*RE1 .+  κSIE*difference[1].*ff
    u0[18X+1:19X] .= κSII*RI1 .+  κSII*difference[3].*ff
    u0[19X+1:20X] .= κSII*RI1 .+  κSII*difference[3].*ff
    if Add_Potassium == 1
        u0[20X+1:21X] = (1/δ)*fK(RE1+RI1,A1,A2,A3) .+ differenceP.*ff
        u0[21X+1:22X] = fηE(K1,B1,B2,B3) .+ differenceηE.*ff
        u0[22X+1:23X] = fηI(K1,C1,C2,C3) .+ differenceηI.*ff
        u0[23X+1:24X] = fκV(K1,D1,D2,D3) .+ differenceDII.*ff
    end
    return u0
end

function find_steady_state(p)
    ΔE, ΔI, κVEE,κVEI, κVIE,κVII,ηE,ηI,τE,τI,
    σEE, σEI,σIE, σII, v, κSEE,κSEI,κSIE, κSII,
    αEE, αEI, αIE, αII, δ,
    A1, A2, A3, B1, B2, B3, C1, C2, C3,
    D1, D2, D3, βK, βηE, βηI,βκV,∇,∇ψ,options=p
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
        steadyStateTest = SteadyStatePotassium(gaps_option)'
        while steadyStateTest[1] < 0 || steadyStateTest[3] < 0
            steadyStateTest =  SteadyStatePotassium(gaps_option)'
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


function guassian_drive(x,t,T,width,spread)
    return exp.(-(x.^2)./width).*exp(-((T-t)^2)/spread)
end



function setup(p,use_rand_init_conds,restard_solve,X)
        ΔE, ΔI, κVEE,κVEI, κVIE,κVII,ηE,ηI,τE,τI,
        σEE, σEI,σIE, σII, v, κSEE,κSEI,κSIE, κSII,
        αEE, αEI, αIE, αII, δ,
        A1, A2, A3, B1, B2, B3, C1, C2, C3,
        D1, D2, D3, βK, βηE, βηI,βκV,∇,∇ψ,options = p
        Add_Potassium,gaps_option,dimension = options
    ssMat = find_steady_state(p) #find steady state(s)
    ssMat = ssMat[sortperm(ssMat[:,1]),:] #arrang from lowest RE to highest RE
    numSS = size(ssMat,1) #number of steady states
    println("",numSS," STEADY STATE(S) FOUND!!!")

    for i = 1:numSS
        RE_ss = ssMat[:,1]
        VE_ss = ssMat[:,2]
        RI_ss = ssMat[:,3]
        VI_ss = ssMat[:,4]
        if Add_Potassium == 0
            println("[",i,"]:RE_ss = ",RE_ss[i]," || VE_ss = ",VE_ss[i]," || RI_ss = ",RI_ss[i]," || VI_ss = ",VI_ss[i])

        elseif Add_Potassium == 1
            K_ss = (1/δ)*fK.(RE_ss+RI_ss,A1,A2,A3)
            ηE_ss = fηE.(K_ss,B1,B2,B3)
            ηI_ss = fηI.(K_ss,C1,C2,C3)
            κV_ss = fκV.(K_ss,D1,D2,D3)
            println("[",i,"]: RE_ss = ",ssMat[i,1]," || VE_ss = ",ssMat[i,2]," || RI_ss = ",ssMat[i,3]," || VI_ss ",ssMat[i,4],
            "\n       K_ss ",K_ss[i]," || ηE_ss =  ",ηE_ss[i]," || ηI_ss ",ηI_ss[i]," || κV_ss = ",κV_ss[i])

        end
    end

    if numSS > 1 && attached_steady_state == 0
    print("\n Choose Steady State to use: (1 - ", numSS,")\n")
    chooseSS = Int(parse(Float64,readline()))
    Steady_State = ssMat[chooseSS,:]
    else
        Steady_State =ssMat[1,:]
    end


    RE_ss = Steady_State[1]
    VE_ss = Steady_State[2]
    RI_ss = Steady_State[3]
    VI_ss = Steady_State[4]
    if Add_Potassium == 1
            K_ss = (1/δ)*fK(RE_ss+RI_ss,A1,A2,A3)
            ηE_ss = fηE(K_ss,B1,B2,B3)
            ηI_ss = fηI(K_ss,C1,C2,C3)
            κV_ss = fκV(K_ss,D1,D2,D3)
    end
    gEE_ss,gEI_ss,gIE_ss,gII_ss = κSEE*RE_ss,κSEI*RI_ss,κSIE*RE_ss,κSII*RI_ss
    ψEE_ss,ψEI_ss,ψIE_ss,ψII_ss = RE_ss,RI_ss,RE_ss,RI_ss

    #Initialise State Vairables

    RE0 = zeros(X)
    VE0=zeros(X)
    RI0 = zeros(X)
    VI0=zeros(X)

    ψEE0 = zeros(X)
    AEE0 = zeros(X)
    ψEI0 = zeros(X)
    AEI0 = zeros(X)
    ψIE0 = zeros(X)
    AIE0 = zeros(X)
    ψII0 = zeros(X)
    AII0 = zeros(X)

    gEE0 = zeros(X)
    pEE0 = zeros(X)
    gEI0 = zeros(X)
    pEI0 = zeros(X)
    gIE0 = zeros(X)
    pIE0 = zeros(X)
    gII0 = zeros(X)
    pII0 = zeros(X)

    K0 = zeros(X)
    ηE0 = zeros(X)
    ηI0 = zeros(X)
    κV0 = zeros(X)

    # # # # # # # # # # # # # # # # # # # #

    # create initial conditions (if perturbing the steady state.. )
    RE0      .= RE_ss .+ perturb
    VE0      .= VE_ss .+ perturb
    RI0      .= RI_ss .+ perturb
    VI0      .= VI_ss .+ perturb

    ψEE0    .= ψEE_ss
    AEE0      .= ψEE_ss/σEE
    ψEI0    .= ψEI_ss
    AEI0      .= ψEI_ss/σEI
    ψIE0    .= ψIE_ss
    AIE0      .= ψIE_ss/σIE
    ψII0    .= ψII_ss
    AII0      .= ψII_ss/σII

    gEE0     .= gEE_ss .+ 1*perturb
    pEE0     .= gEE_ss .+ 1*perturb
    gEI0     .= gEI_ss .+ 1*perturb
    pEI0     .= gEI_ss .+ 1*perturb
    gIE0     .= gIE_ss .+ 1*perturb
    pIE0     .= gIE_ss .+ 1*perturb
    gII0     .= gII_ss .+ 1*perturb
    pII0     .= gII_ss .+ 1*perturb

    if Add_Potassium == 1
        K0       .= K_ss .+ 1*perturb
        ηE0    .= ηE_ss .+ 1*perturb
        ηI0    .= ηI_ss .+ 1*perturb
        κV0    .= κV_ss .+ 1*perturb
    end

    # # # # # # # # # # # # # # # # # # #
    if Add_Potassium == 0
        u0 = zeros(20 * X)
    elseif Add_Potassium == 1
        u0 = zeros(24 * X)
    end

    if restart_solve == 0
    u0[1:X]           = RE0
    u0[X + 1:2 * X]   = VE0
    u0[2*X + 1:3*X]   = RI0
    u0[3*X + 1:4*X]   = VI0
    u0[4*X + 1:5*X]   = ψEE0
    u0[5*X + 1:6*X]   = AEE0
    u0[6*X + 1:7*X]   = ψEI0
    u0[7*X + 1:8*X]   = AEI0
    u0[8*X + 1:9*X]   = ψIE0
    u0[9*X + 1:10*X]  = AIE0
    u0[10*X + 1:11*X] =  ψII0
    u0[11*X + 1:12*X] = AII0
    u0[12*X + 1:13*X] = gEE0
    u0[13*X + 1:14*X] = pEE0
    u0[14*X + 1:15*X] = gEI0
    u0[15*X + 1:16*X] = pEI0
    u0[16*X + 1:17*X] = gIE0
    u0[17*X + 1:18*X] = pIE0
    u0[18*X + 1:19*X] = gII0
    u0[19*X + 1:20*X] = pII0
    if Add_Potassium == 1
        u0[20*X + 1:21*X] = K0
        u0[21*X + 1:22*X] = ηE0
        u0[22*X + 1:23*X] = ηI0
        u0[23*X + 1:24*X] = κV0
    end


    if use_rand_init_conds == 1
        u0 += init_conditions(u0,X,Add_Potassium)
    end


    elseif restart_solve == 1
        if load_data == 1
            u0 = load("/home/james/PhD_Work/Julia_Code/Neural_Field_with_potassium/with_potassium/PDE_allgaps/2_Pop_1D/data/turinghopf_2D.jld","turing_hopf_end")
        else
            u0 = sol[:,end]
        end
    end


    return u0,ssMat

    
end
