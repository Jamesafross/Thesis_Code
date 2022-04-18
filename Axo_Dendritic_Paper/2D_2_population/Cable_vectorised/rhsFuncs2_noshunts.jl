function z1(y, psi, v, V0, W, σ, k, h, type, beta, deltaF, Dx, Dzz,D,Dxx)
    U = y + (W/σ) .* firingrate.(V0 .- h, type, beta).*deltaF


    return  -v*( (1/σ).*U .+ k.*Dx*U .- Dzz * psi .-
            ( W*(1/σ)^2)*(firingrate.(V0 .- h, type, beta).*deltaF) .-
            ((W/σ))*k*Dx*(firingrate.(V0 .- h, type, beta).*deltaF))

end

function z1_noKAPPA(y, psi, v, V0, W, σ, k, h, type, beta, deltaF, Dx, Dzz,tau,D,Dxx)
    U = y + (W/σ) .* firingrate.(V0 .- h, type, beta).*deltaF
    return  -v*((1/σ).*U .- Dzz * psi .-
            (W*(1/σ)^2)*(firingrate.(V0 .- h, type, beta).*deltaF) )

end

function rhsFun_2pop(du, u, p, t)

    g_EXT,type,
    α_EE,α_EI,α_IE,α_II,
    τ_E,τ_I, κ_EE,κ_EI,κ_IE,κ_II,
    v_EE,v_EI, v_IE, v_II,W_EE,W_EI,W_IE,W_II,hE,hI,
    beta,d_EE,d_EI,d_IE,d_II,deltaF_EE,deltaF_EI,
    deltaF_IE,deltaF_II, Dxx, Dx, Drr,linspace_X = p

    Y1_EE = u[1:RX]
    Y1_EI = u[RX+1:2*RX]
    Y1_IE = u[2*RX+1:3*RX]
    Y1_II = u[3*RX+1:4*RX]

    Psi_EE = u[4*RX+1:5*RX]
    Psi_EI = u[5*RX+1:6*RX]
    Psi_IE = u[6*RX+1:7*RX]
    Psi_II = u[7*RX+1:8*RX]

    Y2_EE = u[8*RX+1:9*RX]
    Y2_EI = u[9*RX+1:10*RX]
    Y2_IE = u[10*RX+1:11*RX]
    Y2_II = u[11*RX+1:12*RX]

    g_EE = u[12*RX+1:13*RX]
    g_EI = u[13*RX+1:14*RX]
    g_IE = u[14*RX+1:15*RX]
    g_II = u[15*RX+1:16*RX]

    V_E = u[16*RX+1:17*RX]
    V_I = u[17*RX+1:18*RX]

    VE0 = repeat(V_E[((xzero - 1)*R  + 1):xzero * (R)], X)
    VI0 = repeat(V_I[((xzero - 1)*R  + 1):xzero * (R)], X)

#   Y1_EE
    du[1:RX] .=        z1(Y1_EE, Psi_EE, v_EE, VE0, W_EE, σ_EE, κ_EE, hE, type, beta, deltaF_EE, Dx, Drr,D_E,Dxx)
#   Y1_EI
    du[RX+1:2*RX] .=   z1(Y1_EI, Psi_EI, v_EI, VI0, W_EI, σ_EI, κ_EI, hI, type, beta, deltaF_EI, Dx, Drr,D_I,Dxx)
#   Y1_IE
    du[2*RX+1:3*RX] .= z1(Y1_IE, Psi_IE, v_IE, VE0, W_IE, σ_IE, κ_IE, hE, type, beta, deltaF_IE, Dx, Drr,D_E,Dxx)
#   Y1_II
    du[3*RX+1:4*RX] .= z1(Y1_II, Psi_II, v_II, VI0, W_II, σ_II, κ_II, hI, type, beta, deltaF_II, Dx, Drr,D_I,Dxx)
#   Psi_EE
    du[4*RX + 1:5 * RX] .=  v_EE*(-((1/σ_EE).*Psi_EE .+ κ_EE.*Dx*Psi_EE) + Y1_EE + (W_EE/σ_EE)*firingrate.(VE0 .- hE, type, beta).* deltaF_EE)
#   Psi_EI
    du[5*RX + 1:6 * RX] .=  v_EI*(-((1/σ_EI).*Psi_EI .+ κ_EI.*Dx*Psi_EI) + Y1_EI + (W_EI/σ_EI)*firingrate.(VI0.- hI, type, beta).* deltaF_EI)
#   Psi_IE
    du[6*RX + 1:7 * RX] .=  v_IE*(-((1/σ_IE).*Psi_IE .+ κ_IE.*Dx*Psi_IE) + Y1_IE + (W_IE/σ_IE)*firingrate.(VE0 .- hE, type, beta).* deltaF_IE)
#   Psi_II
    du[7*RX + 1:8 * RX] .=  v_II*(-((1/σ_II).*Psi_II .+ κ_II.*Dx*Psi_II) + Y1_II + (W_II/σ_II)*firingrate.(VI0 .- hI, type, beta).* deltaF_II)
#   Y_EE
    du[8 * RX + 1:9 * RX] .= α_EE * (g_EXT .+  Psi_EE .- Y2_EE)
#   Y_EI
    du[9 * RX + 1:10 * RX] .= α_EI * (g_EXT  .+ Psi_EI .- Y2_EI)
#   Y_IE
    du[10 * RX + 1:11 * RX] .= α_IE * (g_EXT  .+ Psi_IE .- Y2_IE)
#   Y_II
    du[11 * RX + 1:12 * RX] = α_II * (g_EXT .+  Psi_II .- Y2_II)
#   g_EE
    du[12 * RX + 1:13 * RX] .= α_EE * (Y2_EE .- g_EE)
#   g_EI
    du[13 * RX + 1:14 * RX] .= α_EI * (Y2_EI .- g_EI)
#   g_IE
    du[14 * RX + 1:15 * RX] .= α_IE * (Y2_IE .- g_IE)
#   g_II
    du[15 * RX + 1:16 * RX] .= α_II * (Y2_II .- g_II)
#   V_E
    du[16 * RX + 1:17 * RX] .= (-V_E / τ_E) .+  D_E.* Dxx * (V_E) .+  ( V_plus).*g_EE + ( V_neg).*g_EI
#   V_I
    du[17 * RX + 1:18 * RX] .= (-V_I / τ_I) .+  D_I.* Dxx * (V_I) .+  ( V_plus).*g_IE .+ ( V_neg).*g_II
    return du
end
