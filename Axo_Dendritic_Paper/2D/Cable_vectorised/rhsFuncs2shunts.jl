function z1(z, psi, v, V0, W, σ, k, h, type, beta, deltaF, Dx, Dzz)
    y = z + (1/σ)*W .* firingrate.(V0 .- h, type, beta).*deltaF

    return  v*(-((1/σ).*y .+ k.*Dx*y) .+ Dzz * psi .+
         (W * (1/σ)^2)*(firingrate.(V0 .- h, type, beta).*deltaF) .+
         (W*(1/σ))*k*Dx*(firingrate.(V0 .- h, type, beta).*deltaF))

end

function g_ext(d,x,sigma,kappa,W0,dx,type,beta,V_ss)

    if kappa > 0
        return -(exp((d-x)/(kappa*sigma))*W0*heaviside((-d + x)/kappa))/(sigma*kappa)
    else
        return -W0*(firingrate(V_ss,type,beta))*smooth_Delta(-d+x,dx)
    end
end

function smooth_Delta(x,epsilon)

    return (exp.(-(x.^2)/(epsilon^2)))./(sqrt( pi*(epsilon^2)))
end

function diff_delta(x,epsilon)

    return -(2 *x .* exp.(-(x.^2)/(epsilon^2)))./(sqrt(pi)*((epsilon^2)^(3/2)))

end

function firingrate(x,type,steepness)
    #firing rate function (heaviside or sigmoid)
    y = Float64
    if type == 0
        if x >= 0
            y = 1
        elseif x < 0
            y = 0
        end
    elseif type == -1
        if x >= 0
                y = x
        elseif x < 0
            y = 0
        end
    elseif type == 1
        y = 1/(1+exp(-beta*((x))))
    elseif type == 2
        y =  tanh(beta*x)
    end
    return y
end

function W_hat(K,λ,W,d,x)
    y = exp()
    return
end

function eta_tilde()

    return
end




function rhsFun(du, u, p, t)

    du[1:RX] .= z1(u[1:RX], u[RX + 1:2 * RX], v, repeat(u[4 * RX + 1:5 * RX][((xzero - 1)*R  + 1):xzero * (R)], X), W, σ, k, h, type, beta, deltaF, Dx, Drr)
#  Psi
    du[RX + 1:2 * RX] =  v*(-((1/σ).*u[RX+1:2*RX] .+ k.*Dx*u[RX+1:2*RX]) + u[1:RX] + (1/σ)*W*firingrate.(repeat(u[4 * RX + 1:5 * RX][((xzero - 1) * R +  1):xzero * (R)], X) .- h, type, beta).* deltaF)
# y2
    du[2 * RX + 1:3 * RX] = α * (g_c .+  u[RX + 1:2 * RX] - u[2 * RX + 1:3 * RX])
#  g
    du[3 * RX + 1:4 * RX] = α * (u[2 * RX + 1:3 * RX] - u[3 * RX + 1:4 * RX])
#  V
    du[4 * RX + 1:5 * RX] = (-u[4 * RX + 1:5 * RX] / tau) +  D.* Dxx * (u[4 * RX + 1:5 * RX]) +  (V_plus .- u[4 * RX + 1:5 * RX]).*u[3 * RX + 1:4 * RX]
    return du
end
