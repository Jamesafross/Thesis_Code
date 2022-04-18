function z1(y, psi, v, β, V0, W, α, k, h, type, beta, deltaF, Dxx, Dx, Dr1r1,Dr2r2)


    return v*(-(α.*y .+ k.*Dx*y) + (3/2)*(Dr1r1 * psi + Dr2r2 * psi) .+
         (W * α^2).*firingrate.(V0 .- h, type, beta).*deltaF)
end


function rhsFun(du, u, p, t)

#d(u[1])/dt = d((v*w/sigma)*f(V) - y)/dt -- dpsi/dt = z
    du[1:RX] = z1(u[1:RX], u[RX + 1:2 * RX], v, β, repeat(u[4 * RX + 1:5 * RX][((xzero - 1)*R  + 1):xzero * (R)], X), W, α, k, h, type, beta, deltaF, Dxx, Dx, Dr1r1,Dr2r2)
#  Psi
    du[RX + 1:2 * RX] = v*(-(α.*u[RX+1:2*RX] .+ k.*Dx*u[RX+1:2*RX]) + u[1:RX])
# y2
    du[2 * RX + 1:3 * RX] = alfa * (u[RX + 1:2 * RX] - u[2 * RX + 1:3 * RX]  -beta_a*u[5 * RX + 1:6 * RX]) 
#  g
    du[3 * RX + 1:4 * RX] = alfa * (u[2 * RX + 1:3 * RX] - u[3 * RX + 1:4 * RX])
#  V
    du[4 * RX + 1:5 * RX] = (- u[4 * RX + 1:5 * RX] / tau)+ D .* Dxx * (u[4 * RX + 1:5 * RX]) + (V_plus.-u[4*RX + 1:5*RX]) .* u[3 * RX + 1:4 * RX]

    du[5 * RX + 1:6 * RX] = (1/tau_a)*(u[3*RX+1:4*RX] -u[5 * RX + 1:6 * RX] )

    return du
end
