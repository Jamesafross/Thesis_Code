function η_tilde(α,λ)
    return α^2/((α + λ)^2)
end

function epsilon(τ,λ,D)
    return sqrt(complex(((1/τ) + λ)/D))
end

function Afun(σ,λ,v,κ,τ,D,k)
    ϵ = epsilon(τ,λ,D)
    A = ((1/σ) + (λ/v) + κ*ϵ)# abs here??
    return (1/σ)*((A)/((A)^2 + k^2))
end

function G_tilde(τ,λ,D,d)
    return (1/(2*D*epsilon(τ,λ,D)))*exp(-epsilon(τ,λ,D)*d)
end

function F(σ,λ,κ,τ,k,v,D,d,W0)
    return W0*G_tilde(τ,λ,D,d)*Afun(σ,λ,v,κ,τ,D,k)
end

function D_ab(α,γ,σ,λ,κ,τ,k,v,D,d,W0)
    return γ*η_tilde(α,λ)*F(σ,λ,κ,τ,k,v,D,d,W0)
end

function  E(D_matrix,V_plus,V_neg)
    return det(D_matrix*[V_plus 0; 0 V_neg] .- [1 0; 0 1])
end

function f!(F, x)
    F[1] = real(make_E(p,x[1],x[2],k))
    F[2] = imag(make_E(p,x[1],x[2],k))
end


function g!(F, x)
    F[1] = real(make_E(p,x[1],0,k))
end

function make_E(p,reλc,imλc,kc)
    D_EE = D_ab(α_EE,γ_E,σ_EE,reλc + im*imλc,κ_EE,τ_E,kc,v,D_E,d_EE,W0_EE)

    D_EI = D_ab(α_EI,γ_I,σ_EI,reλc + im*imλc,κ_EI,τ_I,kc,v,D_I,d_EI,W0_EI)

    D_IE = D_ab(α_IE,γ_E,σ_IE,reλc + im*imλc,κ_IE,τ_E,kc,v,D_E,d_IE,W0_IE)

    D_II = D_ab(α_II,γ_I,σ_II,reλc + im*imλc,κ_II,τ_I,kc,v,D_I,d_II,W0_II)

    D_mat = [D_EE D_EI; D_IE D_II]

    return E(D_mat,V_plus,V_neg)
end

function plotlines(reL,imL,k)
    y1 = zeros(length(reL))
    y2 = zeros(length(reL))
    for i = 1:length(reL)
        y1[i] = real(make_E(p,reL[i],imL,k))
        y2[i] = imag(make_E(p,reL[i],imL,k))
    end
    return plot(reL,[y1,y2])
end

function make_mov(reL,imL,k)
    @gif for i = 1:length(imL)
        plotlines(reL,imL[i],k)
    end
end
