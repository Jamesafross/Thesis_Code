function firingrate(x,type,beta)
    #firing rate function (heaviside, sigmoid or tanh)
    y = Float64
    if type == "heaviside"
        if x >= beta
            y = 1
        elseif x < beta
            y = 0
        end
    elseif type == "sigmoid"
        y = 1/(1+exp(-beta*((x))))


    elseif type == "tanh"
        y = tanh(beta*x)
    end
    return y
end

function diff_firingrate(x,beta)
    return beta*sech(beta*x)^2
end


function heaviside(x)
    if x > 0
        return 1
    else
        return 0
    end
end


function g_ext(d,x,sigma,kappa,W0,dx,type,beta,V_ss)

    if kappa > 0
        return -(exp((d-x)/(kappa*sigma))*W0*(firingrate(V_ss,type,beta))*heaviside((-d + x)/kappa))/(sigma*kappa)
    else
        return -W0*(firingrate(V_ss,type,beta))*smooth_Delta(-d+x,dx)
    end
end

function W_hat(x,v,kappa,lambda,k,sigma,dx)
    if kappa == 0
        A = (sigma^2)*v*smooth_Delta.(x.-d,dx)
        return real((A/(2*(lambda*sigma + v - im*k*sigma*v))) + (A/(2*(lambda*sigma + v + im*k*sigma*v))))
    end
end
