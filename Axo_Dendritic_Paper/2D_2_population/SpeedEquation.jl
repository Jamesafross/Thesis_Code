using Roots
function get_speed(sigma,alpha1,tau,D,v,h,xi_0)


function nu(lambda,alpha)
    return (alpha^2)/((alpha + lambda)^2)
end

function gamma2(lambda,tau,D)
    return ((1/tau) + lambda)/D
end


function F(gamma2,lambda,tau,D,xi_0)
    return exp(-abs(xi_0)*sqrt(gamma2(lambda,tau,D)))/(2*(D)*sqrt(gamma2(lambda,tau,D)))
end

function f(F,nu,x,gamma2,v,tau,D,xi_0,alpha1,h)
    return F(gamma2,x/(sigma*(1-(x/v))),tau,D,xi_0)*nu(x/(sigma*((1-(x/v)))),alpha1) - 2*h
end



f(x) = f(F,nu,x,gamma2,v,tau,D,xi_0,alpha1,h)
return  find_zero(f, (0, 10), FalsePosition())







end
