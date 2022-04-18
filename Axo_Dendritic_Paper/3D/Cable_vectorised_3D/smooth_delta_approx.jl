function sdelta(ii,d,R,a,dx,min_x)
#Guassian approximation of delta function
    x = (ceil(ii/R)  - d)*dx
    return (1/(abs(a)*sqrt(pi)))*exp(-((x/a)^2))
end
