function sdelta(ii,d,R,a,dx,min_x)
#Guassian approximation of delta function
    if (ceil(ii/R)  - d)*dx == 0
        return 1/dx
    else
        return 0
    end
end
