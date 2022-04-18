function firingrate(x,type,steepness)
    #firing rate function (heaviside or sigmoid)
    y = Float64
    if type == 0
        if x > 0
            y = 1
        elseif x < 0
            y = 0
        else
            y = 0.5
        end
    elseif type == 1
        y = 1/(1+exp(beta*(-(x))))
    end




    return y
end
