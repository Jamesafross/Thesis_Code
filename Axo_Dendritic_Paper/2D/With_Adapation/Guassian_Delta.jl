function smooth_Delta(x,epsilon)

    return (exp.(-(x.^2)/(epsilon^2)))./(sqrt( pi*(epsilon^2)))
end

function diff_delta(x,epsilon)

    return -(2 *x .* exp.(-(x.^2)/(epsilon^2)))./(sqrt(pi)*((epsilon^2)^(3/2)))

end
