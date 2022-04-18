function D2x(R::Int64)
    #second order central difference operator
    A = spzeros(R, R);
    A[1:R + 1:end] .= -2.0
    A[R + 1:R + 1:end] .= 1.0
    A[2:R + 1:end] .= 1.0
    A[1,end] = 1.0
    A[end,1] = 1.0

    return A
end

function D2r1(R1::Int64, R::Int64)
    A = spzeros(R, R);
    A[1:R + 1:end] .= -2.0
    A[R + 1:R + 1:end] .= 1.0
    A[2:R + 1:end] .= 1.0
    A[(R1 - 1) * R + R1 + 1:R1 * (R + 1):end] .= 0.0
    A[R1 * (R + 1):R1 * (R + 1):end] .= 0.0
    A[(R1 - 1) * R + 1:R1 * (R + 1):end] .= 1.0
    A[R1:R1 * (R + 1):end] .= 1.0
    return dropzeros(A)
end

function D2r2(R1::Int64, R::Int64)
    B = spzeros(R, R)
    B[1:R+1:end] .= -2.0
    B[R*R1 + 1:R+1:end] .= 1.0
    B[R1+1:R+1:end - (R1-1)*R] .= 1.0
    B[(R1-1)*R*R1 + 1:R+1:end] .= 1.0
    B[(R-R1) + 1:R+1:R*R1] .= 1.0
    return dropzeros(B)
end


function W_1D(σ,X_space)
        return exp.(-abs.(X_space)/σ)/(2σ)
end

function W_2D(σ,X_space1,X_space2)
    return exp.(-sqrt.(X_space^2 + X_space2^2)/σ)/(2πσ)
end
