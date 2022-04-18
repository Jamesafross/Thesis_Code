function D2r1(R1::Int64, R::Int64, X::Int64, RX::Int64)
    A = spzeros(RX, RX);
    A[1:RX + 1:end] .= -2
    A[RX + 1:RX + 1:end] .= 1
    A[2:RX + 1:end] .= 1
    A[(R1 - 1) * RX + R1 + 1:R1 * (RX + 1):end] .= 0
    A[R1 * (RX + 1):R1 * (RX + 1):end] .= 0
    A[(R1 - 1) * RX + 1:R1 * (RX + 1):end] .= 1
    A[R1:R1 * (RX + 1):end] .= 1

    return A
end
function D2r2(R1::Int64, R2::Int64, R::Int64, X::Int64, RX::Int64)
    A = spzeros(RX, RX)
    B = spzeros(R, R);
    B[1,1] = -2
    B[1,R - R1 + 1] = 1
    B[1,R1 + 1] = 1
    for i = 2:R
        B[i,:] = circshift!(B[i,:],B[1,:], i - 1)
    end

    A[1:R,1:R] = B

    for i = 1:X-1
         A[i * R + 1:(i + 1) * R,i * R + 1:(i + 1) * R] = circshift!(A[i * R + 1:(i + 1) * R,i * R + 1:(i + 1) * R],B, (0, i * R))
    end

    return A
end

function D2x(R::Int64, X::Int64, RX::Int64)
    A = spzeros(RX, RX)
    for i = 1:RX
        A[i,i] = -2
        if i <= RX - R
            A[i,i + R] = 1
        end

        if i < R + 1 || i > RX - R
            A[i,i] = -1
        end
    end

    A[:,:] = Symmetric(A)
    return A


    return A
end

function D1x(R::Int64, X::Int64, RX::Int64)
    A = spzeros(RX, RX);
    A[R * RX + R + 1:RX + 1:end - R * RX] .= 1
    EN = size(A[R * RX + R + 1:RX + 1:end - R * RX], 1)
    A[R + 1:RX + 1:EN * RX] .= -1
    return A
end
