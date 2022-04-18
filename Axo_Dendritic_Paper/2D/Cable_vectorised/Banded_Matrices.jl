

function D2z(R::Int64, X::Int64, RX::Int64)
    A = spzeros(RX, RX);
    A[1:RX + 1:end] .= -2
    A[RX+1:RX+1:end] .= 1
    A[2:RX + 1:end] .= 1
    A[(R-1)*RX+R+1:R*(RX + 1):end] .= 0
    A[R*(RX+1):R*(RX + 1):end] .= 0
    A[(R-1)*RX+1:R*(RX + 1):end] .= 1
    A[R:R*(RX + 1):end] .= 1

    return dropzeros(A)
end

function D2z2(R::Int64, X::Int64, RX::Int64)
    A = spzeros(RX, RX);
    A[1:RX + 1:end] .= -2
    A[RX+1:RX+1:end] .= 1
    A[2:RX + 1:end] .= 1
    A[(R-1)*RX+R+1:R*(RX + 1):end] .= 0
    A[(R-1)*(RX)+R:R*(RX + 1):end] .= -1
    A[1:R*(RX + 1):end] .= -1
    A[R:R*(RX + 1):end] .= 0
    A[(R*RX) + R:R*(RX + 1):end] .= 0

    return dropzeros(A)
end



function D2xB(R::Int64, X::Int64, RX::Int64,d)
    B = spzeros(RX, RX)
    B[1,1] = 1
    B[1,RX-R] = -2
    B[1,RX-2*R] = 1
    for i = 2:X*R
        B[i,:] .= circshift(B[1,:],i-1)
    end

    return B
end

function D2xF(R::Int64, X::Int64, RX::Int64)
    B = spzeros(RX, RX)
    B[1,1] = 1
    B[1,R] = -2
    B[1,2*R] = 1
    for i = 2:X*R
        B[i,:] .= circshift(B[1,:],i-1)
    end
    return B
end

function D1xB(R::Int64, X::Int64, RX::Int64)
    A = spzeros(RX, RX);
    A[R * RX + R + 1:RX + 1:end - R * RX] .= 1
    EN = size(A[R * RX + R + 1:RX + 1:end - R * RX], 1)
    A[R + 1:RX + 1:EN * RX] .= -1
    return A
end

function D1xF(R::Int64, X::Int64, RX::Int64)
    A = spzeros(RX, RX);
    A[R * RX + R + 1:RX + 1:end - R * RX] .= -1
    A[2*(R * RX) + 1 + R:RX + 1:end] .= 1
    #EN = size(A[R * RX + R + 1:RX + 1:end - R * RX], 1)
    #A[R + 1:RX + 1:EN * RX] .= -1
    return A
end


function D2xCent(R::Int64, X::Int64, RX::Int64,xzero)
    A = spzeros(RX,RX)
    for i = 1:RX
        A[i,i] = -2
        if i <= RX-R
        A[i,i+R] = 1
        end

        if i < R+1 || i >RX-R
            A[i,i] = -1
        end
    end
    A[:,:] = Symmetric(A)
    return A
end

function D2xCent2(R::Int64, X::Int64, RX::Int64,xzero)
    B = spzeros(RX, RX)
    B[1,1] = -2
    B[1,RX-R] = 1
    B[1,R] = 1
    for i = 2:X*R
        B[i,:] .= circshift(B[1,:],i-1)
    end
    return B
end
