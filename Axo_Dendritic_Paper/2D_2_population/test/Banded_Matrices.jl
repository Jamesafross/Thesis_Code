

function D2z(R::Int64, X::Int64, RX::Int64)
    A = spzeros(RX, RX);
    A[1:RX + 1:end] .= -2
    A[RX+1:RX+1:end] .= 1
    A[2:RX + 1:end] .= 1
    A[(R-1)*RX+R+1:R*(RX + 1):end] .= 0
    A[R*(RX+1):R*(RX + 1):end] .= 0
    A[(R-1)*RX+1:R*(RX + 1):end] .= 1
    A[R:R*(RX + 1):end] .= 1

    return A
end

function D2x(R::Int64, X::Int64, RX::Int64,xzero)
    A = spzeros(RX,RX)
    B = spzeros(RX, RX);
    B[1,1] = -2
    B[1,RX-R] = 1
    B[1,RX-2*R] = 1
    for i = 2:RX
        B[i,:] .= circshift(A[1,:],i-1)
    end
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
    A[(xzero-6)*R,:] = B[(xzero-4)*R,:]
    A[(xzero-5)*R,:] = B[(xzero-3)*R,:]
    A[(xzero-4)*R,:] = B[(xzero-4)*R,:]
    A[(xzero-3)*R,:] = B[(xzero-3)*R,:]
    A[(xzero-2)*R,:] = B[(xzero-2)*R,:]
    A[(xzero-1)*R,:] = B[(xzero-1)*R,:]


    return A
end

function D2xTest(R::Int64, X::Int64, RX::Int64,xzero)

    B = spzeros(X, X);
    B[1,1] = 1
    B[1,end-1] = 1
    B[1,end] = -2
    for i = 2:X
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
