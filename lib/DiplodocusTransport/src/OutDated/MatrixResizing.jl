function SixDtoThreeD(matrix::Array{Float32,6})

    # Flattens a 6D matrix to a 3D matrix where S[i,j,k,l,m,n] -> S[(j-1)*nump+i,(l-1)*nump+k,(n-1)*nump+m]

    nump3 = size(matrix,1)-1 # ignore overflow bin
    numt3 = size(matrix,2)
    nump1 = size(matrix,3)  
    numt1 = size(matrix,4)
    nump2 = size(matrix,5)
    numt2 = size(matrix,6)

    matrix3D = zeros(Float32,nump3*numt3,nump1*numt1,nump2*numt2)

    for n in axes(matrix,6), m in axes(matrix,5), l in axes(matrix,4), k in axes(matrix,3), j in axes(matrix,2), i in 1:nump3
        matrix3D[(j-1)*nump3+i,(l-1)*nump1+k,(n-1)*nump2+m] = matrix[i,j,k,l,m,n]
    end

    return matrix3D

end

function SixDtoTwoD(matrix::Array{Float32,6})

    # Flattens a 6D matrix to a 3D matrix where S[i,j,k,l,m,n] -> S[(j-1)*nump+i,(l-1)*nump+k,(n-1)*nump+m]

    nump3 = size(matrix,1)-1 # ignore overflow bin
    numt3 = size(matrix,2)
    nump1 = size(matrix,3)  
    numt1 = size(matrix,4)
    nump2 = size(matrix,5)
    numt2 = size(matrix,6)

    matrix3D = zeros(Float32,nump3*numt3,nump1*numt1,nump2*numt2)

    for n in axes(matrix,6), m in axes(matrix,5), l in axes(matrix,4), k in axes(matrix,3), j in axes(matrix,2), i in 1:nump3
        matrix3D[(j-1)*nump3+i,(l-1)*nump1+k,(n-1)*nump2+m] = matrix[i,j,k,l,m,n]
    end

    matrix2D = reshape(matrix3D,nump3*numt3*nump1*numt1,nump2*numt2)

    return matrix2D

end

function FourDtoTwoD(matrix::Array{Float32,4})

    # Flattens a 4D matrix to a 2D matrix where T[i,j,k,l] -> T[(j-1)*nump+i,(l-1)*nump+k]

    nump1 = size(matrix,1)  
    numt1 = size(matrix,2)
    nump2 = size(matrix,3)
    numt2 = size(matrix,4)

    matrix2D = zeros(Float32,nump1*numt1,nump2*numt2)

    for l in axes(matrix,4), k in axes(matrix,3), j in axes(matrix,2), i in axes(matrix,1)
        matrix2D[(j-1)*nump1+i,(l-1)*nump2+k] = matrix[i,j,k,l]
    end

    return matrix2D

end

function u_to_f_list!(f_list::Vector{Vector{Float32}},u::Vector{Float32})

    j = 0
    for i in 1:length(f_list)
        k = length(f_list[i]) # nump_list[i]*numt_list[i]
        f_list[i] .= @view(u[(1+j):(j+k)])
        j += k
    end

end

function f_list_to_u!(u::Vector{Float32},f_list::Vector{Vector{Float32}})

    j = 0
    for i in 1:length(f_list)
        k = length(f_list[i])
        @view(u[(1+j):(j+k)]) .= f_list[i]
        j += k
    end

end

