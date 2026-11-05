mutable struct GainCorrectionStruct <: Function 

    CorrType::Int64 
    Indistinguishable_12::Bool

    E1_d::Vector{Float64}
    E2_d::Vector{Float64}
    E3_d::Vector{Float64}
    E4_d::Vector{Float64}

    high_inds3::Vector{Int64}
    high_inds4::Vector{Int64}

    # for spectrum filters which keeps Float32 precision
    tmpvec3::Vector{Float32}
    tmpvec4::Vector{Float32}

    # for correction which upgrades to Float64 precision
    tmpsortvec3::Vector{Float64}
    tmpsortvec4::Vector{Float64}

    tmpsortvec3_reshape::Array{Float64,3}
    tmpsortvec4_reshape::Array{Float64,3}

    cart_inds3::Vector{CartesianIndex{3}}
    cart_inds4::Vector{CartesianIndex{3}}

    size3::Int64
    size4::Int64

    function GainCorrectionStruct(Parameters::Tuple{String, String, String, String, Float64, Float64, Float64, Float64, Float64, Float64, String, Int64, String, Int64, String, Int64, Float64, Float64, String, Int64, String, Int64, String, Int64, Float64, Float64, String, Int64, String, Int64, String, Int64, Float64, Float64, String, Int64, String, Int64, String, Int64})

        (name1,name2,name3,name4,m1,m2,m3,m4,p1_low,p1_up,p1_grid,p1_num,u1_grid,u1_num,h1_grid,h1_num,p2_low,p2_up,p2_grid,p2_num,u2_grid,u2_num,h2_grid,h2_num,p3_low,p3_up,p3_grid,p3_num,u3_grid,u3_num,h3_grid,h3_num,p4_low,p4_up,p4_grid,p4_num,u4_grid,u4_num,h4_grid,h4_num) = Parameters

        if (name1 == name3) && (name2 == name4) && (name1 != name2)
            CorrType = 1
        else
            CorrType = 2
        end

        Indistinguishable_12::Bool = name1 == name2

        # underflow and overflow bins are taken to have size 0 -> p1_r[1] and p1_r[end] -> 2*p1_r[end] respectively

            p1_r = bounds(p1_low,p1_up,p1_num,p1_grid);
            p1_d = deltaVector(p1_r);

            p2_r = bounds(p2_low,p2_up,p2_num,p2_grid);
            p2_d = deltaVector(p2_r);

            p3_r = bounds(p3_low,p3_up,p3_num,p3_grid);
            p3_r = [0.0 ; p3_r ; 2*p3_r[end]];
            p3_d = deltaVector(p3_r);

            p4_r = bounds(p4_low,p4_up,p4_num,p4_grid);
            p4_r = [0.0 ; p4_r ; 2*p4_r[end]];
            p4_d = deltaVector(p4_r);

            E1_Δ = deltaEVector(p1_r,m1);
            E1_d = E1_Δ ./ p1_d

            E2_Δ = deltaEVector(p2_r,m2);
            E2_d = E2_Δ ./ p2_d

            E3_Δ = deltaEVector(p3_r,m3);
            E3_d = E3_Δ ./ p3_d

            E4_Δ = deltaEVector(p4_r,m4);
            E4_d = E4_Δ ./ p4_d

        size3 = (p3_num+2)*u3_num*h3_num
        size4 = (p4_num+2)*u4_num*h4_num

        high_inds3 = zeros(Int64,size3)
        high_inds4 = zeros(Int64,size4)

        tmpvec3 = zeros(Float32,(p3_num+2))
        tmpvec4 = zeros(Float32,(p4_num+2))

        tmpsortvec3 = zeros(Float64,size3)
        tmpsortvec4 = zeros(Float64,size4)

        tmpsortvec3_reshape = reshape(tmpsortvec3,((p3_num+2),u3_num,h3_num))
        tmpsortvec4_reshape = reshape(tmpsortvec4,((p4_num+2),u4_num,h4_num))

        cart_inds3 = Vector{CartesianIndex{3}}(undef,size3)
        cart_inds4 = Vector{CartesianIndex{3}}(undef,size4)

        # build struct 

        self = new()
        self.CorrType = CorrType
        self.Indistinguishable_12 = Indistinguishable_12
        self.E1_d = E1_d
        self.E2_d = E2_d
        self.E3_d = E3_d
        self.E4_d = E4_d
        self.high_inds3 = high_inds3
        self.high_inds4 = high_inds4
        self.tmpvec3 = tmpvec3
        self.tmpvec4 = tmpvec4

        self.tmpsortvec3 = tmpsortvec3
        self.tmpsortvec4 = tmpsortvec4
        self.tmpsortvec3_reshape = tmpsortvec3_reshape
        self.tmpsortvec4_reshape = tmpsortvec4_reshape

        self.cart_inds3 = cart_inds3
        self.cart_inds4 = cart_inds4

        self.size3 = size3
        self.size4 = size4

        return self

    end

end 

# ================== chunked version ================== #

"""
    GainCorrectionChunk!(Parameters,GainMatrix3,GainMatrix4,LossMatrix1,LossMatrix2)

MC sampling introduces noise that can lead to poor number and energy conservation. `GainCorrection` provides a corrective step to ensure number and energy conservation to numerical precision. If there is no GainMatrix element then the value of the LossMatrix is applied to the same bin as the input state (if they are identical particles); if not identical particles if there is no GainMatrix element then the value of the LossMatrix is set to zero to ensure particle conservation (with good MC sampling this should rarely occur).
This version corrects the input arrays `GainMatrix3`, `GainMatrix4`, and `LossMatrix` in place.
"""
function GainCorrectionChunk!(Parameters::Tuple{String, String, String, String, Float64, Float64, Float64, Float64, Float64, Float64, String, Int64, String, Int64, String, Int64, Float64, Float64, String, Int64, String, Int64, String, Int64, Float64, Float64, String, Int64, String, Int64, String, Int64, Float64, Float64, String, Int64, String, Int64, String, Int64}, GainMatrix3::AbstractArray{Float32, 7}, GainMatrix4::AbstractArray{Float32, 7}, LossMatrix::AbstractArray{Float32, 4},p1::Int64,p2::Int64,GainCorrectionTmp::GainCorrectionStruct)

    
    #= Different possible combinations of identical and different particles that affect how to apply conservation corrections:

        1. n1=n3 != n2=n4 (1 and 3 identical, 2 and 4 identical)
            3 parameters needed:
                alpha: consv. num for 1 and 3
                beta: consv. num for 2 and 4
                gamma: consv. energy for all
            In this setup gain of n3 or n4 could be zero from MC, so have to be careful when applying correction    
        2. all other combinations:
            2 parameters needed:
                alpha: consv. num for all
                beta: consv. energy for all

    =#

    #(name1,name2,name3,name4,m1,m2,m3,m4,p1_low,p1_up,p1_grid,p1_num,u1_grid,u1_num,h1_grid,h1_num,p2_low,p2_up,p2_grid,p2_num,u2_grid,u2_num,h2_grid,h2_num,p3_low,p3_up,p3_grid,p3_num,u3_grid,u3_num,h3_grid,h3_num,p4_low,p4_up,p4_grid,p4_num,u4_grid,u4_num,h4_grid,h4_num) = Parameters

    m1 = Parameters[5]
    m2 = Parameters[6]

    tol = Float64.(eps(Float32)) # tolerance for how low the gain terms can be compared to the loss terms, if less than this tolerance then the gain terms are set to zero and not included in the correction calculation.

    CorrType = GainCorrectionTmp.CorrType
    Indistinguishable_12 = GainCorrectionTmp.Indistinguishable_12

    E1_d = GainCorrectionTmp.E1_d
    E2_d = GainCorrectionTmp.E2_d
    E3_d = GainCorrectionTmp.E3_d
    E4_d = GainCorrectionTmp.E4_d

    high_inds3 = GainCorrectionTmp.high_inds3
    high_inds4 = GainCorrectionTmp.high_inds4
    tmpvec3 = GainCorrectionTmp.tmpvec3
    tmpvec4 = GainCorrectionTmp.tmpvec4
    tmpsortvec3 = GainCorrectionTmp.tmpsortvec3
    tmpsortvec4 = GainCorrectionTmp.tmpsortvec4
    tmpsortvec3_reshape = GainCorrectionTmp.tmpsortvec3_reshape
    tmpsortvec4_reshape = GainCorrectionTmp.tmpsortvec4_reshape
    cart_inds3 = GainCorrectionTmp.cart_inds3
    cart_inds4 = GainCorrectionTmp.cart_inds4
    size3 = GainCorrectionTmp.size3
    size4 = GainCorrectionTmp.size4

    num_wrong::Int64 = 0
    num_right::Int64 = 0
    #cart_inds3 = CartesianIndices(@view(GainMatrix3[:,:,:,1,1,1,1]))[high_inds3]
    #cart_inds4 = CartesianIndices(@view(GainMatrix4[:,:,:,1,1,1,1]))[high_inds4]

    LossSumN1::Float64 = 0.0
    LossSumN2::Float64 = 0.0
    LossSumE1::Float64 = 0.0
    LossSumE2::Float64 = 0.0

    for h2 in axes(GainMatrix3,7), u2 in axes(GainMatrix3,6), h1 in axes(GainMatrix3,5), u1 in axes(GainMatrix3,4)

        # generate filtered gain matrices to smooth out spectrum
        GainMatrix3Filtered = @view(GainMatrix3[:,:,:,u1,h1,u2,h2])
        GainMatrix4Filtered = @view(GainMatrix4[:,:,:,u1,h1,u2,h2])

        # Apply filters
        for u3 in axes(GainMatrix3,2), h3 in axes(GainMatrix3,3) 
            # fixes points that are break monotonically increasing regions
            fix_monotone_center_log!(@view(GainMatrix3Filtered[:,u3,h3]),tmpvec3) 
        end
        for u4 in axes(GainMatrix4,2), h4 in axes(GainMatrix4,3) 
            # fixes points that are break monotonically increasing regions
            fix_monotone_center_log!(@view(GainMatrix4Filtered[:,u4,h4]),tmpvec4) 
        end

        tmpsortvec3_reshape .= E3_d .* GainMatrix3Filtered
        tmpsortvec4_reshape .= E4_d .* GainMatrix4Filtered

        # sort by highest energy contents
        partialsortperm_no_view!(high_inds3,tmpsortvec3, 1:size3; rev=true)
        #cart_inds3 .= CartesianIndices(GainMatrix3Filtered)[high_inds3]
        ci3 = CartesianIndices(GainMatrix3Filtered)
        for i in eachindex(high_inds3)
            cart_inds3[i] = ci3[high_inds3[i]]
        end

        partialsortperm_no_view!(high_inds4,tmpsortvec4, 1:size4; rev=true)
        #cart_inds4 .= CartesianIndices(GainMatrix4Filtered)[high_inds4]
        ci4 = CartesianIndices(GainMatrix4Filtered)
        for i in eachindex(high_inds4)
            cart_inds4[i] = ci4[high_inds4[i]]
        end

        LossSumN1 = Float64(LossMatrix[u1,h1,u2,h2])
        if !Indistinguishable_12
            LossSumN2 = Float64(LossMatrix[u1,h1,u2,h2])
        else
            LossSumN2 = 0e0
        end
        if !Indistinguishable_12 # LossMatrix2 is non-zero
            LossSumE1 = Float64(LossMatrix[u1,h1,u2,h2])*E1_d[p1]
            LossSumE2 = Float64(LossMatrix[u1,h1,u2,h2])*E2_d[p2]
        else
            LossSumE1 = Float64(LossMatrix[u1,h1,u2,h2])*(E1_d[p1]+E2_d[p2])/2e0
            LossSumE2 = 0e0#LossMatrix1[u1,h1,u2,h2]*(E1_d[p1]+E2_d[p2])/2
        end

        if LossSumN1 + LossSumN2 == 0e0
            # no loss term so no interaction at these incoming states

            num_right += 1

            continue
        end
        
        wrong::Bool = true
        nonzero_gain::Bool = true
        max_high_bins::Int64 = 0
        Gain3False::Bool = false
        Gain4False::Bool = false

        p1Big::Bool = false
        p2Big::Bool = false

        p3_offset = length(axes(GainMatrix3,1))
        p4_offset = length(axes(GainMatrix4,1))

        high_3_range::Int64 = 1
        high_4_range::Int64 = 1

        #=if E1_d[p1] > 1e0*E2_d[p2] && m1 > m2
            p1Big = true
        elseif E2_d[p2] > 1e0*E1_d[p1] && m2 > m1
            p2Big = true
        end=#

        #cart_inds3 = CartesianIndices(GainMatrix3Filtered)[(1)]
        #cart_inds4 = CartesianIndices(GainMatrix4Filtered)[(1)]

        GainSumN31::Float64 = 0e0
        GainSumN32::Float64 = 0e0
        GainSumN41::Float64 = 0e0
        GainSumN42::Float64 = 0e0
        GainSumE31::Float64 = 0e0
        GainSumE32::Float64 = 0e0
        GainSumE41::Float64 = 0e0
        GainSumE42::Float64 = 0e0

        GainN1::Float64 = 0e0
        GainN2::Float64 = 0e0
        GainE1::Float64 = 0e0
        GainE2::Float64 = 0e0
        LossN::Float64 = 0e0
        LossE::Float64 = 0e0

        tmpN::Float64 = 0e0
        tmpE::Float64 = 0e0

        a1::Float64 = 0e0
        b1::Float64 = 0e0
        a2::Float64 = 0e0
        b2::Float64 = 0e0
        c1::Float64 = 0e0
        d1::Float64 = 0e0
        c2::Float64 = 0e0
        d2::Float64 = 0e0
        ac1::Float64 = 0e0
        bd1::Float64 = 0e0
        ac2::Float64 = 0e0
        bd2::Float64 = 0e0
        e1::Float64 = 0e0
        e2::Float64 = 0e0
        f::Float64 = 0e0
        L::Float64 = 0e0
        LE::Float64 = 0e0
        alpha1::Float64 = 0e0
        alpha2::Float64 = 0e0
        beta::Float64 = 0e0
        a::Float64 = 0e0
        b::Float64 = 0e0
        c::Float64 = 0e0
        d::Float64 = 0e0
        ac::Float64 = 0e0
        bd::Float64 = 0e0
        e::Float64 = 0e0

        while wrong && nonzero_gain

            GainSumN31 = zero(Float64)
            GainSumN32 = zero(Float64)
            GainSumN41 = zero(Float64)
            GainSumN42 = zero(Float64)
            GainSumE31 = zero(Float64)
            GainSumE32 = zero(Float64)
            GainSumE41 = zero(Float64)
            GainSumE42 = zero(Float64)

            #println("max_high_bins = $max_high_bins")
            
            # search approach is to first find the bin with the highest GainMatrix value 
            high_bins = 0
            searching = true
            #p3_offset = length(axes(GainMatrix3,1))
            
            # sort bins by highest energy contents
            #partialsortperm!(high_inds3,reshape(E3_d .* GainMatrix3Filtered,((p3_num+2)*u3_num*h3_num)), 1:min(max_high_bins+1,size(GainMatrix3,1)); rev=true)

            # get range of indices for the highest bins
            high_3_range = min(max_high_bins+1,size(GainMatrix3,1))
            high_4_range = min(max_high_bins+1,size(GainMatrix4,1))
            #cart_inds3 = CartesianIndices(GainMatrix3Filtered)[high_inds3]

            #println("high_inds3 = $high_inds3")
            #println("cart_inds3 = $cart_inds3")

            for p3 in axes(GainMatrix3,1), u3 in axes(GainMatrix3,2), h3 in axes(GainMatrix3,3) 
                # remove small rates (both in number and energy) to avoid numerical issues and reduce array size in the correction step
                if GainMatrix3Filtered[p3,u3,h3] < tol * (LossSumN1+LossSumN2) && GainMatrix3Filtered[p3,u3,h3]*E3_d[p3] < tol * (LossSumE1+LossSumE2)
                    GainMatrix3Filtered[p3,u3,h3] = zero(eltype(GainMatrix3Filtered)) # remove small values
                    continue
                else 
                    tmpN = Float64(GainMatrix3Filtered[p3,u3,h3])
                    tmpE = Float64(GainMatrix3Filtered[p3,u3,h3])*E3_d[p3]
                end
                if !p2Big && CartesianIndex(p3, u3, h3) in @view(cart_inds3[1:high_3_range])
                    GainSumN32 += tmpN
                    GainSumE32 += tmpE
                else
                    GainSumN31 += tmpN
                    GainSumE31 += tmpE
                end
            end

            #partialsortperm!(high_inds4,reshape(E4_d .* GainMatrix4Filtered,((p4_num+2)*u4_num*h4_num)), 1:min(max_high_bins+1,size(GainMatrix4,1)); rev=true)
            #cart_inds4 = CartesianIndices(GainMatrix4Filtered)[high_inds4]

            for p4 in axes(GainMatrix4,1), u4 in axes(GainMatrix4,2), h4 in axes(GainMatrix4,3) 
                # remove small rates (both in number and energy) to avoid numerical issues and reduce array size in the correction step
                if GainMatrix4Filtered[p4,u4,h4] < tol * (LossSumN1+LossSumN2) && GainMatrix4Filtered[p4,u4,h4]*E4_d[p4] < tol * (LossSumE1+LossSumE2)
                    GainMatrix4Filtered[p4,u4,h4] = zero(eltype(GainMatrix4Filtered)) # remove small values
                    continue
                else
                    tmpN = Float64(GainMatrix4Filtered[p4,u4,h4])
                    tmpE = Float64(GainMatrix4Filtered[p4,u4,h4])*E4_d[p4]
                end
                if !p1Big && CartesianIndex(p4, u4, h4) in @view(cart_inds4[1:high_4_range])
                    GainSumN42 += tmpN
                    GainSumE42 += tmpE
                else
                    GainSumN41 += tmpN
                    GainSumE41 += tmpE
                end
            end

            if CorrType == 1 

                a1 = GainSumN31
                b1 = GainSumN32

                a2 = GainSumN41
                b2 = GainSumN42

                ac1 = GainSumE31
                c1 = ac1 != 0.0 ? GainSumE31/GainSumN31 : 0.0
                bd1 = GainSumE32
                d1 = bd1 != 0.0 ? GainSumE32/GainSumN32 : 0.0

                ac2 = GainSumE41
                c2 = ac2 != 0.0 ? GainSumE41/GainSumN41 : 0.0
                bd2 = GainSumE42
                d2 = bd2 != 0.0 ? GainSumE42/GainSumN42 : 0.0
                
                if (a1+b1 == 0e0 || (max_high_bins==0 && a1==0.0 && b1!=0.0))# There is not gain term produced by MC or only a single bin
                    ##println("here")
                    #c1 = p2Big ? E3_d[p1] : E3_d[p1-1]
                    #d1 = p2Big ? E3_d[p1+1] : E3_d[p1]
                    # underflow bins adds 1 to value of p3 index compared to p1 index for same energy bin.
                    c1 = (E2_d[p2] > E1_d[p1] && m2 > m1)#=p2Big=# ? E3_d[p1+1] : E3_d[p1]
                    d1 = (E2_d[p2] > E1_d[p1] && m2 > m1)#=p2Big=##=p2Big=# ? E3_d[p1+2] : E3_d[p1+1]
                    L = LossSumN1
                    LE = LossSumE1+LossSumE2
                    e2 = a2 + b2 - LossSumN2
                    #c2 = GainSumE41/GainSumN41

                    a1 = ((ac2-c2*e2)+(d1*L-LE))/(d1-c1)
                    b1 = L-a1
                    ac1 = a1*c1
                    bd1 = b1*d1

                    Gain3False = true
                end

                if (a2+b2 == 0e0 || (max_high_bins==0 && a2==0e0 && b2!=0e0))# There is not gain term produced by MC or only a single bin
                    #println("there")
                    #c2 = p1Big ? E4_d[p2] : E4_d[p2-1]
                    #d2 = p1Big ? E4_d[p2+1] : E4_d[p2]
                    # underflow bins adds 1 to value of p4 index compared to p2 index for same energy bin.
                    c2 = (E1_d[p1] > E2_d[p2] && m1 > m2)#=p1Big=# ? E4_d[p2+1] : E4_d[p2]
                    d2 = (E1_d[p1] > E2_d[p2] && m1 > m2)#=p1Big=##=p1Big=# ? E4_d[p2+2] : E4_d[p2+1]
                    L = LossSumN2
                    LE = LossSumE2 + LossSumE1
                    e1 = a1 + b1 - LossSumN1
                    #c1 = GainSumE31/GainSumN31

                    a2 = ((ac1-c1*e1)+(d2*L-LE))/(d2-c2)
                    b2 = L-a2
                    ac2 = a2*c2
                    bd2 = b2*d2
                    Gain4False = true
                end

                e1 = a1 + b1 - LossSumN1
                e2 = a2 + b2 - LossSumN2
                f  = ac1 + bd1 + ac2 + bd2 - LossSumE1 - LossSumE2

                alpha1 = (b2*(d2-c2)*e1+b1*(d1*e1+c2*e2-f))/(a1*(b1*(c1-d1)+b2*(c2-d2))) + 1e0
                alpha2 = (b1*(d1-c1)*e2+b2*(d2*e2+c1*e1-f))/(a2*(b1*(c1-d1)+b2*(c2-d2))) + 1e0
                beta = (f-c1*e1-c2*e2)/(b1*(c1-d1)+b2*(c2-d2)) + 1e0

                if a1 <= 0e0 || a2 <= 0e0 || max_high_bins > max(size(GainMatrix3,1),size(GainMatrix4,1))# loop has not produced any corrected gain terms

                    nonzero_gain = false
                    #println("No valid correction for p1=$p1,p2=$p2, u1=$u1, u2=$u2, h1=$h1, h2=$h2")
                    #println("a1: $a1, b1: $b1, ac1: $ac1, bd1: $bd1, e1: $e1, a2: $a2, b2: $b2, ac2: $ac2, bd2: $bd2, e2: $e2, f: $f, c1: $c1, d1: $d1, c2: $c2, d2: $d2")
                    #println("GainSumN31: $GainSumN31, GainSumN32: $GainSumN32, GainSumN41: $GainSumN41, GainSumN42: $GainSumN42, LossSumN1: $LossSumN1, LossSumN2: $LossSumN2")
                    #println("GainSumE31: $GainSumE31, GainSumE32: $GainSumE32, GainSumE41: $GainSumE41, GainSumE42: $GainSumE42, LossSumE1: $LossSumE1, LossSumE2: $LossSumE2")
                    #println("alpha1: $alpha1, alpha2: $alpha2, beta: $beta, max_high_bins: $max_high_bins")
                    alpha1 = 0e0
                    alpha2 = 0e0
                    beta = 0e0
                    LossMatrix[u1,h1,u2,h2] = 0f0

                    num_wrong += 1

                    continue

                end

                if isnan(alpha1) || isnan(alpha2) || isnan(beta)
                    @warn "NaN encountered for p1=$p1,p2=$p2, u1=$u1, u2=$u2, h1=$h1, h2=$h2"
                    #println("a1: $a1, b1: $b1, ac1: $ac1, bd1: $bd1, e1: $e1, a2: $a2, b2: $b2, ac2: $ac2, bd2: $bd2, e2: $e2, f: $f, c1: $c1, d1: $d1, c2: $c2, d2: $d2")
                    #println("GainSumN31: $GainSumN31, GainSumN32: $GainSumN32, GainSumN41: $GainSumN41, GainSumN42: $GainSumN42")
                    #println("alpha1: $alpha1, alpha2: $alpha2, beta: $beta")
                    nonzero_gain = false

                    # set this bin to zero to be safe
                    alpha1 = 0e0
                    alpha2 = 0e0
                    beta = 0e0
                    LossMatrix[u1,h1,u2,h2] = 0f0

                    num_wrong += 1

                    continue

                end

                if  alpha1 < 0e0 || alpha2 < 0e0 || beta < 0e0

                    #println("p1=$p1,p2=$p2, u1=$u1, u2=$u2, h1=$h1, h2=$h2")
                    #println("alpha1: $alpha1, alpha2: $alpha2, beta: $beta, max_high_bins: $max_high_bins, p3_offset: $p3_offset, p4_offset: $p4_offset")

                    max_high_bins += 1

                    continue # redo calculation
                else
                    #println("p1=$p1,p2=$p2")
                    #println("alpha1: $alpha1, alpha2: $alpha2, beta: $beta, max_high_bins: $max_high_bins, p3_offset: $p3_offset, p4_offset: $p4_offset")
                    #println("new gain1 = $(a1*alpha1+b1*beta), new gain2 = $(a2*alpha2+b2*beta), new gain = $((a1*alpha1+b1*beta)+(a2*alpha2+b2*beta)), Loss = $(LossSumN1+LossSumN2), err = $((a1*alpha1+b1*beta)+(a2*alpha2+b2*beta)- (LossSumN1+LossSumN2)), L1 = $LossSumN1, L2 = $LossSumN2")
                    wrong = false
                    num_right += 1
                end

            elseif CorrType == 2

                GainN1 = GainSumN31 + GainSumN41
                GainE1 = GainSumE31 + GainSumE41
                GainN2 = GainSumN32 + GainSumN42
                GainE2 = GainSumE32 + GainSumE42
                LossN = LossSumN1 + LossSumN2
                LossE = LossSumE1 + LossSumE2

                a = GainN1
                if a == 0e0 || max_high_bins > max(size(GainMatrix3,1),size(GainMatrix4,1))# loop has not produced any gain terms
                    nonzero_gain = false
                    @warn "No valid correction for p1=$p1,p2=$p2, u1=$u1, u2=$u2, h1=$h1, h2=$h2"

                    alpha1 = 0e0
                    alpha2 = 0e0
                    beta = 0e0
                    LossMatrix[u1,h1,u2,h2] = 0f0

                    num_wrong += 1
                    
                    continue
                end

                ac = GainE1
                c = GainE1/GainN1

                b = GainN2
                bd = GainE2
                d = GainE2/GainN2

                e = a+b-LossN
                f = ac+bd-LossE

                if isapprox(c,d)
                    # In this case deltaE is the approx the same for all bins involved, so energy is approx conserved if number is conserved so only one equation and one parameter needed which we take scales b
                    #alpha = -e/(a+b) + 1
                    #beta = alpha
                    beta = -e/b + 1e0
                    alpha = 1e0
                else
                    alpha = (e*d-f)/(ac-a*d) + 1e0
                    beta = (e*c-f)/(bd-b*c) + 1e0
                end

                if isnan(alpha) || isnan(beta)
                    @warn "NaN encountered for p1=$p1,p2=$p2, u1=$u1, u2=$u2, h1=$h1, h2=$h2"
                    #println("a: $a, b: $b, ac: $ac, bd: $bd, e: $e, f: $f, c: $c, d: $d")
                    #println("$(a+b-LossN), $(ac+bd-LossE)")
                    #println("LossE/LossN: $(LossE/LossN), bd/b: $(bd/b), ac/a: $(ac/a)")
                    nonzero_gain = false
                    # set this bin to zero to be safe
                    alpha1 = 0e0
                    alpha2 = 0e0
                    beta = 0e0
                    LossMatrix[u1,h1,u2,h2] = 0f0

                    continue
                end

                if alpha < 0e0 || beta < 0e0

                    #println("alpha: $alpha, beta: $beta, max_high_bins: $max_high_bins, p3_offset: $p3_offset, p4_offset: $p4_offset")

                    max_high_bins += 1

                    continue # redo calculation
                else
                    #println("p1=$p1,p2=$p2")
                    #println("alpha: $alpha, beta: $beta, max_high_bins: $max_high_bins, p3_offset: $p3_offset, p4_offset: $p4_offset")
                    wrong = false
                    num_right += 1
                    alpha1 = alpha
                    alpha2 = alpha
                end

            end # if CorrType

        end # while

        if Gain3False 
            # underflow bins adds 1 to value of p3 index compared to p1 index for same energy bin.
            GainMatrix3[p1+1+1,u1,h1,u1,h1,u2,h2] = Float32(p2Big * b1 * beta)
            GainMatrix3[p1+1,u1,h1,u1,h1,u2,h2] = Float32(p2Big ? a1 * alpha1 : b1 * beta)
            GainMatrix3[p1-1+1,u1,h1,u1,h1,u2,h2] = Float32(!p2Big * a1 * alpha1)
        else
            for p3 in axes(GainMatrix3,1)
                for u3 in axes(GainMatrix3,2), h3 in axes(GainMatrix3,3) 
                    if !p2Big && CartesianIndex(p3, u3, h3) in @view(cart_inds3[1:high_3_range]) 
                        #if GainMatrix3[p3,u3,h3,p1,u1,h1,p2,u2,h2] > tol * LossSumN1
                            GainMatrix3[p3,u3,h3,u1,h1,u2,h2] *= Float32(beta)
                        #else
                        #    CorrectedGainMatrix3[p3,u3,h3,p1,u1,h1,p2,u2,h2] = 0.0
                        #end
                    else
                        #if GainMatrix3[p3,u3,h3,p1,u1,h1,p2,u2,h2] > tol * LossSumN1
                            GainMatrix3[p3,u3,h3,u1,h1,u2,h2] *= Float32(alpha1)
                        #else
                        #    CorrectedGainMatrix3[p3,u3,h3,p1,u1,h1,p2,u2,h2] = 0.0
                        #end
                    end
                end
            end
        end

        if Gain4False 
            # underflow bins adds 1 to value of p4 index compared to p2 index for same energy bin.
            GainMatrix4[p2+1+1,u2,h2,u1,h1,u2,h2] = Float32(p1Big * b2 * beta)
            GainMatrix4[p2+1,u2,h2,u1,h1,u2,h2] = Float32(p1Big ? a2 * alpha2 : b2 * beta)
            GainMatrix4[p2-1+1,u2,h2,u1,h1,u2,h2] = Float32(!p1Big * a2 * alpha2)
        else
            for p4 in axes(GainMatrix4,1)
                for u4 in axes(GainMatrix4,2), h4 in axes(GainMatrix4,3) 
                    if !p1Big && CartesianIndex(p4, u4, h4) in @view(cart_inds4[1:high_4_range])
                        #if GainMatrix4[p4,u4,h4,p1,u1,h1,p2,u2,h2] > tol * LossSumN2
                            GainMatrix4[p4,u4,h4,u1,h1,u2,h2] *= Float32(beta)
                        #else
                        #    CorrectedGainMatrix4[p4,u4,h4,p1,u1,h1,p2,u2,h2] = 0.0
                        #end
                    else
                        #if GainMatrix4[p4,u4,h4,p1,u1,h1,p2,u2,h2] > tol * LossSumN2
                            GainMatrix4[p4,u4,h4,u1,h1,u2,h2] *= Float32(alpha2)
                        #else
                        #    CorrectedGainMatrix4[p4,u4,h4,p1,u1,h1,p2,u2,h2] = 0.0
                        #end
                    end
                end
            end
        end

        # check correction has worked for this bin

        #num_check = sum( @view(GainMatrix3[:,:,:,u1,h1,u2,h2])) + sum(@view(GainMatrix4[:,:,:,u1,h1,u2,h2])) - 2*LossMatrix[u1,h1,u2,h2]

        #eng_check = sum(E3_d .* @view(GainMatrix3[:,:,:,u1,h1,u2,h2])) + sum(E4_d .* @view(GainMatrix4[:,:,:,u1,h1,u2,h2])) - E1_d[p1] * LossMatrix[u1,h1,u2,h2] - E2_d[p2] * LossMatrix[u1,h1,u2,h2]

        #if !isapprox(num_check,0e0,atol=10sqrt(eps(Float32))) || !isapprox(eng_check,0e0,atol=10sqrt(eps(Float32)))
        #    println("Correction failed for p1=$p1,p2=$p2, u1=$u1, u2=$u2, h1=$h1, h2=$h2, Number check: $num_check, Energy check: $eng_check")
        #    println("alpha1: $alpha1, alpha2: $alpha2, beta: $beta, max_high_bins: $max_high_bins")
        #    println("a1: $a1, b1: $b1, ac1: $ac1, bd1: $bd1, e1: $e1, a2: $a2, b2: $b2, ac2: $ac2, bd2: $bd2, e2: $e2, f: $f, c1: $c1, d1: $d1, c2: $c2, d2: $d2")
        #end

        #println("p1=$p1,p2=$p2,u1=$u1,u2=$u2,h1=$h1, h2=$h2")

        #println("$(sum(CorrectedGainMatrix3[:,:,:,p1,u1,h1,p2,u2,h2])) , $(sum(CorrectedGainMatrix4[:,:,:,p1,u1,h1,p2,u2,h2])) , $(CorrectedLossMatrix1[p1,u1,h1,p2,u2,h2]) , $(CorrectedLossMatrix2[p2,u2,h2,p1,u1,h1])")

    end
    #println("Number of bins that couldn't be corrected = $num_wrong")
    #println("Number of bins corrected = $num_right")

    #return nothing #CorrectedGainMatrix3, CorrectedGainMatrix4, CorrectedLossMatrix1, CorrectedLossMatrix2

end


"""
    partialsortperm_no_view!(ix, v, k; by=identity, lt=isless, rev=false)

Like [`partialsortperm!`](@ref), but returns nothing.
"""
function partialsortperm_no_view!(ix::AbstractVector{<:Integer}, v::AbstractVector,k::Union{Integer, OrdinalRange};lt::Function=isless,by::Function=identity,rev::Union{Bool,Nothing}=nothing,order::Base.Ordering=Base.Forward,initialized::Bool=false)

    if axes(ix,1) != axes(v,1)
        throw(ArgumentError("The index vector is used as scratch space and must have the " * "same length/indices as the source vector, $(axes(ix,1)) != $(axes(v,1))"))
    end
    @inbounds for i in eachindex(ix)
        ix[i] = i
    end

    # do partial quicksort
    Base.Sort._sort!(ix, Base.Sort.InitialOptimizations(Base.Sort.ScratchQuickSort(k)), Base.Sort.Perm(Base.Sort.ord(lt, by, rev, order), v), (;))

    return nothing
end


###### LEGACY CODE BELOW, NOT USED ANYMORE BUT KEPT FOR REFERENCE ######

#=

"""
    GainCorrection(Parameters,GainMatrix3,GainMatrix4,LossMatrix1,LossMatrix2)

MC sampling introduces noise that can lead to poor number and energy conservation. `GainCorrection` provides a corrective step to ensure number and energy conservation to numerical precision. If there is no GainMatrix element then the value of the LossMatrix is applied to the same bin as the input state (if they are identical particles); if not identical particles if there is no GainMatrix element then the value of the LossMatrix is set to zero to ensure particle conservation (with good MC sampling this should rarely occur).
"""
function GainCorrection(Parameters::Tuple{String, String, String, String, Float64, Float64, Float64, Float64, Float64, Float64, String, Int64, String, Int64, String, Int64, Float64, Float64, String, Int64, String, Int64, String, Int64, Float64, Float64, String, Int64, String, Int64, String, Int64, Float64, Float64, String, Int64, String, Int64, String, Int64}, GainMatrix3::Array{Float64, 9}, GainMatrix4::Array{Float64, 9}, LossMatrix1::Array{Float64, 6}, LossMatrix2::Array{Float64, 6})

    
    #= Different possible combinations of identical and different particles that affect how to apply conservation corrections:

        1. n1=n3 != n2=n4 (1 and 3 identical, 2 and 4 identical)
            3 parameters needed:
                alpha: consv. num for 1 and 3
                beta: consv. num for 2 and 4
                gamma: consv. energy for all
            In this setup gain of n3 or n4 could be zero from MC, so have to be careful when applying correction    
        2. all other combinations:
            2 parameters needed:
                alpha: consv. num for all
                beta: consv. energy for all

    =#

    (name1,name2,name3,name4,m1,m2,m3,m4,p1_low,p1_up,p1_grid,p1_num,u1_grid,u1_num,h1_grid,h1_num,p2_low,p2_up,p2_grid,p2_num,u2_grid,u2_num,h2_grid,h2_num,p3_low,p3_up,p3_grid,p3_num,u3_grid,u3_num,h3_grid,h3_num,p4_low,p4_up,p4_grid,p4_num,u4_grid,u4_num,h4_grid,h4_num) = Parameters

    tol = sqrt(eps(Float64)) # tolerance for how low the gain terms can be compared to the loss terms, if less than this tolerance then the gain terms are set to zero and not included in the correction calculation.

    if (name1 == name3) && (name2 == name4) && (name1 != name2)
        CorrType = 1
    else
        CorrType = 2
    end

    CorrectedGainMatrix3 = similar(GainMatrix3)
    CorrectedGainMatrix4 = similar(GainMatrix4)
    CorrectedLossMatrix1 = similar(LossMatrix1)
    CorrectedLossMatrix2 = similar(LossMatrix2)
    fill!(CorrectedGainMatrix3,Float64(0))
    fill!(CorrectedGainMatrix4,Float64(0))
    CorrectedLossMatrix1 .= LossMatrix1
    CorrectedLossMatrix2 .= LossMatrix2

    # underflow and overflow bins are taken to have size 0 -> p1_r[1] and p1_r[end] -> 2*p1_r[end] respectively

    p1_r = bounds(p1_low,p1_up,p1_num,p1_grid);
    p1_d = deltaVector(p1_r);
    u1_r = bounds(-1.0,1.0,u1_num,u1_grid);
    u1_d = deltaVector(u1_r);

    p2_r = bounds(p2_low,p2_up,p2_num,p2_grid);
    p2_d = deltaVector(p2_r);
    u2_r = bounds(-1.0,1.0,u2_num,u2_grid);
    u2_d = deltaVector(u2_r);

    p3_r = bounds(p3_low,p3_up,p3_num,p3_grid);
    p3_r = [0.0 ; p3_r ; 2*p3_r[end]];
    p3_d = deltaVector(p3_r);
    u3_r = bounds(-1.0,1.0,u3_num,u3_grid);
    u3_d = deltaVector(u3_r);

    p4_r = bounds(p4_low,p4_up,p4_num,p4_grid);
    p4_r = [0.0 ; p4_r ; 2*p4_r[end]];
    p4_d = deltaVector(p4_r);
    u4_r = bounds(-1.0,1.0,u4_num,u4_grid);
    u4_d = deltaVector(u4_r);

    E1_Δ = deltaEVector(p1_r,m1);
    E1_d = E1_Δ ./ p1_d

    E2_Δ = deltaEVector(p2_r,m2);
    E2_d = E2_Δ ./ p2_d

    E3_Δ = deltaEVector(p3_r,m3);
    E3_d = E3_Δ ./ p3_d

    E4_Δ = deltaEVector(p4_r,m4);
    E4_d = E4_Δ ./ p4_d

    num_wrong = 0
    num_right = 0

    Threads.@threads for idx in CartesianIndices((axes(GainMatrix3,4),axes(GainMatrix3,5),axes(GainMatrix3,6),axes(GainMatrix3,7),axes(GainMatrix3,8),axes(GainMatrix3,9)))

        p1, u1, h1, p2, u2, h2 = Tuple(idx)

        println("Correcting p1=$p1,p2=$p2, u1=$u1, u2=$u2, h1=$h1, h2=$h2 on thread $(Threads.threadid())")

        GainMatrix3Filtered = zeros(size(GainMatrix3)[1:3])
        GainMatrix4Filtered = zeros(size(GainMatrix4)[1:3])

        GainMatrix3Filtered .= @view(GainMatrix3[:,:,:,p1,u1,h1,p2,u2,h2])
        GainMatrix4Filtered .= @view(GainMatrix4[:,:,:,p1,u1,h1,p2,u2,h2])
        
        LossSumN1 = LossMatrix1[p1,u1,h1,p2,u2,h2]
        LossSumN2 = LossMatrix2[p2,u2,h2,p1,u1,h1]
        if name1 != name2 # Indistinguishable_12==false so LossMatrix2 is non-zero
            LossSumE1 = LossMatrix1[p1,u1,h1,p2,u2,h2]*E1_d[p1]
            LossSumE2 = LossMatrix2[p2,u2,h2,p1,u1,h1]*E2_d[p2]
        else
            LossSumE1 = LossMatrix1[p1,u1,h1,p2,u2,h2]*(E1_d[p1]+E2_d[p2])/2
            LossSumE2 = 0.0#LossMatrix1[p1,u1,h1,p2,u2,h2]*(E1_d[p1]+E2_d[p2])/2
        end

        if LossSumN1 + LossSumN2 == 0e0
            # no loss term so no interaction at these incoming states

            num_right += 1

            continue
        end
        
        wrong = true
        nonzero_gain = true
        max_high_bins = 0
        Gain3False = false
        Gain4False = false

        alpha1 = 0.0
        alpha2 = 0.0
        beta = 0.0

        a1 = 0.0
        b1 = 0.0
        a2 = 0.0
        b2 = 0.0

        p1Big = false
        p2Big = false

        p3_offset = length(axes(GainMatrix3,1))
        p4_offset = length(axes(GainMatrix4,1))

        #=if E1_d[p1] > 1e0*E2_d[p2] && m1 > m2
            p1Big = true
        elseif E2_d[p2] > 1e0*E1_d[p1] && m2 > m1
            p2Big = true
        end=#

        cart_inds3 = CartesianIndices(@view(GainMatrix3[:,:,:,p1,u1,h1,p2,u2,h2]))[(1)]
        cart_inds4 = CartesianIndices(@view(GainMatrix4[:,:,:,p1,u1,h1,p2,u2,h2]))[(1)]

        while wrong && nonzero_gain

            GainSumN31 = zero(Float64)
            GainSumN32 = zero(Float64)
            GainSumN41 = zero(Float64)
            GainSumN42 = zero(Float64)
            GainSumE31 = zero(Float64)
            GainSumE32 = zero(Float64)
            GainSumE41 = zero(Float64)
            GainSumE42 = zero(Float64)
            
            # search approach is to first find the bin with the highest GainMatrix value 
            high_bins = 0
            searching = true
            p3_offset = length(axes(GainMatrix3,1))

            # Apply filters
            for u3 in axes(GainMatrix3,2), h3 in axes(GainMatrix3,3) 
                # fixes points that are break monotonically increasing regions
                fix_monotone_center_log!(@view(GainMatrix3Filtered[:,u3,h3])) 
            end
            
            high_inds3 = partialsortperm(reshape(E3_d .* GainMatrix3Filtered,((p3_num+2)*u3_num*h3_num)), 1:min(max_high_bins+1,size(GainMatrix3,1)); rev=true)
            cart_inds3 = CartesianIndices(GainMatrix3Filtered)[high_inds3]

            for p3 in axes(GainMatrix3,1), u3 in axes(GainMatrix3,2), h3 in axes(GainMatrix3,3) 
                #if GainMatrix3[p3,u3,h3,p1,u1,h1,p2,u2,h2] > tol * LossSumN1 # DONT add tolerance as it cuts low gain terms at low energy that might be needed for Graph Laplacian structure
                    tmpN = GainMatrix3Filtered[p3,u3,h3]
                    tmpE = GainMatrix3Filtered[p3,u3,h3]*E3_d[p3]
                #else 
                #    continue
                #end
                if !p2Big && CartesianIndex(p3, u3, h3) in cart_inds3
                    GainSumN32 += tmpN
                    GainSumE32 += tmpE
                else
                    GainSumN31 += tmpN
                    GainSumE31 += tmpE
                end
            end

            # Apply filters
            for u4 in axes(GainMatrix4,2), h4 in axes(GainMatrix4,3) 
                # fixes points that are break monotonically increasing regions
                fix_monotone_center_log!(@view(GainMatrix4Filtered[:,u4,h4])) 
            end

            high_inds4 = partialsortperm(reshape(E4_d .* GainMatrix4Filtered,((p4_num+2)*u4_num*h4_num)), 1:min(max_high_bins+1,size(GainMatrix4,1)); rev=true)
            cart_inds4 = CartesianIndices(GainMatrix4Filtered)[high_inds4]

            for p4 in axes(GainMatrix4,1), u4 in axes(GainMatrix4,2), h4 in axes(GainMatrix4,3) 
                #if GainMatrix4[p4,u4,h4,p1,u1,h1,p2,u2,h2] > tol * LossSumN2 # DONT add tolerance as it cuts low gain terms at low energy that might be needed for Graph Laplacian structure
                    tmpN = GainMatrix4Filtered[p4,u4,h4]
                    tmpE = GainMatrix4Filtered[p4,u4,h4]*E4_d[p4]
                #else
                #    continue
                #end
                if !p1Big && CartesianIndex(p4, u4, h4) in cart_inds4
                    GainSumN42 += tmpN
                    GainSumE42 += tmpE
                else
                    GainSumN41 += tmpN
                    GainSumE41 += tmpE
                end
            end

            if CorrType == 1 

                a1 = GainSumN31
                b1 = GainSumN32

                a2 = GainSumN41
                b2 = GainSumN42

                ac1 = GainSumE31
                c1 = ac1 != 0.0 ? GainSumE31/GainSumN31 : 0.0
                bd1 = GainSumE32
                d1 = bd1 != 0.0 ? GainSumE32/GainSumN32 : 0.0

                ac2 = GainSumE41
                c2 = ac2 != 0.0 ? GainSumE41/GainSumN41 : 0.0
                bd2 = GainSumE42
                d2 = bd2 != 0.0 ? GainSumE42/GainSumN42 : 0.0
                
                if (a1+b1 == 0e0 || (max_high_bins==0 && a1==0.0 && b1!=0.0))# There is not gain term produced by MC or only a single bin
                    #c1 = p2Big ? E3_d[p1] : E3_d[p1-1]
                    #d1 = p2Big ? E3_d[p1+1] : E3_d[p1]
                    # underflow bins adds 1 to value of p3 index compared to p1 index for same energy bin.
                    c1 = p2Big ? E3_d[p1+1] : E3_d[p1]
                    d1 = p2Big ? E3_d[p1+2] : E3_d[p1+1]
                    L = LossSumN1
                    LE = LossSumE1+LossSumE2
                    e2 = a2 + b2 - LossSumN2
                    #c2 = GainSumE41/GainSumN41

                    a1 = ((ac2-c2*e2)+(d1*L-LE))/(d1-c1)
                    b1 = L-a1
                    ac1 = a1*c1
                    bd1 = b1*d1

                    Gain3False = true
                end

                if (a2+b2 == 0e0 || (max_high_bins==0 && a2==0.0 && b2!=0.0))# There is not gain term produced by MC or only a single bin
                    #c2 = p1Big ? E4_d[p2] : E4_d[p2-1]
                    #d2 = p1Big ? E4_d[p2+1] : E4_d[p2]
                    # underflow bins adds 1 to value of p4 index compared to p2 index for same energy bin.
                    c2 = p1Big ? E4_d[p2+1] : E4_d[p2]
                    d2 = p1Big ? E4_d[p2+2] : E4_d[p2+1]
                    L = LossSumN2
                    LE = LossSumE2 + LossSumE1
                    e1 = a1 + b1 - LossSumN1
                    #c1 = GainSumE31/GainSumN31

                    a2 = ((ac1-c1*e1)+(d2*L-LE))/(d2-c2)
                    b2 = L-a2
                    ac2 = a2*c2
                    bd2 = b2*d2
                    Gain4False = true
                end

                e1 = a1 + b1 - LossSumN1
                e2 = a2 + b2 - LossSumN2
                f  = ac1 + bd1 + ac2 + bd2 - LossSumE1 - LossSumE2

                alpha1 = (b2*(d2-c2)*e1+b1*(d1*e1+c2*e2-f))/(a1*(b1*(c1-d1)+b2*(c2-d2))) + 1
                alpha2 = (b1*(d1-c1)*e2+b2*(d2*e2+c1*e1-f))/(a2*(b1*(c1-d1)+b2*(c2-d2))) + 1
                beta = (f-c1*e1-c2*e2)/(b1*(c1-d1)+b2*(c2-d2)) + 1

                if a1 <= 0e0 || a2 <= 0e0 || max_high_bins > max(size(GainMatrix3,1),size(GainMatrix4,1))# loop has not produced any corrected gain terms

                    nonzero_gain = false
                    #println("No valid correction for p1=$p1,p2=$p2, u1=$u1, u2=$u2, h1=$h1, h2=$h2")
                    #println("a1: $a1, b1: $b1, ac1: $ac1, bd1: $bd1, e1: $e1, a2: $a2, b2: $b2, ac2: $ac2, bd2: $bd2, e2: $e2, f: $f, c1: $c1, d1: $d1, c2: $c2, d2: $d2")
                    #println("GainSumN31: $GainSumN31, GainSumN32: $GainSumN32, GainSumN41: $GainSumN41, GainSumN42: $GainSumN42, LossSumN1: $LossSumN1, LossSumN2: $LossSumN2")
                    #println("alpha1: $alpha1, alpha2: $alpha2, beta: $beta, max_high_bins: $max_high_bins")
                    alpha1 = 0.0
                    alpha2 = 0.0
                    beta = 0.0
                    CorrectedLossMatrix1[p1,u1,h1,p2,u2,h2] = 0e0
                    CorrectedLossMatrix2[p2,u2,h2,p1,u1,h1] = 0e0

                    num_wrong += 1

                    continue

                end

                if isnan(alpha1) || isnan(alpha2) || isnan(beta)
                    println("NaN encountered for p1=$p1,p2=$p2, u1=$u1, u2=$u2, h1=$h1, h2=$h2")
                    println("a1: $a1, b1: $b1, ac1: $ac1, bd1: $bd1, e1: $e1, a2: $a2, b2: $b2, ac2: $ac2, bd2: $bd2, e2: $e2, f: $f, c1: $c1, d1: $d1, c2: $c2, d2: $d2")
                    println("GainSumN31: $GainSumN31, GainSumN32: $GainSumN32, GainSumN41: $GainSumN41, GainSumN42: $GainSumN42")
                    #println("$(a+b-LossN), $(ac+bd-LossE)")
                    #println("LossE/LossN: $(LossE/LossN), bd/b: $(bd/b), ac/a: $(ac/a)")
                end

                if  alpha1 < 0e0 || alpha2 < 0e0 || beta < 0e0

                    #println("p1=$p1,p2=$p2, u1=$u1, u2=$u2, h1=$h1, h2=$h2")
                    #println("alpha1: $alpha1, alpha2: $alpha2, beta: $beta, max_high_bins: $max_high_bins, p3_offset: $p3_offset, p4_offset: $p4_offset")

                    max_high_bins += 1

                    continue # redo calculation
                else
                    #println("p1=$p1,p2=$p2")
                    #println("alpha1: $alpha1, alpha2: $alpha2, beta: $beta, max_high_bins: $max_high_bins, p3_offset: $p3_offset, p4_offset: $p4_offset")
                    #println("new gain1 = $(a1*alpha1+b1*beta), new gain2 = $(a2*alpha2+b2*beta), new gain = $((a1*alpha1+b1*beta)+(a2*alpha2+b2*beta)), Loss = $(LossSumN1+LossSumN2), err = $((a1*alpha1+b1*beta)+(a2*alpha2+b2*beta)- (LossSumN1+LossSumN2)), L1 = $LossSumN1, L2 = $LossSumN2")
                    wrong = false
                    num_right += 1
                end

            elseif CorrType == 2

                GainN1 = GainSumN31 + GainSumN41
                GainE1 = GainSumE31 + GainSumE41
                GainN2 = GainSumN32 + GainSumN42
                GainE2 = GainSumE32 + GainSumE42
                LossN = LossSumN1 + LossSumN2
                LossE = LossSumE1 + LossSumE2

                a = GainN1
                if a == 0e0 || max_high_bins > max(size(GainMatrix3,1),size(GainMatrix4,1))# loop has not produced any gain terms
                    nonzero_gain = false
                    println("No valid correction for p1=$p1,p2=$p2, u1=$u1, u2=$u2, h1=$h1, h2=$h2")

                    alpha1 = 0.0
                    alpha2 = 0.0
                    beta = 0.0
                    CorrectedLossMatrix1[p1,u1,h1,p2,u2,h2] = 0e0
                    CorrectedLossMatrix2[p2,u2,h2,p1,u1,h1] = 0e0

                    num_wrong += 1
                    
                    continue
                end

                ac = GainE1
                c = GainE1/GainN1

                b = GainN2
                bd = GainE2
                d = GainE2/GainN2

                e = a+b-LossN
                f = ac+bd-LossE

                if isapprox(c,d)
                    # In this case deltaE is the approx the same for all bins involved, so energy is approx conserved if number is conserved so only one equation and one parameter needed which we take scales b
                    #alpha = -e/(a+b) + 1
                    #beta = alpha
                    beta = -e/b +1
                    alpha = 1.0
                else
                    alpha = (e*d-f)/(ac-a*d)+1
                    beta = (e*c-f)/(bd-b*c)+1
                end

                if isnan(alpha) || isnan(beta)
                    println("NaN encountered for p1=$p1,p2=$p2, u1=$u1, u2=$u2, h1=$h1, h2=$h2")
                    println("a: $a, b: $b, ac: $ac, bd: $bd, e: $e, f: $f, c: $c, d: $d")
                    println("$(a+b-LossN), $(ac+bd-LossE)")
                    println("LossE/LossN: $(LossE/LossN), bd/b: $(bd/b), ac/a: $(ac/a)")
                end

                if alpha < 0e0 || beta < 0e0

                    #println("alpha: $alpha, beta: $beta, max_high_bins: $max_high_bins, p3_offset: $p3_offset, p4_offset: $p4_offset")

                    max_high_bins += 1

                    continue # redo calculation
                else
                    #println("p1=$p1,p2=$p2")
                    #println("alpha: $alpha, beta: $beta, max_high_bins: $max_high_bins, p3_offset: $p3_offset, p4_offset: $p4_offset")
                    wrong = false
                    num_right += 1
                    alpha1 = alpha
                    alpha2 = alpha
                end

            end # if CorrType



        end # while

        if Gain3False 
            # underflow bins adds 1 to value of p3 index compared to p1 index for same energy bin.
            CorrectedGainMatrix3[p1+1+1,u1,h1,p1,u1,h1,p2,u2,h2] = p2Big * b1 * beta
            CorrectedGainMatrix3[p1+1,u1,h1,p1,u1,h1,p2,u2,h2] = p2Big ? a1 * alpha1 : b1 * beta
            CorrectedGainMatrix3[p1-1+1,u1,h1,p1,u1,h1,p2,u2,h2] = !p2Big * a1 * alpha1  
        else
            for p3 in axes(GainMatrix3,1)
                for u3 in axes(GainMatrix3,2), h3 in axes(GainMatrix3,3) 
                    if CartesianIndex(p3, u3, h3) in cart_inds3 && !p2Big
                        #if GainMatrix3[p3,u3,h3,p1,u1,h1,p2,u2,h2] > tol * LossSumN1
                            CorrectedGainMatrix3[p3,u3,h3,p1,u1,h1,p2,u2,h2] = GainMatrix3Filtered[p3,u3,h3] * beta
                        #else
                        #    CorrectedGainMatrix3[p3,u3,h3,p1,u1,h1,p2,u2,h2] = 0.0
                        #end
                    else
                        #if GainMatrix3[p3,u3,h3,p1,u1,h1,p2,u2,h2] > tol * LossSumN1
                            CorrectedGainMatrix3[p3,u3,h3,p1,u1,h1,p2,u2,h2] = GainMatrix3Filtered[p3,u3,h3] * alpha1
                        #else
                        #    CorrectedGainMatrix3[p3,u3,h3,p1,u1,h1,p2,u2,h2] = 0.0
                        #end
                    end
                end
            end
        end

        if Gain4False 
            # underflow bins adds 1 to value of p4 index compared to p2 index for same energy bin.
            CorrectedGainMatrix4[p2+1+1,u2,h2,p1,u1,h1,p2,u2,h2] = p1Big * b2 * beta
            CorrectedGainMatrix4[p2+1,u2,h2,p1,u1,h1,p2,u2,h2] = p1Big ? a2 * alpha2 : b2 * beta
            CorrectedGainMatrix4[p2-1+1,u2,h2,p1,u1,h1,p2,u2,h2] = !p1Big * a2 * alpha2
        else
            for p4 in axes(GainMatrix4,1)
                for u4 in axes(GainMatrix4,2), h4 in axes(GainMatrix4,3) 
                    if CartesianIndex(p4, u4, h4) in cart_inds4 && !p1Big
                        #if GainMatrix4[p4,u4,h4,p1,u1,h1,p2,u2,h2] > tol * LossSumN2
                            CorrectedGainMatrix4[p4,u4,h4,p1,u1,h1,p2,u2,h2] = GainMatrix4Filtered[p4,u4,h4] * beta
                        #else
                        #    CorrectedGainMatrix4[p4,u4,h4,p1,u1,h1,p2,u2,h2] = 0.0
                        #end
                    else
                        #if GainMatrix4[p4,u4,h4,p1,u1,h1,p2,u2,h2] > tol * LossSumN2
                            CorrectedGainMatrix4[p4,u4,h4,p1,u1,h1,p2,u2,h2] = GainMatrix4Filtered[p4,u4,h4] * alpha2
                        #else
                        #    CorrectedGainMatrix4[p4,u4,h4,p1,u1,h1,p2,u2,h2] = 0.0
                        #end
                    end
                end
            end
        end

        # check correction has worked for this bin

        num_check = sum( @view(CorrectedGainMatrix3[:,:,:,p1,u1,h1,p2,u2,h2])) + sum(@view(CorrectedGainMatrix4[:,:,:,p1,u1,h1,p2,u2,h2])) - CorrectedLossMatrix1[p1,u1,h1,p2,u2,h2] - CorrectedLossMatrix2[p2,u2,h2,p1,u1,h1]

        eng_check = sum(E3_d .* @view(CorrectedGainMatrix3[:,:,:,p1,u1,h1,p2,u2,h2])) + sum(E4_d .* @view(CorrectedGainMatrix4[:,:,:,p1,u1,h1,p2,u2,h2])) - E1_d[p1] * CorrectedLossMatrix1[p1,u1,h1,p2,u2,h2] - E2_d[p2] * CorrectedLossMatrix2[p2,u2,h2,p1,u1,h1]

        #=if !isapprox(num_check,0e0,atol=10sqrt(eps(Float64))) || !isapprox(eng_check,0e0,atol=10sqrt(eps(Float64)))
            println("Correction failed for p1=$p1,p2=$p2, u1=$u1, u2=$u2, h1=$h1, h2=$h2, Number check: $num_check, Energy check: $eng_check")
            #println("alpha1: $alpha1, alpha2: $alpha2, beta: $beta, max_high_bins: $max_high_bins")
            #println("a1: $a1, b1: $b1, ac1: $ac1, bd1: $bd1, e1: $e1, a2: $a2, b2: $b2, ac2: $ac2, bd2: $bd2, e2: $e2, f: $f, c1: $c1, d1: $d1, c2: $c2, d2: $d2")
        end=#

        #println("p1=$p1,p2=$p2,u1=$u1,u2=$u2,h1=$h1, h2=$h2")

        #println("$(sum(CorrectedGainMatrix3[:,:,:,p1,u1,h1,p2,u2,h2])) , $(sum(CorrectedGainMatrix4[:,:,:,p1,u1,h1,p2,u2,h2])) , $(CorrectedLossMatrix1[p1,u1,h1,p2,u2,h2]) , $(CorrectedLossMatrix2[p2,u2,h2,p1,u1,h1])")

    end
    println("Number of bins that couldn't be corrected = $num_wrong")
    println("Number of bins corrected = $num_right")

    return CorrectedGainMatrix3, CorrectedGainMatrix4, CorrectedLossMatrix1, CorrectedLossMatrix2

end


=#