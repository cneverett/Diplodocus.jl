function hampel_clip_log!(y::AbstractVector{<:Real};window::Int = 7,n_sigma::Real = 3.0)

    @assert isodd(window) "window size must be odd"

    k = div(window, 2)

    # Work on a copy so earlier corrections don't affect later windows
    y_orig = copy(y)

    for i in (k+1):(length(y)-k)

        # Leave zeros untouched
        y_orig[i] <= 0 && continue

        # Positive values only
        w = @view y_orig[(i-k):(i+k)]
        wpos = filter(>(0), w)

        # Need enough points to compute robust statistics
        length(wpos) < 3 && continue

        logw = log.(wpos)

        med = median(logw)
        mad = median(abs.(logw .- med))

        sigma = 1.4826 * mad

        sigma <= 0 && continue

        lower = med - n_sigma * sigma
        upper = med + n_sigma * sigma

        logyi = log(y_orig[i])

        # Clip in log-space
        logyi_clipped = clamp(logyi, lower, upper)

        y[i] = exp(logyi_clipped)
    end

    return nothing
end



function hampel_clip_log_run!(y::AbstractVector{<:Real};
                          window::Int = 7,
                          n_sigma::Real = 3.0,
                          min_run::Int = 2)

    @assert isodd(window) "window size must be odd"
    @assert min_run >= 1

    n = length(y)
    k = div(window, 2)

    # Work on a copy for decisions, but modify y in place.
    y0 = collect(float.(y))

    # log-space view; zeros are ignored with NaN
    z = Vector{Float64}(undef, n)
    for i in 1:n
        z[i] = y0[i] > 0 ? log(y0[i]) : NaN
    end

    suspect = falses(n)

    # First pass: mark outliers in log-space
    for i in (k+1):(n-k)
        isfinite(z[i]) || continue

        w = @view z[(i-k):(i+k)]
        w = w[isfinite.(w)]
        length(w) < 3 && continue

        med = median(w)
        mad = median(abs.(w .- med))
        sigma = 1.4826 * mad

        sigma <= 0 && continue

        if abs(z[i] - med) > n_sigma * sigma
            suspect[i] = true
        end
    end

    # Second pass: only correct short isolated runs
    i = 1
    while i <= n
        if !suspect[i]
            i += 1
            continue
        end

        j = i
        while j <= n && suspect[j]
            j += 1
        end
        runlen = j - i

        if runlen < min_run
            # Clip each point in the short run using its own local window
            for p in i:(j-1)
                p <= k && continue
                p > n-k && continue
                isfinite(z[p]) || continue

                w = @view z[(p-k):(p+k)]
                w = w[isfinite.(w)]
                length(w) < 3 && continue

                med = median(w)
                mad = median(abs.(w .- med))
                sigma = 1.4826 * mad
                sigma <= 0 && continue

                lower = med - n_sigma * sigma
                upper = med + n_sigma * sigma
                z[p] = clamp(z[p], lower, upper)
            end
        end

        i = j
    end

    # Write back in original scale, preserving zeros
    for i in 1:n
        if isfinite(z[i])
            y[i] = exp(z[i])
        end
    end

    return nothing
end

"""
    localcurve_clip_singletons_log!(y; window=9, degree=2, n_sigma=3.0, floor=1e-300)

In-place robust outlier correction in log-space.

- `y` must be nonnegative.
- Zeros are ignored and left unchanged.
- Each point is compared to a local polynomial fit of degree `degree`
  built from its neighbors in a sliding window.
- Only isolated single-point outliers are corrected.
- Returns `y` and a Bool mask of suspected outliers.

The correction is performed in log-space, then mapped back with `exp`.
"""
function localcurve_clip_singletons_log!(y::AbstractVector{<:Real};
                                         window::Int = 9,
                                         degree::Int = 2,
                                         n_sigma::Real = 3.0,
                                         ridge::Real = 1e-10)

    @assert isodd(window) "window must be odd"
    @assert window >= 3 "window must be at least 3"
    @assert degree >= 1 "degree must be at least 1"
    @assert ridge > 0 "ridge must be positive"

    n = length(y)
    k = div(window, 2)

    # Work on a copy for decisions, but modify y in place
    y0 = collect(float.(y))

    # Log-space data; zeros are ignored
    z = Vector{Float64}(undef, n)
    for i in 1:n
        z[i] = y0[i] > 0 ? log(y0[i]) : NaN
    end

    suspect = falses(n)
    zhat = fill(NaN, n)
    sigma = fill(NaN, n)

    for i in 1:n
        lo = max(1, i - k)
        hi = min(n, i + k)

        idx = Int[]
        for j in lo:hi
            if j != i && isfinite(z[j])
                push!(idx, j)
            end
        end

        # Need enough points to fit at least a linear model
        if length(idx) < degree + 1
            continue
        end

        # Use a lower degree if the window is small
        p = min(degree, length(idx) - 1)

        # Local coordinates centered at i
        x = Float64.(idx .- i)
        xscale = maximum(abs.(x))
        xscale == 0 && continue
        u = x ./ xscale

        # Design matrix [1, u, u^2, ...]
        X = hcat([u .^ d for d in 0:p]...)

        # Positive weights that never vanish completely
        w = exp.(-2.0 .* u.^2)
        sw = sqrt.(w)

        Aw = X .* sw
        bw = z[idx] .* sw

        # Ridge-regularized weighted least squares:
        # (A'A + ridge*I) β = A'b
        G = Aw' * Aw + ridge * I(size(Aw, 2))
        g = Aw' * bw

        β = G \ g

        # Fitted value at the center (u = 0)
        zhat[i] = β[1]

        resid = bw - Aw * β
        medr = median(resid)
        madr = median(abs.(resid .- medr))
        sigma[i] = 1.4826 * madr

        if isfinite(z[i]) && isfinite(zhat[i]) && sigma[i] > 0 &&
           abs(z[i] - zhat[i]) > n_sigma * sigma[i]
            suspect[i] = true
        end
    end

    # Only correct isolated single-point outliers
    for i in 2:(n - 1)
        if suspect[i] && !suspect[i - 1] && !suspect[i + 1] &&
           isfinite(z[i]) && isfinite(zhat[i]) && isfinite(sigma[i]) && sigma[i] > 0

            z[i] = clamp(z[i], zhat[i] - n_sigma * sigma[i], zhat[i] + n_sigma * sigma[i])
        end
    end

    # Write back in original scale, preserving zeros
    for i in 1:n
        if isfinite(z[i])
            y[i] = exp(z[i])
        end
    end

    return y, suspect
end

"""
    detect_and_fix_spikes_log!(y; window=9, degree=2, n_sigma=3.0, ridge=1e-10)

In-place spike detection and correction for nonnegative data.

Method:
- work in log-space: z = log(y), with zeros ignored
- fit a local polynomial to the neighborhood of each point
- flag points whose residual is large relative to a robust local scale
- only correct flagged points that are isolated singletons

Returns:
- `suspect::BitVector`    : all points flagged as suspicious
- `corrected::BitVector`   : points actually modified
- `y`                     : modified in place

Notes:
- `window` must be odd
- `degree=1` or `2` is usually enough
- `ridge` prevents singular fits
"""
function detect_and_fix_spikes_log!(y::AbstractVector{<:Real};
                                         window::Int = 9,
                                         degree::Int = 2,
                                         n_sigma::Real = 3.0,
                                         ridge::Real = 1e-10)

    @assert isodd(window) "window must be odd"
    @assert window >= 3 "window must be at least 3"
    @assert degree >= 1 "degree must be at least 1"

    n = length(y)
    k = div(window, 2)

    # Work on a copy for detection
    y0 = collect(float.(y))

    # Log-space representation
    z = fill(NaN, n)
    for i in 1:n
        if y0[i] > 0
            z[i] = log(y0[i])
        end
    end

    suspect = falses(n)

    #
    # First pass: detect suspicious points
    #
    for i in 1:n

        lo = max(1, i-k)
        hi = min(n, i+k)

        idx = Int[]

        for j in lo:hi
            if j != i && isfinite(z[j])
                push!(idx, j)
            end
        end

        length(idx) < degree + 1 && continue

        p = min(degree, length(idx)-1)

        x = Float64.(idx .- i)

        xmax = maximum(abs.(x))
        xmax == 0 && continue

        x ./= xmax

        X = hcat([x.^d for d in 0:p]...)

        w = exp.(-2 .* x.^2)
        sw = sqrt.(w)

        Aw = X .* sw
        bw = z[idx] .* sw

        G = Aw'Aw + ridge * I(size(Aw,2))
        g = Aw'bw

        β = G \ g

        zpred = β[1]

        resid = bw - Aw*β

        medr = median(resid)
        madr = median(abs.(resid .- medr))

        σ = 1.4826 * madr

        if isfinite(z[i]) && σ > 0 &&
           abs(z[i] - zpred) > n_sigma * σ

            suspect[i] = true
        end
    end

    #
    # Second pass: only correct isolated singleton outliers
    #
    corrected = falses(n)

    for i in 2:(n-1)

        if suspect[i] &&
           !suspect[i-1] &&
           !suspect[i+1]

            lo = max(1, i-k)
            hi = min(n, i+k)

            vals = Float64[]

            for j in lo:hi
                if j != i && isfinite(z[j])
                    push!(vals, z[j])
                end
            end

            if !isempty(vals)
                z[i] = mean(vals)   # mean in log-space
                corrected[i] = true
            end
        end
    end

    #
    # Back to original space
    #
    for i in 1:n
        if isfinite(z[i])
            y[i] = exp(z[i])
        end
    end

    return suspect, corrected
end

function detect_and_fix_spikes_log_mean!(y::AbstractVector{<:Real};
                                         window::Int = 9,
                                         degree::Int = 2,
                                         n_sigma::Real = 3.0,
                                         ridge::Real = 1e-10,
                                         support_radius::Int = 2,
                                         min_support::Int = 2,
                                         support_tol::Real = 0.5)

    @assert isodd(window) "window must be odd"
    @assert window >= 3 "window must be at least 3"
    @assert degree >= 1 "degree must be at least 1"
    @assert support_radius >= 1 "support_radius must be at least 1"
    @assert min_support >= 1 "min_support must be at least 1"

    n = length(y)
    k = div(window, 2)

    y0 = collect(float.(y))

    # Log-space data; zeros are ignored
    z = fill(NaN, n)
    for i in 1:n
        if y0[i] > 0
            z[i] = log(y0[i])
        end
    end

    suspect = falses(n)
    zhat = fill(NaN, n)
    scale = fill(NaN, n)

    # First pass: detect suspicious points from a local robust curve fit
    for i in 1:n
        lo = max(1, i - k)
        hi = min(n, i + k)

        idx = Int[]
        for j in lo:hi
            if j != i && isfinite(z[j])
                push!(idx, j)
            end
        end

        length(idx) < degree + 1 && continue

        p = min(degree, length(idx) - 1)

        x = Float64.(idx .- i)
        xmax = maximum(abs.(x))
        xmax == 0 && continue
        u = x ./ xmax

        X = hcat([u .^ d for d in 0:p]...)
        w = exp.(-2.0 .* u.^2)
        sw = sqrt.(w)

        Aw = X .* sw
        bw = z[idx] .* sw

        G = Aw' * Aw + ridge * I(size(Aw, 2))
        g = Aw' * bw

        β = G \ g

        # Local fitted value at the center, in log space
        zhat[i] = β[1]

        resid = bw - Aw * β
        medr = median(resid)
        madr = median(abs.(resid .- medr))
        scale[i] = 1.4826 * madr

        if isfinite(z[i]) && isfinite(zhat[i]) && scale[i] > 0 &&
           abs(z[i] - zhat[i]) > n_sigma * scale[i]
            suspect[i] = true
        end
    end

    # Helper: decide whether a suspect point is supported by a real local peak
    function peak_supported(i::Int)
        lo = max(1, i - support_radius)
        hi = min(n, i + support_radius)

        vals = Float64[]
        for j in lo:hi
            if j != i && isfinite(z[j])
                push!(vals, z[j])
            end
        end

        length(vals) < 2 && return false

        base = median(vals)

        count_supported = 0
        for v in vals
            if v >= base - support_tol
                count_supported += 1
            end
        end

        return count_supported >= min_support
    end

    # Second pass: only correct isolated singleton suspects that are not peak-supported
    corrected = falses(n)

    for i in 2:(n - 1)
        if suspect[i] && !suspect[i - 1] && !suspect[i + 1] && !peak_supported(i)
            if isfinite(zhat[i])
                z[i] = zhat[i]   # replace with local fitted value in log space
                corrected[i] = true
            end
        end
    end

    # Back to original space
    for i in 1:n
        if isfinite(z[i])
            y[i] = exp(z[i])
        end
    end

    return suspect, corrected
end


function fix_monotone_center_log!(y::AbstractVector{<:Float32},z::Vector{<:Float32};floor::Float32 = 1f-45)

    n = length(y)
    fill!(z,NaN32)

    # Work in log-space; ignore zeros / nonpositive values
    for i in 1:n
        if y[i] > 0
            z[i] = log(max(y[i], floor))
        end
    end

    #corrected = falses(n)

    for i in 3:(n - 2)
        a, b, c, d, e = z[i-2], z[i-1], z[i], z[i+1], z[i+2]

        # Skip if any value is missing
        if any(!isfinite, (a, b, c, d, e))
            continue
        end

        # Four surrounding points are monotone increasing or decreasing
        monotone_increasing = (a < b) && (b < d) && (d < e)
        monotone_decreasing = (a > b) && (b > d) && (d > e)

        if monotone_increasing || monotone_decreasing
            # Center point should lie between its immediate neighbors
            if (c <= min(b, d)) || (c >= max(b, d))
                z[i] = mean((a, b, d, e))
                #corrected[i] = true
            end
        end
    end

    # Back to original scale
    for i in 1:n
        if !isnan(z[i])
            y[i] = exp(z[i])
        end
    end

    return nothing
end