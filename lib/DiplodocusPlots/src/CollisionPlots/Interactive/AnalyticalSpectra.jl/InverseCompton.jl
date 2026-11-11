#=

## TO DO: put paper reference from Marc here

=# 

function omega_min(omega_0, gamma_0)
    p_0 = sqrt(gamma_0^2 - 1)
    return ((gamma_0 - p_0) * omega_0) / (gamma_0 + p_0 + 2 * omega_0)
end

function omega_max(omega_0, gamma_0)
    p_0 = sqrt(gamma_0^2 - 1)
    thigh = gamma_0 + omega_0 - 1
    tlow = omega_c(omega_0, gamma_0)
    val = (omega_0 > (1 + p_0 - gamma_0) / 2) ? thigh  : tlow
    return val
end

function omega_c(omega_0, gamma_0)
    p_0 = sqrt(gamma_0^2 - 1)
    return ((gamma_0 + p_0) * omega_0) / (gamma_0 - p_0 + 2 * omega_0)
end

function omega_bar(omega_0, gamma_0, omega, gamma) 
    p_0 = sqrt(gamma_0^2 - 1)
    p = sqrt(gamma^2 - 1)
    return sqrt(max((omega * omega_0 * (gamma + p)) / (gamma_0 + p_0),0e0))
end

function omega0_bar(omega_0, gamma_0, omega, gamma) 
    p_0 = sqrt(gamma_0^2 - 1)
    p = sqrt(gamma^2 - 1)
    return sqrt(max((omega * omega_0 * (gamma_0 + p_0)) / (gamma + p),0e0))
end

function Lambda_function(omega_0, omega, p_0, p) 
    return (p_0 - p + omega_0 + omega) / 2
end

function Sfun(x)

    threshold = 1e-4
    if x >= threshold
        val = asinh(sqrt(x)) / sqrt(x) 
    else
        val = 1 - x/6 + 3*x^2/40 - 5*x^3/112
    end

    return val
end

function Ffun(x)

    threshold = 1e-4
    if x >= threshold
        val = Sfun(x) - sqrt(1 + x)
    else
        val = -2*x/3 + x^2/5 - 3*x^3/28
    end
    return val
end

function G_fun(w0s, ws, kappa, w0, w, g0)
    p0 = sqrt(g0^2 - 1)
    lambdap = p0^2 + 2*g0*w0 + w0^2 
    lambdam = p0^2 - 2*g0*w + w^2 
    
    term1 = 2 + ((ws - w0s)^2 * (1 + w * w0)) / (w^2 * w0^2)
    term2 = 2 * (Sfun(kappa^2 * lambdap / ws^2) / ws - Sfun(kappa^2 * lambdam / w0s^2) / w0s)
    term3 = (1 + w * w0) * (Ffun(kappa^2 * lambdap / ws^2) / (ws * lambdap) - Ffun(kappa^2 * lambdam / w0s^2) / (w0s * lambdam))
    
    return kappa * (term1 + term2 + term3)
end

function IC_kernel(E_in_lab, bG_el, E_out_lab)

    mec2 = 8.198445204e-07
    w0 = E_in_lab #* 3e8#/ mec2
    w = E_out_lab #* 3e8 #/ mec2
    g0 = sqrt(1 + bG_el^2)
    p0 = bG_el
    p = sqrt(max(p0^2 + 2 * g0 * (w0 - w) + (w0 - w)^2,0e0)) 
    g = sqrt(1 + p^2)
    
    wc = omega_c(w0, g0)
    w0b = omega0_bar(w0, g0, w, g)
    wb = omega_bar(w0, g0, w, g)
    wmin = omega_min(w0, g0)
    wmax = omega_max(w0, g0)
    
    prefac = 3 / (8 * g0 * p0 * w0^2) 
    
    if wc <= w0
        #println("omega_c <= omega_0")
        kappa1 = Lambda_function(w0, w, p0, p) 
        zone1 = prefac * G_fun(w0b, wb, kappa1, w0, w, g0) 
        zone2 = prefac * G_fun(w, w0, p0, w0, w, g0) 
        zone3 = prefac * G_fun(w0, w, p, w0, w, g0) 
        val =  ((w >= wmin && w < wc) * zone1 + (w >= wc && w < w0) * zone2 + (w >= w0 && w < wmax) * zone3)
        return max(val,0e0)
    else
        #println("omega_c > omega_0")
        kappa1 = Lambda_function(w0, w, p0, p) 
        kappa2 = Lambda_function(w, w0, p, p0) 
        zone1 = prefac * G_fun(w0b, wb, kappa1, w0, w, g0) 
        zone2 = prefac * G_fun(wb, w0b, kappa2, w0, w, g0) 
        zone3 = prefac * G_fun(w0, w, p, w0, w, g0) 
        val =  ((w >= wmin && w < w0) * zone1 + (w >= w0 && w < wc) * zone2 + (w >= wc && w < wmax) * zone3)
        return max(val,0e0)
    end
end