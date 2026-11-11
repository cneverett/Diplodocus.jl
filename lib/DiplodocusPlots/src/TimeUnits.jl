""" 
    CodeToSIUnitsTime(t_code)

Returns the time in SI units seconds from the time in code units of seconds * σTc 
"""
function CodeToSIUnitsTime(t_code::T) where T <: AbstractFloat
    # Convert time code to SI units (seconds)
    σT = getfield(DC,Symbol("σT"))
    c = getfield(DC,Symbol("c"))  # σT is the time unit in seconds

    t_SI::T = t_code / (σT*c)  # Convert time code to seconds
    return t_SI
end
function CodeToSIUnitsTime()
    return raw"$[\text{s}]$"
end

""" 
    CodeToCodeUnitsTime(t_code)

Returns the time in code units from the time in code units
"""
function CodeToCodeUnitsTime(t_code::T) where T <: AbstractFloat
    return t_code
end
function CodeToCodeUnitsTime()
    return raw"$[\text{code units}]$"
end

""" 
    SIToCodeUnitsTime(t_SI)

Returns the time in code units seconds from the time in units of seconds
"""
function SIToCodeUnitsTime(t_SI::T) where T <: AbstractFloat
    # Convert time code to SI units (seconds)
    σT = getfield(DC,Symbol("σT"))
    c = getfield(DC,Symbol("c"))  # σT is the time unit in seconds

    t_code::T = t_SI * (σT*c)  # Convert time code to seconds
    return t_code
end

""" 
    SyncToCodeUnitsTime(t_sync,B)

Returns the time in code units from the time in synchrotron units, requires the magnetic field `B` in Tesla to be defined
"""
function SyncToCodeUnitsTime(t_sync::T;B::T=1e-4) where T <: AbstractFloat

    μ0 = getfield(DC,Symbol("μ0"))
    c = getfield(DC,Symbol("c"))
    m = getfield(DC,Symbol("mEle"))

    t_code::T = t_sync * μ0*m*c^2/B^2

    return t_code
end

""" 
    CodeToSyncUnitsTime(t_code,B)

Returns the time in synchrotron units seconds from the time in code units, requires the magnetic field `B` in Tesla to be defined
"""
function CodeToSyncUnitsTime(t_code::T;B::T=1e-4) where T <: AbstractFloat
    # Convert time code to sync units
    μ0 = getfield(DC,Symbol("μ0"))
    c = getfield(DC,Symbol("c"))
    m = getfield(DC,Symbol("mEle"))

    t_sync::T = t_code / (μ0*m*c^2/B^2)
    return t_sync
end
function CodeToSyncUnitsTime()
    return raw"$[t_\text{sync}]$"
end
