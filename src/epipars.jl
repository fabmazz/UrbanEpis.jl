

rate_to_prob(b::Real,t::Real) = 1 - exp(-b*t)

struct CitySimPars1
    R0::AbstractFloat
    ts::AbstractFloat
    a::Real
    b::Real
end

CitySimPars = CitySimPars1