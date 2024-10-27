@with_kw struct AnglesGamma
    θ::Float64
    ϕ::Float64
    γ::Float64
end
#
abstract type DaughterTransformation end

struct HelicityTransformation <: DaughterTransformation
    variables::AnglesGamma
end

function act(T::HelicityTransformation, p)
    @unpack θ, ϕ, γ = T.variables
    p |> Rz(-ϕ) |> Ry(-θ) |> Bz(-γ)
end

function HelicityTransformation(; p, isfirst = true)
    ϕ = azimuthal_angle(p)
    θ = polar_angle(p)
    if pt(p) < 1e-10
        ϕ, θ = 0.0, 0.0
    end
    γ = boost_gamma(p)
    isfirst && return HelicityTransformation(AnglesGamma(; θ, ϕ, γ))
    return HelicityTransformation(AnglesGamma(; θ = π - θ + π, ϕ = mod(ϕ + 2π, 2π) - π, γ))
end
