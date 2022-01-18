struct Rayleigh end

scattering(m, x, ::Rayleigh) = 8 * x^4 / 3 * abs((m^2 - 1) / (m^2 + 2))^2
absorption(m, x, ::Rayleigh) = 4 * x * imag((m^2 - 1) / (m^2 + 2))
extinction(m, x, s::Rayleigh) = scattering(m, x, s) + absorption(m, x, s)
backscattering(m, x, s::Rayleigh) = 3 * scattering(m, x, s) / 2


# s = Scatterer(1.33 + 0.01im, 25., 870.)
# scattering(s, Rayleigh())

@doc raw"""
    rayleigh_scattering(m, x)

Computes Mie efficencies of a spherical particle in the Rayleigh regime (``x = \pi \, d_p / \lambda \ll 1``) given refractive index `m`, `wavelength`, and `diameter`. Uses Rayleigh-regime approximations:
```math
Q_{sca} = \frac{8x^4}{3} \left|{\frac{m^2-1}{m^2+2}}\right|^2
```
```math
Q_{abs} = 4x \: \text{Im}\left\{\frac{m^2-1}{m^2+2}\right\}
```
```math
Q_{ext} = Q_{sca}+ Q_{abs}
```
```math
Q_{back} = \frac{3Q_{sca}}{2}
```
```math
Q_{ratio} = 1.5
```
```math
Q_{pr} = Q_{ext}
```

# Arguments
- `m`: The complex refractive index, with the convention ``m = n+ik``.
- `wavelength`: The wavelength of incident light, in nanometers.
- `diameter`: The diameter of the particle, in nanometers.
- `n_medium`: The refractive index of the surrounding medium. This must be positive, nonzero, and real. Any imaginary part will be discarded.

# Returns
- `qext`, `qsca`, `qabs`, `qback`: The Mie efficencies described above.

# Examples
```julia-repl
julia> rayleigh_scattering(mx(1.33 + 0.01im, 25., 870.)...)

(0.0041753430994240295, 0.00011805645915412197, 0.004057286640269908, 0.00017708468873118297)
```
# Notes
-  different formula on [wiki](https://en.m.wikipedia.org/wiki/Rayleigh_scattering)
"""
function rayleigh_scattering(m::ComplexF64, x::Float64)

    if x ≈ 0
        return 0, 0, 0, 0
    end

    LL = (m^2 - 1) / (m^2 + 2) # Lorentz-Lorenz term

    qsca = 8 * x^4 / 3 * abs(LL)^2  # B&H eq 5.8
    qabs = 4 * x * imag(LL) # B&H eq. 5.11
    qext = qsca + qabs
    qback = 3 * qsca / 2 # B&H eq. 5.9

    # qratio = 1.5
    # qpr = qext

    return qext, qsca, qabs, qback
end


# rayleigh_scattering(mx(1.33 + 0.01im, 0., 870.)...)
# rayleigh_scattering(mx(1.33 + 0.01im, 25., 870.)...)
# #
# lambda = LinRange(600, 1000, 1001)
# qext =  [rayleigh_scattering(mx(1.33 + 0.01im, 25., l)...)[1] for l in lambda]
# plot(lambda, qext)
#


