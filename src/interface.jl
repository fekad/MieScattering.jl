struct Scatterer
    m::ComplexF64
    x::Float64
end

function Scatterer(m::ComplexF64, radius::Float64, wavelength::Float64, n_medium::Float64 = 1.00)
    m /= n_medium
    wavelength /= n_medium  # effective wavelength in the environment

    k = 2π / wavelength  # free space wavenumber in the environment
    x = k .* radius
    return Scatterer(m, x)
end

scattering(s::Scatterer, m) = scattering(s.m, s.x, m)
absorption(s::Scatterer, m) = absorption(s.m, s.x, m)
extinction(s::Scatterer, m) = extinction(s.m, s.x, m)
backscattering(s::Scatterer, m) = backscattering(s.m, s.x, m)



@doc raw"""
mx(m, radius, wavelength, n_medium=1.00)

Conversion of parameters for calculation of light scattering by a small sphere based on Mie theory

# Arguments
- `m`: ``n_{\mathrm{material}}``; refractive index of the material of the sphere
- `radius`: radius of the sphere
- `wavelength`: target vacuum wavelength (range)
- `n_medium`: n_environement; refractive index of the environment
  around the sphere (default: air)

# Return values
- `m`: relative refractive index (``n_\mathrm{material} / n_\mathrm{environment}``)
- `x`: size parameter (wavenumber * radius)
"""
function mx(m, radius, wavelength, n_medium = 1.00)
    m /= n_medium
    wavelength /= n_medium  # effective wavelength in the environment

    k = 2π / wavelength  # free space wavenumber in the environment
    x = k .* radius
    return m, x
end



