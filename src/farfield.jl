export mie_Isff, mie_Isff_EM_n

"""
    mie_Isff_EM_n(nmat, radius, lbd0, nenv, theta, phi, n)

calculate far-feild intensity of scattered field
from n-th order Electric and Magnetic multipoles separately
"""
function mie_Isff_EM_n(nmat, radius, lbd0, nenv, theta, phi, n)
    m = nmat ./ nenv
    lbd = lbd0 ./ nenv  # effective wavelength in the environment
    k = 2π ./ lbd  # free space wavenumber in the environment
    x = k .* radius
    mu = cos.(theta)
    an, bn, _, _ = mie_abcd(m, x, n)
    pin, taun = mie_pt(mu, n)
    S1EDn = (2n + 1) / (n * (n + 1)) .* (an .* pin)
    S2EDn = (2n + 1) / (n * (n + 1)) .* (an .* taun)
    S1MDn = (2n + 1) / (n * (n + 1)) .* (bn .* taun)
    S2MDn = (2n + 1) / (n * (n + 1)) .* (bn .* pin)
    # settings
    E0 = 1  # amplitude of incident electric field
    r = 1  # distance from the sphere to an observation point (1 m)
    # main
    Est_EDn = E0 * exp.(1im .* k * r) ./ (-1im .* k * r) .* cos.(phi) .* S2EDn
    Esp_EDn = -1 * E0 * exp.(1im .* k * r) ./ (-1im .* k * r) .* sin.(phi) .* S1EDn
    Est_MDn = E0 * exp.(1im .* k * r) ./ (-1im .* k * r) .* cos.(phi) .* S2MDn
    Esp_MDn = -1 * E0 * exp.(1im .* k * r) ./ (-1im .* k * r) .* sin.(phi) .* S1MDn
    Is_EDn = abs.(Est_EDn).^2 .+ abs.(Esp_EDn).^2
    Is_MDn = abs.(Est_MDn).^2 .+ abs.(Esp_MDn).^2
    return Is_EDn, Is_MDn
end

"""
    mie_Isff(nmat, radius, lbd0, nenv, theta, phi, nmax::Int=-1)

calculate far-field intensity of scattered light
(as a total from n = 1 to nmax)
"""
function mie_Isff(nmat, radius, lbd0, nenv, theta, phi, nmax::Int=-1)

    if nmax == -1
        nmax = findnmax(mie_mx(nmat, radius, lbd0, nenv)[2])
    end

    Est, Esp = mie_Esff_n(nmat, radius, lbd0, nenv, theta, phi, 1)
    if nmax > 1
        for n = 2:nmax
            Est_n, Esp_n = mie_Esff_n(nmat, radius, lbd0, nenv, theta, phi, n)
            Est = Est .+ Est_n
            Esp = Esp .+ Esp_n
        end
    end
    Ist = abs.(Est).^2
    Isp = abs.(Esp).^2
    Is = Ist .+ Isp
    return Is, Ist, Isp
end

"""
    mie_Esff_n(nmat, radius, lbd0, nenv, theta, phi, n::Int)

calculate scattered electric far-fields E_s,theta and E_s,phi resulting from n-th order resonances
See Sec. 4.4.4 of B&H for derivation

# Arguments
- `nmat`: n_material; refractive index of the material of the sphere
- `radius`: radius of the sphere
- `lbd0`: target vacuum wavelength (range)
- `nenv`: n_environement; refractive index of the environment
      around the sphere
"""
function mie_Esff_n(nmat, radius, lbd0, nenv, theta, phi, n::Int)

    # preparation
    m = nmat ./ nenv
    lbd = lbd0 ./ nenv  # effective wavelength in the environment
    k = 2π ./ lbd  # free space wavenumber in the environment
    x = k .* radius
    mu = cos.(theta)
    S1n, S2n = mie_S1_S2(m, x, mu, n)

    # settings
    E0 = 1  # amplitude of incident electric field
    r = 1  # distance from the sphere to an observation point (1 m)

    # main
    Est_n = E0 * exp.(1im .* k * r) ./ (-1im .* k * r) .* cos.(phi) .* S2n
    Esp_n = -1 * E0 * exp.(1im .* k * r) ./ (-1im .* k * r) .* sin.(phi) .* S1n
    return Est_n, Esp_n
end