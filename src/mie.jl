export mie_mx,mie_scattering, scattering_function



findnmax(x) = round(Int, 2 + x + 4x^(1 / 3))

"""
    mie_mx(m, radius, wavelength, n_medium=1.00)

    Conversion of parameters for calculation of light scattering by a small sphere based on Mie theory

# Arguments
- `m`: n_material; refractive index of the material of the sphere
- `radius`: radius of the sphere
- `wavelength`: target vacuum wavelength (range)
- `n_medium`: n_environement; refractive index of the environment
      around the sphere (default: air)

# Return values
- `m`: relative refractive index (n_material / n_environment)
- `x`: size parameter (wavenumber * radius)
"""
function mie_mx(m, radius, wavelength, n_medium=1.00)
    m /= n_medium
    wavelength /= n_medium  # effective wavelength in the environment

    k = 2π / wavelength  # free space wavenumber in the environment
    x = k .* radius
    return m, x
end

@doc raw"""
    mie_ab(m, x, nmax)

Computes external field coefficients ``a_n`` and ``b_n`` based on inputs of `m` and ``x = \pi \, d_p / \lambda``.

# Arguments
- m: The complex refractive index with the convention ``m = n+ik``.
- x: The size parameter ``x = \pi \, d_p / \lambda``.

# Return values
- an, bn: Arrays of size ``n_max = 2 + x + 4x^{1/3}``

"""
function mie_ab(m, x, nmax)

    mx = m * x
    nmx = max(nmax, round(Int, abs(mx))) + 16

    n = 1:nmax

    nu = n .+ 0.5
    sx = sqrt(0.5 * π * x)

    px = sx * besselj.(nu, x)
    p1x = append!([sin(x)], px[1:end - 1])

    chx = -sx * bessely.(nu, x)
    ch1x = append!([cos(x)], chx[1:end - 1])

    gsx = px - 1im * chx
    gs1x = p1x - 1im * ch1x

    # B&H Equation 4.89
    Dn = zeros(ComplexF64, nmx)
    for i in nmx - 1:-1:1
        Dn[i] = (i / mx) - (1 / (Dn[i + 1] + i / mx))
    end

    D = Dn[2:nmax + 1] # Dn(mx), drop terms beyond nMax

    da = D / m + n / x
    db = m * D + n / x

    an = @. (da * px - p1x) / (da * gsx - gs1x)
    bn = @. (db * px - p1x) / (db * gsx - gs1x)

    return an, bn
end

# using SpecialFunctions
# spherical bessel functions
sbesselj(nu, x) = √(π / 2x) * besselj(nu + .5, x)
sbessely(nu, x) = √(π / 2x) * bessely(nu + .5, x)
shankelh1(nu,x) = sbesselj(nu, x) + 1im * sbessely(nu, x)

"""
    mie_abcd(m, x, n::Int)

Calculate n-th order Mie coefficients a,b,c,d

# Arguments
- `m`: relative refractive index (n_material / n_environment)
- `x`: size parameter (wavenumber * radius)
- `n::Int`
"""
function mie_ab_n(m, x, n::Int)

    # Bessel/Hankel functions
    jx = sbesselj(n, x)
    jmx = sbesselj(n, m * x)
    hx = shankelh1(n, x)

    # derivative
    d1xjx = x * sbesselj(n - 1, x) - n * sbesselj(n, x)
    d1mxjmx = m * x * sbesselj(n - 1, m * x) - n * sbesselj(n, m * x)
    d1xhx = x * shankelh1(n - 1, x) - n * shankelh1(n, x)

    an = (m^2 * jmx * d1xjx - jx * d1mxjmx) / (m^2 * jmx * d1xhx - hx * d1mxjmx)
    bn = (jmx * d1xjx - jx * d1mxjmx) / (jmx * d1xhx - hx * d1mxjmx)

    # cn = (jx * d1xhx - hx * d1xjx) / (jmx * d1xhx - hx * d1mxjmx)
    # dn = (m * jx * d1xhx - m * hx * d1xjx) / (m^2 * jmx * d1xhx - hx * d1mxjmx)

    # return an, bn, cn, dn
    return an, bn
end

#
# m = 1.33 + 0.01im
# x = 3
#
# an, bn = mie_ab(m, x, 3)
# an[1] ≈ mie_ab_n(m, x, 1)[1]
# an[2] ≈ mie_ab_n(m, x, 2)[1]
# an[3] ≈ mie_ab_n(m, x, 3)[1]
# bn[1] ≈ mie_ab_n(m, x, 1)[2]
# bn[2] ≈ mie_ab_n(m, x, 2)[2]
# bn[3] ≈ mie_ab_n(m, x, 3)[2]
#
#
# m = 2 + 1im
# x = 5
#
# an, bn = mie_ab(m, x, 3)
# an[1] ≈ mie_ab_n(m, x, 1)[1]
# an[2] ≈ mie_ab_n(m, x, 2)[1]
# an[3] ≈ mie_ab_n(m, x, 3)[1]
# bn[1] ≈ mie_ab_n(m, x, 1)[2]
# bn[2] ≈ mie_ab_n(m, x, 2)[2]
# bn[3] ≈ mie_ab_n(m, x, 3)[2]


@doc raw"""
    mie_pi_tau(mu, nmax)

Calculates ``\pi_n`` and ``\tau_n``.

This function uses recurrence relations to calculate ``\pi_n`` and ``\tau_n``, beginning with ``\pi_0 = 1``, ``\pi_1 = 3 \mu`` (where ``mu`` is the cosine of the scattering angle), ``\tau_0 = \mu``, and ``\tau_1 = 3 cos(2 cos^{-1}(\mu))``:
```math
\pi_n = \frac{2n-1}{n-1} \mu \pi_{n-1} - \frac{n}{n-1} \pi_{n-2}
```
```math
\tau_n = n \mu\ pi_n - (n+1) \ pi_{n-1}
```

# Arguments
- mu : The cosine of the scattering angle.
- nmax: The number of elements to compute. Typically, ``n_{max} = floor(2+x+4x^{1/3})``, but can be given any integer.

# Return values
- p, t: The ``\pi_n`` and ``\tau_n`` arrays, of length `nmax`.
"""
function mie_pi_tau(mu, nmax)

    p = zeros(nmax)
    t = zeros(nmax)

    p[1] = 1
    p[2] = 3mu

    t[1] = mu
    t[2] = 3 * cos(2 * acos(mu))

    for n in 3:nmax
        p[n] = (2n - 1) / (n - 1) * mu * p[n - 1] - n / (n - 1) * p[n - 2]
        t[n] = n * mu * p[n] - (n + 1) * p[n - 1]
    end

    return p, t
end

"""
    mie_pi_tau(mu, n::Int)

See Sec. 4.3.1 of B&H
ref of pi,tau: [D. Deirmendjian, “Electromagnetic Scattering on Spherical Polydispersions"]
"""
function mie_pi_tau_n(mu, n::Int)
# See Sec. 4.3.1 of B&H
# ref of pi,tau: [D. Deirmendjian, “Electromagnetic Scattering on Spherical Polydispersions"]

    if mu < -1 || mu > 1
        error("mu is cos(theta), -1 <= mu <= 1")
    end

    if n == 0
        pin = 0.  # pi(0, mu) = 0
        taun = 0.  # not used but tau(0, mu) = 0
    elseif n == 1
        pin = 1.  # pi(1, mu) = 1
        taun = mu  # note that pi(0, mu) = 0
    elseif n > 1
        # recurrence, eq.(4.47) of B&H
        pin = (2n - 1) / (n - 1) * mu * mie_pi_tau_n(mu, n - 1)[1] - (n) / (n - 1) * mie_pi_tau_n(mu, n - 2)[1]
        taun = n * mu * pin - (n + 1) * mie_pi_tau_n(mu, n - 1)[1]
    else
        error("n should be larger than 0")
    end
    return pin, taun
end

# # error
# # pin, taun = mie_pi_tau(mu, 1)
#
# mu = cos(0)
# pn, tn = mie_pi_tau(mu, 3)
# pn[1] ≈ mie_pi_tau_n(mu, 1)[1]
# pn[2] ≈ mie_pi_tau_n(mu, 2)[1]
# pn[3] ≈ mie_pi_tau_n(mu, 3)[1]
# tn[1] ≈ mie_pi_tau_n(mu, 1)[2]
# tn[2] ≈ mie_pi_tau_n(mu, 2)[2]
# tn[3] ≈ mie_pi_tau_n(mu, 3)[2]
#
# mu = cos(π / 3)
# pn, tn = mie_pi_tau(mu, 3)
# pn[1] ≈ mie_pi_tau_n(mu, 1)[1]
# pn[2] ≈ mie_pi_tau_n(mu, 2)[1]
# pn[3] ≈ mie_pi_tau_n(mu, 3)[1]
# tn[1] ≈ mie_pi_tau_n(mu, 1)[2]
# tn[2] ≈ mie_pi_tau_n(mu, 2)[2]
# tn[3] ≈ mie_pi_tau_n(mu, 3)[2]
#
# mu = cos(π / 2)
# pn, tn = mie_pi_tau(mu, 3)
# pn[1] ≈ mie_pi_tau_n(mu, 1)[1]
# pn[2] ≈ mie_pi_tau_n(mu, 2)[1]
# pn[3] ≈ mie_pi_tau_n(mu, 3)[1]
# tn[1] ≈ mie_pi_tau_n(mu, 1)[2]
# tn[2] ≈ mie_pi_tau_n(mu, 2)[2]
# tn[3] ≈ mie_pi_tau_n(mu, 3)[2]

@doc raw"""
    mie_scattering(m, x, [nmax])

Computes Mie efficencies *Q* and asymmetry parameter *g* of a single, homogeneous particle. Uses `mie_ab` to calculate ``a_n`` and ``b_n``, and then calculates *Q* via:
```math
Q_{ext}=\frac{2}{x^2}\sum_{n=1}^{n_{max}}(2n+1)\:\text{Re}\left\{a_n+b_n\right\}}
```
```math
Q_{sca}=\frac{2}{x^2}\sum_{n=1}^{n_{max}}(2n+1)(|a_n|^2+|b_n|^2)}
```
```math
Q_{abs}=Q_{ext}-Q_{sca}}
```
```math
Q_{back}=\frac{1}{x^2} \left| \sum_{n=1}^{n_{max}}(2n+1)(-1)^n(a_n-b_n) \right| ^2}
```
```math
Q_{ratio}=\frac{Q_{back}}{Q_{sca}}}
```
```math
g=\frac{4}{Q_{sca}x^2}\left[\sum\limits_{n=1}^{n_{max}}\frac{n(n+2)}{n+1}\text{Re}\left\{a_n a_{n+1}^*+b_n b_{n+1}^*\right\}+\sum\limits_{n=1}^{n_{max}}\frac{2n+1}{n(n+1)}\text{Re}\left\{a_n b_n^*\right\}\right]}
```
```math
Q_{pr}=Q_{ext}-gQ_{sca}}
```
where asterisks denote the complex conjugates.

# Arguments
- m: The complex refractive index, with the convention *m = n+ik*.
- wavelength: The wavelength of incident light, in nanometers.
- diameter: The diameter of the particle, in nanometers.
- n_medium: The refractive index of the surrounding medium. This must be positive, nonzero, and real. Any imaginary part will be discarded.

# Return values
- qext, qsca, qabs, g, qpr, qback, qratio: The Mie efficencies described above.

# Examples

For example, compute the Mie efficencies of a particle ``300 nm`` in diameter with ``m = 1.77+0.63i``, illuminated by ``\lambda = 375 nm``:
```julia-repl
julia> mie_scattering(1.77+0.63im, 375, 300)

(2.858497199156411, 1.3149276685170943, 1.543569530639317, 0.2014551048135256)
```
"""
function mie_scattering(m, x, nmax=findnmax(x))

    if x ≈ 0
        return 0., 0., 0., 0.
    end

    an, bn = mie_ab(m, x, nmax)

    # Note:
    # Q_i = \sigma_i / (\pi a^2)
    # where
    # - Q_i is the efficiency coefficients
    # - \sigma_{i} is the cross section of the respective process

    qsca = 2 / x^2 * sum((2n + 1) * (abs(an[n])^2 + abs(bn[n])^2) for n = 1:nmax)
    qext = 2 / x^2 * sum((2n + 1) * real(an[n] + bn[n]) for n = 1:nmax)
    qabs = qext - qsca
    qback = 1 / x^2 * abs(sum((2n + 1) * (-1)^n * (an[n] - bn[n]) for n = 1:nmax))^2

    # g = 4 / (qsca * x^2) * (sum(n * (n + 2) / (n + 1) * real(an[n]*an[n+1]'+bn[n]*bn[n+1]') for n = 1:nmax-1) + sum((2n + 1) / (n * (n + 1)) * real(an[n]*bn[n]')  for n = 1:nmax-1))
    # qpr = qext - qsca * g
    # qratio = qback / qsca

    return qext, qsca, qabs, qback
end

"""
    mie_scattering_n(m, x, nmax::Int=-1)

Calculate scattering efficiency of a sphere

# Arguments
- `m`: relative refractive index (n_material / n_environment)
- `x`: size parameter (wavenumber * radius)
- `namx`: maximun order of resonances being considered
      (nmax = 1: dipole, 2: dipole+quadrupole, ...,
       -1: sufficiently large value determined by x)

# Return values
- `Qsca`: Scattering efficiency (total up to n-th order resonances)
"""
function mie_scattering_n(m, x, nmax=findnmax(x))

    if x ≈ 0
        return 0.
    end

    qsca, qext = 0., 0.
    for n = 1:nmax
        an, bn = mie_ab_n(m, x, n)
        qsca += (2n + 1) * (abs(an)^2 + abs(bn)^2)
        qext += (2n + 1) * real(an + bn)
    end
    qsca *= 2 / x^2
    qext *= 2 / x^2

    qabs = qext - qsca

    return qext, qsca, qabs, nothing
end

# m = 1.33 + 0.01im
# x = 3
#
# qext1, qsca1, _, _ = mie_scattering(m, x)
# qext2, qsca2, _, _ = mie_scattering_n(m, x)
#
# m = 2 + 1im
# x = 5
#
#
# qext1, qsca1, _, _ = mie_scattering(m, x)
# qext2, qsca2, _, _ = mie_scattering_n(m, x)


@doc raw"""
    mie_S1_S2(m,x,mu)

Calculates ``S_1`` and ``S_2`` at ``\mu = cos(\theta)``, where ``\theta`` is the scattering angle.

Uses `mie_ab` to calculate `a_n` and `b_n`, and `mie_pi_tau` to calculate ``\pi_n`` and ``\tau_n``. ``S_1`` and ``S_2`` are calculated by:

```math
S_1 = \sum\limits_{n=1}^{n_{max}} \frac{2n+1}{n(n+1)} (a_n \pi_n + b_n \tau_n)
```
```math
S_2 = \sum\limits_{n=1}^{n_{max}} \frac{2n+1}{n(n+1)} (a_n \tau_n + b_n \pi_n)
```

# Parameters
- m: The complex refractive index with the convention *m = n+ik*.
- x: The size parameter ``x = \pi \, d_p / \lambda``.
- mu: The cosine of the scattering angle.

# Return values
-S1, S2: The ``S_1`` and ``S_2`` values.

"""
function mie_S1_S2(m, x, mu, nmax=findnmax(x))

    an, bn = mie_ab(m, x, nmax)
    pin, taun = mie_pi_tau(mu, nmax)

    S1 = sum((2n + 1) / (n * (n + 1)) * (an[n] * pin[n] + bn[n] * taun[n]) for n = 1:nmax)
    S2 = sum((2n + 1) / (n * (n + 1)) * (an[n] * taun[n] + bn[n] * pin[n]) for n = 1:nmax)

    return S1, S2
end


"""
    mie_S1_S2(m, x, mu, n::Int)

calculate n-th order S1,S2 component of scattering matrix.
"""
function mie_S1_S2_n(m, x, mu, n)

    a, b = mie_ab_n(m, x, n)
    p, t = mie_pi_tau_n(mu, n)

    S1 = (2n + 1) / (n * (n + 1)) * (a * p + b * t)
    S2 = (2n + 1) / (n * (n + 1)) * (a * t + b * p)

    return S1, S2
end

#
# m = 1.33 + 0.01im
# x = 3
#
# mu = cos(0)
# s1s2 = mie_S1_S2(m, x, mu, 2)
# s1s2 .≈ mie_S1_S2_n(m, x, mu, 1) .+ mie_S1_S2_n(m, x, mu, 2)
#
#
# mu = cos(π / 3)
# s1s2 = mie_S1_S2(m, x, mu, 2)
# s1s2 .≈ mie_S1_S2_n(m, x, mu, 1) .+ mie_S1_S2_n(m, x, mu, 2)
#
# mu = cos(π / 2)
# s1s2 = mie_S1_S2(m, x, mu, 2)
# s1s2 .≈ mie_S1_S2_n(m, x, mu, 1) .+ mie_S1_S2_n(m, x, mu, 2)
#
#
# m = 2 + 1im
# x = 5
#
# mu = cos(0)
# s1s2 = mie_S1_S2(m, x, mu, 2)
# s1s2 .≈ mie_S1_S2_n(m, x, mu, 1) .+ mie_S1_S2_n(m, x, mu, 2)
#
# mu = cos(π / 3)
# s1s2 = mie_S1_S2(m, x, mu, 2)
# s1s2 .≈ mie_S1_S2_n(m, x, mu, 1) .+ mie_S1_S2_n(m, x, mu, 2)
#
# mu = cos(π / 2)
# s1s2 = mie_S1_S2(m, x, mu, 2)
# s1s2 .≈ mie_S1_S2_n(m, x, mu, 1) .+ mie_S1_S2_n(m, x, mu, 2)


@doc raw"""
    scattering_function(m, wavelength, diameter, n_medium=1.0, minAngle=0, maxAngle=180, angularResolution=0.5)

Creates arrays for plotting the angular scattering intensity functions in theta-space with parallel, perpendicular, and unpolarized light. Also includes an array of the angles for each step. This angle can be in either degrees, radians, or gradians for some reason. The angles can either be geometrical angle or the qR vector (see `Sorensen, M. Q-space analysis of scattering by particles: a review. J. Quant. Spectrosc. Radiat. Transfer 2013, 131, 3-12 <http://www.sciencedirect.com/science/article/pii/S0022407313000083>`_). Uses `mie_S1_S2` to compute ``S_1`` and ``S_2``, then computes parallel, perpendicular, and unpolarized intensities by

```math
SL(\theta)=|S_1|^2}
```
```math
SR(\theta)=|S_2|^2}
```
```math
SU(\theta)=\frac{1}{2}(SR + SL)}
```

# Arguments

- theta: An array of the angles used in calculations.
- m: The complex refractive index with the convention *m = n+ik*.
- wavelength: The wavelength of incident light, in nanometers.
- diameter: The diameter of the particle, in nanometers.
- n_medium: The refractive index of the surrounding medium. This must be positive, nonzero, and real. Any imaginary part will be discarded.

# Return values

- SL: An array of the scattered intensity of left-polarized (perpendicular) light. Same size as the **theta** array.
- SR: An array of the scattered intensity of right-polarized (parallel) light. Same size as the **theta** array.
- SU: An array of the scattered intensity of unpolarized light, which is the average of SL and SR. Same size as the **theta** array.
"""
function scattering_function(theta, m, x)

    n = length(theta)

    SL = zeros(n)
    SR = zeros(n)
    SU = zeros(n)

    if x ≈ 0
        return SL, SR, SU
    end

    for i in 1:n
        mu = cos(theta[i])
        s1, s2 = mie_S1_S2(m, x, mu)
        SL[i] = abs(s1^2)
        SR[i] = abs(s2^2)
        SU[i] = (SR[i] + SL[i]) / 2
    end

    return SL, SR, SU
end

function scattering_function_n(theta, m, x)

    n = length(theta)

    SL = zeros(n)
    SR = zeros(n)
    SU = zeros(n)

    if x ≈ 0
        return SL, SR, SU
    end

    for i in 1:n
        mu = cos(theta[i])
        nmax = findnmax(x)
        s1, s2 = 0., 0.
        for n = 1:nmax
            s1n, s2n = mie_S1_S2_n(m, x, mu, n)
            s1 += s1n
            s2 += s2n
        end

        SL[i] = abs(s1^2)
        SR[i] = abs(s2^2)
        SU[i] = (SR[i] + SL[i]) / 2
    end

    return SL, SR, SU
end

# m = 1.33 + 0.01im
# x = 3
#
# theta = 0
# scattering_function(theta, m, x)
# scattering_function_n(theta, m, x)
#
#
# theta = π/3
# scattering_function(theta, m, x)
# scattering_function_n(theta, m, x)
#
# theta = π/2
# scattering_function(theta, m, x)
# scattering_function_n(theta, m, x)
#
#
# m = 2 + 1im
# x = 5
#
# theta = 0
# scattering_function(theta, m, x)
# scattering_function_n(theta, m, x)
#
# theta = π/3
# scattering_function(theta, m, x)
# scattering_function_n(theta, m, x)
#
# theta = π/2
# scattering_function(theta, m, x)
# scattering_function_n(theta, m, x)

