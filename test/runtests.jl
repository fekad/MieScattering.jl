using MieScattering:
    mie_scattering, mie_scattering_n, mie_ab, mie_ab_n, mie_pi_tau, mie_pi_tau_n, mie_S1_S2, mie_S1_S2_n, scattering_function, scattering_function_n
using Test

@testset "scattering values" begin
    @test isapprox(mie_scattering(1.0, 1.0)[2], 0.0, atol = 1e-6)
    @test isapprox(abs(mie_scattering(4.0, 1.0)[2]), 6.062172, atol = 1e-6)
end

@testset "self-consistency" begin

    m = 1.33 + 0.01im
    x = 3

    an, bn = mie_ab(m, x, 3)

    @test isapprox(an[1], mie_ab_n(m, x, 1)[1], atol=1e-6)
    @test isapprox(an[2], mie_ab_n(m, x, 2)[1], atol=1e-6)
    @test isapprox(an[3], mie_ab_n(m, x, 3)[1], atol=1e-6)
    @test isapprox(bn[1], mie_ab_n(m, x, 1)[2], atol=1e-6)
    @test isapprox(bn[2], mie_ab_n(m, x, 2)[2], atol=1e-6)
    @test isapprox(bn[3], mie_ab_n(m, x, 3)[2], atol=1e-6)


    m = 2 + 1im
    x = 5

    an, bn = mie_ab(m, x, 3)
    @test isapprox(an[1], mie_ab_n(m, x, 1)[1])
    @test isapprox(an[2], mie_ab_n(m, x, 2)[1])
    @test isapprox(an[3], mie_ab_n(m, x, 3)[1])
    @test isapprox(bn[1], mie_ab_n(m, x, 1)[2])
    @test isapprox(bn[2], mie_ab_n(m, x, 2)[2])
    @test isapprox(bn[3], mie_ab_n(m, x, 3)[2])


    # error
    # pin, taun = mie_pi_tau(mu, 1)

    mu = cos(0)
    pn, tn = mie_pi_tau(mu, 3)
    @test isapprox(pn[1], mie_pi_tau_n(mu, 1)[1])
    @test isapprox(pn[2], mie_pi_tau_n(mu, 2)[1])
    @test isapprox(pn[3], mie_pi_tau_n(mu, 3)[1])
    @test isapprox(tn[1], mie_pi_tau_n(mu, 1)[2])
    @test isapprox(tn[2], mie_pi_tau_n(mu, 2)[2])
    @test isapprox(tn[3], mie_pi_tau_n(mu, 3)[2])

    mu = cos(π / 3)
    pn, tn = mie_pi_tau(mu, 3)
    @test isapprox(pn[1], mie_pi_tau_n(mu, 1)[1])
    @test isapprox(pn[2], mie_pi_tau_n(mu, 2)[1])
    @test isapprox(pn[3], mie_pi_tau_n(mu, 3)[1])
    @test isapprox(tn[1], mie_pi_tau_n(mu, 1)[2])
    @test isapprox(tn[2], mie_pi_tau_n(mu, 2)[2])
    @test isapprox(tn[3], mie_pi_tau_n(mu, 3)[2])

    mu = cos(π / 2)
    pn, tn = mie_pi_tau(mu, 3)
    @test isapprox(pn[1], mie_pi_tau_n(mu, 1)[1])
    @test isapprox(pn[2], mie_pi_tau_n(mu, 2)[1])
    @test isapprox(pn[3], mie_pi_tau_n(mu, 3)[1])
    @test isapprox(tn[1], mie_pi_tau_n(mu, 1)[2])
    @test isapprox(tn[2], mie_pi_tau_n(mu, 2)[2])
    @test isapprox(tn[3], mie_pi_tau_n(mu, 3)[2])



    m = 1.33 + 0.01im
    x = 3

    qext1, qsca1, _, _ = mie_scattering(m, x)
    qext2, qsca2, _, _ = mie_scattering_n(m, x)

    @test isapprox(qext1, qext2)
    @test isapprox(qsca1, qsca2)

    m = 2 + 1im
    x = 5

    qext1, qsca1, _, _ = mie_scattering(m, x)
    qext2, qsca2, _, _ = mie_scattering_n(m, x)

    @test isapprox(qext1, qext2)
    @test isapprox(qsca1, qsca2)


    m = 1.33 + 0.01im
    x = 3

    mu = cos(0)
    s1s2 = mie_S1_S2(m, x, mu, 2)
    s1s2_n = mie_S1_S2_n(m, x, mu, 1) .+ mie_S1_S2_n(m, x, mu, 2)
    @test all(isapprox.(s1s2, s1s2_n))


    mu = cos(π / 3)
    s1s2 = mie_S1_S2(m, x, mu, 2)
    s1s2_n = mie_S1_S2_n(m, x, mu, 1) .+ mie_S1_S2_n(m, x, mu, 2)
    @test all(isapprox.(s1s2, s1s2_n))

    mu = cos(π / 2)
    s1s2 = mie_S1_S2(m, x, mu, 2)
    s1s2_n = mie_S1_S2_n(m, x, mu, 1) .+ mie_S1_S2_n(m, x, mu, 2)
    @test all(isapprox.(s1s2, s1s2_n))


    m = 2 + 1im
    x = 5

    mu = cos(0)
    s1s2 = mie_S1_S2(m, x, mu, 2)
    s1s2_n = mie_S1_S2_n(m, x, mu, 1) .+ mie_S1_S2_n(m, x, mu, 2)
    @test all(isapprox.(s1s2, s1s2_n))

    mu = cos(π / 3)
    s1s2 = mie_S1_S2(m, x, mu, 2)
    s1s2_n = mie_S1_S2_n(m, x, mu, 1) .+ mie_S1_S2_n(m, x, mu, 2)
    @test all(isapprox.(s1s2, s1s2_n))

    mu = cos(π / 2)
    s1s2 = mie_S1_S2(m, x, mu, 2)
    s1s2_n = mie_S1_S2_n(m, x, mu, 1) .+ mie_S1_S2_n(m, x, mu, 2)
    @test all(isapprox.(s1s2, s1s2_n))


    m = 1.33 + 0.01im
    x = 3

    theta = 0
    r1 = scattering_function(theta, m, x)
    r2 = scattering_function_n(theta, m, x)
    @test all(isapprox.(r1,r2))


    theta = π/3
    r1 = scattering_function(theta, m, x)
    r2 = scattering_function_n(theta, m, x)
    @test all(isapprox.(r1,r2))

    theta = π/2
    r1 = scattering_function(theta, m, x)
    r2 = scattering_function_n(theta, m, x)
    @test all(isapprox.(r1,r2))


    m = 2 + 1im
    x = 5

    theta = 0
    r1 = scattering_function(theta, m, x)
    r2 = scattering_function_n(theta, m, x)
    @test all(isapprox.(r1,r2))

    theta = π/3
    r1 = scattering_function(theta, m, x)
    r2 = scattering_function_n(theta, m, x)
    @test all(isapprox.(r1,r2))

    theta = π/2
    r1 = scattering_function(theta, m, x)
    r2 = scattering_function_n(theta, m, x)
    @test all(isapprox.(r1,r2))


end