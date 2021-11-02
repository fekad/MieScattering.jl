using MieScattering
using Plots

# structure
nmat = 4        # refractive index of the sphere
radius = 100e-9

# environment
nenv = 1.00

# wavelength
lbd0 = range(400e-9, 1000e-9, length=601)

# settings
nmax = -1  # -1: namx large enough (determined from x)

m, x = mie_mx(nmat, radius, lbd0, nenv)
Qsca = mie_Qsca(m, x, nmax)
Qsca_ED = mie_Qsca_a_n(m, x, 1)
Qsca_MD = mie_Qsca_b_n(m, x, 1)
Qsca_EQ = mie_Qsca_a_n(m, x, 2)
Qsca_MQ = mie_Qsca_b_n(m, x, 2)
Qsca_EH = mie_Qsca_a_n(m, x, 3)
Qsca_MH = mie_Qsca_b_n(m, x, 3)

lbdp = lbd0 .* 1e9
plt = plot(lbdp, Qsca,
    xlabel="Wavelength (nm)",
    ylabel="Scattering efficiency",
    # legend=false,
    label="Total",
    # xlims=(-3, 3),
    # ylims=(-3, 3),
    # aspect_ratio=0.5,
    # title="TITLE",
    # linecolor=:blue,
    # linewidth=5,
    # linestyle=:dot,
    # size=(400, 300),
)
plot!(lbdp, Qsca_ED, label="ED", linestyle=:dash)
plot!(lbdp, Qsca_MD, label="MD", linestyle=:dash)
plot!(lbdp, Qsca_EQ, label="EQ", linestyle=:dash)
plot!(lbdp, Qsca_MQ, label="MQ", linestyle=:dash)
plot!(lbdp, Qsca_EH, label="EH", linestyle=:dash)
plot!(lbdp, Qsca_MH, label="MH", linestyle=:dash)
display(plt)

# find peaks
println("maximum of ED is at $(lbdp[findmax(Qsca_ED)[2]])")
println("maximum of MD is at $(lbdp[findmax(Qsca_MD)[2]])")
println("maximum of EQ is at $(lbdp[findmax(Qsca_EQ)[2]])")
println("maximum of MQ is at $(lbdp[findmax(Qsca_MQ)[2]])")
println("maximum of EH is at $(lbdp[findmax(Qsca_EH)[2]])")
println("maximum of MH is at $(lbdp[findmax(Qsca_MH)[2]])")