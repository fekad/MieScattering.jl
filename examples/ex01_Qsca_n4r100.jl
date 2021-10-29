using jlmie
using Plots

# structure
nmat = 4  # refractive index of a sphere
radius = 100e-9

# environment
nenv = 1.00

# wavelength
lbd0 = range(400e-9, 1000e-9, length=601)

# settings
nmax = -1  # -1: namx large enough (determined from x)

m, x = jlmie_mx(nmat, radius, lbd0, nenv)
Qsca = jlmie_Qsca(m, x, nmax)

lbdp = lbd0 .* 1e9
plt = plot(lbdp, Qsca,
    xlabel="Wavelength (nm)",
    ylabel="Scattering efficiency",
    legend=false,
    # label="LABEL",
    # xlims=(-3, 3),
    # ylims=(-3, 3),
    # aspect_ratio=0.5,
    # title="TITLE",
    # linecolor=:blue,
    # linewidth=5,
    # linestyle=:dot,
    # size=(400, 300),
)
display(plt)

