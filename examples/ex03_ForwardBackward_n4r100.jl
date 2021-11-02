using MieScattering
using Plots

# structure
nmat = 4  # refractive index of a sphere
radius = 100e-9

# environment
nenv = 1.00

# wavelength
lbd0 = range(400e-9, 1000e-9, length=601)

# angular range
theta = [0 pi]
phi = [0 0]

# other settings
nmax = -1  # -1: namx large enough (determined from x)

Isff_F, _, _ = mie_Isff(nmat, radius, lbd0, nenv, theta[1], phi[1], nmax)
Isff_B, _, _ = mie_Isff(nmat, radius, lbd0, nenv, theta[2], phi[2], nmax)

lbdp = lbd0 .* 1e9
plt = plot(lbdp,Isff_F,
    xlabel="Wavelength (nm)",
    ylabel="Scattering efficiency",
    # legend=false,
    label="Forward",
    # xlims=(-3, 3),
    # ylims=(-3, 3),
    # aspect_ratio=0.5,
    # title="TITLE",
    # linecolor=:blue,
    # linewidth=5,
    # linestyle=:dot,
    # size=(400, 300)
)
plot!(lbdp, Isff_B, label="Backward")
display(plt)

# find peaks
lbdp_peaks = zeros(2)
lbdp_peaks[1] = lbdp[findmax(Isff_F ./ Isff_B)[2]]
lbdp_peaks[2] = lbdp[200 + findmax(Isff_B[201:end] ./ Isff_F[201:end])[2]]  # avoid high order modes

println("maximum of F/B ratio is at $(lbdp_peaks[1])")
println("maximum of B/F ratio is at $(lbdp_peaks[2])")