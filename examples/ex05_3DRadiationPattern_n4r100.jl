using jlmie
using Plots
plotly()

# settings
# structure
nmat = 4  # refractive index of a sphere
radius = 100e-9

# environment
nenv = 1.00

# wavelength
# lbd0 = 618e-9  # ED
# lbd0 = 831e-9  # MD
lbd0 = 918e-9  # 1nd Kerker
# lbd0 = 756e-9  # 2nd Kerker


# angular range
θ = range(0, π, length=73) # per 2.5 degrees with length 73
ϕ = range(0, 2π, length=145) # per 2.5 degrees with length 145

# settings
nmax = -1  # -1: namx large enough (determined from x)

Isff = zeros(length(θ), length(ϕ))
for i = 1:length(θ)
    for j = 1:length(ϕ)
        Isff[i,j], _, _ = jlmie_Isff(nmat, radius, lbd0, nenv, θ[i], ϕ[j], nmax)
    end
end
Isff = Isff / maximum(Isff)
mat_x = Isff .* (sin.(θ) * cos.(ϕ)')
mat_y = Isff .* (sin.(θ) * sin.(ϕ)')
mat_z = Isff .* (cos.(θ))

# show results
plt = plot(mat_x,mat_y,mat_z,
        fill_z=abs.(Isff),
        st=:surface,
        camera=(45, 30), # azimuth, elevate
        xlabel="x",
        ylabel="y",
        zlabel="z",
        # colorbar_title="Scattering intensity (arb. units)",
        xlims=(-1, 1),
        ylims=(-1, 1),
        zlims=(-1, 1),
);
display(plt)

