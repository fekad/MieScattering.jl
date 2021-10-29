using jlmie
using Plots

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
theta = range(0, 2π, length=73) # per 5 degrees with length 37
phi = 0.0
# phi = pi / 2

# settings
nmax = -1  # -1: namx large enough (determined from x)

Isff = zeros(length(theta))
for i = 1:length(theta)
    Isff[i], _, _ = jlmie_Isff(nmat, radius, lbd0, nenv, theta[i], phi, nmax)
end
Isff = Isff / maximum(Isff)

plt = plot(theta, Isff, proj=:polar, legend=false);
display(plt)

