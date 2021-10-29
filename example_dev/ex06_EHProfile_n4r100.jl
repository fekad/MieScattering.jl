using jlmie
using Plots

# structure
nmat = 4  # refractive index of a sphere
radius = 100e-9

# environment
nenv = 1.00

# wavelength
lbd0 = range(400e-9, 1000e-9, length=601)

# xy plot range
x = range(-150e-9, 150e-9, length=151)
y = 0.0
z = range(-150e-9, 150e-9, length=151)

xp = x * 1e9; yp = y * 1e9; zp = z * 1e9;

# settings
nmax = -1  # -1: namx large enough (determined from x)

Et, Er, Ep = jlmie.jlmie_Enf(nmat, radius, lbd0, nenv, x, y, z, nmax)
