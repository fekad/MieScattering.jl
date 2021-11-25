# Note: this file contains ideas/palns for future interface

# cases:
# sphere : Sphere or Vector(Sphere)
# wavelenght: single or vector

# =====================================================================================================

# Structure
s = Sphere(radius, epsilon)

# epsilon_media: real, positive
sol = solve(s, x, epsilon_media, MieAuto())
# sol = solve(s, wavelenght, epsilon_media, MieAuto()) ??
sol = solve(s, x, Rayleigh())
sol = solve(s, x, MieLarge())
sol = solve(s, x, MieConductive())

qsca = scattering(sol)
qext = extinction(sol)


# =====================================================================================================
# Structure
s = Sphere(radius, epsilon)

# epsilon_media: real, positive
model = MieAuto(s, wavelenght, epsilon_media) # ??

# model = Rayleigh(s, x)
# model = MieLarge(s, x)
# model = MieConductive(s, x)

qsca = scattering(model)
qext = extinction(model)


# =====================================================================================================
# Structure
s = Sphere(radius, epsilon)

# epsilon_media: real, positive
m, x = mx(s, wavelenght, epsilon_media)

model = Mie(m, x)

qsca = scattering(model)
qext = scattering(model)
