# Structure
# epsilon_media: real, positive
s = Sphere(radius, epsilon, [epsilon_media])

sol = solve(s, x, MieAuto())
# sol = solve(s, wavelenght, MieAuto()) ??

qsca = scattering(sol)
qext = scattering(sol)

sol = solve(s, x, [model])
sol = solve(s, x, Rayleigh())
sol = solve(s, x, MieLarge())
sol = solve(s, x, MieConductive())

# =====================================================================================================
# Structure
# epsilon_media: real, positive
s = Sphere(radius, epsilon, [epsilon_media])


model = MieAuto(s, wavelenght) # ??

# model = Rayleigh(s, x)
# model = MieLarge(s, x)
# model = MieConductive(s, x)

qsca = scattering(model)
qext = scattering(model)


# =====================================================================================================
# Structure
# epsilon_media: real, positive
s = Sphere(radius, epsilon, [epsilon_media])

m, x = mx(s, wavelenght)
# cases:
# s : Sphere or Vector(Sphere)
# wavelenght

model = Mie(m, x)

qsca = scattering(model)
qext = scattering(model)
