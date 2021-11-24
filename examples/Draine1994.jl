using Plots
using MieScattering

################################################################
# Fig 1.

# structure
m = 1.33 + .01im   # imag. component to demonstrate absorption

N = 1001
x = range(0.0001, 13, length=N)

Qsca = zeros(N)
Qabs = zeros(N)
for i = 1:N
    local qext, qsca, qabs, qback
    qext, qsca, qabs, qback = mie_scattering(m, x[i])
    Qsca[i] = qsca
    Qabs[i] = qabs
end


plt = plot(x, Qsca,
    xlabel="ka",
    ylabel="Scattering efficiency",
    legend=false,
    label="LABEL",
    title="TITLE",
    yscale=:log10,
    ylims=(0.005, maximum(Qsca))
)
plot!(x, Qabs)
display(plt)

println("max:" , maximum(Qsca), " ind: ",  argmax(Qsca))
println("min:" , minimum(Qsca), " ind: ",  argmin(Qsca))


################################################################
# Fig. 2

# structure
m1 = 2 + 1im   # imag. component to demonstrate absorption

N=1001
x = range(0, 13, length=N)

N = 1001
x = range(0.0001, 13, length=N)

Qsca = zeros(N)
Qabs = zeros(N)
for i = 1:N
    local qext, qsca, qabs, qback
    qext, qsca, qabs, qback = mie_scattering(m, x[i])
    Qsca[i] = qsca
    Qabs[i] = qabs
end


plt = plot(x, Qsca,
    xlabel="ka",
    ylabel="Scattering efficiency",
    legend=false,
    label="LABEL",
    title="TITLE",
    yscale=:log10,
    ylims=(0.005, maximum(Qsca))
)
plot!(x, Qabs)
display(plt)

println("max:" , maximum(Qsca), " ind: ",  argmax(Qsca))
println("min:" , minimum(Qsca), " ind: ",  argmin(Qsca))



################################################################
# Fig. 3.

m = 1.33 + .01im
x = 3

theta = range(0, π, length=2*180+1)
SL, SR, SU = scattering_function(theta, m, x)
SU /= π * x^2

angles = theta*180/π
# TODO: check why sqrt
plot(
    angles, SU,
    label="Unpolarized",
    ylabel="Intensity (\$|S|^2}\$)",
    xlabel="phase angle ϴ",
    yaxis=:log10
)

qext, qsca, qabs, qback = mie_scattering(m, x)
plot!(angles, qsca/4π * ones(length(angles)))


# plot!(theta*180/π, SL, label="left-polarized (perpendicular)")
# plot!(theta*180/π, SR, label="right-polarized (parallel)")

println("max:" , maximum(SU), " ind: ", angles[argmax(SU)])
println("min:" , minimum(SU), " ind: ", angles[argmin(SU)])




m = 1.33+0.01im
x = 7

theta = range(0, π, length=2*180+1)
SL, SR, SU = scattering_function(theta, m, x)
SU /= π * x^2

angles = theta*180/π
# TODO: check why sqrt
plot(
    angles, SU,
    label="Unpolarized",
    ylabel="Intensity (\$|S|^2}\$)",
    xlabel="phase angle ϴ",
    yaxis=:log10
)

qext, qsca, qabs, qback = mie_scattering(m, x)
plot!(angles, qsca/4π * ones(length(angles)))


# plot!(theta*180/π, SL, label="left-polarized (perpendicular)")
# plot!(theta*180/π, SR, label="right-polarized (parallel)")

println("max:" , maximum(SU), " ind: ", angles[argmax(SU)])
println("min:" , minimum(SU), " ind: ", angles[argmin(SU)])


m = 2+1im
x = 7

theta = range(0, π, length=2*180+1)
SL, SR, SU = scattering_function(theta, m, x)
SU /= π * x^2

angles = theta*180/π
# TODO: check why sqrt
plot(
    angles, SU,
    label="Unpolarized",
    ylabel="Intensity (\$|S|^2}\$)",
    xlabel="phase angle ϴ",
    yaxis=:log10
)

qext, qsca, qabs, qback = mie_scattering(m, x)
plot!(angles, qsca/4π * ones(length(angles)))


# plot!(theta*180/π, SL, label="left-polarized (perpendicular)")
# plot!(theta*180/π, SR, label="right-polarized (parallel)")

println("max:" , maximum(SU), " ind: ", angles[argmax(SU)])
println("min:" , minimum(SU), " ind: ", angles[argmin(SU)])

