module MieScattering

    using SpecialFunctions

    export scattering, absorption, extinction, backscattering

    export Scatterer, mx
    include("interface.jl")

    export Mie, mie_scattering, scattering_function
    include("mie.jl")

    export Rayleigh, rayleigh_scattering
    include("rayleigh.jl")


end