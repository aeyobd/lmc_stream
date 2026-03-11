using LilGuys
using Agama


function (@main)(ARGS)
    @assert length(ARGS) == 1 "usage: make_init.jl num_particles"
    N = parse(Int, ARGS[1])

    prof_halo = Agama.Potential(file="./profile.ini")
    pot_lmc_init = Agama.Potential(file="./potential_lmc_init.ini")
    df = Agama.DistributionFunction(type="quasiSpherical", potential=pot_lmc_init._py, density=prof_halo._py)

    gm = Agama.GalaxyModel(pot_lmc_init, df)
    units = Agama.VASILIEV_UNITS
    
    positions, velocities = Agama.sample(gm, N, units)
    snap = LilGuys.Snapshot(positions, velocities, zeros(N))

    LilGuys.write("initial.hdf5", snap)
end
