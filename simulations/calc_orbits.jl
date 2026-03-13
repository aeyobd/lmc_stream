using LilGuys
using Agama

include(joinpath(@__DIR__, "utils.jl"))


function (@main)(ARGS)
    pot = Agama.Potential(file="agama_potential.ini")
    units = Agama.VASILIEV_UNITS

    snap = Snapshot("initial.hdf5")

    lmc_orbit = get_lmc_orbit(".")
    pos0 = lmc_orbit.positions[:, 1]
    vel0 = lmc_orbit.velocities[:, 1]
    time0 = lmc_orbit.times[1]

    coords = [LilGuys.Galactocentric(snap.positions[:, i] .+ pos0, V2KMS * (snap.velocities[:, i] .+ vel0)) for i in 1:length(snap)]

    orbits = LilGuys.agama_orbit(pot, coords; timerange=(time0, 0), N=100, agama_units=units)

    write_final(orbits)
    write_orbits(orbits)
end



"""
    write_orbits(output, orbits; N_max)

Write the first `N_max` orbits to a file "orbits.hdf5" in `output`.
"""
function write_orbits(orbits)
    filename = "orbits.hdf5"

    println("writing orbits")

    positions = cat((orbit.positions for orbit in orbits)..., dims=3)
    velocities = cat((orbit.velocities for orbit in orbits)..., dims=3)

    times = orbits[1].times

    orbits = Orbits(positions, velocities, times)

    LilGuys.write_struct_to_hdf5(filename, orbits)
end



function write_final(orbits)
    positions = hcat([orbit.positions[:, end] for orbit in orbits]...)
    velocities = hcat([orbit.velocities[:, end] for orbit in orbits]...)

    snap = Snapshot(positions, velocities, zeros(length(orbits)))


    LilGuys.write("final_positions.hdf5", snap)
end
