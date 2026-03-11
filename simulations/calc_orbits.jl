using LilGuys
using Agama
import CSV
using DataFrames: DataFrame


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


function get_lmc_orbit(input)
    lmc_file = joinpath(input, "trajlmc.txt")
    df_lmc = lmc_traj = CSV.read(lmc_file, DataFrame, delim=" ", header = [:time, :x, :y, :z, :v_x, :v_y, :v_z], ignorerepeated=true, ntasks=1)

    pos = hcat(df_lmc.x, df_lmc.y, df_lmc.z)'
    vel = hcat(df_lmc.v_x, df_lmc.v_y, df_lmc.v_z)'

    # convert to code units
    t = df_lmc.time .* Agama.time_scale(Agama.VASILIEV_UNITS) 
    vel .*= Agama.velocity_scale(Agama.VASILIEV_UNITS) 

    orbit_lmc = Orbit(times=t, positions=pos, velocities=vel)
end


"""
    write_orbits(output, orbits; N_max)

Write the first `N_max` orbits to a file "orbits.hdf5" in `output`.
"""
function write_orbits(orbits; N_max=100)
    filename = "orbits.hdf5"

    N_max = min(length(orbits), N_max)
    structs = [(string(i) => orbit) for (i, orbit) in enumerate(orbits[1:N_max])]

    LilGuys.write_structs_to_hdf5(filename, structs)
end



function write_final(orbits)
    positions = hcat([orbit.positions[:, end] for orbit in orbits]...)
    velocities = hcat([orbit.velocities[:, end] for orbit in orbits]...)

    snap = Snapshot(positions, velocities, zeros(length(orbits)))


    LilGuys.write("final_positions.hdf5", snap)
end
