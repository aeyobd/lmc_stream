import CSV
using DataFrames: DataFrame
using Agama
using LilGuys

Base.@kwdef struct Orbits
    positions::Array{Float64, 3}
    velocities::Array{Float64, 3}
    times::Vector{Float64}
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


function get_orbit(orbits::Orbits, index::Integer)
    return Orbit(
                 positions = orbits.positions[:, :, index],
                 velocities = orbits.velocities[:, :, index],
                 times = orbits.times
                )
end


function get_positions(orbits::Orbits, index::Integer)
    return orbits.positions[:, index, :]
end

function get_velocities(orbits::Orbits, index::Integer)
    return orbits.velocities[:, index, :]
end
