### A Pluto.jl notebook ###
# v0.20.23

using Markdown
using InteractiveUtils

# ╔═╡ e362baa0-1ef7-11f1-807a-f9861a7276f3
begin
	import Pkg; Pkg.activate()

	using LilGuys
	using CairoMakie
	using Arya
end

# ╔═╡ fa5d0776-fcba-493a-a71b-c3f795975a51
using Agama

# ╔═╡ 63481f41-6a7d-4304-a463-ba699866b1a6
using GeoMakie

# ╔═╡ 3fa62f9b-5430-4f95-af2b-8ea490b0b8df
using Printf

# ╔═╡ 5f47b6e3-6a3e-456f-acfe-14a9ca1fe916
include("plot_utils.jl")

# ╔═╡ fbfdd2f6-a017-4d03-88b7-e69fbca124a3
import LinearAlgebra: ⋅

# ╔═╡ 8ccef492-e3c2-461e-af79-9765db814efc
module Utils
	include("simulations/utils.jl")
end

# ╔═╡ 589f3f59-6687-499b-a7e0-364e59691577
md"""
# Data loading
"""

# ╔═╡ be2c76d4-fb9d-4d1e-a9ea-7370c824e84a
simname = "L2M10first"

# ╔═╡ d8dd1e3c-6587-4282-856a-88ca22575f92
orbits = LilGuys.read_struct_from_hdf5(joinpath("simulations", simname, "orbits.hdf5"), Utils.Orbits)

# ╔═╡ 1fa9621a-5d30-4b84-9d99-c5df16c4e951
lmc_orbit = Utils.get_lmc_orbit(joinpath("simulations", simname))

# ╔═╡ bb44e02d-72e3-43b0-a148-3f27a9b0d944
lmc_orbit_resampled = LilGuys.resample(lmc_orbit, orbits.times)

# ╔═╡ 597099e1-ca84-457e-bbdd-11600d134fcb
times = orbits.times * T2GYR

# ╔═╡ 38797249-6d71-4476-9a72-691f0125f88d
Nt = length(times)

# ╔═╡ c9a10db6-fe27-47b6-aa48-98777b250b70
function dist_to_size(dist)
	return (dist / 50)^-1
end

# ╔═╡ e0ce9870-0fc3-4e3b-add2-ad02bd6efd3b
let 
	fig = Figure()
	ax = GeoAxis(fig[1,1];
	    dest = "+proj=hammer",
	    limits = (0., 360, -90, 90),
	    yticklabelsvisible=false,
	    xticklabelsvisible=false,
	    yticklabelsize=8,
	    xgridwidth=0.3,
	    ygridwidth=0.3,
		xgridcolor = (:black, 0.2),
		ygridcolor = (:black, 0.2),
	    valign=:center,
	             
	)
	xlims!(-180, 180)

	fig

end

# ╔═╡ b56933d3-3ac8-4272-bc2f-a83537ef2ab7
function get_time_label(index)

	time = orbits.times[index] * T2GYR
    label =  @sprintf("today %+2.1f Gyr", time)
    label = replace(label, "-" => "– ") 
    label = replace(label, "+" => "+ ")

	return label

end

# ╔═╡ 81366c68-2390-4475-b59b-eb527b90f217
function add_time!(ax, label; font="Arial", fontsize=7.5)
    text!(ax, 1, 0.0, text=label, color=:grey, align=(:right, :bottom),
        font=font, fontsize=fontsize, space=:relative)
end


# ╔═╡ 6900daeb-ff90-4100-8860-a3607a20a542
function to_icrs(positions)

	gc = [Galactocentric(positions[:, i], zeros(3)) for i in 1:size(positions, 2)]
	icrs = LilGuys.transform.(ICRS, gc)

	return [c.ra for c in icrs], [c.dec for c in icrs], [c.distance for c in icrs]
end

# ╔═╡ c64a7828-5042-4e12-83d3-3d3a2fc91c30
function get_lmc_position(timestep)
	return to_icrs(lmc_orbit_resampled.positions[:, timestep])
end

# ╔═╡ d07cabb4-a3b4-4d38-88de-e8d811e2101c
function get_positions(timestep)
	positions = Utils.get_positions(orbits, timestep)

	return to_icrs(positions)
end

# ╔═╡ bb2334c4-af8a-41a2-9d12-36d5286f77b6
scatter(get_positions(5)[1:2]...)

# ╔═╡ 769a3a7e-531a-48e9-b113-97ac45754440
idxs_wrong = get_positions(length(times))[end] .== 8.122

# ╔═╡ 4ed63804-53b1-40d0-84d9-69711cd8469d
get_positions(1)[end][idxs_wrong]

# ╔═╡ 5b90e09b-6df6-4b2f-8e0f-7f7e45ddeae2
let 
	fig = Figure()
	ax = GeoAxis(fig[1,1];
	    dest = "+proj=hammer",
	    limits = (0., 360, -90, 90),
	    yticklabelsvisible=false,
	    xticklabelsvisible=false,
	    yticklabelsize=8,
	    xgridwidth=0.3,
	    ygridwidth=0.3,
		xgridcolor = (:black, 0.2),
		ygridcolor = (:black, 0.2),
	    valign=:center,
	             
	)
	xlims!(-180, 180)

	lmc_size_scale = 5

	xs, ys, dists = get_positions(1)

	x_obs = Observable(xs)
	y_obs = Observable(ys)
	s_obs = Observable(dist_to_size.(dists))

	x, y, dist = get_lmc_position(1)
	x_lmc = Observable(x)
	y_lmc = Observable(y)
	s_lmc = Observable(lmc_size_scale*dist_to_size.(dist))

	scatter!(x_obs, y_obs, color=:black, markersize=s_obs, alpha=0.05)

	scatter!(x_lmc, y_lmc, color=COLORS[3], marker=:circle, markersize=s_lmc)

	
	time = Observable(get_time_label(1))
	add_time!(ax, time)

	

	record(fig, joinpath("simulations", simname, "sky_animation.mp4"), 1:Nt) do idx
		x_obs[], y_obs[], dists = get_positions(idx)
		s_obs[] = dist_to_size.(dists)

		x_lmc[], y_lmc[], dist = get_lmc_position(idx)
		s_lmc[] = lmc_size_scale*dist_to_size.(dist)
		
		time[] = get_time_label(idx)

	end

	fig

end

# ╔═╡ Cell order:
# ╠═e362baa0-1ef7-11f1-807a-f9861a7276f3
# ╠═fbfdd2f6-a017-4d03-88b7-e69fbca124a3
# ╠═fa5d0776-fcba-493a-a71b-c3f795975a51
# ╠═63481f41-6a7d-4304-a463-ba699866b1a6
# ╠═8ccef492-e3c2-461e-af79-9765db814efc
# ╠═5f47b6e3-6a3e-456f-acfe-14a9ca1fe916
# ╠═3fa62f9b-5430-4f95-af2b-8ea490b0b8df
# ╠═589f3f59-6687-499b-a7e0-364e59691577
# ╠═be2c76d4-fb9d-4d1e-a9ea-7370c824e84a
# ╠═d8dd1e3c-6587-4282-856a-88ca22575f92
# ╠═1fa9621a-5d30-4b84-9d99-c5df16c4e951
# ╠═bb44e02d-72e3-43b0-a148-3f27a9b0d944
# ╠═c64a7828-5042-4e12-83d3-3d3a2fc91c30
# ╠═d07cabb4-a3b4-4d38-88de-e8d811e2101c
# ╠═bb2334c4-af8a-41a2-9d12-36d5286f77b6
# ╠═597099e1-ca84-457e-bbdd-11600d134fcb
# ╠═38797249-6d71-4476-9a72-691f0125f88d
# ╠═c9a10db6-fe27-47b6-aa48-98777b250b70
# ╠═e0ce9870-0fc3-4e3b-add2-ad02bd6efd3b
# ╠═769a3a7e-531a-48e9-b113-97ac45754440
# ╠═4ed63804-53b1-40d0-84d9-69711cd8469d
# ╠═5b90e09b-6df6-4b2f-8e0f-7f7e45ddeae2
# ╠═b56933d3-3ac8-4272-bc2f-a83537ef2ab7
# ╠═81366c68-2390-4475-b59b-eb527b90f217
# ╠═6900daeb-ff90-4100-8860-a3607a20a542
