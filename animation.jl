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
simname = "L3M11"

# ╔═╡ 36c05e90-42ac-4e8c-86e1-306420cc3ef1
Utils.Orbits

# ╔═╡ d8dd1e3c-6587-4282-856a-88ca22575f92
orbits = LilGuys.read_struct_from_hdf5(joinpath("simulations", simname, "orbits.hdf5"), Utils.Orbits)

# ╔═╡ 1fa9621a-5d30-4b84-9d99-c5df16c4e951
lmc_orbit = Utils.get_lmc_orbit(joinpath("simulations", simname))

# ╔═╡ 7b820be9-c4f6-4e60-96bd-eeda903db9e6
X_SUN = [-8.122, 0., 0.]


# ╔═╡ e7e650a8-3015-47b6-b294-645d194c8192
x_vec = [0, 1, 0]

# ╔═╡ 1d84e791-a3e3-4df1-aa06-47849b1d7bb0
y_vec = [0, 0, 1]

# ╔═╡ bb44e02d-72e3-43b0-a148-3f27a9b0d944
lmc_orbit_resampled = LilGuys.resample(lmc_orbit, orbits.times)

# ╔═╡ c64a7828-5042-4e12-83d3-3d3a2fc91c30
function get_lmc_position(timestep)
	xs = lmc_orbit_resampled.positions[:, timestep] ⋅ x_vec 
	ys = lmc_orbit_resampled.positions[:, timestep] ⋅ y_vec

	return xs, ys
end

# ╔═╡ d07cabb4-a3b4-4d38-88de-e8d811e2101c
function get_positions(timestep)
	positions = Utils.get_positions(orbits, timestep)
	xs = vec(x_vec' * positions )
	ys = vec(y_vec' * positions)

	return xs, ys
end

# ╔═╡ bb2334c4-af8a-41a2-9d12-36d5286f77b6
scatter(get_positions(5)...)

# ╔═╡ 72f3c2cc-9ff4-4ed3-bd29-534e57901447
pot = Agama.Potential(file=joinpath("simulations", simname, "potential_mw_init.ini"))

# ╔═╡ c7fa5717-65ed-419a-8ee5-3c430042c983
isodensity_xz = integrate_isodensity_2d(pot, beta=π/2)

# ╔═╡ 79a33e1e-7d06-40af-86c7-71f2ee901c71
function plot_mw_iso!(ax)
	poly!(ax, isodensity_xz..., color=COLORS[2])
end

# ╔═╡ 827a7a67-0b66-4dd0-8185-f20da42ae0cc
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel = "Y", 
		ylabel="Z",
		aspect=DataAspect()
	)

	plot_mw_iso!(ax)

	fig

end

# ╔═╡ 597099e1-ca84-457e-bbdd-11600d134fcb
times = orbits.times * T2GYR

# ╔═╡ 38797249-6d71-4476-9a72-691f0125f88d
Nt = length(times)

# ╔═╡ 4b0ff6a6-76e4-4b35-bdda-1ff5977f60fa
function add_scalebar!(ax, xrange, yrange, scalebar::Real; font="Arial", fontsize=7.5)
    x1 = xrange[1] + 0.05 * (xrange[2] - xrange[1])
    y1 = yrange[1] + 0.05 * (yrange[2] - yrange[1])
    x2 = x1 + scalebar
    y2 = y1

    lines!(ax, [x1, x2], [y1, y2], color=:grey)

    # Format scalebar label
    label = @sprintf("%.0f kpc", scalebar)
    text!(ax, x1, y1, text=label, color=:grey, font=font, fontsize=fontsize)
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
    text!(ax, 0.95, 0.05, text=label, color=:grey, align=(:right, :bottom),
        font=font, fontsize=fontsize, space=:relative)
end


# ╔═╡ 5b90e09b-6df6-4b2f-8e0f-7f7e45ddeae2
let
	xrange = yrange = (-500, 500)
	fig = Figure()
	ax = Axis(fig[1,1], 
			  limits=( xrange, yrange),
			  aspect = DataAspect(),
			  xlabel = "x / kpc",
			  ylabel = "y / kpc",
			 )

	xs, ys = get_positions(1)

	x_obs = Observable(xs)
	y_obs = Observable(ys)

	x, y = get_lmc_position(1)
	x_lmc = Observable(x)
	y_lmc = Observable(y)

	scatter!(x_obs, y_obs, color=:black, markersize=0.3, alpha=1)

	scatter!(x_lmc, y_lmc, color=COLORS[3], marker=:circle)

	
	plot_mw_iso!(ax)

	time = Observable(get_time_label(1))
	add_time!(ax, time)

	
	add_scalebar!(ax, xrange, yrange, 100)

	hidexdecorations!()
	hideydecorations!()

	record(fig, joinpath("simulations", simname, "animation.mp4"), 1:Nt) do idx
		x_obs[], y_obs[] = get_positions(idx)

		x_lmc[], y_lmc[] = get_lmc_position(idx)
		time[] = get_time_label(idx)

	end

	fig

end

# ╔═╡ Cell order:
# ╠═e362baa0-1ef7-11f1-807a-f9861a7276f3
# ╠═fbfdd2f6-a017-4d03-88b7-e69fbca124a3
# ╠═fa5d0776-fcba-493a-a71b-c3f795975a51
# ╠═8ccef492-e3c2-461e-af79-9765db814efc
# ╠═3fa62f9b-5430-4f95-af2b-8ea490b0b8df
# ╠═589f3f59-6687-499b-a7e0-364e59691577
# ╠═be2c76d4-fb9d-4d1e-a9ea-7370c824e84a
# ╠═36c05e90-42ac-4e8c-86e1-306420cc3ef1
# ╠═d8dd1e3c-6587-4282-856a-88ca22575f92
# ╠═1fa9621a-5d30-4b84-9d99-c5df16c4e951
# ╠═7b820be9-c4f6-4e60-96bd-eeda903db9e6
# ╠═e7e650a8-3015-47b6-b294-645d194c8192
# ╠═1d84e791-a3e3-4df1-aa06-47849b1d7bb0
# ╠═bb44e02d-72e3-43b0-a148-3f27a9b0d944
# ╠═c64a7828-5042-4e12-83d3-3d3a2fc91c30
# ╠═d07cabb4-a3b4-4d38-88de-e8d811e2101c
# ╠═bb2334c4-af8a-41a2-9d12-36d5286f77b6
# ╠═5f47b6e3-6a3e-456f-acfe-14a9ca1fe916
# ╠═72f3c2cc-9ff4-4ed3-bd29-534e57901447
# ╠═c7fa5717-65ed-419a-8ee5-3c430042c983
# ╠═827a7a67-0b66-4dd0-8185-f20da42ae0cc
# ╠═79a33e1e-7d06-40af-86c7-71f2ee901c71
# ╠═597099e1-ca84-457e-bbdd-11600d134fcb
# ╠═38797249-6d71-4476-9a72-691f0125f88d
# ╠═5b90e09b-6df6-4b2f-8e0f-7f7e45ddeae2
# ╠═4b0ff6a6-76e4-4b35-bdda-1ff5977f60fa
# ╠═b56933d3-3ac8-4272-bc2f-a83537ef2ab7
# ╠═81366c68-2390-4475-b59b-eb527b90f217
