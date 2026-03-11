### A Pluto.jl notebook ###
# v0.20.23

using Markdown
using InteractiveUtils

# ╔═╡ 82bc2e30-1d7a-11f1-95f6-49e0b164b4f1
begin
	import Pkg; Pkg.activate()

	using LilGuys
	using CairoMakie
	using Arya
end

# ╔═╡ f031667d-e15f-43d7-a832-720e9b73d87c
using GeoMakie

# ╔═╡ 004d2da2-6227-4cd0-ade0-ece085dc06ea
CairoMakie.activate!(type=:png)

# ╔═╡ 936dd128-bfb0-4d5c-b754-8b09a144ff07
simname = "L2M10"

# ╔═╡ adebf069-b4e6-4dd5-bfe1-56b70a06e769
snap = Snapshot(joinpath("simulations", simname, "final_positions.hdf5"))

# ╔═╡ a14c6ff9-e8ef-4cd4-b169-38a3110cfa4b
module Utils
	include("./simulations/calc_orbits.jl")
end

# ╔═╡ 3dc377fe-698d-4715-98ad-20b633b830d2
orbit_lmc = Utils.get_lmc_orbit(joinpath("simulations", simname))

# ╔═╡ a4bcafe3-2f2f-4f43-9223-e70e8188dfc9
scatter(snap.positions)

# ╔═╡ 7ca0c8ee-83ae-4f1f-8c08-b16dc6a4cbc2
plot_rmax = 300

# ╔═╡ eb7ca0a8-7968-40a4-86be-6b0d8a45af55
plot_limits = ((-plot_rmax, plot_rmax), (-plot_rmax, plot_rmax), (-plot_rmax, plot_rmax))

# ╔═╡ f3fff48e-3a11-4dbe-a837-4c8ae36e804c
LilGuys.plot_xyz(snap.positions, plot=:scatter, limits=plot_limits, markersize=1, alpha=0.03, color=:black)

# ╔═╡ 913fa16e-f951-4caf-9359-823760fa47b0
function to_icrs(snap::Union{Snapshot, Orbit})

	gc_coords = [Galactocentric(snap.positions[:, i], snap.velocities[:, i] * V2KMS)
		  for i in 1:length(snap)]
	coords_icrs = LilGuys.transform.(ICRS, gc_coords)
end

# ╔═╡ 07f2c1a1-da0d-4576-9938-80bf5850d206
coords_icrs_lmc = to_icrs(orbit_lmc)

# ╔═╡ e4515a25-08cd-449b-a27d-2be47efe2283
gc_coords = [Galactocentric(snap.positions[:, i], snap.velocities[:, i] * V2KMS)
		  for i in 1:length(snap)]

# ╔═╡ e6b87de5-57bb-488c-bf59-6d853da81893
coords_gsr = LilGuys.transform.(GSR, coords_icrs)

# ╔═╡ 5817a668-73a0-488e-88de-08c62efc5a65
time0 = orbit_lmc.times[1]

# ╔═╡ a54824f1-1ae9-4b05-b508-abaef47456f6
let 
	fig = Figure()
	ax = GeoAxis(fig[1,1];
	    dest = "+proj=hammer",
	    limits = (0., 360, -90, 90),
	    yticklabelsvisible=false,
	    xticklabelsvisible=false,
	    yticklabelsize=8,
	    xgridwidth=0.5,
	    ygridwidth=0.5,
	    valign=:center,
	             
	)
	xlims!(-180, 180)
	
	s = (orbit_lmc.times .- time0)/abs(time0) * 5
	scatter!([c.ra for c in coords_icrs_lmc], [c.dec for c in coords_icrs_lmc], markersize =s, marker=:circle)
	
	scatter!([c.ra for c in coords_icrs], [c.dec for c in coords_icrs], markersize=1, marker=:circle, alpha=0.1)
	
	
	fig

end

# ╔═╡ c6331833-076a-47f5-9ae3-8694969fc806
idx_f_lmc = argmin(abs.(orbit_lmc.times))

# ╔═╡ e1609cc5-b155-48f5-993a-4acba620b742
icrs_lmc = coords_icrs_lmc[idx_f_lmc]

# ╔═╡ dd1adc3f-8ad3-4e1c-a984-24a4b04a33fa
ra0 = icrs_lmc.ra

# ╔═╡ ed0f8515-347c-4024-83f8-c167b2a60226
dec0 = icrs_lmc.dec

# ╔═╡ 98507ade-ed3e-40b0-991f-b0dcba5c232d
gsr_lmc = LilGuys.transform(GSR, icrs_lmc)

# ╔═╡ 13b80ec1-f1b9-4dc7-bd6e-39c6390b1b99
θ_lmc = atand(gsr_lmc.pmra, gsr_lmc.pmdec)

# ╔═╡ a72b3248-9431-4aa6-861c-ed975c6c6725
xi_p_lmc, eta_p_lmc = LilGuys.to_orbit_coords([c.ra for c in coords_icrs_lmc], [c.dec for c in coords_icrs_lmc], ra0, dec0, θ_lmc)

# ╔═╡ ecaa4c7c-ec85-4686-b2f0-c0814ae911b1
xi_p, eta_p = LilGuys.to_orbit_coords([c.ra for c in coords_icrs], [c.dec for c in coords_icrs], ra0, dec0, θ_lmc)

# ╔═╡ 8ecf873e-8874-427e-8e6a-dd2d795ae57e
function stream_plot(; colorbar_label=nothing, kwargs...)
	fig = Figure(size=(6, 3) .* 72)
	ax = Axis(fig[1,1], xlabel=L"\xi' / ^\circ", ylabel=L"\eta'/ ^\circ", aspect=DataAspect())

	s = (orbit_lmc.times .- time0)/abs(time0) * 5
	scatter!(xi_p_lmc, eta_p_lmc, markersize=s, marker=:circle)
	

	
	p = scatter!(xi_p, eta_p, markersize=1; kwargs...)

	if !isnothing(colorbar_label)
		Colorbar(fig[1, 2], p, label=colorbar_label)
	end
	
	fig

end

# ╔═╡ 0664169f-6b2f-421e-a0c7-6fac34f5256b
stream_plot(alpha=0.05, color=:black)

# ╔═╡ a946cea1-fe5b-4522-839b-9b60924394d9
stream_plot(color=[c.distance for c in coords_icrs], colorbar_label="distance / kpc", 
		   colorrange=(0, 200))

# ╔═╡ 9319d3f3-5fc2-4d9a-b37b-311a6c42e086
stream_plot(color=[c.radial_velocity for c in coords_icrs], colorbar_label="v los / km/s", )

# ╔═╡ 0ac64180-eeee-42e3-a765-3726e516c5dc
stream_plot(color=[c.pmra for c in coords_gsr], colorbar_label="pmra gsr / masyr", colorrange=(-1, 1))

# ╔═╡ 038bc3d9-ee0f-454e-88c8-7ece8da76cc6
stream_plot(color=[c.pmdec for c in coords_gsr], colorbar_label="pmdec gsr / masyr", colorrange=(-1, 1))

# ╔═╡ a8b21cd5-e574-43d9-843f-270608a77d2e
coords_icrs = to_icrs(snap)

# ╔═╡ 05662551-c5e1-48f3-b816-141082c59a5a
# ╠═╡ disabled = true
#=╠═╡
coords_icrs = LilGuys.transform.(ICRS, gc_coords)
  ╠═╡ =#

# ╔═╡ Cell order:
# ╠═82bc2e30-1d7a-11f1-95f6-49e0b164b4f1
# ╠═f031667d-e15f-43d7-a832-720e9b73d87c
# ╠═004d2da2-6227-4cd0-ade0-ece085dc06ea
# ╠═936dd128-bfb0-4d5c-b754-8b09a144ff07
# ╠═adebf069-b4e6-4dd5-bfe1-56b70a06e769
# ╠═a14c6ff9-e8ef-4cd4-b169-38a3110cfa4b
# ╠═3dc377fe-698d-4715-98ad-20b633b830d2
# ╠═a4bcafe3-2f2f-4f43-9223-e70e8188dfc9
# ╠═7ca0c8ee-83ae-4f1f-8c08-b16dc6a4cbc2
# ╠═eb7ca0a8-7968-40a4-86be-6b0d8a45af55
# ╠═f3fff48e-3a11-4dbe-a837-4c8ae36e804c
# ╠═913fa16e-f951-4caf-9359-823760fa47b0
# ╠═a8b21cd5-e574-43d9-843f-270608a77d2e
# ╠═07f2c1a1-da0d-4576-9938-80bf5850d206
# ╠═e6b87de5-57bb-488c-bf59-6d853da81893
# ╠═e4515a25-08cd-449b-a27d-2be47efe2283
# ╠═05662551-c5e1-48f3-b816-141082c59a5a
# ╠═5817a668-73a0-488e-88de-08c62efc5a65
# ╠═a54824f1-1ae9-4b05-b508-abaef47456f6
# ╠═13b80ec1-f1b9-4dc7-bd6e-39c6390b1b99
# ╠═dd1adc3f-8ad3-4e1c-a984-24a4b04a33fa
# ╠═ed0f8515-347c-4024-83f8-c167b2a60226
# ╠═e1609cc5-b155-48f5-993a-4acba620b742
# ╠═98507ade-ed3e-40b0-991f-b0dcba5c232d
# ╠═c6331833-076a-47f5-9ae3-8694969fc806
# ╠═a72b3248-9431-4aa6-861c-ed975c6c6725
# ╠═ecaa4c7c-ec85-4686-b2f0-c0814ae911b1
# ╠═8ecf873e-8874-427e-8e6a-dd2d795ae57e
# ╠═0664169f-6b2f-421e-a0c7-6fac34f5256b
# ╠═a946cea1-fe5b-4522-839b-9b60924394d9
# ╠═9319d3f3-5fc2-4d9a-b37b-311a6c42e086
# ╠═0ac64180-eeee-42e3-a765-3726e516c5dc
# ╠═038bc3d9-ee0f-454e-88c8-7ece8da76cc6
