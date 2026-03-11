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

# ╔═╡ 1bb94a81-b4e4-4123-8b36-f68145147471
import CSV

# ╔═╡ 2686ab62-1711-4dcb-97b9-e4ab322497d9
import Agama

# ╔═╡ fce610b5-cfbd-46fe-b570-4d5df7fbf5f0
import DataFrames: DataFrame, rename

# ╔═╡ 004d2da2-6227-4cd0-ade0-ece085dc06ea
CairoMakie.activate!(type=:png)

# ╔═╡ 936dd128-bfb0-4d5c-b754-8b09a144ff07
simname = "L2M10"

# ╔═╡ adebf069-b4e6-4dd5-bfe1-56b70a06e769
snap = Snapshot(joinpath("simulations", simname, "final_positions.hdf5"))

# ╔═╡ dc75f240-e1cd-455e-b431-4c3b7ad5343d
snap_i = Snapshot(joinpath("simulations", simname, "initial.hdf5"))

# ╔═╡ 17e92b56-2d13-4460-a541-80e1e10a30f3
pot_lmc_init = Agama.Potential(file=joinpath("simulations", simname, "potential_lmc_init.ini"))

# ╔═╡ 4c895cf8-4904-40e5-9778-6149999b1e14
r_i = radii(snap_i)

# ╔═╡ cc028c15-ccfe-4c42-b1a7-0e86c36e7e60
snap_i.potential = Agama.potential(pot_lmc_init, snap_i.positions, Agama.VASILIEV_UNITS)

# ╔═╡ a215ed7a-f605-4f43-b772-1e13398e6dbf
ϵ = LilGuys.specific_energy(snap_i)

# ╔═╡ f62020bb-142e-439f-b97d-4612a9ce9ffe
L = LilGuys.angular_momenta(snap_i)

# ╔═╡ a14c6ff9-e8ef-4cd4-b169-38a3110cfa4b
module Utils
	include("./simulations/calc_orbits.jl")
end

# ╔═╡ fab6087b-3f90-4793-a168-7e7c5e54f06d
@assert snap.index == snap_i.index

# ╔═╡ 3dc377fe-698d-4715-98ad-20b633b830d2
orbit_lmc = Utils.get_lmc_orbit(joinpath("simulations", simname))

# ╔═╡ 7ca0c8ee-83ae-4f1f-8c08-b16dc6a4cbc2
plot_rmax = 300

# ╔═╡ eb7ca0a8-7968-40a4-86be-6b0d8a45af55
plot_limits = ((-plot_rmax, plot_rmax), (-plot_rmax, plot_rmax), (-plot_rmax, plot_rmax))

# ╔═╡ 96452520-666e-4ae9-8353-e910a9472e75
hist(r_i, bins=100)

# ╔═╡ f3fff48e-3a11-4dbe-a837-4c8ae36e804c
LilGuys.plot_xyz(snap.positions, plot=:scatter, limits=plot_limits, markersize=1, alpha=0.03, color=:black)

# ╔═╡ 913fa16e-f951-4caf-9359-823760fa47b0
function to_icrs(snap::Union{Snapshot, Orbit})

	gc_coords = [Galactocentric(snap.positions[:, i], snap.velocities[:, i] * V2KMS)
		  for i in 1:length(snap)]
	coords_icrs = LilGuys.transform.(ICRS, gc_coords)
end

# ╔═╡ a8b21cd5-e574-43d9-843f-270608a77d2e
coords_icrs = to_icrs(snap)

# ╔═╡ 07f2c1a1-da0d-4576-9938-80bf5850d206
coords_icrs_lmc = to_icrs(orbit_lmc)

# ╔═╡ e6b87de5-57bb-488c-bf59-6d853da81893
coords_gsr = LilGuys.transform.(GSR, coords_icrs)

# ╔═╡ e4515a25-08cd-449b-a27d-2be47efe2283
gc_coords = [Galactocentric(snap.positions[:, i], snap.velocities[:, i] * V2KMS)
		  for i in 1:length(snap)]

# ╔═╡ 5817a668-73a0-488e-88de-08c62efc5a65
time0 = orbit_lmc.times[1]

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

# ╔═╡ e8b072d4-27e3-4dfc-b2ce-d55565cb2a0d
md"""
# Comparing against observed coordinates
"""

# ╔═╡ a8ef246a-b342-4e79-b794-a61089a79f81


# ╔═╡ 2c205f8f-9321-45dc-90fa-ff03fc4acf60
akshara_df = let
	df = CSV.read("all_MSS_dynamics.csv", DataFrame) 
	df[!, :radial_velocity] = df.rv
	df[!, :radial_velocity_error] = df.rv_error
	df[!, :distance] = df.dist
	df[!, :distance_error] = df.dist_error

	df[!, :xi_p], df[!, :eta_p] = LilGuys.to_orbit_coords(df.ra, df.dec, ra0, dec0, θ_lmc)
	df[!, :eta_p_error] .= 0
	
	df
end

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
	
	scatter!(akshara_df.ra, akshara_df.dec)
	fig

end

# ╔═╡ e5a3b299-bac3-4723-8f72-65bf58d21463
sim_df = let
	df = LilGuys.to_frame(coords_icrs)
	df[!, :xi_p], df[!, :eta_p] = LilGuys.to_orbit_coords(df.ra, df.dec, ra0, dec0, θ_lmc)

	df[!, :r_0] = r_i
	df[!, :energy_0] = ϵ
	df[!, :L_x] = L[1, :]
	df[!, :L_y] = L[2, :]
	df[!, :L_z] = L[3, :]
	df[!, :L_0] = radii(L)
	
	df
end

# ╔═╡ 7edd507b-2ca9-4502-94f6-567b75c2a8ef
df_lmc = let
	df = LilGuys.to_frame(coords_icrs_lmc)
	df[!, :xi_p], df[!, :eta_p] = LilGuys.to_orbit_coords(df.ra, df.dec, ra0, dec0, θ_lmc)

	df
end

# ╔═╡ ee3f0069-31eb-4472-bec3-598c98ef6f6b
akshara_icrs = LilGuys.coords_from_df(akshara_df)

# ╔═╡ b1587ebb-4353-4998-a221-97009ffe2a88
akshara_xi_p, akshara_eta_p = LilGuys.to_orbit_coords([c.ra for c in akshara_icrs], [c.dec for c in akshara_icrs], ra0, dec0, θ_lmc)

# ╔═╡ 8ecf873e-8874-427e-8e6a-dd2d795ae57e
function stream_plot(; colorbar_label=nothing, kwargs...)
	fig = Figure(size=(6, 3) .* 72)
	ax = Axis(fig[1,1], xlabel=L"\xi' / ^\circ", ylabel=L"\eta'/ ^\circ", aspect=DataAspect())

	s = (orbit_lmc.times .- time0)/abs(time0) * 5
	scatter!(xi_p_lmc, eta_p_lmc, markersize=s, marker=:circle)
	

	
	p = scatter!(xi_p, eta_p, markersize=1; kwargs...)

	scatter!(akshara_xi_p, akshara_eta_p)

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

# ╔═╡ cd631337-61fe-440b-a38e-1012aa1dd7c7
stream_plot(color=sim_df.r_0, colorbar_label="r_i", colorrange=(0, 50))

# ╔═╡ 49f5b518-1dec-4e98-8b49-3d89e704533a
function observable_vs_stream(ykey; ylims=(nothing, nothing),
							  colorbar_label=nothing, kwargs...)
	fig = Figure(size=(6, 3) .* 72)
	ax = Axis(fig[1,1], xlabel=L"\xi' / ^\circ", ylabel=string(ykey), 
			  limits=(-180, 180, ylims[1], ylims[2])
			 )

	s = (orbit_lmc.times .- time0)/abs(time0) * 5
	scatter!(df_lmc.xi_p, df_lmc[!, ykey], markersize=s, marker=:circle)
	

	
	p = scatter!(sim_df.xi_p, sim_df[!, ykey], markersize=1, color=:black, alpha=0.1)


	errorscatter!(akshara_df.xi_p, akshara_df[!, ykey], 
				  yerror=akshara_df[!, "$(ykey)_error"], color=COLORS[4])

	if !isnothing(colorbar_label)
		Colorbar(fig[1, 2], p, label=colorbar_label)
	end
	
	fig

end

# ╔═╡ 447b6b59-ca3d-40ff-9037-272b3ef5b0e1
observable_vs_stream(:eta_p, ylims=(-40, 40))

# ╔═╡ 922830a7-736d-4d49-9a83-247d3722ad95
hist(sim_df.xi_p, bins=300, normalization=:pdf, axis=(
	xlabel = L"\xi'", 
	ylabel="density"
))

# ╔═╡ 9d6d21e9-35ee-4fee-b1bc-2964ec935843
let
	fig = Figure()
	ax = Axis(fig[1,1],  ylabel="L")

	scatter!(sim_df.xi_p, sim_df.L_x, markersize=1, alpha=0.05)


	ax = Axis(fig[2,1], ylabel="L")

	scatter!(sim_df.xi_p, sim_df.L_y, markersize=1, alpha=0.05)


	ax = Axis(fig[3,1], xlabel="xi'", ylabel="L")

	scatter!(sim_df.xi_p, sim_df.L_z, markersize=1, alpha=0.05)

	fig
end

# ╔═╡ d4c2bb50-77f3-44ef-8f0b-8126964a2f70
observable_vs_stream(:pmra, ylims=(-1, 4))

# ╔═╡ c0dd1c42-8c6c-42a4-9ec6-316bb5652d4a
observable_vs_stream(:pmdec, ylims=(-5, 5))

# ╔═╡ 05a3d78d-571c-4bcd-ae17-aa7f6202203f
observable_vs_stream(:radial_velocity)

# ╔═╡ f71fd567-59fa-4564-9ee0-1e6bd9d580e8
observable_vs_stream(:distance, ylims=(0, 200))

# ╔═╡ fec39bf2-fa3c-43ca-ba89-0437791e08e5
md"""
# Getting closest particles in phase space
"""

# ╔═╡ e78858a9-0181-4e2e-a0b6-180ea728f8c6
function phase_plot(df, akshara_df)
	fig = Figure(size=(6, 3) .* 72)
	ax = Axis(fig[1,1], xlabel="ra", ylabel="dec")
	scatter!(df.ra, df.dec)
	scatter!(akshara_df.ra, akshara_df.dec)


	ax = Axis(fig[1,2], xlabel="pmra", ylabel="pmdec")
	scatter!(df.pmra, df.pmdec)
	scatter!(akshara_df.pmra, akshara_df.pmdec)


	ax = Axis(fig[1,3], xlabel="distance", ylabel="RV")
	scatter!(df.distance, df.radial_velocity)
	scatter!(akshara_df.distance, akshara_df.radial_velocity)

	fig

end

# ╔═╡ 4831c1d8-4782-4a62-85a4-9785b530789a
r_nearby = 1 # degree

# ╔═╡ a59b71d9-6e7b-4915-83dd-056db44b5f15
function get_coord_dists(sim_df, obs)
	χ2 = @. (sim_df.ra - obs.ra)^2 / r_nearby^2
	χ2 .+= @. (sim_df.dec - obs.dec)^2 / r_nearby^2
	χ2 .+= @. (sim_df.pmra - obs.pmra)^2 / obs.pmra_error^2
	χ2 .+= @. (sim_df.pmdec - obs.pmdec)^2 / obs.pmdec_error^2
	χ2 .+= @. (sim_df.distance - obs.distance)^2 / obs.distance_error^2
	χ2 .+= @. (sim_df.radial_velocity - obs.radial_velocity)^2 / obs.radial_velocity_error^2

	return χ2
end

# ╔═╡ 08a354b6-4305-4938-93af-af9157474ee3
function get_closest_particles(sim_df, obs, N=100)
	return sortperm(get_coord_dists(sim_df, obs))[1:N]
end

# ╔═╡ c5cc4baf-7be3-49c9-a4f1-3955625d4a8a
sort(get_coord_dists(sim_df, akshara_df[1, :]))

# ╔═╡ ab4b578f-37db-4958-859b-807ff0025670
closest_particles = [sim_df[get_closest_particles(sim_df, akshara_df[i, :]), :] for i in 1:size(akshara_df, 1)]

# ╔═╡ 05b8a2d5-eaf2-458b-aef2-e7d1bec1f3ad
Nstream = size(akshara_df, 1)

# ╔═╡ bf89fb53-e16f-4f84-adfd-75f6586f2795
let
	fig = Figure(size=(6, 24) .* 72)

	for i in 1:Nstream
		ax = Axis(fig[i,1], xlabel="ra", ylabel="dec")
		df = closest_particles[i]
		obs_df = akshara_df[i:i, :]
		chi2 = get_coord_dists(df, obs_df)
		
		scatter!(df.ra, df.dec, color=chi2)
		scatter!(obs_df.ra, obs_df.dec, color=:black, markersize=10)
	
	
		ax = Axis(fig[i,2], xlabel="pmra", ylabel="pmdec")
		scatter!(df.pmra, df.pmdec, color=chi2)
		
		errorscatter!(obs_df.pmra, obs_df.pmdec, 
					  xerror = obs_df.pmra_error, 
					  yerror=obs_df.pmdec_error,
				 color=:black, markersize=10)
	
	
		ax = Axis(fig[i,3], xlabel="distance", ylabel="RV")
		scatter!(df.distance, df.radial_velocity, color=chi2)
		errorscatter!(obs_df.distance, obs_df.radial_velocity, 
					  xerror = obs_df.distance_error, 
					  yerror=obs_df.radial_velocity_error,
				 color=:black, markersize=10)
	end

	
	fig
end
	
	

# ╔═╡ a42ef1a5-3305-4ed6-a843-c58b213af617
let 
	fig = Figure()
	ax = Axis(fig[1,1], xlabel="initial distance from LMC", ylabel="density")
	all_closest = vcat(closest_particles...)
		# stephist!(closest_particles[i].r_0, normalization=:pdf)
	# end
	stephist!(all_closest.r_0, bins=100, normalization=:pdf, label="near stream stars")

	stephist!(sim_df.r_0, bins=100, normalization=:pdf, color=:black, label="ALL")

	axislegend()
	fig
end

# ╔═╡ 24c627f0-bab4-409d-b1e0-b6c08ea1dee1
let 
	fig = Figure()
	ax = Axis(fig[1,1], xlabel="binding energy", ylabel="density")
	all_closest = vcat(closest_particles...)
		# stephist!(closest_particles[i].r_0, normalization=:pdf)
	# end
	stephist!(all_closest.energy_0, bins=100, normalization=:pdf, label="near stream stars")

	stephist!(sim_df.energy_0, bins=100, normalization=:pdf, color=:black, label="ALL")

	axislegend()
	fig
end

# ╔═╡ 40261052-1f77-45f4-b492-091db0d39740
let 
	fig = Figure()
	ax = Axis(fig[1,1], xlabel="initial distance from LMC", ylabel="initial angular momentum"
			)

	scatter!(sim_df.energy_0, sim_df.L_0, markersize=1, color=:black, alpha=0.05)

	for i in eachindex(closest_particles)
		scatter!(closest_particles[i].energy_0, closest_particles[i].L_0, markersize=3)
	end


	fig
end

# ╔═╡ eeb6adcf-8672-4053-9ccb-99569216caf2
begin
	function L_max(ϵ)
	    r = r_circ(ϵ)
	    v = v_circ(r)
	
	    return r*v
	end


	function r_circ(ϵ)
	    result = NaN
	    if nextfloat(Φ_lmc(0), 5) <= ϵ
	        result = LilGuys.find_zero(
	        r -> Φ_lmc(r) + 1/2 * M_lmc(r) / r - ϵ,
	        (1e-10, 10000.))
	    end
	    return result
	end


	Φ_lmc(r) = Agama.potential(pot_lmc_init, [1,0,0] * r, Agama.VASILIEV_UNITS)
	M_lmc(r) = Agama.enclosed_mass(pot_lmc_init, r, Agama.VASILIEV_UNITS)

	v_circ(r) = sqrt(M_lmc(r) / r)
end

# ╔═╡ 840629ea-b961-4cf0-9470-48786db2c5c7
L_max(-0.5)

# ╔═╡ a5735a10-72b9-4fc9-83ac-78bb5c6af259
L_max_interp = let
	x = LinRange(Φ_lmc(0)+0.001, -0.001, 1000)
	L = L_max.(x)

	LilGuys.lerp(x, L)
end

# ╔═╡ f9cf4d59-2baf-42bb-91e0-a6b94a62e633
L_max(Φ_lmc(0))

# ╔═╡ 447f6f98-4a25-4314-b320-e9bc4b8d2b3b
let 
	fig = Figure()
	ax = Axis(fig[1,1], xlabel="initial distance from LMC",
			  xscale=log10,
			)
	scatter!(sim_df.energy_0, sim_df.L_0 ./ L_max_interp.(-sim_df.energy_0), markersize=1, color=:black, alpha=0.05)
	
	for i in eachindex(closest_particles)
		scatter!(closest_particles[i].energy_0, closest_particles[i].L_0 ./ L_max_interp.(-closest_particles[i].energy_0), markersize=3)
	end


	fig
end

# ╔═╡ 17b01e69-4bea-4139-b152-fac1b5f80133
closest_particles[1].L_0

# ╔═╡ fb76a263-05ee-4030-add5-f33b5fd27084
let
	fig  = Figure()
	ax = Axis(fig[1,1])

	df = sim_df[sim_df.energy_0 .< 0.4, :]

	scatter!(df.xi_p, df.eta_p, alpha=0.1, color=:black, markersize=1,)
	fig

end

# ╔═╡ 24a8f484-e305-44c2-940d-b9d7de384216
scatter(sim_df.r_0, sim_df.energy_0)

# ╔═╡ 975d5a0a-f8b3-467d-ba37-93bbce5fbfc0
let
	x = LinRange(0, 100, 100)
	v = v_circ.(x)
	phi = Φ_lmc.(x)
	e = @. 1/2*v^2 + phi
	lines(x, -e)

end

# ╔═╡ Cell order:
# ╠═82bc2e30-1d7a-11f1-95f6-49e0b164b4f1
# ╠═1bb94a81-b4e4-4123-8b36-f68145147471
# ╠═2686ab62-1711-4dcb-97b9-e4ab322497d9
# ╠═fce610b5-cfbd-46fe-b570-4d5df7fbf5f0
# ╠═f031667d-e15f-43d7-a832-720e9b73d87c
# ╠═004d2da2-6227-4cd0-ade0-ece085dc06ea
# ╠═936dd128-bfb0-4d5c-b754-8b09a144ff07
# ╠═adebf069-b4e6-4dd5-bfe1-56b70a06e769
# ╠═dc75f240-e1cd-455e-b431-4c3b7ad5343d
# ╠═17e92b56-2d13-4460-a541-80e1e10a30f3
# ╠═4c895cf8-4904-40e5-9778-6149999b1e14
# ╠═cc028c15-ccfe-4c42-b1a7-0e86c36e7e60
# ╠═a215ed7a-f605-4f43-b772-1e13398e6dbf
# ╠═f62020bb-142e-439f-b97d-4612a9ce9ffe
# ╠═a14c6ff9-e8ef-4cd4-b169-38a3110cfa4b
# ╠═fab6087b-3f90-4793-a168-7e7c5e54f06d
# ╠═3dc377fe-698d-4715-98ad-20b633b830d2
# ╠═7ca0c8ee-83ae-4f1f-8c08-b16dc6a4cbc2
# ╠═eb7ca0a8-7968-40a4-86be-6b0d8a45af55
# ╠═96452520-666e-4ae9-8353-e910a9472e75
# ╠═f3fff48e-3a11-4dbe-a837-4c8ae36e804c
# ╠═913fa16e-f951-4caf-9359-823760fa47b0
# ╠═a8b21cd5-e574-43d9-843f-270608a77d2e
# ╠═07f2c1a1-da0d-4576-9938-80bf5850d206
# ╠═e6b87de5-57bb-488c-bf59-6d853da81893
# ╠═e4515a25-08cd-449b-a27d-2be47efe2283
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
# ╠═cd631337-61fe-440b-a38e-1012aa1dd7c7
# ╠═e8b072d4-27e3-4dfc-b2ce-d55565cb2a0d
# ╠═a8ef246a-b342-4e79-b794-a61089a79f81
# ╠═2c205f8f-9321-45dc-90fa-ff03fc4acf60
# ╠═e5a3b299-bac3-4723-8f72-65bf58d21463
# ╠═7edd507b-2ca9-4502-94f6-567b75c2a8ef
# ╠═ee3f0069-31eb-4472-bec3-598c98ef6f6b
# ╠═b1587ebb-4353-4998-a221-97009ffe2a88
# ╠═49f5b518-1dec-4e98-8b49-3d89e704533a
# ╠═447b6b59-ca3d-40ff-9037-272b3ef5b0e1
# ╠═922830a7-736d-4d49-9a83-247d3722ad95
# ╠═9d6d21e9-35ee-4fee-b1bc-2964ec935843
# ╠═d4c2bb50-77f3-44ef-8f0b-8126964a2f70
# ╠═c0dd1c42-8c6c-42a4-9ec6-316bb5652d4a
# ╠═05a3d78d-571c-4bcd-ae17-aa7f6202203f
# ╠═f71fd567-59fa-4564-9ee0-1e6bd9d580e8
# ╠═fec39bf2-fa3c-43ca-ba89-0437791e08e5
# ╠═e78858a9-0181-4e2e-a0b6-180ea728f8c6
# ╠═4831c1d8-4782-4a62-85a4-9785b530789a
# ╠═a59b71d9-6e7b-4915-83dd-056db44b5f15
# ╠═08a354b6-4305-4938-93af-af9157474ee3
# ╠═c5cc4baf-7be3-49c9-a4f1-3955625d4a8a
# ╠═ab4b578f-37db-4958-859b-807ff0025670
# ╠═05b8a2d5-eaf2-458b-aef2-e7d1bec1f3ad
# ╠═bf89fb53-e16f-4f84-adfd-75f6586f2795
# ╠═a42ef1a5-3305-4ed6-a843-c58b213af617
# ╠═24c627f0-bab4-409d-b1e0-b6c08ea1dee1
# ╠═40261052-1f77-45f4-b492-091db0d39740
# ╠═eeb6adcf-8672-4053-9ccb-99569216caf2
# ╠═840629ea-b961-4cf0-9470-48786db2c5c7
# ╠═a5735a10-72b9-4fc9-83ac-78bb5c6af259
# ╠═f9cf4d59-2baf-42bb-91e0-a6b94a62e633
# ╠═447f6f98-4a25-4314-b320-e9bc4b8d2b3b
# ╠═17b01e69-4bea-4139-b152-fac1b5f80133
# ╠═fb76a263-05ee-4030-add5-f33b5fd27084
# ╠═24a8f484-e305-44c2-940d-b9d7de384216
# ╠═975d5a0a-f8b3-467d-ba37-93bbce5fbfc0
