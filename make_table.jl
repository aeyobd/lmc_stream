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

# ╔═╡ c26401e9-87f3-4570-a65a-834f14c3de00
using PyFITS

# ╔═╡ 1bb94a81-b4e4-4123-8b36-f68145147471
import CSV

# ╔═╡ 2686ab62-1711-4dcb-97b9-e4ab322497d9
import Agama

# ╔═╡ fce610b5-cfbd-46fe-b570-4d5df7fbf5f0
import DataFrames: DataFrame, rename

# ╔═╡ 004d2da2-6227-4cd0-ade0-ece085dc06ea
CairoMakie.activate!(type=:png)

# ╔═╡ 936dd128-bfb0-4d5c-b754-8b09a144ff07
simname = "L2M11"

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
ϵ = LilGuys.specific_energy(snap_i) * V2KMS^2

# ╔═╡ f62020bb-142e-439f-b97d-4612a9ce9ffe
L = LilGuys.angular_momenta(snap_i) * V2KMS 

# ╔═╡ a14c6ff9-e8ef-4cd4-b169-38a3110cfa4b
module Utils
	include("./simulations/calc_orbits.jl")
end

# ╔═╡ fab6087b-3f90-4793-a168-7e7c5e54f06d
@assert snap.index == snap_i.index

# ╔═╡ 3dc377fe-698d-4715-98ad-20b633b830d2
orbit_lmc = Utils.get_lmc_orbit(joinpath("simulations", simname))

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

# ╔═╡ e5a3b299-bac3-4723-8f72-65bf58d21463
sim_df = let
	df = LilGuys.to_frame(coords_icrs)
	df_gsr = LilGuys.to_frame(coords_gsr)
	
	df[!, :index] = snap.index

	df[!, :pmra_gsr] = df_gsr.pmra
	df[!, :pmdec_gsr] = df_gsr.pmdec
	df[!, :radial_velocity_gsr] = df_gsr.radial_velocity

	df[!, :r_lmc_i] = r_i
	df[!, :energy_lmc_i] = ϵ
	df[!, :L_x_lmc_i] = L[1, :]
	df[!, :L_y_lmc_i] = L[2, :]
	df[!, :L_z_lmc_i] = L[3, :]
	df[!, :L_lmc_i] = radii(L)

	df[!, :x_lmc_i] = snap_i.positions[1, :]
	df[!, :y_lmc_i] = snap_i.positions[2, :]
	df[!, :z_lmc_i] = snap_i.positions[3, :]
	df[!, :v_x_lmc_i] = snap_i.velocities[1, :] * V2KMS
	df[!, :v_y_lmc_i] = snap_i.velocities[2, :] * V2KMS
	df[!, :v_z_lmc_i] = snap_i.velocities[3, :] * V2KMS

	df
end

# ╔═╡ 56ce33d6-8a1f-4606-b4b3-9ade161e4ae6
write_fits(joinpath("simulations", simname, "icrs_final_coords.fits"), sim_df)

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
# ╠═913fa16e-f951-4caf-9359-823760fa47b0
# ╠═a8b21cd5-e574-43d9-843f-270608a77d2e
# ╠═07f2c1a1-da0d-4576-9938-80bf5850d206
# ╠═e6b87de5-57bb-488c-bf59-6d853da81893
# ╠═e4515a25-08cd-449b-a27d-2be47efe2283
# ╠═5817a668-73a0-488e-88de-08c62efc5a65
# ╠═e5a3b299-bac3-4723-8f72-65bf58d21463
# ╠═c26401e9-87f3-4570-a65a-834f14c3de00
# ╠═56ce33d6-8a1f-4606-b4b3-9ade161e4ae6
