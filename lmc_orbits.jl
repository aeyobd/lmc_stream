### A Pluto.jl notebook ###
# v0.20.23

using Markdown
using InteractiveUtils

# ╔═╡ 8ada3564-1f17-11f1-a1b4-915da9683951
begin
	import Pkg
	Pkg.activate()
	using CairoMakie
	using Arya
	using OrderedCollections
end

# ╔═╡ 2e6379ee-53be-427d-a2fd-ed707713cc72
using LilGuys; FIGDIR = "figures"

# ╔═╡ c66b07d2-728e-4057-955c-d40c0dfc94b9
module Utils
	include("simulations/utils.jl")
end

# ╔═╡ 67b65512-8128-4be7-8d6c-1177d00dea4f
function get_lmc_orbit(name)
	return Utils.get_lmc_orbit(joinpath("simulations", name))

end

# ╔═╡ ecaf3a85-a928-4cf1-bf4e-ccf45418571c
lmc_orbits = OrderedDict(
	"V+21" => get_lmc_orbit("V21"),
	"L2M10first" => get_lmc_orbit("L2M10first"),
	"L2M10" => get_lmc_orbit("L2M10"),
	"L2M11" => get_lmc_orbit("L2M11"),
	"L3M10" => get_lmc_orbit("L3M10"),
	"L3M11" => get_lmc_orbit("L3M11"),
)

# ╔═╡ 70722339-2133-4843-ac54-55ead0c41f42
function plot_r_t!(orbit; kwargs...)
	lines!(orbit.times * T2GYR, radii(orbit); kwargs...)
end

# ╔═╡ 9dca24d5-87b6-449d-91fb-916aaf8316fe
let
	fig = Figure(size=(5, 3) .* 72)
	ax = Axis(fig[1,1],
			 xlabel = "time / Gyr",
			 ylabel = "galcen radius / kpc",)


	for (label, orbit) in lmc_orbits
		plot_r_t!(orbit, label=label)
	end

	Legend(fig[1,2], ax)

	@savefig "lmc_orbits"

	fig
end

# ╔═╡ Cell order:
# ╠═8ada3564-1f17-11f1-a1b4-915da9683951
# ╠═2e6379ee-53be-427d-a2fd-ed707713cc72
# ╠═c66b07d2-728e-4057-955c-d40c0dfc94b9
# ╠═67b65512-8128-4be7-8d6c-1177d00dea4f
# ╠═ecaf3a85-a928-4cf1-bf4e-ccf45418571c
# ╠═70722339-2133-4843-ac54-55ead0c41f42
# ╠═9dca24d5-87b6-449d-91fb-916aaf8316fe
