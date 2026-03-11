### A Pluto.jl notebook ###
# v0.20.23

using Markdown
using InteractiveUtils

# ╔═╡ 647c0e28-1d94-11f1-912f-0f8c68233ab1
begin
	import Pkg; Pkg.activate()

	using WGLMakie
	using LilGuys
end

# ╔═╡ fe8591ed-47d7-4c4f-a77c-05b4b7f11f9a
simname = "L2M10"

# ╔═╡ 50b04e3a-f834-44c9-823b-3f91e278d99a
snap = Snapshot(joinpath("simulations", simname, "final_positions.hdf5"))

# ╔═╡ ac14c72c-11b1-4e0b-b75e-43d7bf8e654b
let
	scatter(snap.positions)
end

# ╔═╡ Cell order:
# ╠═647c0e28-1d94-11f1-912f-0f8c68233ab1
# ╠═fe8591ed-47d7-4c4f-a77c-05b4b7f11f9a
# ╠═50b04e3a-f834-44c9-823b-3f91e278d99a
# ╠═ac14c72c-11b1-4e0b-b75e-43d7bf8e654b
