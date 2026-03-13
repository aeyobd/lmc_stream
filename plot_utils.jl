
X_SUN = [-8.122, 0., 0.]

function integrate_isodensity_2d(pot, initial=[-X_SUN[1], 0.]; s_scale=0.0003, h_scale=0.0001, h0=0.0001, kwargs...)
	x = initial[1]
	y = initial[2]

	xs = [x]
	ys = [y]
	h = h0
    s = h0

	θ = atan(y, x)
	x_vec = [1, 0]

	y_vec = [0, 1]
	ρ(x) = pot._py.projectedDensity(x; kwargs...) |> py2f
	
    ρ_0 = ρ(x_vec * x .+ y_vec * y)
    dlρ_max = 0

    x0 = [x, y]
    dx = (ρ(x0 .+ x_vec * h) - ρ(x0))/h
    dy = (ρ(x0 .+ y_vec * h) - ρ(x0))/h

    s = s_scale * (sqrt(x^2 + y^2) / sqrt(dx^2 + dy^2) )
    h = h_scale * s

	for i in 1:100000
		x0 = [x, y]
		dx = (ρ(x0 .+ x_vec * h) - ρ(x0))/h
		dy = (ρ(x0 .+ y_vec * h) - ρ(x0))/h
		x += s .* dy
		y += -s .* dx

		push!(xs, x)
		push!(ys, y)

        s = s_scale * (sqrt(x^2 + y^2) / sqrt(dx^2 + dy^2) )
        h = h_scale * s

		θ_new = atan(y, x)
        dlρ_max = max(dlρ_max, abs(log10(ρ(x0)) - log10(ρ_0)))

		if θ_new > 0 && (θ < 0)
			break
		end

		θ = θ_new
	end
    @info "max log rel error = $dlρ_max"

    return xs, ys
end
