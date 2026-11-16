#!/usr/bin/env -S julia --project=. --threads auto
using Random
using Statistics
using QuadGK
using Plots

"""
    f(x::Float64) -> Float64

Integrando del ejercicio 1: cos(πx / 2).
"""
f(x::Float64)::Float64 = cos(pi * 0.5 * x)

"""
    rho_base(x::Float64) -> Float64

Densidad de probabilidad no normalizada: 1 - 1.2x³.
"""
rho_base(x::Float64)::Float64 = 1.0 - 1.2 * x^3

"""
    _single_uniform(seed::Int, n_samples::Int, limit_inf::Float64, limit_sup::Float64, integrand::Function)

Ejecuta una única realización vectorizada de Monte Carlo uniforme en [limit_inf, limit_sup].
"""
function _single_uniform(
	seed::Int,
	n_samples::Int,
	limit_inf::Float64,
	limit_sup::Float64,
	integrand::Function
)
	rng = Xoshiro(seed)
	ancho = limit_sup - limit_inf
	xs = limit_inf .+ ancho .* rand(rng, n_samples)
	vals = integrand.(xs)
	return (
		integral = ancho * mean(vals),
		std_error = ancho * std(vals) / sqrt(n_samples)
	)
end

"""
    estimacion_uniforme(seed::Int=1, n_samples::Int=1_000_000, repeticiones_experimento::Int=1;
                        limit_inf::Float64=0.0, limit_sup::Float64=0.8, integrand::Function=f)

Estima la integral mediante Monte Carlo uniforme en [limit_inf, limit_sup] a lo largo de
`repeticiones_experimento` realizaciones independientes (semillas `seed, seed+1, ...`).
"""
function estimacion_uniforme(
	seed::Int = 1,
	n_samples::Int = 1_000_000,
	repeticiones_experimento::Int = 1;
	limit_inf::Float64 = 0.0,
	limit_sup::Float64 = 0.8,
	integrand::Function = f
)
	@assert limit_sup > limit_inf && n_samples > 1 && repeticiones_experimento >= 1

	results = [_single_uniform(seed + r - 1, n_samples, limit_inf, limit_sup, integrand)
	           for r in 1:repeticiones_experimento]
	integrals = [r.integral for r in results]
	std_errors = [r.std_error for r in results]

	return (
		integral = mean(integrals),
		std_error = mean(std_errors),
		std_between = repeticiones_experimento > 1 ? std(integrals) : NaN,
		integrals = integrals,
		std_errors = std_errors
	)
end

"""
    _single_non_uniform(seed::Int, n_samples::Int, limit_inf::Float64, limit_sup::Float64,
                        M::Float64, integrand::Function, rho_unnorm::Function)

Ejecuta una única realización vectorizada de Monte Carlo por aceptación-rechazo.
"""
function _single_non_uniform(
	seed::Int,
	n_samples::Int,
	limit_inf::Float64,
	limit_sup::Float64,
	M::Float64,
	integrand::Function,
	rho_unnorm::Function
)
	rng = Xoshiro(seed)
	ancho = limit_sup - limit_inf
	norm_const = 1.0 / quadgk(rho_unnorm, limit_inf, limit_sup)[1]
	techo = M / ancho

	x0 = limit_inf .+ ancho .* rand(rng, n_samples)
	u = rand(rng, n_samples)
	rho_vals = norm_const .* rho_unnorm.(x0)
	@assert all(0.0 .<= rho_vals .<= techo) "Envolvente M insuficiente"

	mask = u .<= (rho_vals ./ techo)
	num_accepted = count(mask)
	@assert num_accepted > 1 "Se requieren al menos 2 muestras aceptadas"

	vals = integrand.(x0[mask]) ./ rho_vals[mask]
	return (
		integral = mean(vals),
		std_error = std(vals) / sqrt(num_accepted),
		accepted = num_accepted,
		efficiency = num_accepted / n_samples
	)
end

"""
    estimacion_no_uniforme(seed::Int=1, n_samples::Int=1_000_000, repeticiones_experimento::Int=1;
                           limit_inf::Float64=0.0, limit_sup::Float64=0.8, M::Float64=1.5,
                           integrand::Function=f, rho_unnorm::Function=rho_base)

Estima la integral mediante muestreo por aceptación-rechazo en [limit_inf, limit_sup] a lo largo de
`repeticiones_experimento` realizaciones independientes (semillas `seed, seed+1, ...`).
"""
function estimacion_no_uniforme(
	seed::Int = 1,
	n_samples::Int = 1_000_000,
	repeticiones_experimento::Int = 1;
	limit_inf::Float64 = 0.0,
	limit_sup::Float64 = 0.8,
	M::Float64 = 1.5,
	integrand::Function = f,
	rho_unnorm::Function = rho_base
)
	@assert limit_sup > limit_inf && n_samples > 1 && repeticiones_experimento >= 1

	results = [_single_non_uniform(seed + r - 1, n_samples, limit_inf, limit_sup, M, integrand, rho_unnorm)
	           for r in 1:repeticiones_experimento]
	integrals = [r.integral for r in results]
	std_errors = [r.std_error for r in results]
	accepted_vec = [r.accepted for r in results]

	return (
		integral = mean(integrals),
		std_error = mean(std_errors),
		std_between = repeticiones_experimento > 1 ? std(integrals) : NaN,
		accepted = round(Int, mean(accepted_vec)),
		efficiency = mean([r.efficiency for r in results]),
		integrals = integrals,
		std_errors = std_errors
	)
end

"""
    main()

Ejecuta los apartados a), b) y c) y guarda las gráficas `apartado_b.png` y `apartado_c.png`.
"""
function main()
	valor_exacto = (2.0 / pi) * sin(pi * 0.5 * 0.8)

	# Apartado a): Estimación con ambos métodos (repeticiones_experimento = 1)
	res_a_unif = estimacion_uniforme(1, 1_000_000, 1)
	res_a_no = estimacion_no_uniforme(1, 1_000_000, 1)
	println("=== Apartado a) (N = 1_000_000, 1 repetición) ===")
	println("Uniforme:            $(res_a_unif.integral) +/- $(res_a_unif.std_error)")
	println("Aceptación-Rechazo:  $(res_a_no.integral) +/- $(res_a_no.std_error) (eficiencia = $(round(res_a_no.efficiency, digits=4)))")
	println("Valor exacto:        $valor_exacto\n")

	# Apartado b): Comparación de precisión con repeticiones_experimento = 50
	K = 50
	res_b_unif = estimacion_uniforme(1, 100_000, K)
	res_b_no = estimacion_no_uniforme(1, 100_000, K)
	mejora_se = res_b_unif.std_error / res_b_no.std_error
	mejora_emp = res_b_unif.std_between / res_b_no.std_between

	println("=== Apartado b) (N = 100_000, $K repeticiones) ===")
	println("Uniforme:            std_error = $(res_b_unif.std_error), std_between = $(res_b_unif.std_between)")
	println("Aceptación-Rechazo:  std_error = $(res_b_no.std_error), std_between = $(res_b_no.std_between)")
	println("Factor de mejora:    $(round(mejora_se, digits=2))x (por error estándar), $(round(mejora_emp, digits=2))x (empírico)\n")

	default(
		fontfamily = "Computer Modern",
		titlefontsize = 11,
		guidefontsize = 10,
		tickfontsize = 9,
		legendfontsize = 8,
		framestyle = :box,
		grid = true,
		gridalpha = 0.25,
		minorgrid = true,
		minorgridalpha = 0.1,
		size = (520, 360),
		dpi = 300
	)
	c_unif = "#0072B2"
	c_no = "#D55E00"

	p_b = histogram(
		res_b_unif.integrals,
		bins = 15,
		normalize = :pdf,
		alpha = 0.45,
		color = c_unif,
		linecolor = c_unif,
		linewidth = 0.8,
		label = "Uniforme",
		xlabel = "Estimación de la integral",
		ylabel = "Densidad de probabilidad",
		title = "Distribucion de estimaciones (N = 10^5, K = $K)",
		legend = :topright
	)
	histogram!(
		p_b,
		res_b_no.integrals,
		bins = 15,
		normalize = :pdf,
		alpha = 0.55,
		color = c_no,
		linecolor = c_no,
		linewidth = 0.8,
		label = "Aceptación-Rechazo"
	)
	vline!(p_b, [valor_exacto], color = :black, linestyle = :dash, linewidth = 1.6, label = "Valor exacto")
	savefig(p_b, "tarea2_apartado_b_figura.pdf")

	# Apartado c): Variación del error con N (escala log-log, pendiente -1/2)
	Ns = [1_000, 3_000, 10_000, 30_000, 100_000, 300_000, 1_000_000]
	err_unif = [estimacion_uniforme(1, n, 1).std_error for n in Ns]
	err_no = [estimacion_no_uniforme(1, n, 1).std_error for n in Ns]
	ref_inv_sqrt = err_unif[1] .* sqrt.(Ns[1] ./ Ns)

	println("=== Apartado c) (Error estándar vs N) ===")
	for (n, eu, en) in zip(Ns, err_unif, err_no)
		println("N = $(lpad(n, 7)):  Uniforme = $(round(eu, sigdigits=4)),  Acep-Rech = $(round(en, sigdigits=4))")
	end

	p_c = plot(
		Ns, err_unif,
		xscale = :log10, yscale = :log10,
		color = c_unif,
		marker = :circle, markersize = 4.5, markerstrokewidth = 0.6,
		linewidth = 1.8,
		label = "Uniforme",
		xlabel = "Número de propuestas N",
		ylabel = "Error estándar",
		title = "Convergencia Monte Carlo",
		legend = :bottomleft
	)
	plot!(
		p_c, Ns, err_no,
		color = c_no,
		marker = :square, markersize = 4.5, markerstrokewidth = 0.6,
		linewidth = 1.8,
		label = "Aceptación-Rechazo"
	)
	plot!(
		p_c, Ns, ref_inv_sqrt,
		linestyle = :dash, color = :black, linewidth = 1.4,
		label = "Referencia ~ N^(-1/2)"
	)
	savefig(p_c, "tarea2_apartado_c_figura.pdf")
	println("\nGráficas guardadas: tarea2_apartado_b_figura.pdf, tarea2_apartado_c_figura.pdf")
end

main()
