using Random

"""
    f(x::Float64) -> Float64

Integrando del ejercicio 1: cos(πx / 2).
"""
f(x::Float64)::Float64 = cos(pi * 0.5 * x)

"""
    welford_step(mean_val::Float64, m2::Float64, x::Float64, count::Int) -> Tuple{Float64, Float64}

Actualiza media y suma de cuadrados centrados (algoritmo de Welford) en una sola pasada.
"""
@inline function welford_step(mean_val::Float64, m2::Float64, x::Float64, count::Int)
	delta = x - mean_val
	mean_val += delta / count
	m2 += delta * (x - mean_val)
	return mean_val, m2
end

"""
    uno_ai(; seed::Int=1, n_samples::Int=1_000_000, limit_inf::Float64=0.0, limit_sup::Float64=0.8)

Ejercicio 1a(i): Estima ∫ f(x) dx mediante Monte Carlo crudo con muestreo uniforme en [limit_inf, limit_sup].
"""
function uno_ai(
	seed::Int = 1,
	n_samples::Int = 1_000_000,
	limit_inf::Float64 = 0.0,
	limit_sup::Float64 = 0.8
)
	@assert limit_sup > limit_inf && n_samples > 1
	rng = Xoshiro(seed)
	ancho = limit_sup - limit_inf

	mean_f = 0.0
	m2 = 0.0

	for i in 1:n_samples
		x_i = limit_inf + ancho * rand(rng)
		mean_f, m2 = welford_step(mean_f, m2, f(x_i), i)
	end

	int_uniform = ancho * mean_f
	std_error = ancho * sqrt((m2 / (n_samples - 1)) / n_samples)

	println("Ejercicio 1a(i)")
	println("Valor de la integral (Uniforme): $int_uniform ± $std_error")
	return (integral = int_uniform, std_error = std_error)
end

"""
    uno_aii(; seed::Int=1, n_samples::Int=1_000_000, limit_inf::Float64=0.0, limit_sup::Float64=0.8)

Ejercicio 1a(ii): Genera muestras de ρ(x) ∝ 1 - 1.2x³ en [0.0, 0.8] mediante aceptación-rechazo
con propuesta uniforme sobre el mismo soporte y estima ∫ f(x) dx = E_ρ[f(X)/ρ(X)].
"""
function uno_aii(
	seed::Int = 1,
	n_samples::Int = 1_000_000,
	limit_inf::Float64 = 0.0,
	limit_sup::Float64 = 0.8
)
	@assert limit_sup > limit_inf && n_samples > 1
	rng = Xoshiro(seed)
	ancho = limit_sup - limit_inf

	# Constante de normalización exacta en [0, 0.8]: ∫(1 - 1.2x³)dx = 0.67712
	A = 1.0 / 0.67712
	rho_nu(x) = A * (1.0 - 1.2 * x^3)

	# Propuesta uniforme g(x) = 1 / ancho; M_opt = ancho * A ≈ 1.1815
	g_val = 1.0 / ancho
	M = 1.2
	techo = M * g_val

	num_samples_accepted = 0
	mean_val = 0.0
	m2 = 0.0

	for _ in 1:n_samples
		x0 = limit_inf + ancho * rand(rng)
		rho_val = rho_nu(x0)
		@assert 0.0 <= rho_val <= techo "Envolvente M insuficiente en x = $x0"

		if rand(rng) <= (rho_val / techo)
			num_samples_accepted += 1
			mean_val, m2 = welford_step(mean_val, m2, f(x0) / rho_val, num_samples_accepted)
		end
	end

	@assert num_samples_accepted > 1 "Se requieren al menos 2 muestras aceptadas"
	std_error = sqrt((m2 / (num_samples_accepted - 1)) / num_samples_accepted)
	efficiency = num_samples_accepted / n_samples

	println("Ejercicio 1a(ii)")
	println("Valor de la integral (Aceptación-Rechazo): $mean_val ± $std_error")
	println("Número de muestras aceptadas: $num_samples_accepted")
	println("Eficiencia: $efficiency")
	return (
		integral = mean_val,
		std_error = std_error,
		accepted = num_samples_accepted,
		efficiency = efficiency
	)
end

res_i = uno_ai()
println("----------------------------------------")
res_ii = uno_aii()
println("----------------------------------------")
valor_exacto = (2.0 / pi) * sin(pi * 0.5 * 0.8)
println("Valor analítico exacto: $valor_exacto")
