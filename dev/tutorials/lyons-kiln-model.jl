### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 6deb7420-f864-11ef-1354-d78734aa776a
begin
    @info("Initializing toolbox...")
    using Pkg

    open("pluto_init.log", "w") do logs
        root = "D:/Kompanion/bin/pkgs/AuChimiste.jl"
        Pkg.activate(root; io=logs)
        Pkg.instantiate(; io=logs)
    end

    # Pkg.add("SciMLBaseMLStyleExt")
    # Pkg.status()
    
    push!(LOAD_PATH, @__DIR__)

    using PlutoLinks
    using PlutoUI: TableOfContents
	import PlutoUI
	
    TableOfContents()
end

# ╔═╡ a2c6f9df-2029-4115-8ae8-a4cfe633bfe1
begin
    @info("Required tools...")

	using CairoMakie
	using DifferentialEquations
	using ModelingToolkit
end

# ╔═╡ 385546eb-0c61-4aab-a2ba-40b691a6c4e4
@independent_variables Θ

# ╔═╡ b6ef75af-66ba-4b7d-b2d1-5eb86faa0f39
D = Differential(Θ)

# ╔═╡ 283fdf3f-4d6a-42c9-b50a-e7ed01820ed8
# Not sure, check: https://en.wikipedia.org/wiki/Gas_constant
RGAS = 1.985875279009

# ╔═╡ 214a4851-a9a2-4a61-a60a-c47b62ff6624
arrhenius(T, A, E) = A * exp(-E / (RGAS * T))

# ╔═╡ 3a8076fc-a09c-4bb8-bd90-448060bd427d
@variables (begin
	L(Θ)
	X(Θ)
	Y(Θ)
	C(Θ)
	S(Θ)
	C₁(Θ)
	W(Θ)
	Gₛ(Θ)
	μₛ(Θ)
	c′(Θ)
	W′(Θ)
	Gᵧ(Θ)
	μᵧ(Θ)
	Tₛ(Θ)
	Tᵧ(Θ)
	Tᵤ(Θ)
	Tᵥ(Θ)
	h₁(Θ)
	h₂(Θ)
	h₃(Θ)
	h₄(Θ)
	h₅(Θ)
end)

# ╔═╡ 17821e33-2e34-4dbf-b6a1-8c1dad2c3900
# Subscripts 
# 1 = gas to wall 
# 2 = gas to solid 
# 3 = inside wall to solid 
# 4 = inside wall to outside wall 
# 5 = outside wall to ambient
@parameters (begin
	Mx
	My
	Mc
	Ms
	Mc′
	
	Aw = 7.080e+07
	Ac = 1.393e+27
	Ax = 5.245e+22
	Ay = 2.780e+08
	
	Ew = 18071.0
	Ec = 261046.0
	Ex = 82644.0
	Ey = 110192.0
	
	ΔHw = 970.0
	ΔHc = 1275.0
	ΔHx = 381.0
	ΔHy = 11.0
	
	ΔHf = 21050.0
	
	Gᵦ
	Gₙ
	B
	N

	RF
	RW
	Tₐ
	ϕ
	Cpₛ
	Cpᵧ
	A₁
	A₂
	A₃
	A₄
	A₅

	# XXX: f's depend on position actually? How?
	f₁
	f₂
	f₃
	f₄ = 3.18
	f₅
	ϵᵧ = 0.382  # 10% of this near discharge?
	ϵₛ = 0.716
	ϵᵤ = 1.0
	ϵᵥ = 0.078
end)

# ╔═╡ cbdb545d-505a-4b00-b110-47755da99cd3
begin
	kw = arrhenius(Aw, Ew, Tₛ)
	kc = arrhenius(Ac, Ec, Tₛ)
	kx = arrhenius(Ax, Ex, Tₛ)
	ky = arrhenius(Ay, Ey, Tₛ)
	
	r_xc = Mx / Mc
	r_yc = My / Mc
	r_sc = Ms / Mc
	r_cc = Mc′ / Mc
	r_βn = Gᵦ / Gₙ
	
	rr_1 = kx * C^2 * S
	rr_2 = ky * C * S
	rr_3 = ky * C * X
	rr_4 = kx * C^2 * S
	rr_5 = ky * C * X
	rr_6 = kc * C₁
	
	Fₛ = 1 + r_cc * C₁ + W
	Fᵧ = 1 + ϕ + W′ + c′

	hA₁ = h₁ * A₁
	hA₂ = h₂ * A₂
	hA₃ = h₃ * A₃
	hA₄ = h₄ * A₄
	hA₅ = h₅ * A₅

	sᵤ = hA₁ + hA₃ + hA₄
	sᵥ = hA₄ + hA₅
	
	qc = ΔHc * B * rr_6
	qw = ΔHw * B * RW
	qx = ΔHx * B * rr_1
	qy = ΔHy * B * rr_3 # Unreadable: g or y?

	qG = B * Cpₛ * (Tₛ - Tᵧ) * D(Gₛ)
	qF = N * ΔHf * RF
	
	q₂ₛ = hA₂ * (Tᵧ - Tₛ)
	q₃ₛ = hA₃ * (Tᵤ - Tₛ)
	
	q₁ᵧ = hA₁ * (Tᵤ - Tᵧ)
	q₂ᵧ = -q₂ₛ
	
	qₛ = q₂ₛ + q₃ₛ - qc - qw + qx + qy
	qᵧ = q₁ᵧ + q₂ᵧ + qG - qF

	σ = 3.4203e-09
	
	eqs = [
		D(L)  ~ Gᵦ / B
		D(X)  ~ r_xc * (rr_1 / 2 - rr_2)
		D(Y)  ~ r_yc * rr_3
		D(C)  ~ -rr_3 - rr_5 + rr_6
		D(S)  ~ r_sc * rr_1  / 2
		D(C₁) ~ -rr_6
		D(W)  ~ -RW
		D(c′) ~ r_βn * r_cc * rr_6
		D(W′) ~ r_βn * RW
		D(Tₛ) ~ qₛ / (Cpₛ * μₛ)
		D(Tᵧ) ~ r_βn * (N / B) * qᵧ / (Cpᵧ * μᵧ)
		Gₛ    ~ Gᵦ * Fₛ
		μₛ    ~ B * Fₛ
		Gᵧ    ~ Gₙ * Fᵧ
		μᵧ    ~ N * Fᵧ
		Tᵤ    ~ (hA₁ * Tᵧ + hA₃ * Tₛ + hA₄ * Tᵥ) / sᵤ
		Tᵥ    ~ (hA₄ * Tᵤ + hA₅ * Tₐ) / sᵥ
		h₁    ~ f₁ + σ * ϵᵧ * ϵᵤ * (Tᵧ^3 + Tᵤ^3)
		h₂    ~ f₂ + σ * ϵᵧ * ϵₛ * (Tᵧ^3 + Tₛ^3)
		h₃    ~ f₃ + σ * ϵᵤ * ϵₛ * (Tᵤ^3 + Tₛ^3) * (A₂ / A₃)
		h₄    ~ f₄
		h₅    ~ f₅ + σ * ϵᵥ * (Tᵤ^3 + Tₐ^3) * 100 # Not Tᵥ?
	]
end

# ╔═╡ Cell order:
# ╟─6deb7420-f864-11ef-1354-d78734aa776a
# ╟─a2c6f9df-2029-4115-8ae8-a4cfe633bfe1
# ╠═385546eb-0c61-4aab-a2ba-40b691a6c4e4
# ╠═b6ef75af-66ba-4b7d-b2d1-5eb86faa0f39
# ╠═283fdf3f-4d6a-42c9-b50a-e7ed01820ed8
# ╠═214a4851-a9a2-4a61-a60a-c47b62ff6624
# ╠═3a8076fc-a09c-4bb8-bd90-448060bd427d
# ╠═17821e33-2e34-4dbf-b6a1-8c1dad2c3900
# ╠═cbdb545d-505a-4b00-b110-47755da99cd3
