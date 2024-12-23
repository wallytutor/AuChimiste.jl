### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 3aeadbf4-c14d-11ef-3fd8-09cedaa2b25d
begin
    @info("Initializing toolbox...")
    using Pkg

    open("pluto_init.log", "w") do logs
        Pkg.activate(dirname(@__DIR__); io=logs)
        Pkg.instantiate(; io=logs)
    end

    # If you need to install something else:
    # Pkg.add("Latexify")

    push!(LOAD_PATH, @__DIR__)

    using PlutoLinks
    using PlutoUI: TableOfContents

    TableOfContents()
end

# ╔═╡ 72358599-fbdc-49c5-af59-b69c9ce3d0dd
begin
    @info("Local toolbox...")
    @revise using AuChimiste
end

# ╔═╡ fd0a3399-8589-428c-9036-0d8a57ae6a92
md"""
# AuChimiste Sandbox
"""

# ╔═╡ d1b76b85-5d4f-4c95-9228-9b4527b16155
begin
	A = component(:stoichiometry; Al=2, O=3)
	3quantity(A, 1.0)	
end

# ╔═╡ 7099477c-8ec3-46eb-940f-cc491a6a24d3
let
	A = component(:stoichiometry; Al=2, O=3)
	C = component(:stoichiometry; Ca=1, O=1)

	ma = quantity(A, A.molar_mass)
	mc = quantity(C, C.molar_mass)
	ma + mc
end

# ╔═╡ 3db1f45f-7da5-4281-a9f3-9408b7f64ba0
begin
	ma = quantity(:stoichiometry, 1.0; Al=2, O=3)
	mc = quantity(:stoichiometry, 1.0; Ca=1, O=1)
	ma + mc
end

# ╔═╡ Cell order:
# ╟─fd0a3399-8589-428c-9036-0d8a57ae6a92
# ╟─3aeadbf4-c14d-11ef-3fd8-09cedaa2b25d
# ╟─72358599-fbdc-49c5-af59-b69c9ce3d0dd
# ╠═d1b76b85-5d4f-4c95-9228-9b4527b16155
# ╠═7099477c-8ec3-46eb-940f-cc491a6a24d3
# ╠═3db1f45f-7da5-4281-a9f3-9408b7f64ba0
