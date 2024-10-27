### A Pluto.jl notebook ###
# v0.19.45

using Markdown
using InteractiveUtils

# ╔═╡ 40db7730-5ffd-414b-bacd-39f1ebfab4c7
# ╠═╡ show_logs = false
begin
	using Pkg
	Pkg.activate(mktempdir())
	Pkg.add([
		Pkg.PackageSpec(url="https://github.com/mmikhasenko/DecayAngles.jl"),
		Pkg.PackageSpec(url="https://github.com/JuliaHEP/LorentzVectorBase.jl"),
		Pkg.PackageSpec(url="https://github.com/mmikhasenko/FourVectors.jl"),
		Pkg.PackageSpec("DataFrames"),
		Pkg.PackageSpec("LinearAlgebra"),
		Pkg.PackageSpec("Parameters"),
	])
	using DecayAngles
	using DecayAngles.AbstractTrees
	using LinearAlgebra
	using DataFrames
	using FourVectors
	using Parameters
end

# ╔═╡ a07d250e-9caa-49f4-ad35-c4946f61906c
md"""
## Test with real-world example
"""

# ╔═╡ e2202a7f-607c-4d71-9601-aae272e16d9f
begin
	const masses = Dict(
		"pi" => 0.13957018,
		"eta" => 0.547862,
		"proton" => 0.938272046
	)
	const labels = ["pi+", "pi-_1", "pi-_2", "eta"];
	const all_labels = ["pi+", "pi-_1", "pi-_2", "eta", "pi-_0", "proton_t"];
	const masses_id = Dict(all_labels .=> ["pi", "pi", "pi", "eta", "pi", "proton"]);
end;

# ╔═╡ 1ab67f6d-b937-4513-95e5-e132dce627b1
momenta = Dict(
    "pi-_1" => [ -0.24818327141877264, -0.20112688257643338, 27.190349612337307 ],
    "pi+" => [ 0.5604238188629799, 0.28189575586139787, 30.538753073867316 ],
    "pi-_2"=> [ 0.12005509350839914, -0.0885967539890001, 35.286454346616935 ],
    "eta" => [ -0.7657807716990657, -0.012833700231363973, 97.08805230456898 ],
	"pi-_0" => [0.00558236529759858, 0.005463061735694731, 190.63513271129372],
	"proton_t" => [0.0, 0.0, 0.0]
);

# ╔═╡ 6eafffb6-a602-4e30-9c23-2db2d2119a39
four_momenta = [
	k => FourVector(v...; M = masses[masses_id[k]])
		for (k, v) in momenta];

# ╔═╡ 2592df4e-fd0a-461f-9bf0-a25e129286df
four_momenta_gj = let
	_dmomenta = Dict(four_momenta)
	pR = sum(_dmomenta[l] for l in labels)
	# in cms
	pv′ = transform_to_cmf.(last.(four_momenta), pR |> Ref)
	# 
	_dmomenta′ = Dict(first.(four_momenta) .=> pv′)
	which_z = _dmomenta′["pi-_0"]
	which_xplus = -_dmomenta′["proton_t"]
	pv′′ = rotate_to_plane.(pv′, which_z |> Ref, which_xplus |> Ref)
	# 
	first.(four_momenta) .=> pv′′
end;

# ╔═╡ c7c99cc5-0d54-4335-a791-aa0b84963482
md"""
## Test with four-vectors
"""

# ╔═╡ f63650f3-77d1-4e18-a229-74a70096efef
labeled_topology = ((("pi-_1", "eta"), "pi+"), "pi-_2")

# ╔═╡ 62f63bb9-03ed-4075-af1e-c9f370862259
begin
	dn0 = DecayNode(labeled_topology);
	dn1 = add_indices_order(dn0);
	dn2 = add_transform_through(HelicityTransformation, dn1, Dict(four_momenta_gj));
end

# ╔═╡ ce433384-4277-468a-8da1-b8c5ad558417
dn2 |> to_table |> DataFrame

# ╔═╡ 8f180a10-ca23-40bc-80ef-df2e9d4c60d7
decay_angles(dn2) |> DataFrame

# ╔═╡ Cell order:
# ╠═40db7730-5ffd-414b-bacd-39f1ebfab4c7
# ╟─a07d250e-9caa-49f4-ad35-c4946f61906c
# ╠═e2202a7f-607c-4d71-9601-aae272e16d9f
# ╠═1ab67f6d-b937-4513-95e5-e132dce627b1
# ╠═6eafffb6-a602-4e30-9c23-2db2d2119a39
# ╠═2592df4e-fd0a-461f-9bf0-a25e129286df
# ╟─c7c99cc5-0d54-4335-a791-aa0b84963482
# ╠═f63650f3-77d1-4e18-a229-74a70096efef
# ╠═62f63bb9-03ed-4075-af1e-c9f370862259
# ╠═ce433384-4277-468a-8da1-b8c5ad558417
# ╠═8f180a10-ca23-40bc-80ef-df2e9d4c60d7
