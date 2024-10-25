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
		Pkg.PackageSpec("Setfield"),
		Pkg.PackageSpec("LinearAlgebra"),
		Pkg.PackageSpec("Parameters"),
	])
	using DecayAngles
	using DecayAngles.AbstractTrees
	using Setfield
	using LinearAlgebra
	using DataFrames
	using FourVectors
	using Parameters
end

# ╔═╡ ce97968a-0a5a-4908-b87b-410a51639cb8
let
	topology = DecayNode( (1,(2,3)); apply=x->(; tuple=x))
	to_table(topology)
end |> DataFrame

# ╔═╡ 7ec9dacc-c703-45ca-80e3-375c1682acac
function tuple_if_collection(n)
	leaves = flatten_sort_nested_tuple(n)
	length(leaves)==1 && return leaves[1]
	leaves
end

# ╔═╡ 496da1e2-0f8a-419f-afab-9130504483d6
flatten_sort_nested_tuple

# ╔═╡ eb0f7f93-3d43-445c-80e9-f437a799441c
md"""
## Test
"""

# ╔═╡ 935df519-3061-4e43-9105-6d72144a8b8c
topology = (((4,2),3),1)

# ╔═╡ 10e771dc-6ac7-4c56-a678-a02b92d07b31
tupled_tree = DecayNode(topology; apply = tuple_if_collection)

# ╔═╡ ccbcdcd3-db74-42af-92e4-40f12226f51b
indexed_tree = let
	index = 0
	t = map_tree(tupled_tree) do v,n
		index += 1
		braket = v
		(; braket, index)
	end
end

# ╔═╡ 129f69da-460f-46b0-a118-ee03f69065b6
to_table(indexed_tree) |> DataFrame

# ╔═╡ a07d250e-9caa-49f4-ad35-c4946f61906c
md"""
## Test with real-world example
"""

# ╔═╡ e2202a7f-607c-4d71-9601-aae272e16d9f
begin
	const masses = Dict(
		"pi" => 0.13957018,
		"eta" => 0.547862,
	)
	const labels = ["pi+", "pi-_1", "pi-_2", "eta"];
	const masses_id = Dict(labels .=> ["pi", "pi", "pi", "eta"]);
	const names = Dict(1:length(labels) .=> labels)
end

# ╔═╡ 1ab67f6d-b937-4513-95e5-e132dce627b1
momenta = Dict(
    "pi-_1" => [ -0.24818327141877264, -0.20112688257643338, 27.190349612337307 ],
    "pi+" => [ 0.5604238188629799, 0.28189575586139787, 30.538753073867316],
    "pi-_2"=> [ 0.12005509350839914, -0.0885967539890001, 35.286454346616935 ],
    "eta" => [ -0.7657807716990657, -0.012833700231363973, 97.08805230456898 ]
);

# ╔═╡ 6eafffb6-a602-4e30-9c23-2db2d2119a39
four_momenta = Dict(
	k => FourVector(v...; M = masses[masses_id[k]]) for (k, v) in momenta);

# ╔═╡ c7c99cc5-0d54-4335-a791-aa0b84963482
md"""
## Test with four-vectors
"""

# ╔═╡ 62f63bb9-03ed-4075-af1e-c9f370862259
dn = let
	labeled_topology = ((("pi-_1", "eta"),"pi+"),"pi-_2")
	_dn = DecayNode(labeled_topology; apply = flatten_sort_nested_tuple)
	map_tree(_dn) do v, n
		isfirst = true
		(; isfirst, indices = v)
	end
end

# ╔═╡ 84f16934-39ce-4bb0-a916-a8b677817982
with_angles = 
	map_with_parent(dn, (; _momenta=four_momenta)) do child_value, parent_value
		# 
		@unpack indices = child_value
		@unpack _momenta = parent_value
		#
		this_momenta = Dict(k=>_momenta[k] for k in indices)
		p_this = sum(values(this_momenta))
		@unpack cosθ, ϕ = spherical_coordinates(p_this)
		θ = acos(cosθ)
		γ = boost_gamma(p_this)
		#
		this_momenta_transformed = Dict(
			k=>_momenta[k] |> Rz(-ϕ) |> Ry(-θ) |> Bz(-γ)
				for k in indices)	
		# 
		(; child_value..., m=mass(p_this), θ = acos(cosθ), ϕ, _momenta=this_momenta_transformed, )
	end

# ╔═╡ 122264b4-3e96-4562-ac73-9e2c0a6b170f
to_table(with_angles) |> DataFrame

# ╔═╡ Cell order:
# ╠═40db7730-5ffd-414b-bacd-39f1ebfab4c7
# ╠═ce97968a-0a5a-4908-b87b-410a51639cb8
# ╠═7ec9dacc-c703-45ca-80e3-375c1682acac
# ╠═496da1e2-0f8a-419f-afab-9130504483d6
# ╟─eb0f7f93-3d43-445c-80e9-f437a799441c
# ╠═935df519-3061-4e43-9105-6d72144a8b8c
# ╠═10e771dc-6ac7-4c56-a678-a02b92d07b31
# ╠═ccbcdcd3-db74-42af-92e4-40f12226f51b
# ╠═129f69da-460f-46b0-a118-ee03f69065b6
# ╟─a07d250e-9caa-49f4-ad35-c4946f61906c
# ╠═e2202a7f-607c-4d71-9601-aae272e16d9f
# ╠═1ab67f6d-b937-4513-95e5-e132dce627b1
# ╠═6eafffb6-a602-4e30-9c23-2db2d2119a39
# ╟─c7c99cc5-0d54-4335-a791-aa0b84963482
# ╠═62f63bb9-03ed-4075-af1e-c9f370862259
# ╠═84f16934-39ce-4bb0-a916-a8b677817982
# ╠═122264b4-3e96-4562-ac73-9e2c0a6b170f
