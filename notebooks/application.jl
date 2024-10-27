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

# ╔═╡ 7ec9dacc-c703-45ca-80e3-375c1682acac
function tuple_if_collection(n)
	leaves = flatten_sort_nested_tuple(n)
	length(leaves)==1 && return leaves[1]
	leaves
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
end

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

# ╔═╡ 90f1e513-008d-47da-98d4-c1f90ed128d0
dn0 = DecayNode(labeled_topology);

# ╔═╡ a4dd6d7d-3c6f-433d-aed4-b2b1ef50bbc5
function add_indices_order(basic_tree)
	# add indices
	_dn = map_tree(basic_tree) do original_info, node
		(; names = original_info, indices = flatten_nested_tuple(original_info))
	end
	# add isfirst
	_dn = map_with_parent_node(_dn) do value, parent
		_ch = children(parent)
		isfirst = true
		if (length(_ch) >= 2)
			first_child = first(_ch)
			isfirst = !(nodevalue(first_child) != value)
		end
		(; value..., isfirst)
	end
end

# ╔═╡ 62f63bb9-03ed-4075-af1e-c9f370862259
dn1 = add_indices_order(dn0);

# ╔═╡ 800775d8-13ca-4938-ae9f-6242cfb1eae6
let
	@assert add_indices_order(DecayNode(((1,2)))) |> to_table |> DataFrame == 
		DataFrame(
			names = [(1,2), 1, 2],
			indices=[(1,2),(1,),(2,)],
			isfirst=[true, true, false])
	@assert add_indices_order(DecayNode(("ah",("bh","cx")))) |> to_table |> DataFrame == 
		DataFrame(
			names = [("ah",("bh","cx")), "ah", ("bh","cx"), "bh", "cx"],
			indices = [("ah","bh","cx"), ("ah",), ("bh", "cx"), ("bh",), ("cx",)],
			isfirst = [true, true, false, true, false]
		)
end

# ╔═╡ e7cca695-7db3-4357-94e7-5840f66e0814
begin
	@with_kw struct AnglesGamma
		θ::Float64
		ϕ::Float64
		γ::Float64
	end
	# 
	abstract type DaughterTransformation end
	struct HelicityTransformation <: DaughterTransformation
		variables::AnglesGamma
	end
	function act(T::HelicityTransformation, p)
		@unpack θ, ϕ, γ = T.variables
		p |> Rz(-ϕ) |> Ry(-θ) |> Bz(-γ)
	end
	function HelicityTransformation(p_this; isfirst=true)
		ϕ = azimuthal_angle(p_this)
		θ = polar_angle(p_this)
		if pt(p_this) < 1e-10
			ϕ, θ = 0.0, 0.0
		end
		γ = boost_gamma(p_this)
		#
		isfirst && return HelicityTransformation(AnglesGamma(; θ, ϕ, γ))
		return HelicityTransformation(AnglesGamma(; θ=π-θ+π, ϕ=mod(ϕ+2π, 2π)-π,γ))
	end
end

# ╔═╡ 51b16d87-1d80-4092-84fe-9b128b42bbae
function add_transform_through(::Type{X}, tree_with_indices_order, momenta::T) where {X <: DaughterTransformation, T}
	fake_parent = DecayNode(
		(; _momenta=T(momenta)), [tree_with_indices_order])
	# the map loops the tree recursively given a function
	filled_tree = map_with_parent_node(tree_with_indices_order, fake_parent) do child_value, parent
		@unpack indices, isfirst = child_value
		@unpack _momenta = nodevalue(parent)
		# using indices compute the daugther four-vector 
		this_momenta = T(k=>_momenta[k] for k in indices)
		p_this = sum(values(this_momenta))
		# make transformation and apply it to all vectors
		h = X(p_this; isfirst)
		this_momenta_transformed = T(k=>act(h, _momenta[k]) for k in indices)	
		#
		(; child_value...,
			m=mass(p_this), h.variables.θ, h.variables.ϕ,
			_momenta=this_momenta_transformed)	
	end
	return filled_tree
end

# ╔═╡ cdc31627-f5a9-4945-a804-4584f9836524
dn2 = add_transform_through(HelicityTransformation, dn1, Dict(four_momenta_gj));

# ╔═╡ ce433384-4277-468a-8da1-b8c5ad558417
dn2 |> to_table |> DataFrame

# ╔═╡ fa908cda-ba06-46ad-aa70-01737b1eb951
function decay_angles(tree_with_momenta)
	all_nodes = collect(PreOrderDFS(tree_with_momenta))
	all_nodes_with_children = filter(all_nodes) do xi
		length(children(xi)) > 1
	end
	map(all_nodes_with_children) do node
		ch1 = first(children(node))
		ch1_value = nodevalue(ch1)
		value = nodevalue(node)
		child1_θ, child1_ϕ = ch1_value.θ, ch1_value.ϕ
		(; value.names, value.m, child1_θ, child1_ϕ)
	end |> DataFrame
end

# ╔═╡ 8f180a10-ca23-40bc-80ef-df2e9d4c60d7
decay_angles(dn2)

# ╔═╡ Cell order:
# ╠═40db7730-5ffd-414b-bacd-39f1ebfab4c7
# ╠═7ec9dacc-c703-45ca-80e3-375c1682acac
# ╟─a07d250e-9caa-49f4-ad35-c4946f61906c
# ╠═e2202a7f-607c-4d71-9601-aae272e16d9f
# ╠═1ab67f6d-b937-4513-95e5-e132dce627b1
# ╠═6eafffb6-a602-4e30-9c23-2db2d2119a39
# ╠═2592df4e-fd0a-461f-9bf0-a25e129286df
# ╟─c7c99cc5-0d54-4335-a791-aa0b84963482
# ╠═f63650f3-77d1-4e18-a229-74a70096efef
# ╠═90f1e513-008d-47da-98d4-c1f90ed128d0
# ╠═62f63bb9-03ed-4075-af1e-c9f370862259
# ╠═cdc31627-f5a9-4945-a804-4584f9836524
# ╠═ce433384-4277-468a-8da1-b8c5ad558417
# ╠═8f180a10-ca23-40bc-80ef-df2e9d4c60d7
# ╠═800775d8-13ca-4938-ae9f-6242cfb1eae6
# ╠═a4dd6d7d-3c6f-433d-aed4-b2b1ef50bbc5
# ╠═e7cca695-7db3-4357-94e7-5840f66e0814
# ╠═51b16d87-1d80-4092-84fe-9b128b42bbae
# ╠═fa908cda-ba06-46ad-aa70-01737b1eb951
