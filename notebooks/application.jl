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

# ╔═╡ 21d1f0a5-3825-48dc-8264-3584ceb02ca2


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

# ╔═╡ 62f63bb9-03ed-4075-af1e-c9f370862259
dn = let
	# constructor from labeled_topology
	_dn = DecayNode(labeled_topology; apply = flatten_sort_nested_tuple)
	# add indices
	_dn = map_tree(_dn) do v, n
		(; indices = v)
	end
	# add isfirst
	_dn = map_with_parent_node(_dn) do v, p
		_ch = children(p)
		isfirst = true
		if (length(_ch) >= 2)
			x = first(children(_ch))
			isfirst = !(nodevalue(x) != v)
		end
		(; v..., isfirst)
	end
end;

# ╔═╡ 84f16934-39ce-4bb0-a916-a8b677817982
with_angles = 
	map_with_parent_node(dn, DecayNode((; _momenta=Dict(four_momenta_gj)), [dn])) do child_value, parent
		# 
		@unpack indices, isfirst = child_value
		@unpack _momenta = nodevalue(parent)
		#
		this_momenta = Dict(k=>_momenta[k] for k in indices)
		p_this = sum(values(this_momenta))
		# 
		whose_p = isfirst ? p_this : -p_this
		ϕ = azimuthal_angle(whose_p)
		θ = polar_angle(whose_p)
		if pt(p_this) < 1e-10
			ϕ, θ = 0.0, 0.0
		end
		#
		γ = boost_gamma(p_this)
		extra = isfirst ? identity : Ry(-π) 
		#
		this_momenta_transformed = Dict(
			k=>_momenta[k] |> Rz(-ϕ) |> Ry(-θ) |> extra |> Bz(-γ)
				for k in indices)	
		# 
		(; child_value...,
			m=mass(p_this), θ, ϕ,
			_momenta=this_momenta_transformed)
	end;

# ╔═╡ 122264b4-3e96-4562-ac73-9e2c0a6b170f
to_table(with_angles) |> DataFrame

# ╔═╡ fa908cda-ba06-46ad-aa70-01737b1eb951
let
	all_nodes = collect(PreOrderDFS(with_angles))
	all_nodes_with_children = filter(all_nodes) do xi
		length(children(xi)) > 1
	end
	map(all_nodes_with_children) do node
		ch1 = first(children(node))
		ch1_value = nodevalue(ch1)
		value = nodevalue(node)
		(; value.indices, value.m, ch1_value.θ, ch1_value.ϕ)
	end |> DataFrame
end

# ╔═╡ Cell order:
# ╠═40db7730-5ffd-414b-bacd-39f1ebfab4c7
# ╠═21d1f0a5-3825-48dc-8264-3584ceb02ca2
# ╠═7ec9dacc-c703-45ca-80e3-375c1682acac
# ╟─a07d250e-9caa-49f4-ad35-c4946f61906c
# ╠═e2202a7f-607c-4d71-9601-aae272e16d9f
# ╠═1ab67f6d-b937-4513-95e5-e132dce627b1
# ╠═6eafffb6-a602-4e30-9c23-2db2d2119a39
# ╠═2592df4e-fd0a-461f-9bf0-a25e129286df
# ╟─c7c99cc5-0d54-4335-a791-aa0b84963482
# ╠═f63650f3-77d1-4e18-a229-74a70096efef
# ╠═62f63bb9-03ed-4075-af1e-c9f370862259
# ╠═84f16934-39ce-4bb0-a916-a8b677817982
# ╠═122264b4-3e96-4562-ac73-9e2c0a6b170f
# ╠═fa908cda-ba06-46ad-aa70-01737b1eb951
