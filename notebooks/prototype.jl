### A Pluto.jl notebook ###
# v0.19.45

using Markdown
using InteractiveUtils

# ╔═╡ 40db7730-5ffd-414b-bacd-39f1ebfab4c7
begin
	using AbstractTrees
	using Setfield
	using LinearAlgebra
	using DataFrames
	using Tables
end

# ╔═╡ 414c1fc8-73e9-4afd-93f8-f25c96e34dc5


# ╔═╡ 88491cba-acae-44f9-8828-536277e7422a


# ╔═╡ d66a0a22-5c61-4bc0-833d-f3d2c2ab478c


# ╔═╡ 3a651de8-9222-11ef-245d-016a66f04cd7
begin
	struct DecayNode{T}
	    value::T
	    children::Vector{}
	end
	
	AbstractTrees.children(node::DecayNode) = node.children
	AbstractTrees.nodevalue(node::DecayNode) = node.value
	
	"""
		DecayNode(node::Tuple, f0=identity)

	Wrap a bracket topology notation into a DecayNode structure

	## Example
	```julia
	julia> map_tree(DecayNode((1,(2,3)))) do value, node
		tuple_representation = value
		(; name=string(value), tuple_representation)
	end
	```
	"""
	function DecayNode(node; apply=identity)
	    new_value = apply(node)
	    new_children = map(child -> DecayNode(child; apply), children(node))
	    return DecayNode(new_value, collect(new_children))
	end
	#
	"""
		map_tree(f, node::DecayNode)

	Iterate over the three and apply a function f(value, node) to every nodes, storing the result in the value field.

	## Example
	```julia
	julia> map_tree(build_tree((1,(2,3)))) do value, node
		tuple_representation = value
		(; name=string(value), tuple_representation)
	end
	```
	"""
	function map_tree(f, node::DecayNode)
	    new_value = f(nodevalue(node), node)
	    new_children = map(child -> map_tree(f, child), node.children)
	    return DecayNode(new_value, new_children)
	end
	# 
	"""
		map_with_parent(f, node::DecayNode)

	Porpagates the values from parent to all children and apply operation on it.
	The signature of the function is f(value, parent_value)

	## Example
	```julia
	julia> map_with_parent(build_tree((1,(2,3)))) do value, parent_value
		parent_value
		(; name=string(value), tuple_representation)
	end
	```
	"""	
	function map_with_parent(f, node::DecayNode, parent_value=nothing)
	    new_value = f(nodevalue(node), parent_value)
	    new_children = map(children(node)) do child
			map_with_parent(f, child, new_value)
		end
	    return DecayNode(new_value, new_children)
	end
	# 
	function totable(tree::DecayNode)
	    all_nodes = collect(PreOrderDFS(tree))
	    return map(all_nodes) do node
			(; node.value...)
		end
	end
end

# ╔═╡ ce97968a-0a5a-4908-b87b-410a51639cb8
let
	topology = DecayNode( (1,(2,3)); apply=x->(; tuple=x))
	totable(topology)
end |> DataFrame

# ╔═╡ 7ec9dacc-c703-45ca-80e3-375c1682acac
function tuple_if_collection(n)
	leaves = Tuple(sort(collect(Leaves(n))))
	length(leaves)==1 && return leaves[1]
	leaves
end

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
totable(indexed_tree) |> DataFrame

# ╔═╡ a07d250e-9caa-49f4-ad35-c4946f61906c
md"""
## Test with real-world example
"""

# ╔═╡ e2202a7f-607c-4d71-9601-aae272e16d9f
begin
	const labels = ["pi+", "pi-_1", "pi-_2", "eta"];
	const masses_id = Dict(labels .=> ["pi", "pi", "pi", "eta"]);
	const names = Dict(1:length(labels) .=> labels)
end

# ╔═╡ bc9f5338-78f5-41cf-b163-e0e655c48760
function chain_name(tupled_tree, labels)
	# replace indices with labels
	named_tree = map_tree(tupled_tree) do v,n
		string_dict = Dict(string.(1:length(labels)) .=> labels)
		name = string(v)
		replace(name, string_dict..., ","=>"")
	end
	# zero decay particle has name X
	new_named_tree = @set named_tree.value = "X"
	# 
	formed_tree = map_tree(new_named_tree) do v,n
		ch_names = getproperty.(children(n), :value)
		v*"->"*join(ch_names, ",")
	end
	# combine strings
	leaves = collect(PreOrderDFS(formed_tree))
	filter!(leaves) do n
		length(children(n)) != 0
	end
	list_of_names = getproperty.(leaves, :value)
	# join all with ;
	join(list_of_names, ";")
end	

# ╔═╡ 7cd365cb-ee40-49dc-b6c1-c5d75d3f73eb
chain_name(tupled_tree, labels)

# ╔═╡ c7c99cc5-0d54-4335-a791-aa0b84963482
md"""
## Test with four-vectors
"""

# ╔═╡ 1ab67f6d-b937-4513-95e5-e132dce627b1
momenta = Dict(
    "pi-_1" => [ -0.24818327141877264, -0.20112688257643338, 27.190349612337307 ],
    "pi+" => [ 0.5604238188629799, 0.28189575586139787, 30.538753073867316],
    "pi-_2"=> [ 0.12005509350839914, -0.0885967539890001, 35.286454346616935 ],
    "eta" => [ -0.7657807716990657, -0.012833700231363973, 97.08805230456898 ]
)

# ╔═╡ 84f16934-39ce-4bb0-a916-a8b677817982
map_with_parent(tupled_tree, momenta) do child_indices, _momenta
	@show child_indices
	Dict([names[i]=>_momenta[names[i]] for i in child_indices])
end

# ╔═╡ 4a9593f9-9762-40c9-b3e9-fb5d91930c58
test_tree = map_tree(tupled_tree) do v,n
	total = sum(v) do ind
		name = names[ind]
		momenta[name]
	end
	(; node_name = v, p_norm = norm(total))
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
AbstractTrees = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Setfield = "efcf1570-3423-57d1-acb7-fd33fddbac46"
Tables = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"

[compat]
AbstractTrees = "~0.4.5"
DataFrames = "~1.7.0"
Setfield = "~1.1.1"
Tables = "~1.12.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.1"
manifest_format = "2.0"
project_hash = "7b3d89e8d7e1f8d5b3d7480cafbd433274dcb1a4"

[[deps.AbstractTrees]]
git-tree-sha1 = "2d9c9a55f9c93e8887ad391fbae72f8ef55e1177"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.5"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "8ae8d32e09f0dcf42a36b90d4e17f5dd2e4c4215"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.16.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.ConstructionBase]]
git-tree-sha1 = "76219f1ed5771adbb096743bff43fb5fdd4c1157"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.8"

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseLinearAlgebraExt = "LinearAlgebra"
    ConstructionBaseStaticArraysExt = "StaticArrays"

    [deps.ConstructionBase.weakdeps]
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "DataStructures", "Future", "InlineStrings", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrecompileTools", "PrettyTables", "Printf", "Random", "Reexport", "SentinelArrays", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "fb61b4812c49343d7ef0b533ba982c46021938a6"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.7.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "1d0a14036acb104d9e89698bd408f63ab58cdc82"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.20"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"
version = "1.11.0"

[[deps.InlineStrings]]
git-tree-sha1 = "45521d31238e87ee9f9732561bfee12d4eebd52d"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.4.2"

    [deps.InlineStrings.extensions]
    ArrowTypesExt = "ArrowTypes"
    ParsersExt = "Parsers"

    [deps.InlineStrings.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"
    Parsers = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.InvertedIndices]]
git-tree-sha1 = "0dc7b50b8d436461be01300fd8cd45aa0274b038"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.3.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "2fa9ee3e63fd3a4f7a9a4f4744a52f4856de82df"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.13"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.OrderedCollections]]
git-tree-sha1 = "dfdf5519f235516220579f949664f1bf44e741c5"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.3"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "36d8b4b899628fb92c2749eb488d884a926614d3"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.3"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.PrettyTables]]
deps = ["Crayons", "LaTeXStrings", "Markdown", "PrecompileTools", "Printf", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "1101cd475833706e4d0e7b122218257178f48f34"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "2.4.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "ff11acffdb082493657550959d4feb4b6149e73a"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.4.5"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "66e0a8e672a0bdfca2c3f5937efb8538b9ddc085"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.1"

[[deps.StaticArraysCore]]
git-tree-sha1 = "192954ef1208c7019899fbf8049e717f92959682"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.3"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

    [deps.Statistics.weakdeps]
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.StringManipulation]]
deps = ["PrecompileTools"]
git-tree-sha1 = "a6b1675a536c5ad1a60e5a5153e1fee12eb146e3"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.4.0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "598cd7c1f68d1e205689b1c2fe65a9f85846f297"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.12.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"
"""

# ╔═╡ Cell order:
# ╠═40db7730-5ffd-414b-bacd-39f1ebfab4c7
# ╠═414c1fc8-73e9-4afd-93f8-f25c96e34dc5
# ╠═ce97968a-0a5a-4908-b87b-410a51639cb8
# ╠═88491cba-acae-44f9-8828-536277e7422a
# ╠═d66a0a22-5c61-4bc0-833d-f3d2c2ab478c
# ╠═3a651de8-9222-11ef-245d-016a66f04cd7
# ╠═7ec9dacc-c703-45ca-80e3-375c1682acac
# ╟─eb0f7f93-3d43-445c-80e9-f437a799441c
# ╠═935df519-3061-4e43-9105-6d72144a8b8c
# ╠═10e771dc-6ac7-4c56-a678-a02b92d07b31
# ╠═ccbcdcd3-db74-42af-92e4-40f12226f51b
# ╠═129f69da-460f-46b0-a118-ee03f69065b6
# ╟─a07d250e-9caa-49f4-ad35-c4946f61906c
# ╠═e2202a7f-607c-4d71-9601-aae272e16d9f
# ╠═bc9f5338-78f5-41cf-b163-e0e655c48760
# ╠═7cd365cb-ee40-49dc-b6c1-c5d75d3f73eb
# ╟─c7c99cc5-0d54-4335-a791-aa0b84963482
# ╠═1ab67f6d-b937-4513-95e5-e132dce627b1
# ╠═84f16934-39ce-4bb0-a916-a8b677817982
# ╠═4a9593f9-9762-40c9-b3e9-fb5d91930c58
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
