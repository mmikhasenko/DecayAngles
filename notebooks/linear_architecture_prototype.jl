### A Pluto.jl notebook ###
# v0.19.45

using Markdown
using InteractiveUtils

# ╔═╡ 75b6b4d3-5a64-4f1e-a40d-1b476f8b6b7d
# ╠═╡ show_logs = false
begin
	import Pkg
	Pkg.activate(mktempdir())
	Pkg.add([
		Pkg.PackageSpec(url="https://github.com/mmikhasenko/DecayTreeDataFrames.jl"),
		Pkg.PackageSpec("DataFrames"),
		Pkg.PackageSpec("Parameters")
	])
    using DecayTreeDataFrames
    using DataFrames
    using Parameters
end

# ╔═╡ a04ac952-4289-42a0-a693-e2251473b225
md"""
# Prototype for linear architecture

Idea was to store the tree nodes in an array where the children are registered with array index.
"""

# ╔═╡ b33f4073-434b-4f90-92ad-6d512cdc0913
md"""
## Application
"""

# ╔═╡ b4e8e8ea-b416-4876-84e3-35243487575b
tree_structure = ((:e, (:r, :u)), (:h, :g))

# ╔═╡ 162685ef-ec60-4bcf-aa93-7b1b6f6fa57c
logistics = build_tree(tree_structure)

# ╔═╡ 7d2e372d-2d98-4051-82a5-0a91e9611f04
df = DataFrame(; logistics)

# ╔═╡ a4a94f22-6b93-4241-b898-13ab1451236e
let
    df.notation = getproperty.(df.logistics, :info)
    transform_from_children!(df, :notation) do logistics
        @unpack left, right = logistics
        i1, i2 = df[left, :notation], df[right, :notation]
		"($i1, $i2)"
    end
    df
end

# ╔═╡ 3a630e49-45bd-4260-8bf4-a9db42d90ccf
let
    df.subsystem = getproperty.(df.logistics, :info) .|> vcat
    transform_from_children!(df, :subsystem) do logistics
        @unpack left, right = logistics
        i1, i2 = df[left, :subsystem], df[right, :subsystem]
		vcat(i1, i2)
    end
    df
end

# ╔═╡ 69810976-7abd-45cf-96f3-7a20c248d1b1
let
	df.chain .= "0"
    transform_from_parent!(df, :chain) do logistics
        @unpack parent, me, left, right = logistics
        df[parent, :chain] * " -> " * string(df[me, :logistics].info)
    end
    df
end

# ╔═╡ 97e82d88-ba47-4a9a-8cc7-9a0a589f10eb
md"""
## Lorentz Transformations
"""

# ╔═╡ b73b65f7-e747-4c99-849f-3f5482f59e21
let
	df.transformation .= "I"
    transform_from_parent!(df, :transformation) do logistics
        @unpack parent, me, left, right = logistics
		# 
		l_me = df[me, :logistics].info
		l_pe = df[parent, :logistics].info
		label = l_me * " -> " * l_pe
		# determine transformation: `particle-one` or `particle-two`
		which = df[parent, :logistics].left == me ? "Tf" : "Ts"
		my_transform = which * "⁻¹($label)"
        my_transform * " ∘ " * df[parent, :transformation]
    end
    df
end

# ╔═╡ 870aadb7-2054-493b-8fe1-2c98a454e365
let
    transform!(df, :transformation => ByRow() do t
		t * " [p1,p2, ...]"
	end => :four_vectors)
end

# ╔═╡ 432f56a5-475b-43a5-b3d4-404362ce85d9
let
	df.angles .= ""
	transform_from_children!(df, :angles) do logistics
		@unpack left = logistics
		s = df[left, :subsystem]
		"angles of p̄_$s in [p̄1,p̄2, ...]"
	end
	df
end

# ╔═╡ Cell order:
# ╟─a04ac952-4289-42a0-a693-e2251473b225
# ╠═75b6b4d3-5a64-4f1e-a40d-1b476f8b6b7d
# ╟─b33f4073-434b-4f90-92ad-6d512cdc0913
# ╠═b4e8e8ea-b416-4876-84e3-35243487575b
# ╠═162685ef-ec60-4bcf-aa93-7b1b6f6fa57c
# ╠═7d2e372d-2d98-4051-82a5-0a91e9611f04
# ╠═a4a94f22-6b93-4241-b898-13ab1451236e
# ╠═3a630e49-45bd-4260-8bf4-a9db42d90ccf
# ╠═69810976-7abd-45cf-96f3-7a20c248d1b1
# ╟─97e82d88-ba47-4a9a-8cc7-9a0a589f10eb
# ╠═b73b65f7-e747-4c99-849f-3f5482f59e21
# ╠═870aadb7-2054-493b-8fe1-2c98a454e365
# ╠═432f56a5-475b-43a5-b3d4-404362ce85d9
