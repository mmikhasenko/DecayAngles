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
julia> topology = DecayNode( (1,(2,3)) )
julia> map_tree(topology) do value, node
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

"""
	map_with_parent_node(f, node::DecayNode, parent_node = DecayNode(nothing, [node]))

Propagates the values from parent to all children and apply operation on it.
The signature of the function is f(value, parent_node)

## Example
```julia
julia> map_with_parent(DecayNode((1,(2,3))), "0") do value, parent_value
    string(parent_value) * " => " * string(value)
end
```
"""
function map_with_parent_node(f,
	node::DecayNode,
	parent_node = DecayNode(nothing, [node]))
	#
    new_value = f(nodevalue(node), parent_node)
	new_node = @set node.value = new_value
    new_children = map(children(node)) do child
        map_with_parent_node(f, child, new_node)
    end
    return DecayNode(new_value, new_children)
end

"""
    totable(tree::DecayNode)

Iterate over the tree and table-like collection of named tuples from the node values.

## Example
```julia
julia> topology = DecayNode( (1,(2,3)); apply=x->(; tuple=x))
julia> totable(topology) |> DataFrame
```
"""
function to_table(tree::DecayNode)
    all_nodes = collect(PreOrderDFS(tree))
    return map(all_nodes) do node
        (; nodevalue(node)...)
    end
end
