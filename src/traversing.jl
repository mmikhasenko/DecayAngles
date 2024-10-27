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

function add_transform_through(::Type{X}, tree_with_indices_order, momenta::T) where {X <: DaughterTransformation, T}
    fake_parent = DecayNode(
        (; _momenta = T(momenta)), [tree_with_indices_order])
    # the map loops the tree recursively given a function
    filled_tree = map_with_parent_node(tree_with_indices_order, fake_parent) do child_value, parent
        @unpack indices, isfirst = child_value
        @unpack _momenta = nodevalue(parent)
        # using indices compute the daugther four-vector
        this_momenta = T(k => _momenta[k] for k in indices)
        p_this = sum(values(this_momenta))
        # make transformation and apply it to all vectors
        h = X(; p = p_this, isfirst)
        this_momenta_transformed = T(k => act(h, _momenta[k]) for k in indices)
        #
        (; child_value...,
            m = mass(p_this), h.variables.θ, h.variables.ϕ,
            _momenta = this_momenta_transformed)
    end
    return filled_tree
end

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
