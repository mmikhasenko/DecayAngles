var documenterSearchIndex = {"docs":
[{"location":"95-reference/#reference","page":"Reference","title":"Reference","text":"","category":"section"},{"location":"95-reference/#Contents","page":"Reference","title":"Contents","text":"","category":"section"},{"location":"95-reference/","page":"Reference","title":"Reference","text":"Pages = [\"95-reference.md\"]","category":"page"},{"location":"95-reference/#Index","page":"Reference","title":"Index","text":"","category":"section"},{"location":"95-reference/","page":"Reference","title":"Reference","text":"Pages = [\"95-reference.md\"]","category":"page"},{"location":"95-reference/","page":"Reference","title":"Reference","text":"Modules = [DecayAngles]","category":"page"},{"location":"95-reference/#DecayAngles.DecayNode-Tuple{Any}","page":"Reference","title":"DecayAngles.DecayNode","text":"DecayNode(node::Tuple, f0=identity)\n\nWrap a bracket topology notation into a DecayNode structure\n\nExample\n\njulia> topology = DecayNode( (1,(2,3)) )\njulia> map_tree(topology) do value, node\n\ttuple_representation = value\n\t(; name=string(value), tuple_representation)\nend\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#DecayAngles.flatten_nested_tuple-Tuple{Any}","page":"Reference","title":"DecayAngles.flatten_nested_tuple","text":"flatten_nested_tuple(n)\n\nFlatten a nested tuple into a single tuple\n\nExample\n\njulia julia> flatten_nested_tuple(((4,1),(2,3))) (4,1,2,3)`\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#DecayAngles.flatten_sort_nested_tuple-Tuple{Any}","page":"Reference","title":"DecayAngles.flatten_sort_nested_tuple","text":"flatten_sort_nested_tuple(n)\n\nFlatten a nested tuple into a single tuple which is also sorted.\n\nExample\n\njulia julia> flatten_sort_nested_tuple(((4,1),(2,3))) (1,2,3,4)`\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#DecayAngles.map_tree-Tuple{Any, DecayNode}","page":"Reference","title":"DecayAngles.map_tree","text":"map_tree(f, node::DecayNode)\n\nIterate over the three and apply a function f(value, node) to every nodes, storing the result in the value field.\n\nExample\n\njulia> map_tree(build_tree((1,(2,3)))) do value, node\n\ttuple_representation = value\n\t(; name=string(value), tuple_representation)\nend\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#DecayAngles.map_with_parent_node","page":"Reference","title":"DecayAngles.map_with_parent_node","text":"map_with_parent_node(f, node::DecayNode, parent_node = DecayNode(nothing, [node]))\n\nPropagates the values from parent to all children and apply operation on it. The signature of the function is f(value, parent_node)\n\nExample\n\njulia> map_with_parent(DecayNode((1,(2,3))), \"0\") do value, parent_value\n    string(parent_value) * \" => \" * string(value)\nend\n\n\n\n\n\n","category":"function"},{"location":"95-reference/#DecayAngles.to_table-Tuple{DecayNode}","page":"Reference","title":"DecayAngles.to_table","text":"totable(tree::DecayNode)\n\nIterate over the tree and table-like collection of named tuples from the node values.\n\nExample\n\njulia> topology = DecayNode( (1,(2,3)); apply=x->(; tuple=x))\njulia> totable(topology) |> DataFrame\n\n\n\n\n\n","category":"method"},{"location":"","page":"DecayAngles","title":"DecayAngles","text":"CurrentModule = DecayAngles","category":"page"},{"location":"#DecayAngles","page":"DecayAngles","title":"DecayAngles","text":"","category":"section"},{"location":"","page":"DecayAngles","title":"DecayAngles","text":"Documentation for DecayAngles.","category":"page"}]
}
