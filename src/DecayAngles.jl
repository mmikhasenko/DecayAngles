module DecayAngles

using AbstractTrees
using Setfield
using Parameters
using FourVectors

export DecayNode
export map_tree
export map_with_parent_node
export to_table
export nodevalue
export children
include("decaynode.jl")

export flatten_nested_tuple
export flatten_sort_nested_tuple
include("utils.jl")

export AnglesGamma
export HelicityTransformation
export DaughterTransformation
include("transformations.jl")

export add_transform_through
export add_indices_order
export decay_angles
include("traversing.jl")

end # module DecayAngles
