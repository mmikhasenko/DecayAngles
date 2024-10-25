# DecayAngles

<!-- [![Stable Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://mmikhasenko.github.io/DecayAngles.jl/stable) -->
[![In development documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://mmikhasenko.github.io/DecayAngles.jl/dev)
[![Build Status](https://github.com/mmikhasenko/DecayAngles.jl/workflows/Test/badge.svg)](https://github.com/mmikhasenko/DecayAngles.jl/actions)
[![Test workflow status](https://github.com/mmikhasenko/DecayAngles.jl/actions/workflows/Test.yml/badge.svg?branch=main)](https://github.com/mmikhasenko/DecayAngles.jl/actions/workflows/Test.yml?query=branch%3Amain)
[![Lint workflow Status](https://github.com/mmikhasenko/DecayAngles.jl/actions/workflows/Lint.yml/badge.svg?branch=main)](https://github.com/mmikhasenko/DecayAngles.jl/actions/workflows/Lint.yml?query=branch%3Amain)
[![Docs workflow Status](https://github.com/mmikhasenko/DecayAngles.jl/actions/workflows/Docs.yml/badge.svg?branch=main)](https://github.com/mmikhasenko/DecayAngles.jl/actions/workflows/Docs.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/mmikhasenko/DecayAngles.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/mmikhasenko/DecayAngles.jl)
[![BestieTemplate](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/JuliaBesties/BestieTemplate.jl/main/docs/src/assets/badge.json)](https://github.com/JuliaBesties/BestieTemplate.jl)

**DecayAngles.jl** is a Julia package designed to simplify the calculation of angular observables and transformations for decay processes. The package provides a tree structure that can represent decay chains and facilitates the extraction and manipulation of angular information at each node. The design is flexible and efficient, enabling the computation of complex decay observables using recursive tree traversal and parent-child relationships between nodes.


## Features

- **Tree Representation of Decay Chains**: Represents decay chains as hierarchical tree structures, where each node holds the decay parameters and observables.
- **Recursive Traversal**: Allows recursive operations on the tree, providing access to both the parent and child nodes during calculations.
- **Table Interface**: values stored at the node are easy to convert into Table.

## Installation

To install **DecayAngles.jl**, you can add the package via the Julia package manager:

```julia
] add https://github.com/yourusername/DecayAngles.jl
```

## Usage

### Basic Usage

To start using **DecayAngles.jl**, you can create a tree structure to represent a decay chain and then perform operations on each node in the tree.

```julia
using DecayAngles

# Define a decay topology as a nested tuple structure
topology = (1, (2, 3))

# a special constructor from a tuple
dn = DecayNode(topology)
```

In this example, `dn` is a `DecayNode` that represents a decay topology where particle `1` and the pair `(2, 3)` form the branching structure.

#### Mapping over the Tree

The `map_tree` function allows you to apply an operation to each node in the tree. Here, we create a representation of each node by converting its value to a string.

```julia
name_repr = map_tree(dn) do value, node
    tuple_representation = value
    (; name=string(value), tuple_representation)
end
```

After this operation, `name_repr` is a new tree where each node contains a `NamedTuple` with a `name` (a string representation of the original node’s value) and a `tuple_representation`.

#### Mapping with Parent Information

Using `map_with_parent`, you can apply a function that has access to both the current node and its parent’s value, which can be useful for calculations that depend on the hierarchical structure of the decay chain.

```julia
with_parent = map_with_parent(dn, "X") do value, parent_value
    string(parent_value) * " => " * string(value)
end
```

In this example, each node in `with_parent` contains a string showing the path from the root node to the current node in the form of `"parent => child"`. For the root node, the initial `parent_value` is set to `"X"`, which you can customize.

### Summary of Key Functions

- **`DecayNode`**: Represents nodes in a tree structure for decay chains, storing kinematic or angular data.
- **`map_tree`**: Traverses the tree and applies a function to each node, generating a new tree with transformed node values.
- **`map_with_parent`**: Similar to `map_tree`, but provides access to the parent’s value, allowing calculations that depend on both the parent and child nodes.

This setup allows you to represent decay chains as trees, traverse and manipulate them flexibly, and perform calculations at each node or across parent-child relationships.

---

## Development and Contribution

We welcome contributions and feedback! Feel free to submit pull requests or issues on the [GitHub repository](https://github.com/mmikhasenko/DecayAngles.jl).

---

### Acknowledgments

This package was inspired by the work on [decayangle](https://github.com/KaiHabermann/decayangle), co-developed with Kai Habermann, and extends the concept of handling decay angular observables using a more flexible, tree-based approach.

---

### License

This package is licensed under the MIT License. See the [LICENSE](link_to_license) file for details.
