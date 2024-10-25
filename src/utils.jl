"""
    flatten_nested_tuple(n)

Flatten a nested tuple into a single tuple

## Example
```julia
julia> flatten_nested_tuple(((4,1),(2,3)))
(4,1,2,3)
````
"""
flatten_nested_tuple(n) = Tuple(collect(Leaves(n)))

"""
    flatten_sort_nested_tuple(n)

Flatten a nested tuple into a single tuple which is also sorted.

## Example
```julia
julia> flatten_sort_nested_tuple(((4,1),(2,3)))
(1,2,3,4)
````
"""
flatten_sort_nested_tuple(n) = Tuple(sort(collect(Leaves(n))))
