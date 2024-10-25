using DecayAngles
using Test

topology = (1,(2,3))
dn = DecayNode(topology)

@testset "Structure" begin
    @test dn isa DecayNode
end


name_repr = map_tree(dn) do value, node
    tuple_representation = value
    (; name=string(value), tuple_representation)
end

@testset "Map tree" begin
    @test nodevalue(name_repr).name == "(1, (2, 3))"
    @test nodevalue(children(name_repr)[1]).name == "1"
    @test nodevalue(children(name_repr)[2]).name == "(2, 3)"
end


with_parent = map_with_parent(dn, "X") do value, parent_value
	string(parent_value) * " => " * string(value)
end

@testset "Map with parent" begin
    @test nodevalue(children(with_parent)[1]) == ("X => (1, (2, 3)) => 1")
    @test nodevalue(children(with_parent)[2]) == ("X => (1, (2, 3)) => (2, 3)")
end

@testset "Flatten sorted tuple" begin
    topology = ((4,1),(2,3))
    @test flatten_nested_tuple(topology) == (4,1,2,3)
    @test flatten_sort_nested_tuple(topology) == (1,2,3,4)
end
