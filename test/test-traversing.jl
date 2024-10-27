using DecayAngles
using Test
using DataFrames

@testset "Add Indices Order" begin
    @test add_indices_order(DecayNode(((1, 2)))) |> to_table |> DataFrame ==
          DataFrame(
        names = [(1, 2), 1, 2],
        indices = [(1, 2), (1,), (2,)],
        isfirst = [true, true, false])

    @test add_indices_order(DecayNode(("ah", ("bh", "cx")))) |> to_table |> DataFrame ==
          DataFrame(
        names = [("ah", ("bh", "cx")), "ah", ("bh", "cx"), "bh", "cx"],
        indices = [("ah", "bh", "cx"), ("ah",), ("bh", "cx"), ("bh",), ("cx",)],
        isfirst = [true, true, false, true, false])
end
