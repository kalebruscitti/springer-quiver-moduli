using Test
using Oscar
include(joinpath(@__DIR__, "..", "src", "TwoRowSpringer.jl"))
using .TwoRowSpringer

# helpers to write small integer matrices
col(v...) = reshape(collect(v), length(v), 1)   # column vector as n×1 matrix
row(v...) = reshape(collect(v), 1, length(v))   # row vector as 1×n matrix
m11(x) = reshape([x], 1, 1)                      # 1×1 matrix

@testset "dimension vectors (eq:2rowdv)" begin
    @test QuiverData(4, 2).v == [1, 2, 1]
    @test QuiverData(4, 2).d == [0, 2, 0]
    @test QuiverData(4, 1).v == [1, 1, 1]
    @test QuiverData(4, 1).d == [1, 0, 1]
    @test QuiverData(5, 2).v == [1, 2, 2, 1]
    @test QuiverData(5, 2).d == [0, 1, 1, 0]
    @test QuiverData(6, 2).v == [1, 2, 2, 2, 1]
    @test QuiverData(6, 2).d == [0, 1, 0, 1, 0]
    @test_throws ArgumentError QuiverData(3, 2)   # k > n-k
end

@testset "symbolic stable locus is nonempty machinery" begin
    qd = QuiverData(5, 2)
    rep = symbolic_rep(qd)
    I = adhm_ideal(rep)
    @test I isa Oscar.Ideal
    @test length(gens(I)) == sum(getV(qd, i)^2 for i in 1:(qd.n-1))  # v_i×v_i per vertex
    sd = stability_data(rep)
    @test length(sd) == qd.n - 1
    @test all(t -> ncols(t.matrix) >= t.target, sd)   # enough columns to be stable
end

@testset "cup diagram enumeration" begin
    @test length(cup_diagrams(4, 2)) == 2
    @test length(cup_diagrams(4, 1)) == 3
    @test length(cup_diagrams(6, 3)) == 5    # Catalan C_3
    @test length(cup_diagrams(6, 2)) == 9
    # no ray strictly inside a cup ⇒ cups have odd length
    for D in cup_diagrams(8, 3), (i, j) in D.cups
        @test isodd(j - i)
    end
end

@testset "paper example (2,2)  (lines 2233–2250)" begin
    qd = QuiverData(4, 2)
    rep = numeric_rep(qd;
        A = Dict(1 => col(1, 0), 2 => row(0, 1)),
        B = Dict(1 => row(0, 1), 2 => col(1, 0)),
        Gamma = Dict(2 => [1 0; 0 1]))          # Δ defaults to 0
    @test is_admissible(rep)
    @test is_stable(rep)
    @test in_stable_locus(rep)

    diags = cup_diagrams(4, 2)
    d_nested   = first(filter(D -> (1, 4) in D.cups, diags))   # {(1,4),(2,3)}
    d_adjacent = first(filter(D -> (1, 2) in D.cups, diags))   # {(1,2),(3,4)}
    @test in_component(rep, d_nested)
    @test !in_component(rep, d_adjacent)
    # lands in exactly one component
    @test count(D -> in_component(rep, D), diags) == 1
end

@testset "paper example (3,1)  (lines 2317–2330)" begin
    qd = QuiverData(4, 1)
    rep = numeric_rep(qd;
        A = Dict(1 => m11(1), 2 => m11(0)),
        B = Dict(1 => m11(0), 2 => m11(1)),
        Gamma = Dict(1 => m11(1), 3 => m11(1)))
    @test is_admissible(rep)
    @test is_stable(rep)

    diags = cup_diagrams(4, 1)
    d2 = first(filter(D -> (2, 3) in D.cups, diags))   # component a2: ker B_1 = ker A_2
    @test in_component(rep, d2)
    @test count(D -> in_component(rep, D), diags) == 1
end

@testset "gauge invariance of the equations (def:GL action)" begin
    qd = QuiverData(4, 2)
    rep = numeric_rep(qd;
        A = Dict(1 => col(1, 0), 2 => row(0, 1)),
        B = Dict(1 => row(0, 1), 2 => col(1, 0)),
        Gamma = Dict(2 => [1 0; 0 1]))
    g = Dict(1 => m11(3), 2 => [1 1; 0 1], 3 => m11(2))   # invertible per vertex
    rep2 = act(g, rep)
    @test is_admissible(rep2)
    @test is_stable(rep2)
    d_nested = first(filter(D -> (1, 4) in D.cups, cup_diagrams(4, 2)))
    @test in_component(rep2, d_nested)   # component membership is gauge-invariant
end
