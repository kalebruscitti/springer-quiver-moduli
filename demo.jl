# Print the quiver, stable-locus equations, and Springer-component relations for
# a two-row shape.  Usage:  julia demo.jl [n] [k]   (default n=5, k=2 → (3,2)).

include(joinpath(@__DIR__, "src", "TwoRowSpringer.jl"))
using .TwoRowSpringer
using Oscar

n = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 5
k = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 2

qd  = QuiverData(n, k)
rep = symbolic_rep(qd)

println("="^70)
println(qd)
println("\n" * "="^70)
println("STABLE LOCUS  Λ⁺(d,v)")
println("="^70)

println("\nADHM / admissibility equations  (eq:L1):  B_iA_i − A_{i-1}B_{i-1} − Γ_iΔ_i = 0")
I = adhm_ideal(rep)
println("  polynomial ring in ", ngens(rep.R), " variables; ",
        length(gens(I)), " generators. Sample:")
for g in gens(I)[1:min(6, end)]
    println("    ", g)
end

println("\nStability  (eq:L2):  each M_i = [ A_{i-1} | Γ_{j→i}, j≥i ] has full row rank v_i")
for t in stability_data(rep)
    println("  vertex $(t.vertex):  M_$(t.vertex) is $(nrows(t.matrix))×$(ncols(t.matrix)), ",
            "need rank $(t.target);  #maximal minors = $(length(t.maximal_minors))")
end

println("\n" * "="^70)
println("SPRINGER FIBRE  (Δ = 0)  and IRREDUCIBLE COMPONENTS")
println("="^70)
srep = springer_restrict(rep)
diags = cup_diagrams(qd)
println("\n$(length(diags)) cup diagram(s) ↔ irreducible components:\n")
for (idx, D) in enumerate(diags)
    println("  component $idx:  cups = $(D.cups),  rays = $(D.rays)")
    for r in cup_relations(srep, D)
        println("     cup $(r.cup):  ker($(nrows(r.M))×$(ncols(r.M)) B-composite) = ",
                "ker($(nrows(r.N))×$(ncols(r.N)) A-composite)")
    end
    for r in ray_relations(srep, D)
        eq = r.kind === :BA ? "B_$(r.ray) A_$(r.ray) = 0" : "Γ_{$(n-k)→$(r.ray)} = 0"
        println("     ray $(r.ray):  $eq   ($(nrows(r.matrix))×$(ncols(r.matrix)))")
    end
end
println()
