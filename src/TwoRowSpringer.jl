"""
    TwoRowSpringer

Computer-algebra models of two-row Springer fibres `B_{(n-k,k)}` as Nakajima
quiver varieties, following Im–Lai–Wang (arXiv:2009.08778, arXiv:1910.03010).

Given a shape `(n-k, k)` this provides:
  * the framed double quiver of type `A_{n-1}` and its dimension vectors — `QuiverData`;
  * the generic representation over a polynomial ring — `symbolic_rep`;
  * the equations of the stable locus `Λ⁺(d,v)`: ADHM (`adhm_ideal`) and Maffei
    stability (`stability_data`);
  * the Springer restriction `Δ = 0` (`springer_restrict`) and the cup/ray kernel
    relations of each irreducible component (`cup_diagrams`, `cup_relations`,
    `ray_relations`, `in_component`).
"""
module TwoRowSpringer

using Oscar

include("quiver.jl")
include("representation.jl")
include("equations.jl")
include("springer.jl")
include("gauge.jl")

export QuiverData, dimV, dimD, getV
export Rep, symbolic_rep, numeric_rep, point_hom, getA, getB, getGamma, getDelta
export A_comp, B_comp, composite_Gamma
export adhm_residuals, adhm_ideal, is_admissible
export stability_matrix, stability_minors, stability_data,
       is_stable, in_stable_locus
export springer_restrict, cup_diagrams, cup_relations, ray_relations,
       ray_ideal_gens, satisfies_cups, satisfies_rays, in_component,
       component_ideal, component_coordinate_ring
export act

end # module
