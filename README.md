# TwoRowSpringer

Julia/[Oscar](https://www.oscar-system.org/) code that realises **two-row Springer
fibres** `B_{(n-k,k)}` as **Nakajima quiver varieties** and produces, for a given
shape `(n-k, k)`:

1. the framed double quiver of type `A_{n-1}` and its dimension vectors `d, v`;
2. the **equations cutting out the stable locus** `Λ⁺(d,v)` — ADHM/admissibility
   (polynomial equations) plus the Maffei stability condition (an open full-rank
   condition);
3. the **Springer restriction** `Δ = 0` and the **cup/ray kernel relations** that
   cut out each irreducible component of the fibre.

Following Im–Lai–Wang: `arXiv:2009.08778` (type-A quiver part) and `arXiv:1910.03010`
(original combined manuscript). Equation labels below refer to `ILWquiver20210905.tex`.

## Construction (formula ↔ code)

Type `A_{n-1}` framed double quiver: spaces `V_i` (vertices `1..n-1`), framing `D_i`;
maps `A_i:V_i→V_{i+1}`, `B_i:V_{i+1}→V_i`, `Γ_j:D_j→V_j`, `Δ_j:V_j→D_j`.

| Object | Paper | Code |
|---|---|---|
| dims `dim V_i`, `dim D_i` | `eq:2rowdv` | `QuiverData(n,k)` → `.v`, `.d` |
| generic representation | `R(d,v)` | `symbolic_rep(qd)` |
| composite `Γ_{j→i}` | `def:Gaij` | `composite_Gamma(rep,j,i)` |
| ADHM `B_iA_i=A_{i-1}B_{i-1}+Γ_iΔ_i` | `eq:L1` | `adhm_ideal(rep)`, `adhm_residuals(rep)` |
| stability `Im A_{i-1}+ΣIm Γ_{j→i}=V_i` | `eq:L2` | `stability_data(rep)`, `is_stable` |
| Springer locus `Δ=0` | `prop:Spr` | `springer_restrict(rep)` |
| components ↔ cup diagrams | `defi:…comp…` | `cup_diagrams(n,k)` |
| cup relation `ker B_{s→i-1}=ker A_{s→j}` | `defi:…comp…` | `cup_relations(rep,D)` |
| ray relation `B_iA_i=0` / `Γ_{n-k→i}=0` | `prop:ray` | `ray_relations(rep,D)` |
| component ideal / coordinate ring | — | `component_ideal`, `component_coordinate_ring` |
| `GL(V)` action | `def:GL action` | `act(g,rep)` |

**Stability is an open condition.** `M_i = [A_{i-1} | Γ_{j→i} (j≥i)]` must have full
row rank `v_i`; the stable locus is where, for every `i`, **not all** maximal minors
of `M_i` vanish (`stability_minors(rep,i)`). `is_stable`/`in_stable_locus` check this
numerically by rank.

**Cup relations are kernel equalities** (`ker M = ker N`), a locally closed
condition — not a single vanishing ideal. `cup_relations` returns the map pair
`(M, N, stacked=[M;N])` (impose `rank[M;N]=rank M=rank N` symbolically);
`satisfies_cups` checks equality exactly on numeric points. **Ray relations** are
genuine polynomial equations: `ray_relations` / `ray_ideal_gens`.

## Usage

```julia
include("src/TwoRowSpringer.jl"); using .TwoRowSpringer, Oscar

qd  = QuiverData(5, 2)          # shape (3,2):  dim V = [1,2,2,1], dim D = [0,1,1,0]
rep = symbolic_rep(qd)          # generic rep over a 24-variable polynomial ring

I   = adhm_ideal(rep)           # the ADHM polynomial equations (an Oscar ideal)
sd  = stability_data(rep)       # per-vertex stacked matrix M_i + its maximal minors

srep = springer_restrict(rep)   # Δ = 0
for D in cup_diagrams(qd)       # irreducible components
    cup_relations(srep, D)      # kernel-equality data per cup
    ray_relations(srep, D)      # B_iA_i=0 or Γ_{n-k→i}=0 per ray
end
```

Build a concrete point and test membership with `numeric_rep` + `in_component`
(defaults `Δ=0`); see `test/runtests.jl` for the paper's `(2,2)` and `(3,1)` points.

### Coordinate ring of a component (for Oscar)

```julia
Q, proj = component_coordinate_ring(rep, D)   # D a cup diagram
krull_dim(Q)                                  # dimension of the component locus
I = component_ideal(rep, D)
is_prime(I)                                   # irreducibility
```

`component_ideal(rep, D)` returns the ideal in `rep.R` cutting out (the closure
of) the stable component locus `Λ^a`: admissible (`eq:L1`), `Δ=0`, ray relations,
and the cup relations `ker M = ker N`. Because a cup relation is a **kernel
equality**, not a vanishing condition, each inclusion is encoded with an
auxiliary matrix (`ker M ⊆ ker N ⟺ ∃C: N = CM`) and those variables are then
**eliminated**, so the result lives in the original representation coordinates.
**Stability** is an open condition, imposed by **saturating** against the maximal
minors of each `M_i` (removing unstable components). Keywords:
`delta_zero`, `include_rays`, `saturate_stability` (set the last to `false` for a
much cheaper, unsaturated ideal). E.g. for `(3,1)` the component `cup (2,3)`
comes out as the prime ideal `(Δ, B_1, A_2)`, `krull_dim = 4` (dim 1 modulo `GL(V)`).

Pretty-print everything for a shape:

```
julia demo.jl 5 2      # (n-k,k) = (3,2)
```

## Files

```
src/TwoRowSpringer.jl   module + exports
src/quiver.jl           QuiverData, dimension vectors (eq:2rowdv)
src/representation.jl   Rep; symbolic_rep, numeric_rep
src/equations.jl        composites, adhm_ideal, stability_data (eq:L1, eq:L2)
src/springer.jl         Δ=0, cup_diagrams, cup/ray relations, component_ideal
src/gauge.jl            GL(V) action
test/runtests.jl        verified against the paper's (2,2) & (3,1) examples
demo.jl                 CLI pretty-printer
```

## Running the tests

```
julia test/runtests.jl
```

(uses an environment with Oscar available; or `julia --project=. -e 'using Pkg;
Pkg.instantiate()'` first to use the bundled `Project.toml`).

## Scope / not yet implemented

- The **flag realization** map `M_{(n-k,k)} → S̃_{(n-k,k)}` (Thm `thm:Fi`,
  `F_i = ker[A_{1→i}Γ_{→1} | … | Γ_{→i}]`) is deferred; the composite-map plumbing
  (`composite_Gamma`, `A_comp`, `B_comp`) is in place to add it.
- Base ring defaults to `QQ`; pass `base = GF(p)` / a number field to `symbolic_rep`
  / `numeric_rep` for other fields.
- Equal-size case `n = 2k` uses a single 2-dimensional framing at vertex `k`.
