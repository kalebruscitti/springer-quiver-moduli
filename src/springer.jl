# The Springer fibre inside the quiver variety (Δ = 0, prop:Spr) and the
# cup/ray kernel relations cutting out its irreducible components
# (defi:Springer_component_in_quiver line 3007, prop:ray line 3409).

# ---------------------------------------------------------------------------
# Springer-fibre restriction:  Δ = 0.
# ---------------------------------------------------------------------------

"`springer_restrict(rep)` — the same representation with `Δ = 0` (prop:Spr)."
springer_restrict(rep::Rep) = Rep(rep.qd, rep.R, rep.A, rep.B, rep.Gamma, Dict{Int,Any}())

# ---------------------------------------------------------------------------
# Cup diagrams: k cups, n-2k rays; cups are non-crossing with no ray strictly
# inside a cup (this forces odd cup length automatically).
# A diagram is a NamedTuple (cups = Vector{Tuple{Int,Int}}, rays = Vector{Int}).
# ---------------------------------------------------------------------------

function _diagrams(lo::Int, hi::Int, allow_rays::Bool)
    lo > hi && return [(Tuple{Int,Int}[], Int[])]
    res = Tuple{Vector{Tuple{Int,Int}},Vector{Int}}[]
    for m in (lo+1):hi                                   # cup (lo, m)
        for (icups, irays) in _diagrams(lo + 1, m - 1, false)   # inside: pure cups
            for (rcups, rrays) in _diagrams(m + 1, hi, allow_rays)
                push!(res, (vcat([(lo, m)], icups, rcups), vcat(irays, rrays)))
            end
        end
    end
    if allow_rays                                        # lo is a ray
        for (rcups, rrays) in _diagrams(lo + 1, hi, true)
            push!(res, (rcups, vcat([lo], rrays)))
        end
    end
    return res
end

"""
    cup_diagrams(n, k)

All standard cup diagrams on `n` vertices with `k` cups and `n-2k` rays; these
index the irreducible components of the two-row Springer fibre `B_{(n-k,k)}`.
"""
function cup_diagrams(n::Integer, k::Integer)
    return [(cups = sort(c), rays = sort(r))
            for (c, r) in _diagrams(1, Int(n), true) if length(c) == k]
end

cup_diagrams(qd::QuiverData) = cup_diagrams(qd.n, qd.k)

# ---------------------------------------------------------------------------
# Cup relations (defi:Springer_component_in_quiver):  for a cup {i, j},
#   ker B_{s → i-1} = ker A_{s → j},   s = (i+j-1)/2   (common source V_s).
# Encoded symbolically as the map pair (M, N) and their stacking; verified
# numerically by exact kernel comparison.
# ---------------------------------------------------------------------------

"""
    cup_relations(rep, diagram)

For each cup `{i,j}` return `(cup, M, N, stacked)` where `M = B_{s→i-1}`,
`N = A_{s→j}` (`s = (i+j-1)/2`) and `stacked = [M; N]`. The defining relation is
`ker M = ker N`, i.e. `rank(stacked) = rank(M) = rank(N)`.
"""
function cup_relations(rep::Rep, diagram)
    out = []
    for (i, j) in diagram.cups
        s = (i + j - 1) ÷ 2
        M = B_comp(rep, s, i - 1)      # B_{s → i-1} : V_s → V_{i-1}
        N = A_comp(rep, s, j)          # A_{s → j}   : V_s → V_j
        push!(out, (cup = (i, j), M = M, N = N, stacked = vcat(M, N)))
    end
    return out
end

"`_ker_equal(M, N)` — exact test `ker M = ker N` (right kernels, shared source)."
function _ker_equal(M, N)
    KM = kernel(M; side = :right)
    KN = kernel(N; side = :right)
    return iszero(N * KM) && iszero(M * KN)
end

"`satisfies_cups(rep, diagram)` — numeric: every cup relation `ker M = ker N` holds."
satisfies_cups(rep::Rep, diagram) =
    all(r -> _ker_equal(r.M, r.N), cup_relations(rep, diagram))

# ---------------------------------------------------------------------------
# Ray relations (prop:ray):  at a ray i with c(i) = #cups fully left of i,
#   B_i A_i = 0            if c(i) ≥ 1,
#   Γ_{n-k → i} = 0        if c(i) = 0.
# ---------------------------------------------------------------------------

"`_cups_left(diagram, i)` — number of cups `(a,b)` with `b < i`."
_cups_left(diagram, i) = count(((a, b),) -> b < i, diagram.cups)

"""
    ray_relations(rep, diagram)

For each ray `i` return `(ray, kind, matrix)`: `kind = :BA` with matrix `B_iA_i`
(must be 0) when `c(i) ≥ 1`, or `kind = :Gamma` with matrix `Γ_{n-k→i}` (must be
0) when `c(i) = 0`. These are genuine polynomial equations (all entries vanish).
"""
function ray_relations(rep::Rep, diagram)
    nk = rep.qd.n - rep.qd.k
    out = []
    for i in diagram.rays
        if _cups_left(diagram, i) >= 1
            push!(out, (ray = i, kind = :BA, matrix = getB(rep, i) * getA(rep, i)))
        else
            push!(out, (ray = i, kind = :Gamma, matrix = composite_Gamma(rep, nk, i)))
        end
    end
    return out
end

"`ray_ideal_gens(rep, diagram)` — all entrywise generators of the ray relations."
function ray_ideal_gens(rep::Rep, diagram)
    gens = elem_type(rep.R)[]
    for r in ray_relations(rep, diagram)
        M = r.matrix
        for a in 1:nrows(M), b in 1:ncols(M)
            push!(gens, M[a, b])
        end
    end
    return gens
end

"`satisfies_rays(rep, diagram)` — numeric: every ray relation matrix is zero."
satisfies_rays(rep::Rep, diagram) =
    all(r -> iszero(r.matrix), ray_relations(rep, diagram))

"""
    in_component(rep, diagram)

Numeric membership test for the irreducible component `Λ^a` indexed by `diagram`:
Springer locus (`Δ = 0`), admissible, stable, and all cup + ray relations hold.
"""
function in_component(rep::Rep, diagram)
    all(iszero, values(rep.Delta)) || return false
    is_admissible(rep) && is_stable(rep) &&
        satisfies_cups(rep, diagram) && satisfies_rays(rep, diagram)
end

# ---------------------------------------------------------------------------
# The defining ideal / coordinate ring of a component locus, for Oscar.
#
# The locus is  { admissible (eq:L1),  Δ = 0,  ray relations,  ker M = ker N
#                 for every cup,  stable (eq:L2) }  ⊆ R(d,v).
#
# A cup relation `ker M = ker N` is not a set of vanishing polynomials directly;
# it is the double row-space inclusion, encoded with auxiliary matrices
#     N = C·M   (⟺ ker M ⊆ ker N)   and   M = D·N   (⟺ ker N ⊆ ker M),
# whose variables are then eliminated. Stability is the open condition "not all
# maximal minors of each M_i vanish", imposed by saturating against those minors.
# ---------------------------------------------------------------------------

"""
    component_ideal(rep, diagram;
                    delta_zero = true, include_rays = true, saturate_stability = true)

Ideal in `rep.R` cutting out (the Zariski closure of) the stable component locus
`Λ^a` for `diagram`. Built from a **symbolic** `Rep`.

- `delta_zero`         : impose `Δ = 0` (the Springer fibre, prop:Spr).
- `include_rays`       : add the ray relations (prop:ray; implied on the stable
                         locus, but they sharpen the ideal before saturation).
- `saturate_stability` : saturate against the stability minors so that
                         unstable components are removed (eq:L2). Set `false` to
                         keep the raw (unsaturated) ideal — much cheaper.

Cup relations are encoded with auxiliary matrices that are then eliminated, so
the result lives in the original representation coordinates and can be handed to
Oscar (`dim`, `is_prime`, `radical`, Hilbert series, …).
"""
function component_ideal(rep::Rep, diagram;
                         delta_zero::Bool = true,
                         include_rays::Bool = true,
                         saturate_stability::Bool = true)
    R = rep.R
    base = base_ring(R)
    srep = delta_zero ? springer_restrict(rep) : rep
    cups = cup_relations(rep, diagram)

    # auxiliary variables: C_ci (b×a) and D_ci (a×b) for each cup with a,b > 0
    auxnames = String[]
    for (ci, r) in enumerate(cups)
        a, b = nrows(r.M), nrows(r.N)
        if a > 0 && b > 0
            append!(auxnames, ["cC$(ci)_$(x)_$(y)" for x in 1:b for y in 1:a])
            append!(auxnames, ["cD$(ci)_$(x)_$(y)" for x in 1:a for y in 1:b])
        end
    end

    orignames = string.(gens(R))
    Rext, extg = polynomial_ring(base, vcat(orignames, auxnames))
    phi = hom(R, Rext, extg[1:ngens(R)])
    pushcol(gs, M) = for x in 1:nrows(M), y in 1:ncols(M); push!(gs, M[x, y]); end

    gs = elem_type(Rext)[]
    for Mres in adhm_residuals(srep)                 # admissibility (eq:L1)
        pushcol(gs, map_entries(phi, Mres))
    end
    if delta_zero                                    # Δ = 0
        for j in rep.qd.framing
            pushcol(gs, map_entries(phi, getDelta(rep, j)))
        end
    end
    if include_rays                                  # ray relations (prop:ray)
        for g in ray_ideal_gens(srep, diagram)
            push!(gs, phi(g))
        end
    end

    off = ngens(R)                                   # cup relations (kernel equality)
    for (ci, r) in enumerate(cups)
        a, b = nrows(r.M), nrows(r.N)
        Mext, Next = map_entries(phi, r.M), map_entries(phi, r.N)
        if a > 0 && b > 0
            C = matrix(Rext, b, a, [extg[off + (x-1)*a + y] for x in 1:b for y in 1:a]); off += a*b
            D = matrix(Rext, a, b, [extg[off + (x-1)*b + y] for x in 1:a for y in 1:b]); off += a*b
            pushcol(gs, Next - C * Mext)             # N = C M  ⟺ ker M ⊆ ker N
            pushcol(gs, Mext - D * Next)             # M = D N  ⟺ ker N ⊆ ker M
        elseif a == 0 && b > 0
            pushcol(gs, Next)                        # ker M = V_s ⇒ N = 0
        elseif b == 0 && a > 0
            pushcol(gs, Mext)                        # ker N = V_s ⇒ M = 0
        end
    end

    Iext = ideal(Rext, gs)
    auxvars = extg[(ngens(R)+1):end]
    Iel = isempty(auxvars) ? Iext : eliminate(Iext, auxvars)

    psi = hom(Rext, R, vcat(gens(R), fill(R(0), length(auxnames))))
    IR = ideal(R, [psi(g) for g in gens(Iel)])

    if saturate_stability
        for i in 1:(rep.qd.n - 1)
            mins = stability_minors(rep, i)
            isempty(mins) && continue
            IR = saturation(IR, ideal(R, mins))
        end
    end
    return IR
end

"""
    component_coordinate_ring(rep, diagram; kwargs...) -> (Q, proj)

The coordinate ring `Q = rep.R / component_ideal(rep, diagram; …)` as an Oscar
quotient ring, together with the projection. Pass to Oscar to check properties,
e.g. `krull_dim(Q)`, `is_prime(modulus(Q))`. Keyword arguments are forwarded to
[`component_ideal`](@ref).
"""
function component_coordinate_ring(rep::Rep, diagram; kwargs...)
    return quo(rep.R, component_ideal(rep, diagram; kwargs...))
end
