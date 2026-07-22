# Equations cutting out the stable locus Λ⁺(d,v):
#   * ADHM / admissibility  (eq:L1)  — polynomial equations,
#   * Maffei stability      (eq:L2)  — an open full-rank condition.
# Plus the composite maps Γ_{j→i}, A_{p→q}, B_{q→p} they are built from.

# ---------------------------------------------------------------------------
# Composite maps (def:Gaij line 1170, eq:to).
# ---------------------------------------------------------------------------

"`A_comp(rep, p, q)` = A_{p→q} = A_{q-1}⋯A_p : V_p → V_q  (p ≤ q)."
function A_comp(rep::Rep, p::Integer, q::Integer)
    M = identity_matrix(rep.R, getV(rep.qd, p))
    for t in p:(q-1)
        M = getA(rep, t) * M
    end
    return M
end

"`B_comp(rep, q, p)` = B_{q→p} = B_p⋯B_{q-1} : V_q → V_p  (q ≥ p)."
function B_comp(rep::Rep, q::Integer, p::Integer)
    M = identity_matrix(rep.R, getV(rep.qd, p))
    for t in p:(q-1)
        M = M * getB(rep, t)
    end
    return M
end

"""
    composite_Gamma(rep, j, i)

The composite `Γ_{j→i} : D_j → V_i` (def:Gaij): `B_i⋯B_{j-1} Γ_j` if `j ≥ i`,
`A_{i-1}⋯A_j Γ_j` if `j ≤ i`.
"""
function composite_Gamma(rep::Rep, j::Integer, i::Integer)
    if j >= i
        return B_comp(rep, j, i) * getGamma(rep, j)
    else
        return A_comp(rep, j, i) * getGamma(rep, j)
    end
end

# ---------------------------------------------------------------------------
# ADHM / admissibility  (eq:L1):  B_i A_i = A_{i-1} B_{i-1} + Γ_i Δ_i.
# ---------------------------------------------------------------------------

"`adhm_residuals(rep)` — the matrices `B_iA_i − A_{i-1}B_{i-1} − Γ_iΔ_i`, `i = 1..n-1`."
function adhm_residuals(rep::Rep)
    return [getB(rep, i) * getA(rep, i) -
            getA(rep, i - 1) * getB(rep, i - 1) -
            getGamma(rep, i) * getDelta(rep, i)
            for i in 1:(rep.qd.n - 1)]
end

"""
    adhm_ideal(rep)

The admissibility ideal: all entries of the ADHM residuals vanish. Returns an
Oscar ideal in the polynomial ring `rep.R` (use on a symbolic `Rep`).
"""
function adhm_ideal(rep::Rep)
    gens = elem_type(rep.R)[]
    for M in adhm_residuals(rep)
        for a in 1:nrows(M), b in 1:ncols(M)
            push!(gens, M[a, b])
        end
    end
    return ideal(rep.R, gens)
end

"`is_admissible(rep)` — numeric check that every ADHM residual vanishes."
is_admissible(rep::Rep) = all(iszero, adhm_residuals(rep))

# ---------------------------------------------------------------------------
# Stability  (eq:L2):  Im A_{i-1} + Σ_{j≥i} Im Γ_{j→i} = V_i  ⇔  M_i has full row rank v_i.
# ---------------------------------------------------------------------------

"""
    stability_matrix(rep, i)

The stacked map `M_i = [ A_{i-1} | Γ_{j→i} (framing j ≥ i) ] : V_{i-1} ⊕ ⨁ D_j → V_i`.
Stability at vertex `i` is exactly: `M_i` has full row rank `v_i`.
"""
function stability_matrix(rep::Rep, i::Integer)
    blocks = Any[getA(rep, i - 1)]
    for j in rep.qd.framing
        j >= i && push!(blocks, composite_Gamma(rep, j, i))
    end
    return reduce(hcat, blocks)
end

"""
    stability_minors(rep, i)

The maximal (`v_i × v_i`) minors of `M_i`. The stable locus is where, for every
`i`, **not all** of these vanish (an open condition). Returns `[]` if `M_i` has
fewer than `v_i` columns (then vertex `i` can never be stable).
"""
function stability_minors(rep::Rep, i::Integer)
    M = stability_matrix(rep, i)
    r = getV(rep.qd, i)
    ncols(M) < r && return elem_type(rep.R)[]
    return minors(M, r)
end

"`is_stable(rep)` — numeric check: every `M_i` has full row rank `v_i` (eq:L2)."
function is_stable(rep::Rep)
    for i in 1:(rep.qd.n - 1)
        M = stability_matrix(rep, i)
        rank(M) == getV(rep.qd, i) || return false
    end
    return true
end

"`in_stable_locus(rep)` — numeric: admissible (eq:L1) **and** stable (eq:L2)."
in_stable_locus(rep::Rep) = is_admissible(rep) && is_stable(rep)

"""
    stability_data(rep)

Human-readable summary: for each vertex `i`, the stacked matrix `M_i`, the target
dimension `v_i`, and the list of maximal minors whose common vanishing is the
*unstable* locus at `i`.
"""
function stability_data(rep::Rep)
    return [(vertex = i,
             target = getV(rep.qd, i),
             matrix = stability_matrix(rep, i),
             maximal_minors = stability_minors(rep, i))
            for i in 1:(rep.qd.n - 1)]
end
