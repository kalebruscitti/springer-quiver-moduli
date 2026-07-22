# A quiver representation of the framed double quiver: the quadruple (A, B, Γ, Δ).
#
#   A_i : V_i → V_{i+1}   (matrix  v_{i+1} × v_i)
#   B_i : V_{i+1} → V_i   (matrix  v_i × v_{i+1})
#   Γ_j : D_j → V_j       (matrix  v_j × d_j)   at framing vertices j
#   Δ_j : V_j → D_j       (matrix  d_j × v_j)   at framing vertices j
#
# The same struct carries either symbolic entries (over a polynomial ring, one
# indeterminate per matrix entry) or numeric entries (over QQ, GF(p), …); the
# equation routines in equations.jl / springer.jl are generic over the ring.

"""
    Rep

A representation of the framed double quiver for a `QuiverData`. Access maps
through `getA`, `getB`, `getGamma`, `getDelta`, which return correctly-shaped
zero matrices for out-of-range / absent indices (the `rmk:0` conventions).
"""
struct Rep
    qd::QuiverData
    R                       # base ring (poly ring for symbolic, field for numeric)
    A::Vector               # A[i] = matrix of A_i, i = 1..n-1
    B::Vector               # B[i] = matrix of B_i, i = 1..n-1
    Gamma::Dict{Int,Any}    # framing vertex j => Γ_j
    Delta::Dict{Int,Any}    # framing vertex j => Δ_j
end

getA(rep::Rep, i::Integer) =
    (1 <= i <= rep.qd.n - 1) ? rep.A[i] :
    zero_matrix(rep.R, getV(rep.qd, i + 1), getV(rep.qd, i))

getB(rep::Rep, i::Integer) =
    (1 <= i <= rep.qd.n - 1) ? rep.B[i] :
    zero_matrix(rep.R, getV(rep.qd, i), getV(rep.qd, i + 1))

getGamma(rep::Rep, j::Integer) =
    haskey(rep.Gamma, j) ? rep.Gamma[j] :
    zero_matrix(rep.R, getV(rep.qd, j), dimD(rep.qd, j))

getDelta(rep::Rep, j::Integer) =
    haskey(rep.Delta, j) ? rep.Delta[j] :
    zero_matrix(rep.R, dimD(rep.qd, j), getV(rep.qd, j))

# ---------------------------------------------------------------------------
# Symbolic representation: one indeterminate per matrix entry.
# ---------------------------------------------------------------------------

"""
    symbolic_rep(qd; base = QQ) -> Rep

Build the generic representation over a polynomial ring `base[…]` with one
indeterminate per entry of every `A_i, B_i, Γ_j, Δ_j`. Variable names are
`A{i}_{r}_{c}`, `B{i}_{r}_{c}`, `G{j}_{r}_{c}`, `D{j}_{r}_{c}`.
"""
function symbolic_rep(qd::QuiverData; base = QQ)
    n = qd.n
    # specs in a fixed order; only blocks with positive size get variables
    specs = Tuple{Symbol,Int,Int,Int}[]
    for i in 1:(n-1)
        r, c = getV(qd, i + 1), getV(qd, i)
        (r > 0 && c > 0) && push!(specs, (:A, i, r, c))
    end
    for i in 1:(n-1)
        r, c = getV(qd, i), getV(qd, i + 1)
        (r > 0 && c > 0) && push!(specs, (:B, i, r, c))
    end
    for j in qd.framing
        r, c = getV(qd, j), dimD(qd, j)
        (r > 0 && c > 0) && push!(specs, (:G, j, r, c))
        (r > 0 && c > 0) && push!(specs, (:D, j, c, r))
    end

    prefix = Dict(:A => "A", :B => "B", :G => "G", :D => "D")
    names = String[]
    for (tag, i, r, c) in specs, a in 1:r, b in 1:c
        push!(names, "$(prefix[tag])$(i)_$(a)_$(b)")
    end

    R, g = polynomial_ring(base, names)

    A = Any[zero_matrix(R, getV(qd, i + 1), getV(qd, i)) for i in 1:(n-1)]
    B = Any[zero_matrix(R, getV(qd, i), getV(qd, i + 1)) for i in 1:(n-1)]
    Gamma = Dict{Int,Any}()
    Delta = Dict{Int,Any}()

    off = 0
    for (tag, i, r, c) in specs
        M = matrix(R, r, c, [g[off + (a - 1) * c + b] for a in 1:r for b in 1:c])
        off += r * c
        tag === :A && (A[i] = M)
        tag === :B && (B[i] = M)
        tag === :G && (Gamma[i] = M)
        tag === :D && (Delta[i] = M)
    end

    return Rep(qd, R, A, B, Gamma, Delta)
end

# ---------------------------------------------------------------------------
# Numeric representation: fill in explicit matrices (for tests / examples).
# ---------------------------------------------------------------------------

"""
    point_hom(sym, num; base = num.R)

Ring homomorphism `sym.R → base` sending each indeterminate of the symbolic
representation `sym` to its value in the numeric representation `num` (built for
the same `QuiverData`). Use it to evaluate any symbolic ideal / polynomial at a
concrete point, e.g. `all(iszero, phi.(gens(I)))`.
"""
function point_hom(sym::Rep, num::Rep; base = num.R)
    R = sym.R
    vals = elem_type(base)[]
    for g in gens(R)
        m = match(r"^([A-Z])(\d+)_(\d+)_(\d+)$", string(g))
        m === nothing && error("cannot parse variable name $(string(g))")
        p = m.captures[1][1]
        i = parse(Int, m.captures[2])
        a = parse(Int, m.captures[3])
        b = parse(Int, m.captures[4])
        M = p == 'A' ? getA(num, i) :
            p == 'B' ? getB(num, i) :
            p == 'G' ? getGamma(num, i) : getDelta(num, i)
        push!(vals, base(M[a, b]))
    end
    return hom(R, base, vals)
end

"""
    numeric_rep(qd; A, B, Gamma, Delta = nothing, base = QQ) -> Rep

Build a representation from explicit matrices. `A`, `B`, `Gamma`, `Delta` are
`Dict{Int, <:AbstractMatrix}` (plain Julia arrays are fine): `A[i]` is `A_i`,
`Gamma[j]` is `Γ_j` at a framing vertex `j`. Missing `Δ` defaults to zero (the
Springer-fibre locus). Sizes are validated against `qd`.
"""
function numeric_rep(qd::QuiverData;
                     A::AbstractDict = Dict{Int,Any}(),
                     B::AbstractDict = Dict{Int,Any}(),
                     Gamma::AbstractDict = Dict{Int,Any}(),
                     Delta = nothing,
                     base = QQ)
    n = qd.n
    toM(arr, r, c, name) = begin
        size(arr, 1) == r && size(arr, 2) == c ||
            throw(ArgumentError("$name has size $(size(arr)), expected ($r, $c)"))
        matrix(base, r, c, [base(arr[a, b]) for a in 1:r for b in 1:c])
    end

    Amat = Any[zero_matrix(base, getV(qd, i + 1), getV(qd, i)) for i in 1:(n-1)]
    Bmat = Any[zero_matrix(base, getV(qd, i), getV(qd, i + 1)) for i in 1:(n-1)]
    for (i, arr) in A
        Amat[i] = toM(arr, getV(qd, i + 1), getV(qd, i), "A[$i]")
    end
    for (i, arr) in B
        Bmat[i] = toM(arr, getV(qd, i), getV(qd, i + 1), "B[$i]")
    end
    Gam = Dict{Int,Any}()
    Del = Dict{Int,Any}()
    for (j, arr) in Gamma
        Gam[j] = toM(arr, getV(qd, j), dimD(qd, j), "Gamma[$j]")
    end
    for j in qd.framing
        if Delta !== nothing && haskey(Delta, j)
            Del[j] = toM(Delta[j], dimD(qd, j), getV(qd, j), "Delta[$j]")
        else
            Del[j] = zero_matrix(base, dimD(qd, j), getV(qd, j))
        end
    end
    return Rep(qd, base, Amat, Bmat, Gam, Del)
end
