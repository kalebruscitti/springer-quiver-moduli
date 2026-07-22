# Quiver data for a two-row shape (n-k, k).
#
# Framed double quiver of type A_{n-1} (ILWquiver §2, eq:2rowdv, lines 976-1002).
# Vertices 1..n-1 carry spaces V_i; framing vertices carry D_i.
#   dim V_i = i        (i <= k)
#           = k        (k <= i <= n-k)
#           = n-i      (i >= n-k)
#   dim D_i = 1        at i = k and i = n-k          (k < n-k)
#           = 2        at i = k                       (n = 2k, single framing)
# Out-of-range spaces are zero (rmk:0): V_0 = V_n = 0.

"""
    QuiverData(n, k)

Combinatorial data of the framed double quiver of type `A_{n-1}` attached to the
two-row Jordan type `λ = (n-k, k)` with `1 ≤ k ≤ n-k`.

Fields:
- `n`, `k`               : the shape is `(n-k, k)`; `∑ = n`.
- `v :: Vector{Int}`     : dimension vector `(dim V_1, …, dim V_{n-1})`.
- `d :: Vector{Int}`     : framing dimension vector `(dim D_1, …, dim D_{n-1})`.
- `framing :: Vector{Int}`   : the vertices carrying a nonzero `D_i` (sorted).
- `equal :: Bool`        : whether `n = 2k` (single dim-2 framing at vertex `k`).
"""
struct QuiverData
    n::Int
    k::Int
    v::Vector{Int}
    d::Vector{Int}
    framing::Vector{Int}
    equal::Bool
end

function QuiverData(n::Integer, k::Integer)
    n = Int(n); k = Int(k)
    1 <= k || throw(ArgumentError("need k ≥ 1, got k=$k"))
    2k <= n || throw(ArgumentError("need k ≤ n-k, got (n-k,k)=($(n-k),$k)"))

    v = Int[dimV(n, k, i) for i in 1:(n-1)]
    d = zeros(Int, n - 1)
    equal = (n == 2k)
    if equal
        d[k] = 2
        framing = [k]
    else
        d[k] = 1
        d[n-k] = 1
        framing = [k, n - k]
    end
    return QuiverData(n, k, v, d, framing, equal)
end

"`dimV(n, k, i)` — dimension of `V_i`; 0 if `i` is out of range `1..n-1`."
function dimV(n::Integer, k::Integer, i::Integer)
    (1 <= i <= n - 1) || return 0
    if i <= k
        return i
    elseif i <= n - k
        return k
    else
        return n - i
    end
end

"`dimD(qd, i)` — dimension of the framing space `D_i` (0 if none)."
dimD(qd::QuiverData, i::Integer) = (1 <= i <= qd.n - 1) ? qd.d[i] : 0

"`getV(qd, i)` — dimension of `V_i`, honouring the zero conventions `V_0 = V_n = 0`."
getV(qd::QuiverData, i::Integer) = dimV(qd.n, qd.k, i)

function Base.show(io::IO, qd::QuiverData)
    println(io, "QuiverData for two-row shape (n-k, k) = ($(qd.n - qd.k), $(qd.k))")
    println(io, "  type A_{$(qd.n - 1)} framed double quiver")
    println(io, "  dim V = ", qd.v)
    println(io, "  dim D = ", qd.d, "   (framing at vertices ", qd.framing, ")")
    print(io,   "  ", qd.equal ? "equal-size case n = 2k (single 2-dim framing)" :
                                 "distinct-rows case k < n-k")
end
