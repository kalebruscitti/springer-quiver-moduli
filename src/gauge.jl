# The GL(V) = ∏_i GL(V_i) action (def:GL action line 1244), used to move numeric
# points to normal form and to check gauge-invariance of the equations.
#
#   g·(A,B,Γ,Δ) = ( g_{i+1} A_i g_i^{-1},  g_i B_i g_{i+1}^{-1},  g_i Γ_i,  Δ_i g_i^{-1} ).

"""
    act(g, rep) -> Rep

Act by `g = Dict{Int, <:AbstractMatrix}` (or a length `n-1` vector), where `g[i]`
is an invertible `v_i × v_i` matrix at vertex `i`. Returns a new numeric `Rep`.
"""
function act(g, rep::Rep)
    qd = rep.qd
    R = rep.R
    gg = Dict{Int,Any}()
    for i in 1:(qd.n - 1)
        gi = _gget(g, i)
        gi === nothing && continue
        M = matrix(R, getV(qd, i), getV(qd, i),
                   [R(gi[a, b]) for a in 1:getV(qd, i) for b in 1:getV(qd, i)])
        gg[i] = M
    end
    getg(i) = (1 <= i <= qd.n - 1 && haskey(gg, i)) ? gg[i] :
              identity_matrix(R, getV(qd, i))

    A = Any[getg(i + 1) * getA(rep, i) * inv(getg(i)) for i in 1:(qd.n - 1)]
    B = Any[getg(i) * getB(rep, i) * inv(getg(i + 1)) for i in 1:(qd.n - 1)]
    Gamma = Dict{Int,Any}(j => getg(j) * getGamma(rep, j) for j in qd.framing)
    Delta = Dict{Int,Any}(j => getDelta(rep, j) * inv(getg(j)) for j in qd.framing)
    return Rep(qd, R, A, B, Gamma, Delta)
end

_gget(g::AbstractDict, i) = get(g, i, nothing)
_gget(g::AbstractVector, i) = (1 <= i <= length(g)) ? g[i] : nothing
