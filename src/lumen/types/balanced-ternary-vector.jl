
# Note: this was ported from a version for Julia v0.6
# Source: https://rosettacode.org/wiki/Balanced_ternary#Julia
struct BalancedTernaryVector <: Signed
    digits::Vector{Int8}
    BalancedTernaryVector() = zero(BalancedTernaryVector)
    BalancedTernaryVector(n) = convert(BalancedTernaryVector, n)
end

const sgn2chr = Dict{Int8,Char}(-1 => '-', 0 => '0', +1 => '+')
Base.show(io::IO, bt::BalancedTernaryVector) =
    print(io, join(sgn2chr[x] for x in reverse(bt.digits)))
Base.copy(bt::BalancedTernaryVector) = BalancedTernaryVector(copy(bt.digits))
Base.zero(::Type{BalancedTernaryVector}) = BalancedTernaryVector(Int8[0])
Base.iszero(bt::BalancedTernaryVector) = bt.digits == Int8[0]
Base.convert(::Type{T}, bt::BalancedTernaryVector) where {T<:Integer} =
    sum(3^T(ex - 1) * s for (ex, s) in enumerate(bt.digits))
Int(bt::BalancedTernaryVector) = Base.convert(Int, bt)
function Base.convert(::Type{BalancedTernaryVector}, n::Signed)
    r = BalancedTernaryVector(Int8[])
    if iszero(n)
        push!(r.digits, 0)
    end
    while n != 0
        if mod(n, 3) == 0
            push!(r.digits, 0)
            n = fld(n, 3)
        elseif mod(n, 3) == 1
            push!(r.digits, 1)
            n = fld(n, 3)
        else
            push!(r.digits, -1)
            n = fld(n + 1, 3)
        end
    end
    return r
end
const chr2sgn = Dict{Char,Int8}('-' => -1, '0' => 0, '+' => 1)
function Base.convert(::Type{BalancedTernaryVector}, s::AbstractString)
    return BalancedTernaryVector(map(x -> getindex(chr2sgn, x), collect(reverse(s))))
end

macro bt_str(s)
    convert(BalancedTernaryVector, s)
end

const table = NTuple{2,Int8}[(0, -1), (1, -1), (-1, 0), (0, 0), (1, 0), (-1, 1), (0, 1)]
function _add(a::Vector{Int8}, b::Vector{Int8}, c::Int8 = Int8(0))
    if isempty(a) || isempty(b)
        if c == 0
            return isempty(a) ? b : a
        end
        return _add([c], isempty(a) ? b : a)
    else
        d, c = table[4+(isempty(a) ? 0 : a[1])+(isempty(b) ? 0 : b[1])+c]
        r = _add(a[2:end], b[2:end], c)
        if !isempty(r) || d != 0
            return pushfirst!(r, d)
        else
            return r
        end
    end
end
function Base.:+(a::BalancedTernaryVector, b::BalancedTernaryVector)
    v = _add(a.digits, b.digits)
    return isempty(v) ? BalancedTernaryVector(0) : BalancedTernaryVector(v)
end
Base.:-(bt::BalancedTernaryVector) = BalancedTernaryVector(-bt.digits)
Base.:-(a::BalancedTernaryVector, b::BalancedTernaryVector) = a + (-b)
function _mul(a::Vector{Int8}, b::Vector{Int8})
    if isempty(a) || isempty(b)
        return Int8[]
    else
        if a[1] == -1
            x = (-BalancedTernaryVector(b)).digits
        elseif a[1] == 0
            x = Int8[]
        elseif a[1] == 1
            x = b
        end
        y = append!(Int8[0], _mul(a[2:end], b))
        return _add(x, y)
    end
end
function Base.:*(a::BalancedTernaryVector, b::BalancedTernaryVector)
    v = _mul(a.digits, b.digits)
    return isempty(v) ? BalancedTernaryVector(0) : BalancedTernaryVector(v)
end
