# This Trit implementation was created with the help of Claude (claude.ai)
struct TritVector
    # Store all values in a single bit vector
    # Each trit takes 2 consecutive bits
    # Encoding: 00 -> 0, 01 -> 1, 10 -> -1
    bits::BitVector
    length::Int64

    function TritVector(len::Integer)
        # We need 2 bits per value
        new(falses(2 * len), len)
    end
end

# Consistency checking function
function is_consistent(arr::TritVector)::Tuple{Bool,Vector{Int}}
    invalid_indices = Int[]

    for i = 1:arr.length
        base_idx = 2(i - 1) + 1
        # Check if both bits are set (invalid 11 state)
        if arr.bits[base_idx] && arr.bits[base_idx+1]
            push!(invalid_indices, i)
        end
    end

    return isempty(invalid_indices), invalid_indices
end

# Enhanced constructor with immediate consistency check
function TritVector(bits::BitVector)
    if length(bits) % 2 != 0
        throw(ArgumentError("Bit vector length must be even"))
    end

    arr = new(bits, div(length(bits), 2))
    consistent, invalid_indices = is_consistent(arr)
    if !consistent
        throw(ArgumentError("Invalid bit patterns at indices: $invalid_indices"))
    end
    return arr
end

# Verify consistency after any operation that modifies the vector
function verify_consistency!(arr::TritVector)
    consistent, invalid_indices = is_consistent(arr)
    if !consistent
        throw(ErrorException("Vector consistency violated at indices: $invalid_indices"))
    end
end

# Get the value at a specific index
function Base.getindex(arr::TritVector, i::Integer)::Int
    if i < 1 || i > arr.length
        throw(BoundsError(arr, i))
    end

    # Calculate the bit positions for this trit
    base_idx = 2(i - 1) + 1

    # Read both bits at once
    bits = (arr.bits[base_idx] << 1) | arr.bits[base_idx+1]

    if bits == 0b00
        return 0
    elseif bits == 0b01
        return 1
    elseif bits == 0b10
        return -1
    else
        error("Invalid bit pattern detected at index $i! Both bits are set!")
    end
end

# Set the value at a specific index
function Base.setindex!(arr::TritVector, value::Integer, i::Integer)
    if i < 1 || i > arr.length
        throw(BoundsError(arr, i))
    end

    if !(value in -1:1)
        throw(ArgumentError("Value must be -1, 0, or 1"))
    end

    # Calculate the bit positions for this trit
    base_idx = 2(i - 1) + 1

    if value == 0
        arr.bits[base_idx] = false
        arr.bits[base_idx+1] = false
    elseif value == 1
        arr.bits[base_idx] = false
        arr.bits[base_idx+1] = true
    else  # value == -1
        arr.bits[base_idx] = true
        arr.bits[base_idx+1] = false
    end
end

# Bulk operations for better performance
function Base.fill!(arr::TritVector, value::Integer)
    if !(value in -1:1)
        throw(ArgumentError("Value must be -1, 0, or 1"))
    end

    if value == 0
        fill!(arr.bits, false)
    elseif value == 1
        for i = 1:arr.length
            base_idx = 2(i - 1) + 1
            arr.bits[base_idx] = false
            arr.bits[base_idx+1] = true
        end
    else  # value == -1
        for i = 1:arr.length
            base_idx = 2(i - 1) + 1
            arr.bits[base_idx] = true
            arr.bits[base_idx+1] = false
        end
    end
    verify_consistency!(arr)
    return arr
end

# String representation
function Base.show(io::IO, arr::TritVector)
    consistent, invalid_indices = is_consistent(arr)
    if !consistent
        print(io, "WARNING: Inconsistent TritVector at indices $invalid_indices: [")
    else
        print(io, "TritVector([")
    end
    for i = 1:arr.length
        i > 1 && print(io, ", ")
        if i in invalid_indices
            print(io, "ERROR")
        else
            print(io, arr[i])
        end
    end
    print(io, "])")
end


# Length of the vector
Base.length(arr::TritVector) = arr.length

# Memory size in bytes
Base.sizeof(arr::TritVector) = ceil(Int, length(arr.bits) / 8)

# Iterator interface
Base.iterate(arr::TritVector, state = 1) =
    state > arr.length ? nothing : (arr[state], state + 1)
Base.eltype(::Type{TritVector}) = Int


# Copy constructor
function Base.copy(arr::TritVector)
    new_arr = TritVector(arr.length)
    new_arr.bits .= arr.bits
    verify_consistency!(new_arr)
    return new_arr
end

# Arithmetic operations
function Base.:(+)(a::TritVector, b::TritVector)
    if length(a) != length(b)
        throw(DimensionMismatch("Vectors have different lengths"))
    end
    result = TritVector(length(a))
    for i = 1:length(a)
        # Clamp sum to [-1, 1]
        sum = clamp(a[i] + b[i], -1, 1)
        result[i] = sum
    end
    verify_consistency!(result)
    return result
end

# Multiplication
function Base.:(*)(a::TritVector, b::TritVector)
    if length(a) != length(b)
        throw(DimensionMismatch("Vectors have different lengths"))
    end
    result = TritVector(length(a))
    for i = 1:length(a)
        result[i] = clamp(a[i] * b[i], -1, 1)
    end
    verify_consistency!(result)
    return result
end

# Negation
function Base.:(-)(arr::TritVector)
    result = TritVector(length(arr))
    for i = 1:length(arr)
        result[i] = -arr[i]
    end
    return result
end

# Aggiungere questo dopo le altre definizioni di Base.
import Base: isless

function Base.isless(a::TritVector, b::TritVector)
    if length(a) != length(b)
        throw(DimensionMismatch("Vectors have different lengths"))
    end
    
    # Prima confronta per numero di elementi non-zero
    a_nonzero = count(x -> x != 0, a)
    b_nonzero = count(x -> x != 0, b)
    
    if a_nonzero != b_nonzero
        return a_nonzero < b_nonzero
    end
    
    # Se hanno lo stesso numero di elementi non-zero, confronta elemento per elemento
    for i in 1:length(a)
        if a[i] != b[i]
            # Definisce un ordine: -1 < 0 < 1
            return a[i] < b[i]
        end
    end
    
    # Se tutti gli elementi sono uguali, i vettori sono uguali
    return false
end

# Per completezza, aggiungiamo anche l'operatore ==
function Base.:(==)(a::TritVector, b::TritVector)
    length(a) == length(b) || return false
    for i in 1:length(a)
        a[i] == b[i] || return false
    end
    return true
end