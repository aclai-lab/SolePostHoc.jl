struct BalancedTritVector
    # We'll pack 3 trits into every 5 bits
    bits::BitVector
    length::Int64

    function BalancedTritVector(len::Integer)
        # Calculate number of bits needed:
        # Each group of 3 trits needs 5 bits
        # Round up to nearest group of 3 for the total length
        num_groups = ceil(Int, len / 3)
        num_bits = num_groups * 5
        new(falses(num_bits), len)
    end
end

# Convert 3 trits to their 5-bit compressed representation
function compress_trits(t1::Int, t2::Int, t3::Int)::UInt8
    # Convert from {-1,0,1} to {0,1,2} for easier arithmetic
    v1 = t1 + 1
    v2 = t2 + 1
    v3 = t3 + 1

    # Compute base-3 number: t1*3² + t2*3¹ + t3*3⁰
    value = v1 * 9 + v2 * 3 + v3

    return UInt8(value)
end

# Convert 5-bit compressed representation back to 3 trits
function decompress_bits(bits::Integer)::Tuple{Int,Int,Int}
    # Convert to integer value
    value = Int(bits)

    # Extract each trit using base-3 arithmetic
    t3 = value % 3
    value ÷= 3
    t2 = value % 3
    t1 = value ÷ 3

    # Convert back to {-1,0,1} representation
    return (t1 - 1, t2 - 1, t3 - 1)
end
# Get the value at a specific index
function Base.getindex(arr::BalancedTritVector, i::Integer)::Int
    if i < 1 || i > arr.length
        throw(BoundsError(arr, i))
    end

    # Find which group of 3 trits contains our value
    group_idx = ceil(Int, i / 3)
    position_in_group = ((i - 1) % 3) + 1

    # Extract 5 bits for this group
    base_idx = (group_idx - 1) * 5 + 1
    group_bits = UInt8(0)
    for j = 0:4
        if base_idx + j <= length(arr.bits)
            group_bits |= (arr.bits[base_idx+j] << j)
        end
    end

    # Decompress and return the requested trit
    t1, t2, t3 = decompress_bits(group_bits)
    return [t1, t2, t3][position_in_group]
end

# Set the value at a specific index
function Base.setindex!(arr::BalancedTritVector, value::Integer, i::Integer)
    if i < 1 || i > arr.length
        throw(BoundsError(arr, i))
    end

    if !(value in -1:1)
        throw(ArgumentError("Value must be -1, 0, or 1"))
    end

    # Find which group of 3 trits contains our value
    group_idx = ceil(Int, i / 3)
    position_in_group = ((i - 1) % 3) + 1

    # Get current values for this group
    base_idx = (group_idx - 1) * 5 + 1
    group_bits = UInt8(0)
    for j = 0:4
        if base_idx + j <= length(arr.bits)
            group_bits |= (arr.bits[base_idx+j] << j)
        end
    end

    # Get current trits
    t1, t2, t3 = decompress_bits(group_bits)
    trits = [t1, t2, t3]

    # Update the requested trit
    trits[position_in_group] = value

    # Compress back to bits
    new_bits = compress_trits(trits[1], trits[2], trits[3])

    # Update the bits in the array
    for j = 0:4
        if base_idx + j <= length(arr.bits)
            arr.bits[base_idx+j] = ((new_bits >> j) & 1) == 1
        end
    end
end

# Memory size in bytes (actual size used)
function Base.sizeof(arr::BalancedTritVector)
    num_groups = ceil(Int, arr.length / 3)
    num_bits = num_groups * 5
    return ceil(Int, num_bits / 8)
end

# String representation
function Base.show(io::IO, arr::BalancedTritVector)
    print(io, "BalancedTritVector([")
    for i = 1:arr.length
        i > 1 && print(io, ", ")
        print(io, arr[i])
    end
    print(io, "])")
end

# Iterator interface
Base.iterate(arr::BalancedTritVector, state = 1) =
    state > arr.length ? nothing : (arr[state], state + 1)
Base.eltype(::Type{BalancedTritVector}) = Int
Base.length(arr::BalancedTritVector) = arr.length

# Copy constructor
function Base.copy(arr::BalancedTritVector)
    new_arr = BalancedTritVector(arr.length)
    new_arr.bits .= arr.bits
    return new_arr
end

# Arithmetic operations
function Base.:(+)(a::BalancedTritVector, b::BalancedTritVector)
    if length(a) != length(b)
        throw(DimensionMismatch("Arrays have different lengths"))
    end
    result = BalancedTritVector(length(a))
    for i = 1:length(a)
        result[i] = clamp(a[i] + b[i], -1, 1)
    end
    return result
end

function Base.:(*)(a::BalancedTritVector, b::BalancedTritVector)
    if length(a) != length(b)
        throw(DimensionMismatch("Arrays have different lengths"))
    end
    result = BalancedTritVector(length(a))
    for i = 1:length(a)
        result[i] = clamp(a[i] * b[i], -1, 1)
    end
    return result
end

function Base.:(-)(arr::BalancedTritVector)
    result = BalancedTritVector(length(arr))
    for i = 1:length(arr)
        result[i] = -arr[i]
    end
    return result
end
