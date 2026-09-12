# function apply(
#     f::LumenEnsemble{R,T},
#     d::Matrix{T}
# ) where {R<:Unsigned,T<:AbstractFloat}
#     n = size(d, 1)
#     preds = Vector{R}(undef, n)

#     @inbounds for i in 1:n
#         fill!(counts, 0)
#         for r in f.roots
#             while f.feat[r] > 0
#                 r = d[i, f.feat[r]] < f.thr[r] ? f.left[r] : f.right[r]
#             end
#             counts[f.leaf[r]] += 1
#         end
#         preds[i] = argmax(counts)
#     end

#     return preds
# end
