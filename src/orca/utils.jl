# ─────────────────────────────────────────────────────────────────────────────
# PRUNING HELPERS  
# ─────────────────────────────────────────────────────────────────────────────
"""
    tree_depth(leaf)

Base case: return 1 as the depth of the node.
Used by the recursive `tree_depth(node)` to get the depth of the entire tree
"""

function tree_depth(leaf::ConstantModel)
    return 1
end

"""
    tree_depth(node)

Recursively visits every child node starting from the starting node's children. 
"""

function tree_depth(node::Branch)
    return 1 + max(tree_depth(posconsequent(node)), tree_depth(negconsequent(node)))
end

"""
    get_majority_class(leaf)

Base case: return the class values stored in a leaf node.
Used by the recursive `get_majority_class(node)` to collect all leaf values
in a subtree.
"""
function get_majority_class(leaf::ConstantModel)
    return [outcome(leaf)] 
end

"""
    get_majority_class(node)

Recursively collect all leaf values under an internal node and return the most
frequent class (majority vote).

This is used when a subtree needs to be replaced by a single synthetic leaf:
the new leaf will predict the class that was most common among all leaves below
the cut point.
"""
function get_majority_class(node::Branch)
    left_val  = get_majority_class(posconsequent(node))
    right_val = get_majority_class(negconsequent(node))
    return vcat(left_val, right_val)
end

"""
    prune_tree(leaf, curr_depth, max_depth)

Base case for tree pruning: a leaf is already a terminal node, so it is
returned unchanged regardless of depth.
"""
function prune_tree(leaf::ConstantModel, curr_depth::Int, max_depth::Int)
    # Leaves have no children, nothing to prune
    return leaf
end

"""
    prune_tree(node::Branch, curr_depth, max_depth)

Recursively prune an internal node so that no node beyond `max_depth` has
children. When the recursion reaches depth `max_depth`, the node is replaced
by a synthetic leaf labelled with the majority class of all leaves below it.

# Arguments
- `node`        : The current `DecisionTree.Node` being visited.
- `curr_depth`  : Depth of `node` in the original tree (root = 1).
- `max_depth`   : Maximum allowed depth; subtrees below this level are cut off.

# Returns
Either the original node (possibly with pruned children) or a new `leaf::ConstantModel`
that collapses the subtree into a single prediction.
"""
function prune_tree(node::Branch, curr_depth::Int, max_depth::Int)

    if curr_depth >= max_depth
        # We have reached the pruning boundary.
        # Find the most common class in all leaves below this node ...
        all_classes = get_majority_class(node)
        majority_value = mode(all_classes)

        # ... and replace the entire subtree with a single leaf.
        return ConstantModel(majority_value)
    else
        # Still above the cut depth: recurse into both children,
        # incrementing the depth counter at each level.
        new_left  = prune_tree(posconsequent(node),  curr_depth + 1, max_depth)
        new_right = prune_tree(negconsequent(node), curr_depth + 1, max_depth)

        # Reconstruct the current node with (potentially pruned) children.
        # The split condition (feature id + threshold) is preserved unchanged.
        return Branch(antecedent(node), new_left, new_right, info(node))
    end
end

#-------------------------------------------------------------------

function extract_feature_value(formula)
    
    atom = tree(formula)
    cond_obj = atom.value
    
    feat_id = cond_obj.metacond.feature 
    
    feat_val = cond_obj.threshold    
    
    return feat_id, feat_val
end


function update_formula_value(formula, new_value)
    atom = tree(formula)
    cond_obj = atom.value

    new_condition = SoleData.ScalarCondition(cond_obj.metacond, new_value)

    return typeof(atom)(new_condition)
end



"""
    collect_alphabet!(leaf::ConstantModel, alphabet_table)

Base case: return nothing.
Used by the recursive `collect_alphabet!(node::Branch, alphabet_table)`
"""
function collect_alphabet!(leaf::ConstantModel, alphabet_table)

    # If the passed node is a leaf i have nothing to collect
    return nothing
end

"""
    collect_alphabet!(node::Branch, alphabet_table)

Recursively collects the features and values contained in each node
of the tree passed tree.

# Arguments
- `node`            : The current `DecisionTree.Node` being visited.
- `alphabet_table`  : A dictionary that maps each unique feature to a
                      list of all the values they operates with.
"""
function collect_alphabet!(node::Branch, alphabet_table)
    # I identify the current node's feature and value
    form = antecedent(node)
    feat_id, feat_val = extract_feature_value(form)

    # I check if the current feature is already present
    # in the alphabet table
    if !haskey(alphabet_table, feat_id)

        # If the feature is not already present
        # I insert it in the table, creating a new row
        alphabet_table[feat_id] = Float64[]
    end

    # I insert the value of the current node in its
    # feature's node.
    push!(alphabet_table[feat_id], feat_val)

    # I recursively do the same thing with the
    # current's node children
    collect_alphabet!(posconsequent(node), alphabet_table)
    collect_alphabet!(negconsequent(node), alphabet_table)

end

"""
    modify_alphabet(leaf::ConstantModel, alphabet_table)

Base case: a leaf cannot be modified so it is returned as is.
Used by the recursive `modify_alphabet(node::Branch, alphabet_table)`
"""
function modify_alphabet(leaf::ConstantModel, alphabet_table)

    # If the passed node is a leaf there is nothing
    # to modify, so I return the leaf.
    return leaf
end

"""
    modify_alphabet(node::Branch, alphabet_table)

Recursively

# Arguments
- `node`            : The current `DecisionTree.Node` being visited.
- `alphabet_table`  : A dictionary that maps each unique feature to a
                      list of all the values they operates with.

# Returns
Either a modified version of the original tree, where one or more nodes
has been changed, wheter in its feature or the value the feature operates with
"""
function modify_alphabet(node::Branch, alphabet_table)

    # I identify the current node's feature and value
    form = antecedent(node)
    feat_id, feat_val = extract_feature_value(form)

    # Case 1: the feature i'm looking for is not in the table
    # so i prune the tree to this level
    if !haskey(alphabet_table, feat_id)
        all_classes = get_majority_class(node)
        majority_value = mode(all_classes)

        return ConstantModel(majority_value)
    end

    # If i'm here it means that the feature i'm looking for is
    # in the table

    # This contains all the values the analized feature operates with
    allowed_values = alphabet_table[feat_id]

    #= Case 2 and 3:

        Case 2) The value the current feature is operating with is NOT
        in the feature's row in the alphabet table.

        I replace the value in the node with the closest value i have
        in the feature's row in the alphabet table
        --------------------------------------------------------------
        Case 3) The value the current feature is operating with is in
        the feature's row in the alphabet table.

        I do nothing
    =#

        #=
         I subtract the node's value from all values in allowed_values,
         to have it become an array with the differences beetwen the
         node's value and all the previous values.

         Then i take the absolute value to avoid any shenanigan.

         Finally i use the argmin() function to select the index
         of the minimum element in the array, that is the index
         of the element with the least difference from the node's value.

         That means that if this is a case 2, than we take the closest
         element and replace the original one with it.

         If this is a case 3 then the index will point to the actual
         element in the allowed_values, due to it having difference = 0
        =#
        closest_idx = argmin(abs.(allowed_values .- feat_val))
        new_feat_val = allowed_values[closest_idx]

    new_form = update_formula_value(form, new_feat_val)

    #  I recursively do the same thing with the children of this node
    left_child = modify_alphabet(posconsequent(node), alphabet_table)
    right_child = modify_alphabet(negconsequent(node), alphabet_table)

    # I recreate the node with the new value i found in the table
    return Branch(new_form, left_child, right_child, info(node))
end