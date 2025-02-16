using DelimitedFiles

function load_data(name)
    data_path = "/home/lele7x/results/evolutionary-rule-extraction/datasets/"

    if name == "breast"
        breast = DelimitedFiles.readdlm(joinpath(data_path, "breast-cancer-wisconsin.csv"), ',')
        (X, idx_deleted) = begin
            IX = breast[:, 1:10]

            X = nothing
            idx_deleted = []
            for i in 1:size(IX,1)
                row = breast[i,:]
                if !any(row .== "?")
                    X = isnothing(X) ? row' : [X; row']
                else
                    push!(idx_deleted,i)
                end
            end

            (X,idx_deleted)
        end
        Y = breast[setdiff(1:size(breast,1),idx_deleted), 11]
        return X, Y
    end

    error("Unknown dataset")
end

X, Y = load_data("breast")
