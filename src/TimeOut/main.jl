module TimeOut

"""
    get_using_sph(algo::Symbol)::String

Generate a Julia `using` statement for importing a specific algorithm from SolePostHoc.

This function creates a properly formatted import statement that will be injected into 
the `script_file` string when spawning a new Julia subprocess.

# Arguments
- `algo::Symbol`: The name of the algorithm/function to import from SolePostHoc

# Returns
- `String`: A formatted using statement in the form "using SolePostHoc: algo\\n"
"""
get_using_sph(algo::Symbol)::String = "using SolePostHoc: $algo\n"

"""
This constructor is specific for the [`SolePostHoc.jl`](https://github.com/aclai-lab/SolePostHoc.jl) package.

The Julia environment in which the process will be launched must have
the dependent packages pre-installed.
It is assumed that whoever uses this module will already have a properly configured
environment with all the necessary dependencies.

# Necessary dependencies:
- [`JLD2.jl`](https://github.com/JuliaIO/JLD2.jl)
- [`SolePostHoc.jl`](https://github.com/aclai-lab/SolePostHoc.jl)
- [`DecisionTree.jl`](https://github.com/JuliaAI/DecisionTree.jl)
"""
function sph_script_builder()
    base_script = "using JLD2\n"
end

macro run_with_timeout(expr, model, timeout_sec, kwargs)
    return quote
        local sph_expr  = $(string(expr))
        local sph_model = $(esc(model))
        local timeout   = $(esc(timeout_sec))
        local kwargs    = $(esc(kwargs))

        # create tmp directory and files
        local current_dir = pwd()
        local tmp_dir     = joinpath(current_dir, "tmp")
        isdir(tmp_dir) || mkdir(tmp_dir)
        
        local model_file  = joinpath(tmp_dir, basename(tempname()) * ".jld2")
        local result_file = joinpath(tmp_dir, basename(tempname()) * ".jld2")
        
        local result    = nothing
        local timed_out = false
        
        try
            # save the model to tmp file
            JLD2.jldsave(model_file; model=sph_model)
            
            # subprocess script
            local script_file = """
                using JLD2, DecisionTree
                using SolePostHoc: Lumen
                
                # load model
                model = JLD2.load("$model_file", "model")
                
                # execute PostHoc algo
                result = $sph_expr(model; $kwargs...)
                
                # save result
                JLD2.jldsave("$result_file"; result=result)
            """
            
            # start the subprocess
            local proc = run(`julia -e $script_file`, wait=false)
            # for debug only
            # local proc = run(pipeline(`julia -e $script_file`, stdout=stdout, stderr=stderr), wait=false)
            
            # create timeout timer
            local timer = Timer(timeout) do t
                if process_running(proc)
                    @warn "Process exceeded timeout of $timeout seconds, killing..."
                    kill(proc, Base.SIGTERM)
                    sleep(0.5)
                    process_running(proc) && kill(proc, Base.SIGKILL)
                end
            end
            
            # wait for process to complete
            wait(proc)
            close(timer)
            
            # check if result file exists and load it
            if isfile(result_file)
                try
                    result = JLD2.load(result_file, "result")
                catch e
                    @warn "Error loading result: $e"
                    timed_out = true
                end
            else
                @warn "Result file not found, something went wrong in subprocess."
                timed_out = true
            end
            
        catch e
            @warn "Error in @run_with_timeout: $e"
            timed_out = true
        finally
            # cleanup temporary files
            rm(model_file,  force=true)
            rm(result_file, force=true)
        end
        
        (result, timed_out)
    end
end

end