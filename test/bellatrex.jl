pwd()
Pkg.activate()
Pkg.add("PyCall")

using PyCall
@pyimport bellatrex
using SolePostHoc


const req_py_pkgs = ["bellatrex"]
const fs = PyNULL()

function __init__()
    pypkgs = getindex.(PyCall.Conda.parseconda(`list`, PyCall.Conda.ROOTENV), "name")
    needinstall = !all(p -> in(p, pypkgs), req_py_pkgs)

    if (needinstall)
        PyCall.Conda.pip_interop(true, PyCall.Conda.ROOTENV)
        PyCall.Conda.add(req_py_pkgs)
    end

    copy!(fs, pyimport_conda("bellatrex.BellatrexExplain", "bellatrex"))
end
