using Test
using SolePostHoc

# get_using_sph
@test get_using_sph(:lumen) == "using SolePostHoc: lumen\n"

# mk_tmp_dir
current_dir = pwd()
tmp_dir = mk_tmp_dir()
@test tmp_dir == current_dir * "/tmp"
@test mk_tmp_dir() == current_dir * "/tmp" # return the path to the already created dir
tmp_dir = mk_tmp_dir(:symbol_tmp)
@test tmp_dir == current_dir * "/symbol_tmp"

# rm_tmp_dir
@test rm_tmp_dir("/invalid/invalid") == false
@test rm_tmp_dir(:invalid) == false
@test rm_tmp_dir("tmp") == true
@test rm_tmp_dir(:symbol_tmp) == true
