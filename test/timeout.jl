using Test
using SolePostHoc

# get_using_sph
@test get_using_sph(:lumen) == "using SolePostHoc: lumen\n"

# create_tmp_dir
current_dir = pwd()
tmp_dir = create_tmp_dir()
@test tmp_dir == current_dir * "/tmp"
@test create_tmp_dir() == current_dir * "/tmp" # return the path to the already created dir
tmp_dir = create_tmp_dir(:symbol_tmp)
@test tmp_dir == current_dir * "/symbol_tmp"

