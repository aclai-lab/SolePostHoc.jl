# test su concatenazione stringhe

script_test = """
    using JLD2, DecisionTree
"""
test = """
    using SolePostHoc: Lumen
"""
@btime script_test * test
# 25.041 ns (1 allocation: 80 bytes)

using StaticStrings
sscript_test = static"""
    using JLD2, DecisionTree
"""
stest = static"""
    using SolePostHoc: Lumen
"""
@btime sscript_test * stest
# 314.873 ns (5 allocations: 288 bytes)

# https://discourse.julialang.org/t/avoid-allocations-for-string-concatenation/104666/2
usage = StaticString((Tuple(sscript_test)..., Tuple(stest)...))
@btime StaticString((Tuple(sscript_test)..., Tuple(stest)...))
# 1.635 μs (4 allocations: 256 bytes)

using InlineStrings
iscript_test = String31("""
    using JLD2, DecisionTree
""")
itest = String31("""
    using SolePostHoc: Lumen
""")
@btime iscript_test * itest
# 28.817 ns (1 allocation: 80 bytes)

using StaticTools
sscript_test = MallocString("""
    using JLD2, DecisionTree
""")
stest = MallocString("""
    using SolePostHoc: Lumen
""")
@btime sscript_test * stest
# 37.024 ns (1 allocation: 32 bytes)