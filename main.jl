include("RipleysK1DMod.jl")

using .RipleysK1D

T = [3, 4, 10]
my_test = [1, 6, 12, 15, 17, 24]
K̂ = k1d_univ(my_test, T)

