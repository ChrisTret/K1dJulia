include("k1d_functionsMod.jl")
using K1dUniv
using K1dBiv

# univariate
# hand-calculated 4/3, 8/3. 32/3
my_test = [1, 6, 12, 15, 17, 24]
k1d_univ(my_test, [3, 4, 10])


# bivariate
test1 = [1, 7, 8, 14, 15, 21, 27, 28]
test2 = [4, 5, 10, 19, 23, 24, 25, 30]
# K_12 = K_21 = 2.8125 6.5625 15.0000
println(k1d_biv(test1, test2, [3, 4, 10]))
println(k1d_biv(test2, test1, [3, 4, 10]))