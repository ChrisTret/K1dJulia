module RipleysK1D

using Plots

include("k1d_univ_Mod.jl")
include("k1d_bivariateMod.jl")
include("k1d_speed_tests_Mod.jl")
include("k_plotsMod.jl")
using .K1dUniv, .K1dBiv, .K1dSpeedTests, .K1dPlots

export k1d_univ, k1d_biv, speed_comparison_univ, k_plot

end # end of module
