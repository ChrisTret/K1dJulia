module RipleysK1D

using CSV, DataFrames, Plots

include("data_wranglingMod.jl")
include("k1d_univMod.jl")
include("k1d_bivariateMod.jl")
include("k1d_speed_testsMod.jl")
include("k_plotsMod.jl")
using .DataWrangling, .K1dUniv, .K1dBiv, .K1dSpeedTests, .K1dPlots

export process_genome_data, k1d_univ, k1d_biv, k1d_all_comparisons, speed_comparison_univ, k_plot

end # end of module
