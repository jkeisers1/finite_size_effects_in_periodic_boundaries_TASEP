using Plots
using JLD
using Pkg
using Statistics

include("AnalyticalFunctions.jl")
include("interactive_plotting_helpers.jl")
#include("package_dependencies.jl")


@load "data/c_data_250.jld"
@load "data/J_data_250.jld"
@load "data/ρ_data_250.jld"
@load "data/t_data_250.jld"

J250 = J 

@load "data/c_data_500.jld"
@load "data/J_data_500.jld"
@load "data/ρ_data_500.jld"
@load "data/t_data_500.jld"

J500 = J 

@load "data/c_data_750.jld"
@load "data/J_data_750.jld"
@load "data/ρ_data_750.jld"
@load "data/t_data_750.jld"

data = J

J750 = J

J[1,1]
ρ_list = collect(range(0.01, stop= 0.99, length = 50))
k₊_list = [10^-6, 10^-5, 10^-4, 10^-3, 10^-2, 10^-1, 10^0, 10]
k₋_list = [10^-6, 10^-5, 10^-4, 10^-3, 10^-2, 10^-1, 10^0, 10]
γ = 1

data
k₋_list, k₊_list
colors = palette(:tab20)
L = 250

titlefontsize=14 
labelfontsize=16 
tickfontsize=12
markersize=6.5
linewidth=5
legendfontsize= 10
xlims=nothing
ylims=(0,0.1)
shapes = [:circle, :square, :diamond, :star5, :cross, :utriangle]

color_palette = [
    RGBA(0/255, 0/255, 255/255, 1),      # Dark Blue
    RGBA(255/255, 128/255, 0/255, 1),    # Dark Orange
    RGBA(0/255, 128/255, 0/255, 1),      # Dark Green
    RGBA(255/255, 255/255, 0/255, 1),    # Yellow
    RGBA(255/255, 0/255, 255/255, 1),    # Pink
    RGBA(0/255, 255/255, 255/255, 1),    # Light Blue
    RGBA(128/255, 0/255, 0/255, 1)       # Dark Red
]

plot_current_whoosh_cluster(J500, ρ_list, 500, error = true)

