module UrbanEpis


using Random

using DataFrames

using Graphs
#using StatsBase

export draw_people_same_tile, extract_tile_data

include("city.jl")
include("plot_utils.jl")

export draw_misinformed, find_fraction_infected, draw_misinformed_tile

end
