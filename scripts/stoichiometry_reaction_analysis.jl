"""
    stoichiometry_reaction_analysis.jl

Analyzes the stoichiometry of the denitrification reaction.
It calculates the observed differences in Nitrate, DIC, and DOC concentrations between outflow and inflow (or control),
and computes stoichiometric ratios (DIC/NO3, DOC/NO3, DIC/DOC) to characterize the metabolic pathway (e.g., autotrophic vs heterotrophic).

# Outputs:
- `data/processed_results/stoichiometry_dic_doc.csv`: Table of calculated stoichiometric differences and ratios.
"""

using CairoMakie
using DataFrames, Statistics
using Dates
using JLD2
using CSV
include("plot_theme.jl")
set_theme!(figtheme)
# Load the outflow data and calculate the average nitrate decrease the average DOC difference,
# and the average DIC difference.

# calculate the average outflow - inflow difference for NO3-
include("model_data_structures.jl")

@load "data/processed_results/outflow_data.jld2"
@load "data/processed_results/inflow_data.jld2"
L = 0.08 #m (8 cm)  # Spatial locations

# include results from lagrangian analysis
df = CSV.read("data/processed_results/no3_rates.csv", DataFrame)

# Column 4 data: DIC, DOC data for comparison
dic4 = all_ds[4].dic.conc
doc4= all_ds[4].doc.conc
meddic = median(dic4)
meddoc = median(doc4)

# df of results
df_results = DataFrame(
    "no3_diff" => [],
    "dic_diff" => [],
    "doc_diff" => [],
    "dic_no3_ratio" => [],
    "doc_no3_ratio" => [],
    "dic_doc_ratio" => [],
)

for c in 1:3
    dic = all_ds[c].dic.conc
    doc = all_ds[c].doc.conc
    # dic concentration difference
    dic_diff = abs(median(dic) - meddic) * 1e-3 # convert to mol/L
    # doc concentration difference
    doc_diff = abs(median(doc) - meddoc) * 1e-3 # convert to mol/L
    # get the nitrate concentration difference
    no3_diff = abs(df[!, :concentration_diff][c])

    #stoichiometric ratio
    dic_no3 = dic_diff / no3_diff
    doc_no3 = doc_diff / no3_diff
    dic_doc = dic_diff / doc_diff

    push!(df_results, [no3_diff*1e3, dic_diff*1e3, doc_diff*1e3, dic_no3, doc_no3, dic_doc])


end

display(df_results)

CSV.write("data/processed_results/stoichiometry_dic_doc.csv", df_results)