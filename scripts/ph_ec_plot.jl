"""
    ph_ec_plot.jl

This script generates a plot of the pH and EC (Electrical Conductivity) over time for all columns.

# Outputs:
- `figs/ph_ec_plot.png`: The generated plot.
- `figs/ph_ec_plot.pdf`: The generated plot in PDF format.
"""

using CairoMakie
using DataFrames, Statistics
using Dates
using JLD2
include("plot_theme.jl")
set_theme!(figtheme)
# Data structures used in the model building and simulation
include("model_data_structures.jl")

@load "data/processed_results/outflow_data.jld2"

fig = Figure(size = (238, 350))
axph = Axis(fig[1, 1], title = "a. pH",
    xlabel = "Time (days)", ylabel = "pH",
    xticks = 5:5:28
    )
axec = Axis(fig[2, 1], title = "b. EC",
    xlabel = "Time (days)", ylabel = "EC [μS/cm]",
    xticks = 5:5:28
    )

# loop into all columns and plot the data in each figure/axis
for c in 1:4
    col = all_ds[c]

    if c == 4
        label = "Column $c (control)"
    else
        label = "Column $c"
    end

    scatter!(axph, col.pH.t ./ (24*60*60), col.pH.conc, label = label, color = colors[c], markersize = 8, marker = markers[c])
    lines!(axph, col.pH.t ./ (24*60*60), col.pH.conc, color = colors[c], linestyle = :dash)
    
    scatter!(axec, col.ec.t ./ (24*60*60), col.ec.conc, label = label, color = colors[c], markersize = 8, marker = markers[c])
    lines!(axec, col.ec.t ./ (24*60*60), col.ec.conc, color = colors[c], linestyle = :dash)
end

xlims!(axph, (4.0, 27.0))
xlims!(axec, (4.0, 27.0))

Legend(fig[3, :], axph, framevisible=false, merge=true, orientation = :horizontal, nbanks=2)

save("figs/ph_ec_plot.png", fig, px_per_unit = 8.33)
save("figs/ph_ec_plot.pdf", fig, pt_per_unit = 1)

println("pH and EC plot generated.")
