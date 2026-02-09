using CairoMakie
using DataFrames, Statistics
using Dates
using JLD2
include("plot_theme.jl")
set_theme!(figtheme)
# Data structures used in the model building and simulation
include("model_data_structures.jl")

@load "data/processed_results/outflow_data.jld2"
@load "data/processed_results/inflow_data.jld2"
L = 0.08 #m (8 cm)  # Spatial locations

# First figure overview of data outflow
fig_height = 300*3
fig = Figure(size = (1000, fig_height))
axn = Axis(fig[1:2, 1], title = "a. Nitrate",
    xlabel = "Time (days)", ylabel = "(NO₃⁻) [mmol L⁻¹]",
    yticks = 0:5e-1:2.1,
    xticks = 5:5:28
    )
axdoc = Axis(fig[1, 2], title = "b. DOC",
    xlabel = "Time (days)", ylabel = "DOC [mmol L⁻¹]",
    #yticks = 0:50:300
    )
axdic = Axis(fig[2, 2], title = "c. DIC",
    xlabel = "Time (days)", ylabel = "DIC [mmol L⁻¹]",
    #yticks = 0:50:300
    )
axso4 = Axis(fig[3, 2], title = "e. Sulphate",
    xlabel = "Time (days)", ylabel = "(SO₄²⁻) [mmol L⁻¹]",
    #yticks = 0:2e-1:1.1
    )
axno2 = Axis(fig[3, 1], title = "d. Nitrite",
    xlabel = "Time (days)", ylabel = "(NO₂⁻) [mmol L⁻¹]",
    #yticks = 0:2e-1:1.1
    )


# loop into all columns and plot the data in each figure/axis
for c in 1:4
    # check the outflow data
    col = all_ds[c]

    # Plot results
    if c == 4
        label = "Column $c (control)"
    else
        label = "Column $c"
    end

    # plotting outflow data.
    scatter!(axn, col.no3.t ./ (24*60*60), col.no3.conc, label = label, color = colors[c], markersize = 8, marker = markers[c])
    scatter!(axdoc, col.doc.t ./ (24*60*60), col.doc.conc, label = label, color = colors[c], markersize = 8, marker = markers[c])
    lines!(axdoc, col.doc.t ./ (24*60*60), col.doc.conc, label = label, color = colors[c], linestyle = :dash)
    scatter!(axdic, col.dic.t ./ (24*60*60), col.dic.conc, label = label, color = colors[c], markersize = 8, marker = markers[c])
    lines!(axdic, col.dic.t ./ (24*60*60), col.dic.conc, label = label, color = colors[c], linestyle = :dash)
    scatter!(axso4, col.so4.t ./ (24*60*60), col.so4.conc, label = label, color = colors[c], markersize = 8, marker = markers[c])
    lines!(axso4, col.so4.t ./ (24*60*60), col.so4.conc, label = label, color = colors[c], linestyle = :dash)
    scatter!(axno2, col.no2.t ./ (24*60*60), col.no2.conc, label = label, color = colors[c], markersize = 8, marker = markers[c])
    lines!(axno2, col.no2.t ./ (24*60*60), col.no2.conc, label = label, color = colors[c], linestyle = :dash)
end


# Finally, we Plot the inflow concentration
c_in_plot = []
c_indata = c_ins[1] # using column 1 inflow data (all are very similar)

plot_t = 0:0.1:27
for t in plot_t
    # find the last switch time before t
    idx = findlast(c_indata.t_in .<= t*24*60*60)
    if idx === nothing
        c_in_loc = c_indata.c_in[1][1]*1e3
    else
        c_in_loc = c_indata.c_in[idx+1][1]*1e3
    end
    push!(c_in_plot, c_in_loc)
end
lines!(axn, plot_t, c_in_plot, label = "Inflow concentration", color = :black, linestyle = :dash)
lines!(axn, plot_t, c_in_plot.-1, label = "Inflow concentration - 1 mM",
    color = :gray, linestyle = :dot, linewidth = 2)
# add some vertical lines to indicate the -1 mM shifts
idx = 1
for t_shift in [6.0, 16.0, 22.0, 26.0]
    lines!(axn, [t_shift, t_shift], [c_indata.c_in[idx+1][1]*1e3-1, c_indata.c_in[idx+1][1]*1e3],
        color = :gray, linestyle = :dot, linewidth = 2)
    # add a text next to the line
    text!(axn, t_shift, c_indata.c_in[idx+1][1]*1e3-1+0.3, text = "-1 mM", align = (:left, :bottom),
        fontsize = 12, rotation = π/2)
    global idx += 1
end

# inflow concentration for the remaining plots
lines!(axdoc, plot_t, zeros(length(plot_t)), label = "Inflow concentration", color = :black, linestyle = :dash)
lines!(axdic, plot_t, ones(length(plot_t)).*30/12, label = "Inflow concentration", color = :black, linestyle = :dash)
lines!(axso4, plot_t, ones(length(plot_t)).*0.74/96, label = "Inflow concentration", color = :black, linestyle = :dash)
lines!(axno2, plot_t, zeros(length(plot_t)), label = "Inflow concentration", color = :black, linestyle = :dash)
xlims!(axn, (4.0, 27.0))
xlims!(axdoc, (4.0, 27.0))
xlims!(axdic, (4.0, 27.0))
xlims!(axso4, (4.0, 27.0))
xlims!(axno2, (4.0, 27.0))

Legend(fig[4, :], axn, framevisible=false, merge=true, orientation = :horizontal)
# linkxaxes!(axn, axf, axs)
resize_to_layout!(fig)
fig
save("figs/outflow_plot.png", fig, px_per_unit = 2.0)
