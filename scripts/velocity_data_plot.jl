"""
    velocity_data_plot.jl

Generates a plot of the time-variable flow velocity and discharge rates for all columns.
It visualizes the piecewise constant velocity profiles derived from flow rate measurements.

# Outputs:
- `figs/flow_velocity.png`: The generated plot.
- `figs/flow_velocity.pdf`: The generated plot in PDF format.
"""

using CairoMakie
using DataFrames
using Statistics
using Dates
using JLD2
include("plot_theme.jl")
set_theme!(figtheme)

include("model_data_structures.jl")
@load "data/processed_results/flow_velocity_data.jld2"

# Re-loading tracer params to get phi for back-calc is safest if we want exact v->Q relation used in model
tracer_params = load("data/processed_results/tracer_params.jld2")["tracer_params"]
D = 3.5*1e-2 
A = π * D^2 / 4
ϕ_avg = mean([tracer_params[k][1] for k in 1:3])
k_conv = (ϕ_avg * A * 1e6 / 24)


# Second figure: overview of velocity and discharge data.
fig2 = Figure(size = (400, 300))
axv2 = Axis(fig2[1, 1], title = "Flow velocity",
    xlabel = "Time (days)", ylabel = "Velocity [m day⁻¹]",
    xticks = 5:5:28,
    xgridvisible = false,
    ygridvisible = false,
    )
axq2 = Axis(fig2[1, 1],
    ylabel = "Flow rate [mL hr⁻¹]",
    yaxisposition = :right,
    xgridvisible = false,
    ygridvisible = false,
    )
hidespines!(axq2, :l, :t, :b)
hidexdecorations!(axq2)

for c in 1:4
     # Plot results
    if c == 4
        label = "Column $c (control)"
    else
        label = "Column $c"
    end
    #flow times
    plot_flowt = 0:0.0001:27
    # discharge rate times (middle of the interval)
    avg_times = (v_da[c].start_times .+ v_da[c].end_times) ./ 2
    # flow velocity and discharge datasets
    v_c = v_da[c].v

    end_times = v_da[c].end_times
    function v_func(t)
        @inbounds for i in eachindex(v_c)
            if t <= end_times[i]
                return v_c[i]
            end
        end
        return v_c[end]
    end
    
    v_plot = [v_func(t*24*60*60) * 24*60*60 for t in plot_flowt]
    lines!(axv2, plot_flowt, v_plot, label = label, color = colors[c])
    scatter!(axv2, avg_times ./ (24*60*60), v_c .* (24*60*60), label = label, color = colors[c],
        marker=markers[c], markersize = 6)
end

linkxaxes!(axv2, axq2)
linkyaxes!(axv2, axq2)

# Set yticks for axq2 (flow rate) based on axv2 (velocity) using average porosity
q_ticks = 0:0.5:3.0
v_tick_pos = q_ticks ./ k_conv
axq2.yticks = (v_tick_pos, string.(q_ticks))

Legend(fig2[2, :], axv2, merge=true, orientation = :horizontal)

resize_to_layout!(fig2)
fig2
save("figs/flow_velocity.png", fig2, px_per_unit = 2.0)
save("figs/flow_velocity.pdf", fig2, pt_per_unit = 1)