using CairoMakie
using DataFrames
using Dates
using JLD2
include("plot_theme.jl")
set_theme!(figtheme)
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
hidespines!(axq2)
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
    # discharge rate times
    avg_times = v_da[c].end_times .- diff([0.0; v_da[c].start_times])./2
    # flow velocity and discharge datasets
    v_c = v_da[c].v
    q_c = q_disch[c].(avg_times)
    v_plot = [v_interp[c](t*24*60*60) * 24*60*60 for t in plot_flowt]
    lines!(axv2, plot_flowt, v_plot, label = label, color = colors[c])
    scatter!(axv2, avg_times ./ (24*60*60), v_c .* (24*60*60), label = label, color = colors[c], markersize = 8)
    scatter!(axq2, avg_times ./ (24*60*60), q_c .* (1e6*3600), label = label, color = colors[c], markersize = 8)
end

linkxaxes!(axv2, axq2)
Legend(fig2[2, :], axn2, merge=true, orientation = :horizontal)

resize_to_layout!(fig2)
fig2
save("figs/flow_velocity.png", fig2, px_per_unit = 2.0)