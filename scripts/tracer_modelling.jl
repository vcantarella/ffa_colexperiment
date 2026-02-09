"""
    tracer_modelling.jl

Fits the Ogata-Banks analytical solution to the experimental bromide breakthrough curves 
to estimate transport parameters: porosity (ϕ) and longitudinal dispersivity (αₗ).

# Outputs:
- `figs/breakthrough_fit.png`: Plot of the model fit vs experimental data.
- `figs/breakthrough_fit.pdf`: Plot of the model fit vs experimental data in PDF format.
- `data/processed_results/tracer_params.jld2`: Saved optimized parameters (ϕ, αₗ).
"""

using JLD2
using CairoMakie
using Statistics
using SpecialFunctions
using ForwardDiff
using nonlinearlstr
include("plot_theme.jl")
set_theme!(figtheme)
function ogata_banks(t, p, q, De, L, C₀)
    ϕ, αₗ = p
    v = q/ϕ
    D = De + αₗ*v
    ξ = v*t/L
    η = D/(v*L)
    C = C₀/2*erfc((1-ξ)/(2√(ξ*η)))
    return C
end

function make_problem(c_data, t_data, q, De, L, C₀)

    function residuals(p)
        c_model = ogata_banks.(t_data, Ref(p), Ref(q), Ref(De), Ref(L), Ref(C₀))
        res = c_data .- c_model
        return res
    end

    jac = (p) -> ForwardDiff.jacobian(residuals, p)

    return residuals, jac
end


# Define parametersL = 0.08 #m (8 cm)  # Spatial locations
D = 3.5*1e-2 # diameter cm to m
A = π * D^2 / 4 # Cross-sectional area
L = 0.08 #m column length
t = range(0.1, stop=24*3600, length=100)  # Time locations
# -- starting values
ϕ = 0.3 #porosity
αₗ = 8e-5 #longitudinal dispersivity (m)
De = 1e-9 #dispersion coefficient (m2/s)

# load the data and plot the results
bromide_ds = load("data/processed_results/bromide_breakthrough_data.jld2")
flow_rate_ds = load("data/processed_results/bromide_flow_rate_data.jld2")
Br_dic = bromide_ds["Br"]
time_exp = bromide_ds["avg_time"]

ps = Dict(
    1 => [ϕ, αₗ],
    2 => [ϕ, αₗ],
    3 => [ϕ, αₗ],
    # 4 => [ϕ, αₗ],
)

large_fig = Figure()
large_ax = Axis(large_fig[1, 1],
    xlabel = "time from injection (hr)",
    ylabel = "bromide concentration (mmol L⁻¹)",
    title = "tracer breakthrough fit")
colors = [:crimson, :steelblue, :forestgreen, :darkorange]

meanQ = mean(flow_rate_ds["Qs"][1])
sdQ = std(flow_rate_ds["Qs"][1])
for column in [1, 2, 3]
    global ϕ, αₗ  # Declare as global to access/modify global variables
    #column = 1
    Q_1 = flow_rate_ds["Qs"][column]
    meanQ = mean(Q_1)
    meanQ = meanQ/1e6
    # skip missing
    q = meanQ / A  # specific flow rate in m^3/s
    # v_est = q / ϕ  # Initial velocity in m/s
    c_in = 1e-3
    p_tracer = [ϕ, αₗ]
    c_data = Br_dic[column]*1e-3
    t_data = time_exp[column]
    f, jac = make_problem(c_data, t_data, q, De, L, c_in)

    sol = nonlinearlstr.lm_trust_region(f, jac, p_tracer)

    p_sol = sol[1]
    c_plot = ogata_banks.(t, Ref(p_sol), Ref(q), Ref(De), Ref(L), Ref(c_in))
    ps[column] = p_sol
    #     title = "Bromide breakthrough fit")
    lines!(large_ax, t ./ 3600, c_plot*1e3, color = colors[column], label = "Column $column Model Fit")
    scatter!(large_ax, t_data ./ 3600, c_data*1e3, color = colors[column], marker = markers[column], label = "Column $column Data")
end
axislegend(large_ax, position = :rb, "Legend")
resize_to_layout!(large_fig)
large_fig


save("figs/breakthrough_fit.png", large_fig)
save("figs/breakthrough_fit.pdf", large_fig)

# Save the parameters
save("data/processed_results/tracer_params.jld2", "tracer_params", ps)