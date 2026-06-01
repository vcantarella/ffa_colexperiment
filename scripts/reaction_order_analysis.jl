"""
    reaction_order_analysis.jl

Supporting information: reaction order discrimination for denitrification
in plug-flow sediment column experiments.

ΔC/C_in is the fractional NO₃⁻ removal: (C_in - C_out) / C_in.
Diagnostic plots follow standard PFR kinetic analysis:
  a. ΔC/C_in vs τ/C_in          (linear through origin → zero order)
  b. -ln(1 - ΔC/C_in) vs τ      (linear through origin → first order)
  c. (ΔC/C_in)/(1 - ΔC/C_in) vs τ·C_in (linear through origin → second order)

# Outputs:
- figs/reaction_order_diagnostic.pdf/png
"""

using CairoMakie
using DataInterpolations
using QuadGK
using Statistics
using JLD2
include("plot_theme.jl")
set_theme!(figtheme)

include("model_data_structures.jl")
@load "data/processed_results/outflow_data.jld2"
@load "data/processed_results/inflow_data.jld2"
@load "data/processed_results/flow_velocity_data.jld2"

function load_tracer_parameters(path::String)
    params = load(path)["tracer_params"]
    params[4] = [mean([params[k][1] for k in 1:3]),
                 mean([params[k][2] for k in 1:3])]
    return params
end

tracer_params = load_tracer_parameters("data/processed_results/tracer_params.jld2")
velocity_data, inflow_data, experimental_data = v_da, c_ins, all_ds

function make_v_func(v_da::VDataA)
    function v_inst(t)
        for i in eachindex(v_da.end_times)
            if t ≤ v_da.end_times[i]
                return v_da.v[i]
            end
        end
        return v_da.v[end]
    end
    return v_inst
end

function make_c_in_func(c_in::CinData)
    function cin(t)
        for i in eachindex(c_in.t_in)
            if t ≤ c_in.t_in[i]
                return c_in.c_in[i][1]
            end
        end
        return c_in.c_in[end][1]
    end
    return cin
end

const L_COLUMN = 0.08

# ── Collect (τ, C_in, C_out) from all columns ────────────────────────────────

all_tau  = Float64[]  # residence time [days]
all_cin  = Float64[]  # inflow concentration [mmol L⁻¹]
all_cout = Float64[]  # outflow concentration [mmol L⁻¹]
col_ids  = Int[]

for c in 1:3
    v_inst   = make_v_func(velocity_data[c])
    c_in_fun = make_c_in_func(inflow_data[c])

    X_traj(t) = quadgk(v_inst, 0, t)[1]
    dense_t = 1:(3*3600):(27*24*60*60)
    dense_x = [X_traj(t) for t in dense_t]
    T_interp = DataInterpolations.LinearInterpolation(dense_t, dense_x)
    mint = T_interp(L_COLUMN)
    τ_func(t) = t - T_interp(X_traj(t) - L_COLUMN)

    no3_exp   = experimental_data[c].no3
    valid_idx = findall(t -> t > mint, no3_exp.t)

    for i in valid_idx
        t_val    = no3_exp.t[i]
        tau_s    = τ_func(t_val)
        cin_mol  = c_in_fun(t_val - tau_s)     # mol/L
        cout_mol = no3_exp.conc[i] * 1e-3      # mmol/L → mol/L

        push!(all_tau,  tau_s / 86400.0)        # → days
        push!(all_cin,  cin_mol * 1e3)          # → mmol/L
        push!(all_cout, cout_mol * 1e3)         # → mmol/L
        push!(col_ids, c)
    end
end

N = length(all_cout)
println("Total data points: $N")

# ── Diagnostic quantities ─────────────────────────────────────────────────────
# ΔC/C_in: fractional NO₃⁻ removal

frac_removal = (all_cin .- all_cout) ./ all_cin

x_z = all_tau ./ all_cin               # zero order  [days L mmol⁻¹]
y_z = frac_removal

x_f = all_tau                          # first order [days]
y_f = -log.(1 .- frac_removal)

x_s = all_tau .* all_cin               # second order [days mmol L⁻¹]
y_s = frac_removal ./ (1 .- frac_removal)

# ── Linear regression through origin: y = a·x ────────────────────────────────

function fit_origin(x, y)
    a = sum(x .* y) / sum(x .^ 2)
    ŷ = a .* x
    SS_res = sum((y .- ŷ) .^ 2)
    SS_tot = sum((y .- mean(y)) .^ 2)
    R² = 1 - SS_res / SS_tot
    return a, R²
end

k0, R2_0 = fit_origin(x_z, y_z)
k1, R2_1 = fit_origin(x_f, y_f)
k2, R2_2 = fit_origin(x_s, y_s)

println("\n=== Diagnostic plot linear fits (through origin) ===")
println("Zero order:   k = $(k0) mmol L⁻¹ day⁻¹,  R² = $(round(R2_0, digits=4))")
println("First order:  k = $(k1) day⁻¹,            R² = $(round(R2_1, digits=4))")
println("Second order: k = $(k2) L mmol⁻¹ day⁻¹,   R² = $(round(R2_2, digits=4))")

# ── Diagnostic figure ─────────────────────────────────────────────────────────

fig = Figure(size = (500, 700))

ax_a = Axis(fig[1, 1],
    title = "a. Zero-order test",
    xlabel = "τ / Cᵢₙ [days L mmol⁻¹]",
    ylabel = "ΔC / Cᵢₙ [-]",
)

ax_b = Axis(fig[2, 1],
    title = "b. First-order test",
    xlabel = "τ [days]",
    ylabel = "-ln(1 - ΔC/Cᵢₙ) [-]",
)

ax_c = Axis(fig[3, 1],
    title = "c. Second-order test",
    xlabel = "τ ⋅ Cᵢₙ [days mmol L⁻¹]",
    ylabel = "(ΔC/Cᵢₙ) / (1 - ΔC/Cᵢₙ) [-]",
)

for c in 1:3
    mask = col_ids .== c
    scatter!(ax_a, x_z[mask], y_z[mask],
        color = colors[c], marker = markers[c], markersize = 8, label = "Column $c")
    scatter!(ax_b, x_f[mask], y_f[mask],
        color = colors[c], marker = markers[c], markersize = 8, label = "Column $c")
    scatter!(ax_c, x_s[mask], y_s[mask],
        color = colors[c], marker = markers[c], markersize = 8, label = "Column $c")
end

xr_z = collect(range(0, maximum(x_z) * 1.05, length = 100))
lines!(ax_a, xr_z, k0 .* xr_z,
    color = :black, linestyle = :dash, label = "Linear fit")
text!(ax_a, 0.05, 0.92, text = "R² = $(round(R2_0, digits=3))",
    space = :relative, fontsize = 11)

xr_f = collect(range(0, maximum(x_f) * 1.05, length = 100))
lines!(ax_b, xr_f, k1 .* xr_f,
    color = :black, linestyle = :dash, label = "Linear fit")
text!(ax_b, 0.05, 0.92, text = "R² = $(round(R2_1, digits=3))",
    space = :relative, fontsize = 11)

xr_s = collect(range(0, maximum(x_s) * 1.05, length = 100))
lines!(ax_c, xr_s, k2 .* xr_s,
    color = :black, linestyle = :dash, label = "Linear fit")
text!(ax_c, 0.05, 0.92, text = "R² = $(round(R2_2, digits=3))",
    space = :relative, fontsize = 11)

Legend(fig[4, 1], ax_a, merge = true, orientation = :horizontal, nbanks = 2)

save("figs/reaction_order_diagnostic.png", fig)
save("figs/reaction_order_diagnostic.pdf", fig)
println("\nDiagnostic figure saved.")
