"""
    comparison_full_profile.jl

This script generates a full profile comparison plot of denitrification rates, 
organic carbon (Corg), and total sulphur (Stotal) versus depth. 
It compares the current project data with literature values from Eschenbach et al. (2013, 2015) 
and Weymann et al. (2010).

# Outputs:
- `figs/comparison_full_profile.png`: The generated plot.
- `figs/comparison_full_profile.pdf`: The generated plot in PDF format.
- `data/processed_results/no3_rates_final.csv`: Updated project data with 10°C rate estimates.
"""

using CSV
using DataFrames
using Statistics
using CairoMakie
include("plot_theme.jl")
set_theme!(figtheme)
# --- Constants ---
const m_N = 14.007  # g/mol
const m_S = 32.065  # g/mol
const m_C = 12.011  # g/mol

# --- Helpers ---
function parse_depth(d_str)
    ismissing(d_str) && return missing
    try
        parts = split(string(d_str), "-")
        return length(parts) == 2 ? (parse(Float64, parts[1]) + parse(Float64, parts[2])) / 2 : parse(Float64, d_str)
    catch
        return missing
    end
end

function normalize_zone(z_raw, ref_dataset)
    z = lowercase(string(z_raw))
    if occursin("weymann", lowercase(ref_dataset))
        return occursin("autotrophic", z) ? "Sulphidic" : occursin("heterotrophic", z) ? "Non-Sulphidic" : titlecase(z)
    end
    return occursin("non-sulphidic", z) || occursin("non sulphidic", z) ? "Non-Sulphidic" : 
           occursin("transition", z) ? "Transition Zone" : 
           occursin("sulphidic", z) ? "Sulphidic" : titlecase(z)
end

function get_col_robust(df, prefix)
    for n in names(df)
        if occursin(lowercase(prefix), lowercase(n))
            return n
        end
    end
    return nothing
end

# --- Data Containers ---
rates_rows = []
conc_rows = []

# --- 1. Load Project Rates ---
df_project = CSV.read("data/processed_results/no3_rates.csv", DataFrame)
for row in eachrow(df_project)
    push!(rates_rows, (
        Reference = "Current Project",
        Zone = "Sulphidic",
        Depth_Mid = 17.5,
        Rate_Min = row.rate_mol_kg_d / 3.8, # Conservative Q10
        Rate_Max = row.rate_mol_kg_d / 1.4  # Optimistic Q10
    ))
end

# --- 1b. Load Project Concentrations ---
df_soli = CSV.read("data/processed_results/soliphase_results.csv", DataFrame)
for row in eachrow(df_soli)
    push!(conc_rows, (
        Reference = "Current Project",
        Zone = "Sulphidic",
        Depth_Mid = 17.5,
        C_org = row["TOC [mol/kg]"],
        S_total = row["total S [mol/kg]"]
    ))
end

# --- 2. Load Eschenbach et al. (2015) ---
df_esch = CSV.read("data/external/eschenbach_2013_2015.csv", DataFrame)
col_dcum = get_col_robust(df_esch, "Dcum")
col_insitu = get_col_robust(df_esch, "Dr in situ")
col_corg_esch = get_col_robust(df_esch, "Corg")
col_stot_esch = get_col_robust(df_esch, "Total-S")

for row in eachrow(df_esch)
    !occursin("FFA", string(row["Sample Location"])) && continue
    d_mid = parse_depth(row["Depth interval (m)"])
    ismissing(d_mid) && continue
    
    z = normalize_zone(row["Sediment Group"], "Eschenbach")
    ref = "Eschenbach et al.\n (2013,2015)"
    
    # Rates: Collect ANY available rate
    available_rates = Float64[]
    
    if !isnothing(col_insitu) && !ismissing(row[col_insitu])
        push!(available_rates, (row[col_insitu] * 1e-6) / m_N)
    end
    if !isnothing(col_dcum) && !ismissing(row[col_dcum])
        push!(available_rates, (row[col_dcum] * 1e-3) / m_N / 365.0)
    end
    
    if !isempty(available_rates)
        push!(rates_rows, (
            Reference=ref, 
            Zone=z, 
            Depth_Mid=d_mid, 
            Rate_Min=minimum(available_rates), 
            Rate_Max=maximum(available_rates)
        ))
    end
    
    # Concentrations
    c_val = (!isnothing(col_corg_esch) && !ismissing(row[col_corg_esch])) ? (row[col_corg_esch] * 1e-3) / m_C : missing
    s_val = (!isnothing(col_stot_esch) && !ismissing(row[col_stot_esch])) ? (row[col_stot_esch] * 1e-3) / m_S : missing
    
    if !ismissing(c_val) || !ismissing(s_val)
        push!(conc_rows, (Reference=ref, Zone=z, Depth_Mid=d_mid, C_org=c_val, S_total=s_val))
    end
end

# --- 3. Load Weymann et al. (2010) ---
df_wey = CSV.read("data/external/weymann_et_al_2010.csv", DataFrame)
# Strict column finding for Weymann to avoid partial matches
cols_di = [n for n in names(df_wey) if startswith(n, "Di") || startswith(n, "Di\u200b")]
col_di = isempty(cols_di) ? nothing : first(cols_di)

cols_dmax = [n for n in names(df_wey) if startswith(n, "Dmax") || startswith(n, "Dmax\u200b")]
col_dmax = isempty(cols_dmax) ? nothing : first(cols_dmax)

col_corg_wey = get_col_robust(df_wey, "Org C")
col_stot_wey = get_col_robust(df_wey, "Total S")

for row in eachrow(df_wey)
    d_mid = parse_depth(row["Depth (m)"])
    ismissing(d_mid) && continue
    
    z = normalize_zone(row["Zone"], "Weymann")
    ref = "Weymann et al. (2010)"
    
    # Rates
    available_rates = Float64[]
    
    if !isnothing(col_di) && !ismissing(row[col_di])
        push!(available_rates, (row[col_di] * 1e-6) / m_N)
    end
    if !isnothing(col_dmax) && !ismissing(row[col_dmax])
        push!(available_rates, (row[col_dmax] * 1e-6) / m_N)
    end
    
    if !isempty(available_rates)
        push!(rates_rows, (
            Reference=ref, 
            Zone=z, 
            Depth_Mid=d_mid, 
            Rate_Min=minimum(available_rates), 
            Rate_Max=maximum(available_rates)
        ))
    end
    
    # Concentrations
    c_val = (!isnothing(col_corg_wey) && !ismissing(row[col_corg_wey])) ? (row[col_corg_wey] * 1e-3) / m_C : missing
    s_val = (!isnothing(col_stot_wey) && !ismissing(row[col_stot_wey])) ? (row[col_stot_wey] * 1e-3) / m_S : missing
    
    if !ismissing(c_val) || !ismissing(s_val)
        push!(conc_rows, (Reference=ref, Zone=z, Depth_Mid=d_mid, C_org=c_val, S_total=s_val))
    end
end

df_rates = DataFrame(rates_rows)
df_conc = DataFrame(conc_rows)

# --- Plotting ---
f = Figure(size = (1200, 700))

# Layout: 3 Columns [Rates | Corg | S_total] + Legend
ax_rates = Axis(f[1, 1], xlabel = "Rate (mol N kg⁻¹ d⁻¹)",
 ylabel = "Depth (m)", title = "a. Denitrification Rates",
 yreversed = true, xscale = log10, titlealign= :left)
ax_c = Axis(f[1, 2], xlabel = "Corg (mol C kg⁻¹)",
 title = "b. Organic Carbon",
 yreversed = true, xscale = log10,
 titlealign = :left)
ax_s = Axis(f[1, 3], xlabel = "Stotal (mol S kg⁻¹)",
 title = "c. Total Sulphur",
  yreversed = true, xscale = log10,
  titlealign = :left)

linkyaxes!(ax_rates, ax_c, ax_s)
hideydecorations!(ax_c, grid = false)
hideydecorations!(ax_s, grid = false)

# Styling
zone_colors = Dict("Non-Sulphidic" => :darkgoldenrod4,
                                    "Sulphidic" => :gold1,
                                    "Transition Zone" => :darkkhaki)
ref_markers = Dict("Current Project" => :circle,
                                    "Eschenbach et al.\n (2013,2015)" => :rect,
                                    "Weymann et al. (2010)" => :diamond)
ref_markersisze = Dict("Current Project" => 19,
                                    "Eschenbach et al.\n (2013,2015)" => 16,
                                    "Weymann et al. (2010)" => 16)

linewidth = 5
strokewidth = 1.0
strokecolor = :grey21 
# 1. Rates Plot
for row in eachrow(df_rates)
    c = get(zone_colors, row.Zone, :gray)
    m = get(ref_markers, row.Reference, :circle)
    msize = get(ref_markersisze, row.Reference, 16)
    # Range Bar (only if Min != Max)
    if row.Rate_Min != row.Rate_Max
        rangebars!(ax_rates, [row.Depth_Mid], [row.Rate_Min], [row.Rate_Max], direction = :x, color = (c, 0.8), linewidth = linewidth)
    end
    # Markers
    scatter!(ax_rates, [row.Rate_Min, row.Rate_Max], [row.Depth_Mid, row.Depth_Mid], color = c, marker = m, markersize = msize,
     strokewidth = strokewidth, strokecolor = strokecolor)
end

# 2. Corg Plot
for row in eachrow(filter(r -> !ismissing(r.C_org), df_conc))
    c = get(zone_colors, row.Zone, :gray)
    m = get(ref_markers, row.Reference, :circle)
    msize = get(ref_markersisze, row.Reference, 16)
    scatter!(ax_c, row.C_org, row.Depth_Mid, color = c, marker = m, markersize = msize,
     strokewidth = strokewidth, strokecolor = strokecolor)
end

# 3. Stotal Plot
for row in eachrow(filter(r -> !ismissing(r.S_total), df_conc))
    c = get(zone_colors, row.Zone, :gray)
    m = get(ref_markers, row.Reference, :circle)
    msize = get(ref_markersisze, row.Reference, 16)
    scatter!(ax_s, row.S_total, row.Depth_Mid, color = c, marker = m, markersize = msize,
     strokewidth = strokewidth, strokecolor = strokecolor)
end

# --- Highlight Current Project Data ---
# We add a red dashed box to highlight the location of the project data (depth ~17.5m)
p_depth_min, p_depth_max = 17.1, 17.9

# Highlight in Rates plot
p_rates = filter(r -> r.Reference == "Current Project", df_rates)
if !isempty(p_rates)
    r_min = minimum([p_rates.Rate_Min; p_rates.Rate_Max])
    r_max = maximum([p_rates.Rate_Min; p_rates.Rate_Max])
    # Increased padding for visibility on log scale (0.5x to 2.0x)
    poly!(ax_rates, Rect2f(r_min * 0.5, p_depth_min, r_max * 2.0 - r_min * 0.5, p_depth_max - p_depth_min),
          color = :transparent, strokecolor = :red, strokewidth = 2.5, linestyle = :dash)
end

# Highlight in Corg plot
p_c = filter(r -> r.Reference == "Current Project" && !ismissing(r.C_org), df_conc)
if !isempty(p_c)
    c_min, c_max = minimum(p_c.C_org), maximum(p_c.C_org)
    poly!(ax_c, Rect2f(c_min * 0.8, p_depth_min, c_max * 1.2 - c_min * 0.8, p_depth_max - p_depth_min),
          color = :transparent, strokecolor = :red, strokewidth = 2.5, linestyle = :dash)
end

# Highlight in Stotal plot
p_s = filter(r -> r.Reference == "Current Project" && !ismissing(r.S_total), df_conc)
if !isempty(p_s)
    s_min, s_max = minimum(p_s.S_total), maximum(p_s.S_total)
    # Adjusted padding for better fit on log scale (0.7x to 1.4x)
    poly!(ax_s, Rect2f(s_min * 0.7, p_depth_min, s_max * 1.4 - s_min * 0.7, p_depth_max - p_depth_min),
          color = :transparent, strokecolor = :red, strokewidth = 2.5, linestyle = :dash)
end

# Legends
dataset_elements = Any[]
dataset_labels = String[]
# Define the order explicitly for the legend
for name in ["Current Project", "Eschenbach et al.\n (2013,2015)", "Weymann et al. (2010)"]
    m = MarkerElement(marker = ref_markers[name], color = :black, markersize = 12)
    if name == "Current Project"
        # Combine marker and highlight box for the Current Project
        push!(dataset_elements, [m, PolyElement(color = :transparent, strokecolor = :red, strokewidth = 2, linestyle = :dash)])
    else
        push!(dataset_elements, m)
    end
    push!(dataset_labels, name)
end

zone_leg = [MarkerElement(marker = :circle, color = zone_colors[z], markersize = 12) => z for z in sort(collect(keys(zone_colors)))]

Legend(f[1, 4], 
    [dataset_elements, first.(zone_leg)], 
    [dataset_labels, last.(zone_leg)], 
    ["Dataset", "Zone"])

save("figs/comparison_full_profile.png", f)
save("figs/comparison_full_profile.pdf", f)
println("Full profile plot updated: figs/comparison_full_profile.png and pdf")

# --- Export Project Data with 10C Rates ---
println("\n--- Exporting Updated Project Data ---")
# Calculate 10C estimates
df_project.rate_min_10C = df_project.rate_mol_kg_d ./ 3.8
df_project.rate_max_10C = df_project.rate_mol_kg_d ./ 1.4

# Write back to CSV
CSV.write("data/processed_results/no3_rates_final.csv", df_project)
println("Updated 'data/no3_rates.csv' with 10°C rate estimates.")

# --- Linear Regression: Total S vs TOC (All Data) ---
println("\n--- Performing S-C Stoichiometry Analysis (All Data) ---")
df_reg_all = filter(r -> !ismissing(r.C_org) && !ismissing(r.S_total), df_conc)

# Extract values
x_all = df_reg_all.C_org
y_all = df_reg_all.S_total

if length(x_all) > 1
    m_all = cov(x_all, y_all) / var(x_all)
    b_all = mean(y_all) - m_all * mean(x_all)

    println("Equation (Combined Data): Total S = $(round(m_all, digits=4)) * TOC + $(round(b_all, digits=4))")

    # --- Plotting ---
    fig_reg = Figure(size = (600, 500))
    ax_reg = Axis(fig_reg[1, 1], 
        xlabel = "TOC (mol C kg⁻¹)", 
        ylabel = "Total S (mol S kg⁻¹)",
        title = "S-C Stoichiometry (Current Project + Literature)")

    # Data points colored by Reference
    for ref in unique(df_reg_all.Reference)
        df_sub = filter(r -> r.Reference == ref, df_reg_all)
        scatter!(ax_reg, df_sub.C_org, df_sub.S_total, 
            marker = ref_markers[ref], 
            markersize = 12, 
            strokewidth = 1, 
            strokecolor = :black,
            label = ref)
    end

    # Regression line
    x_range = range(minimum(x_all)*0.8, maximum(x_all)*1.2, length=100)
    y_line = m_all .* x_range .+ b_all
    lines!(ax_reg, x_range, y_line, color = :red, linewidth = 2, label = "Regression (All Data)")

    axislegend(ax_reg, position = :rt, labelsize = 10)

    save("figs/s_toc_regression.png", fig_reg)
    println("Regression plot saved to figs/s_toc_regression.png")
else
    println("Not enough data points for combined regression.")
end