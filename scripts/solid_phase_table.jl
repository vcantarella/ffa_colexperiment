"""
    solid_phase_table.jl

Aggregates solid-phase analysis results (Iron and TOC/TIC) for the sediment columns.
It reads raw Excel files, cleans and unifies the data, converts percentages to molar concentrations,
calculates mean values per column/depth, and exports the results to a CSV file.

# Outputs:
- `data/processed_results/soliphase_results.csv`: Table of solid phase concentrations (Fe, TOC, TIC).
"""

using XLSX
using DataFrames
using CairoMakie
using Statistics
include("plot_theme.jl")
set_theme!(figtheme)

# Helper to parse sample names
function parse_sample_name(name)
    # Check for Fe pattern: "C1 0-4"
    m_fe = match(r"^C(\d)\s+(\d-\d)", name)
    if !isnothing(m_fe)
        return (column="column " * m_fe.captures[1], depth=m_fe.captures[2])
    end
    
    # Check for TOC pattern: "C1-0-4"
    m_toc = match(r"^C(\d)-(\d-\d)", name)
    if !isnothing(m_toc)
        return (column="column " * m_toc.captures[1], depth=m_toc.captures[2])
    end

    # Check for the background sample (Br 17-18)
    m_br = match(r"^(Br)\s{1}(17-18)", name)
    if !isnothing(m_br)
        return (column=m_br.captures[1], depth=m_br.captures[2])
    end
    
    return nothing
end

println("--- Loading and Unifying Data ---")

# --- Load Fe Data ---
file_path_fe = "data/raw_lab_data/fe_analysis.xlsx"
sheet_name_fe = "fe_plate_1"
df_fe = DataFrame(XLSX.readtable(file_path_fe, sheet_name_fe))

# Select relevant columns and clean
df_fe_clean = select(df_fe, "name", "fe2+ [mol/kg]", "fe3+ [mol/kg]")
parsed_fe = parse_sample_name.(convert.(String, df_fe.name))
df_fe_clean.Column = [p !== nothing ? p.column : missing for p in parsed_fe]
df_fe_clean.Depth = [p !== nothing ? p.depth : missing for p in parsed_fe]
# Filter for valid C-columns
filter!(row -> !ismissing(row.Column), df_fe_clean)

# --- Load TOC Data ---
file_path_toc = "data/raw_lab_data/TOC_Analysis_Fuhrberg.xlsx"
sheet_name_toc = "Tabelle1"
df_toc = DataFrame(XLSX.readtable(file_path_toc, sheet_name_toc))

# Select relevant columns and clean
df_toc_clean = select(df_toc, "Name", "TOC%", "TIC%")
parsed_toc = parse_sample_name.(convert.(String, df_toc.Name))
df_toc_clean.Column = [p !== nothing ? p.column : missing for p in parsed_toc]
df_toc_clean.Depth = [p !== nothing ? p.depth : missing for p in parsed_toc]
# Filter for valid C-columns
filter!(row -> !ismissing(row.Column), df_toc_clean)

# Select data from CHNS
file_name = "data/raw_lab_data/Fuhrberg CHNS.xlsx"
sheet_name = "260212_Johann_KS"
df_s = DataFrame(XLSX.readtable(file_name, sheet_name))
df_s_clean = select(df_s, "Name", "S  [%]")
rename!(df_s_clean, "S  [%]"=> "S%")
parsed_s = parse_sample_name.(convert.(String, df_s.Name))
df_s_clean.Column = [p !== nothing ? p.column : missing for p in parsed_s]
df_s_clean.Depth = [p !== nothing ? p.depth : missing for p in parsed_s]
# Filter for valid C-columns
filter!(row -> !ismissing(row.Column), df_s_clean)


# --- Unify ---
# Join on Column and Depth
df_joined = outerjoin(df_fe_clean, df_toc_clean, on=[:Column, :Depth], makeunique=true)
df_joined = outerjoin(df_joined, df_s_clean, on=[:Column, :Depth], makeunique=true)
# --- Process Data ---
# 1. Convert TOC/TIC/Carbon % to mol/kg
# 1% = 10 g/kg
# Molar mass of C = 12.011 g/mol
const M_C = 12.011
const M_S = 32.065

function percent_to_mol_kg(percent, M)
    if ismissing(percent)
        return missing
    end
    # Ensure input is Float64
    p_float = try 
        Float64(percent) 
    catch 
        missing 
    end
    if ismissing(p_float) return missing end
    
    return (p_float * 10) / M
end

df_joined[!, "TOC [mol/kg]"] = percent_to_mol_kg.(df_joined[!, "TOC%"], M_C)
df_joined[!, "TIC [mol/kg]"] = percent_to_mol_kg.(df_joined[!, "TIC%"], M_C)
df_joined[!, "total S [mol/kg]"] = percent_to_mol_kg.(df_joined[!, "S%"], M_S)

# 2. Group by Column and Depth and calculate means
# Select columns to aggregate
cols_to_mean = ["fe2+ [mol/kg]", "fe3+ [mol/kg]", "TOC [mol/kg]", "TIC [mol/kg]", "total S [mol/kg]"]

# Group and combine
df_final = combine(groupby(df_joined, [:Column, :Depth]), 
    cols_to_mean .=> (x -> isempty(skipmissing(x)) ? missing : mean(skipmissing(x))) .=> cols_to_mean
)

# Sort
sort!(df_final, [:Column, :Depth])

println("Processed Unified Data (Means, mol/kg):")
println(df_final)

# rename background
df_final[df_final[!, :Column] .== "Br",:Depth] .= ""
df_final[df_final[!, :Column] .== "Br",:Column] .= "background"

df_final[!,:name] = df_final[!,:Column].*" - ".*df_final[!,:Depth]
df_final[!,:name] .= strip.(df_final[!,:name],' ')
df_final[!,:name] .= strip.(df_final[!,:name],'-')
df_final[!,:name] .= strip.(df_final[!,:name],' ')
df_final[!,:name] .= ifelse.(isnothing.(match.(r"(\d)$", df_final[!,:name])),df_final[!,:name], df_final[!,:name].*" cm")
names = df_final[!, :name]
using CSV
CSV.write("data/processed_results/soliphase_results.csv", df_final)

# calculate ratios
df_ratios = copy(df_final)
df_ratios[!, "C/S"] = df_ratios[!,"TOC [mol/kg]"] ./ df_ratios[!,"total S [mol/kg]"]
df_ratios[!, "Fe/S"] = df_ratios[!,"fe2+ [mol/kg]"] ./ df_ratios[!,"total S [mol/kg]"]

# mean values
meanC_S = mean(df_ratios[!, "C/S"])
meanFe_S = mean(df_ratios[!, "Fe/S"])

# Plot results in a Figure.
fig_height = 600
fig = Figure(size = (700, fig_height))
axtoc = Axis(fig[1, 1], title = "a. TOC",
    ylabel = "TOC [mol kg⁻¹]",
    )
axs = Axis(fig[1, 2], title = "b. total Sulphur",
    ylabel = "total S [mol kg⁻¹]",
    #yticks = 0:50:300
    )
axtic = Axis(fig[1, 3], title = "c. TIC",
    ylabel = "TIC [mol kg⁻¹]",
    #yticks = 0:50:300
    )
axfe2 = Axis(fig[2, 1], title = "d. Fe²⁺",
    ylabel = "Fe²⁺ [mol kg⁻¹]",
    #yticks = 0:2e-1:1.1
    )
axfe3 = Axis(fig[2, 2], title = "e. Fe³⁺",
    ylabel = "Fe³⁺ [mol kg⁻¹]",
    #yticks = 0:2e-1:1.1
    )

function get_color_ind(name)
    if contains(name, "column 1")
        return 1,"Column 1"
    elseif contains(name, "column 2")
        return 2,"Column 2"
    elseif contains(name, "column 3")
        return 3,"Column 3"
    elseif contains(name, "column 4")
        return 4,"Column 4"
    else
        return 5,"Background"
    end
end
for (i,name) in enumerate(names)
    color_ind, label_name = get_color_ind(name)
    barplot!(axtoc, i, df_final[df_final[!,:name].== name, "TOC [mol/kg]"][1], color = colors[color_ind], label = label_name)
    barplot!(axtic, i, df_final[df_final[!,:name].== name, "TIC [mol/kg]"][1], color = colors[color_ind])
    barplot!(axs, i, df_final[df_final[!,:name].== name, "total S [mol/kg]"][1], color = colors[color_ind])
    barplot!(axfe2, i, df_final[df_final[!,:name].== name, "fe2+ [mol/kg]"][1], color = colors[color_ind])
    barplot!(axfe3, i, df_final[df_final[!,:name].== name, "fe3+ [mol/kg]"][1], color = colors[color_ind])

end
# 3. Replace numbers with names on the X-axis
# We apply this to all axes in your figure
all_axes = [axtoc, axtic, axs, axfe2, axfe3]

for ax in [axfe2, axfe3]
    # Set ticks at positions 1 to N, and label them with sample_names
    ax.xticks = (1:length(names), names)
    
    # Rotate labels 45 degrees (π/4) so they don't overlap
    ax.xticklabelrotation = π/4
    ax.xticklabelalign = (:right, :center) # Aligns rotated text nicely
end

for ax in [axtoc, axtic, axs]
    # Set ticks at positions 1 to N, and label them with sample_names
    ax.xticks = (1:length(names), names)
    ax.xticklabelsvisible=false
end

Legend(fig[2,3], axtoc, "Sample", merge=true, tellwidth = false, halign = :left)

# 4. Final layout tweak
autolimits!.(all_axes) # Ensures the view matches the data
display(fig)

fig
save("figs/solid_phase_barplots.png", fig, px_per_unit = 2)
save("figs/solid_phase_barplots.pdf", fig, pt_per_unit = 1)