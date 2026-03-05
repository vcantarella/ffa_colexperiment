using JLD2
using DataFrames
using XLSX
using Dates

# Include the data structures so JLD2 can reconstruct the objects
include("model_data_structures.jl")

println("Starting data export to Excel...")

output_file = "data/processed_results/experimental_processed_dataset.xlsx"

# 1. Load Outflow Data (Main Chemistry)
load_outflow = load("data/processed_results/outflow_data.jld2")
all_ds = load_outflow["all_ds"]

df_outflow = DataFrame()
for col_id in sort(collect(keys(all_ds)))
    d = all_ds[col_id]
    
    # Helper to extract conc_ds into a temp dataframe
    function extract_conc(cds, name, cid)
        return DataFrame(
            Column = cid,
            Time_s = cds.t,
            Variable = name,
            Value_mM = cds.conc
        )
    end
    
    append!(df_outflow, extract_conc(d.no3, "NO3", col_id))
    append!(df_outflow, extract_conc(d.doc, "DOC", col_id))
    append!(df_outflow, extract_conc(d.dic, "DIC", col_id))
    append!(df_outflow, extract_conc(d.no2, "NO2", col_id))
    append!(df_outflow, extract_conc(d.so4, "SO4", col_id))
    append!(df_outflow, extract_conc(d.pH, "pH", col_id))
    append!(df_outflow, extract_conc(d.ec, "EC", col_id))
end

# 2. Load Bromide Tracer Data
load_br = load("data/processed_results/bromide_breakthrough_data.jld2")
br_conc = load_br["Br"]
br_time = load_br["avg_time"]

df_bromide = DataFrame()
for col_id in sort(collect(keys(br_conc)))
    append!(df_bromide, DataFrame(
        Column = col_id,
        Time_s = br_time[col_id],
        Br_mM = br_conc[col_id]
    ))
end

# 3. Load Flow Rate Data
load_flow = load("data/processed_results/bromide_flow_rate_data.jld2")
qs = load_flow["Qs"]
st = load_flow["start_times_dict"]
et = load_flow["end_times_dict"]

df_flow = DataFrame()
for col_id in sort(collect(keys(qs)))
    append!(df_flow, DataFrame(
        Column = col_id,
        Start_Time_s = st[col_id],
        End_Time_s = et[col_id],
        Flow_Rate_cm3_s = qs[col_id]
    ))
end

# 4. Create Metadata/Comments Sheet
metadata = [
    "Sheet" "Description" "Units / Comments";
    "Outflow_Concentrations" "Processed concentrations of NO3, NO2, SO4, DOC, DIC, pH, EC" "Time in seconds since t0; Concentrations in mM (corrected for dead volume); pH (unitless); EC (μS/cm)";
    "Bromide_Tracer" "Breakthrough curve data for Bromide" "Time in seconds since t0; Concentration in mM";
    "Flow_Rates" "Measured flow rates during the experiment" "Flow rate in cm³/s; Times in seconds since t0";
    "General_Note" "This dataset contains the processed and corrected results from the column experiments." "t0 is the start of the experiment for each column (corrected for inflow dead volume)."
]

# Write to Excel
XLSX.openxlsx(output_file, mode="w") do xf
    # Outflow Concentrations
    sheet1 = xf[1]
    XLSX.rename!(sheet1, "Outflow_Concentrations")
    XLSX.writetable!(sheet1, df_outflow)
    
    # Bromide
    XLSX.addsheet!(xf, "Bromide_Tracer")
    XLSX.writetable!(xf["Bromide_Tracer"], df_bromide)
    
    # Flow
    XLSX.addsheet!(xf, "Flow_Rates")
    XLSX.writetable!(xf["Flow_Rates"], df_flow)
    
    # Metadata
    XLSX.addsheet!(xf, "Metadata")
    sh_meta = xf["Metadata"]
    for i in 1:size(metadata, 1)
        for j in 1:size(metadata, 2)
            sh_meta[XLSX.CellRef(i, j)] = metadata[i, j]
        end
    end
end

println("Export complete: $output_file")
