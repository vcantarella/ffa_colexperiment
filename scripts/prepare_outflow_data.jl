using DataFrames
using CSV
using XLSX
using Dates

include("prepare_flow_velocity_cinflow.jl")
# Dead volumes per column. This refer only to the dead volume of the outflow tube
tube_diam = 1.52e-3 # m
# load processed DIC data
df_cp = CSV.read("data/raw_lab_data/DIC_DOC_raw_data.csv", DataFrame)
df_cp[!, "t"] .= NaN
# process time
for i in 1:4
    index = df_cp[!, "column"] .== i
    end_times = df_cp[index, "end_time"]
    t_uncorrected = Dates.value.(Dates.Second.(df_cp[index, "avg_time"] - t0s[i]))
    df_cp[index, "t"] = t_uncorrected .-
        dvs[i]/100*π*tube_diam^2/4 ./ disch_function.(t_uncorrected, Ref(disch_ds[i].Q), Ref(disch_ds[i].end_times))
end



start_times_dict = Dict(1 => [Dates.value(Dates.Second(t - t0s[1])) for t in start_time[column .== 1]],
                        2 => [Dates.value(Dates.Second(t - t0s[2])) for t in start_time[column .== 2]],
                        3 => [Dates.value(Dates.Second(t - t0s[3])) for t in start_time[column .== 3]],
                        4 => [Dates.value(Dates.Second(t - t0s[4])) for t in start_time[column .== 4]])
end_times_dict = Dict(1 => [Dates.value(Dates.Second(t - t0s[1])) for t in end_time[column .== 1]],
                      2 => [Dates.value(Dates.Second(t - t0s[2])) for t in end_time[column .== 2]],
                      3 => [Dates.value(Dates.Second(t - t0s[3])) for t in end_time[column .== 3]],
                      4 => [Dates.value(Dates.Second(t - t0s[4])) for t in end_time[column .== 4]])

avg_times_dict = Dict(1 => (start_times_dict[1] .+ end_times_dict[1]) ./ 2 .-
    dvs[1]/100*π*tube_diam^2/4 ./ disch_function.(end_times_dict[1], Ref(disch_ds[1].Q), Ref(disch_ds[1].end_times)),
    2 => (start_times_dict[2] .+ end_times_dict[2]) ./ 2 .-
    dvs[2]/100*π*tube_diam^2/4 ./ disch_function.(end_times_dict[2], Ref(disch_ds[2].Q), Ref(disch_ds[2].end_times)),
    3 => (start_times_dict[3] .+ end_times_dict[3]) ./ 2 .- 
    dvs[3]/100*π*tube_diam^2/4 ./ disch_function.(end_times_dict[3], Ref(disch_ds[3].Q), Ref(disch_ds[3].end_times)),
    4 => (start_times_dict[4] .+ end_times_dict[4]) ./ 2 .-
    dvs[4]/100*π*tube_diam^2/4 ./ disch_function.(end_times_dict[4], Ref(disch_ds[4].Q), Ref(disch_ds[4].end_times)),
)

# Now we have the average times for each column, we can create the datasets per column
no3 = vcat(df[!,"no3- [mgN/L]"],df_tr[!,"no3- [mgN/L]"]) # concentration in mM
no3_ic = vcat(df[!,"no3_mg_L"], repeat([missing], size(df_tr[!,"no3- [mgN/L]"], 1))) # concentration in mM
std_no3 = vcat(df[!,"CI_no3"], df_tr[!,"CI_no3"]) # standard deviation in mM
no2 = vcat(df[!,"no2_mmol_L"], repeat([missing], size(df_tr[!,"no3- [mgN/L]"], 1))) # concentration in mM
so4 = vcat(df[!,"so4_mg_L"], repeat([missing], size(df_tr[!,"no3- [mgN/L]"], 1))) # concentration in mM

#fe = df[!,"Fe2+"] # concentration in mM
pH = vcat(df[!,"pH"], repeat([missing], size(df_tr[!,"no3- [mgN/L]"], 1))) # pH values
ec = vcat(df[!,"EC"], repeat([missing], size(df_tr[!,"no3- [mgN/L]"], 1))) # electrical conductivity in μS/cm


# For each column, we create a dataset with the concentrations and the average times
all_ds = Dict{Int, ds_m2}()
for i in 1:4
    # only the values where there are no missing values in the concentration
    id1 = column .== i
    # no3 dataset is a mix of cuvette and IC measurements. We prioritize IC measurements when available.
    t = avg_times_dict[i][.!ismissing.(no3_ic[id1])]
    t_cuv = avg_times_dict[i][.!ismissing.(no3[id1]) .& ismissing.(no3_ic[id1])]
    no3_v = vcat(convert.(Float64, no3_ic[id1][.!ismissing.(no3_ic[id1])])./62,
                convert.(Float64, no3[id1][.!ismissing.(no3[id1]) .& ismissing.(no3_ic[id1])])./14)
    no3_v = convert.(Float64, no3_v)
    avg_times = vcat(t, t_cuv)
    argsort_idx = sortperm(avg_times)
    avg_times_sorted = avg_times[argsort_idx]
    no3_sorted = no3_v[argsort_idx]
    no3_ds = conc_ds(no3_sorted, avg_times_sorted)
    cp_index = df_cp[!,:column] .== i
    doc = df_cp[cp_index, "DOC"]./12
    dic = df_cp[cp_index, "DIC"]./12
    doc_ds = conc_ds(convert.(Float64, doc[.!ismissing.(doc)]), df_cp[cp_index,"t"][.!ismissing.(doc)])
    dic_ds = conc_ds(convert.(Float64, dic[.!ismissing.(dic)]), df_cp[cp_index,"t"][.!ismissing.(dic)])
    no2_col = conc_ds(convert.(Float64, no2[id1][.!ismissing.(no2[id1])]), avg_times_dict[i][.!ismissing.(no2[id1])])
    so4_col = conc_ds(convert.(Float64, so4[id1][.!ismissing.(so4[id1])]./96), avg_times_dict[i][.!ismissing.(so4[id1])])
    # create a dataset for the column
    col_ds = ds_m2(i, no3_ds, doc_ds, dic_ds, no2_col, so4_col)
    all_ds[i] = col_ds
    # save the dataset to a file
    
end

@save "data/processed_results/outflow_data.jld2" all_ds

println("Data preparation complete. Datasets for each column are ready.")