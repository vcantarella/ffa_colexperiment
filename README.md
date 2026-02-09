# Sediment Column Experiment Analysis

This repository contains the Julia scripts and data for analyzing sediment column experiments investigating nitrate transport and reduction (denitrification) under time-variable flow conditions.

## Project Overview

The project involves four sediment columns (3 active, 1 control) subjected to varying flow rates and inflow concentrations. The analysis pipeline includes:
- Processing raw laboratory data (flow rates, concentrations of NO3, NO2, SO4, DIC, DOC, etc.).
- Estimating transport parameters (porosity, dispersivity) using bromide tracer tests.
- Modeling nitrate transport using a Lagrangian approach with zero-order reaction kinetics.
- Calculating denitrification rates and comparing them with solid-phase characteristics and literature data.
- Analyzing reaction stoichiometry.

## Repository Structure

```
.
├── data/
│   ├── external/              # Literature data (Eschenbach et al., Weymann et al.)
│   ├── raw_lab_data/          # Raw experimental data (Excel, CSV)
│   └── processed_results/     # Generated intermediate data files (.jld2, .csv)
├── figs/                      # Generated figures (PNG, PDF)
├── scripts/                   # Julia analysis scripts
├── Project.toml               # Julia project dependencies
└── README.md                  # This file
```

## Requirements

The project uses **Julia**. The dependencies are listed in `Project.toml`. 

**Key Packages:**
- `CairoMakie`: For plotting.
- `DataFrames`, `CSV`, `XLSX`: For data manipulation.
- `JLD2`: For saving/loading Julia data structures.
- `DataInterpolations` / `QuadGK`: For modeling math.
- `nonlinearlstr`: For parameter estimation: (https://github.com/vcantarella/nonlinearlstr.jl)

## Usage

To run the full analysis pipeline and generate all figures (including PDF and EPS formats), you can use the provided bash script:

```bash
./run_analysis.sh
```

Alternatively, to run manually, ensure you have instantiated the Julia environment:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

### Recommended Execution Order

The scripts should generally be run in the following order to ensure necessary data dependencies are met:

1.  **Data Preparation:**
    - `scripts/preparation_tracer_modelling.jl`: Processes bromide tracer data.
    - `scripts/prepare_flow_velocity_cinflow.jl`: Processes flow velocity and inflow schedules.
    - `scripts/prepare_outflow_data.jl`: Processes outflow concentrations.
    - `scripts/solid_phase_table.jl`: Aggregates solid-phase analysis (Fe, TOC).

2.  **Modeling & Analysis:**
    - `scripts/tracer_modelling.jl`: Fits the Ogata-Banks model to tracer data to find transport parameters.
    - `scripts/lagrangian.jl`: Runs the Lagrangian nitrate transport model and estimates reaction rates.
    - `scripts/stoichiometry_reaction_analysis.jl`: Analyzes metabolic stoichiometry (DIC/NO3, DOC/NO3).

3.  **Visualization & Comparison:**
    - `scripts/outflow_data_overview.jl`: Plots overview of all outflow concentrations.
    - `scripts/velocity_data_plot.jl`: Plots flow velocity and discharge profiles.
    - `scripts/comparison_full_profile.jl`: Compares calculated rates with literature profiles.

## Script Descriptions

| Script | Description | Outputs (Data/Figs) |
| :--- | :--- | :--- |
| `model_data_structures.jl` | Defines common data structures (`QData`, `ds`, etc.). | N/A |
| `plot_theme.jl` | Sets the visual style for Makie plots. | N/A |
| `preparation_tracer_modelling.jl` | Prepares bromide tracer data (cleaning, dead volume correction). | `bromide_breakthrough_data.jld2`, `figs/bromide_breakthrough_prep.*` |
| `prepare_flow_velocity_cinflow.jl` | Constructs velocity profiles and inflow concentration schedules. | `flow_velocity_data.jld2`, `inflow_data.jld2` |
| `prepare_outflow_data.jl` | Cleans and aligns outflow concentration data (NO3, DIC, DOC, etc.). | `outflow_data.jld2` |
| `solid_phase_table.jl` | Aggregates solid-phase Fe and TOC/TIC data. | `soliphase_results.csv` |
| `tracer_modelling.jl` | Fits transport model to estimate porosity and dispersivity. | `tracer_params.jld2`, `figs/breakthrough_fit.*` |
| `lagrangian.jl` | Main simulation: calculates residence times and reaction rates. | `no3_rates.csv`, `figs/nitrate_and_rates.*` |
| `stoichiometry_reaction_analysis.jl` | Calculates stoichiometric ratios for reaction characterization. | `stoichiometry_dic_doc.csv` |
| `outflow_data_overview.jl` | Plots comprehensive time-series of outflow data. | `figs/outflow_plot.*` |
| `velocity_data_plot.jl` | Plots flow velocity and discharge over time. | `figs/flow_velocity.*` |
| `comparison_full_profile.jl` | Compares project rates with Eschenbach & Weymann datasets. | `no3_rates_final.csv`, `figs/comparison_full_profile.*` |

## Data Sources

- **Lab Data:** Located in `data/raw_lab_data/`, includes Excel files for tracer curves, general samples, and solid phase analysis.
- **External Data:** `data/external/` contains CSVs from *Eschenbach et al. (2015)* and *Weymann et al. (2010)* used for comparison.

