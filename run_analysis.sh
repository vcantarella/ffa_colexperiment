#!/bin/bash

# Exit on any error
set -e

echo "--- Starting Sediment Column Experiment Analysis ---"

# 1. Julia Environment Setup
echo "Instantiating Julia environment..."
julia --project=. -e 'using Pkg; Pkg.instantiate()'

# 2. Run Analysis Scripts in Order
echo "Running data preparation..."
julia --project=. scripts/preparation_tracer_modelling.jl
julia --project=. scripts/prepare_flow_velocity_cinflow.jl
julia --project=. scripts/prepare_outflow_data.jl
julia --project=. scripts/solid_phase_table.jl

echo "Running modeling and stoichiometry analysis..."
julia --project=. scripts/tracer_modelling.jl
julia --project=. scripts/lagrangian.jl
julia --project=. scripts/stoichiometry_reaction_analysis.jl

echo "Generating overview and comparison plots..."
julia --project=. scripts/outflow_data_overview.jl
julia --project=. scripts/velocity_data_plot.jl
julia --project=. scripts/comparison_full_profile.jl

# 3. Convert PDF figures to EPS
echo "Converting PDF figures to EPS using pdftops..."
if command -v pdftops >/dev/null 2>&1; then
    for pdf_file in figs/*.pdf; do
        # Check if file exists to handle empty glob
        [ -e "$pdf_file" ] || continue
        
        eps_file="${pdf_file%.pdf}.eps"
        echo "Converting $pdf_file -> $eps_file"
        pdftops -eps -level3 "$pdf_file" "$eps_file"
    done
    echo "Conversion complete."
else
    echo "Warning: pdftops (poppler) not found. Skipping EPS conversion."
fi

echo "--- Analysis pipeline completed successfully ---"
