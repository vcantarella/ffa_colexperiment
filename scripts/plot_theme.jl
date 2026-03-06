"""
    plot_theme.jl

Defines the Makie plot theme used across the project for consistent visualization.
Sets fonts, axis styles, and legend properties.
"""

using CairoMakie

# marker for columns 1:4
markers = [:circle, :rect, :diamond, :pentagon]

# color for columns 1:4(+background)
colors = [:crimson, :steelblue, :forestgreen, :darkorange, :darkgrey]

figtheme = Theme(
    Figure = (
        fonts = :Helvetica,
    ),
    Axis = (
        titlealign = :left,
        titlesize = 10,
        xlabelsize = 9,
        ylabelsize = 9,
        xticklabelsize = 8,
        yticklabelsize = 8,
        xticksize = 7,
        yticksize = 7,
        topspinevisible = false,
        rightspinevisible = false,
        spinewidth = 1.2,
    ),
    Legend = (
        titlesize = 8,
        titlehalign = :left,
        labelsize = 9,
        framevisible = false,
        rowgap = 1,
        colgap = 8,             # Space between items in horizontal layout
        patchlabelgap = 2,      # Space between marker and text
        patchsize = (10, 10),   # Smaller markers/lines in legend
        padding = (2, 2, 2, 2), # Minimal internal padding
    )

)