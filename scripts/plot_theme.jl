"""
    plot_theme.jl

Defines the Makie plot theme used across the project for consistent visualization.
Sets fonts, axis styles, and legend properties.
"""

using CairoMakie

# marker for columns 1:4
markers = [:circle, :rect, :diamond, :pentagon]

# color for columns 1:4
colors = [:crimson, :steelblue, :forestgreen, :darkorange]

figtheme = Theme(
    Figure = (
        fonts = :Helvetica,
    ),
    Axis = (
        titlealign = :left,
        titlesize = 13,
        xlabelsize = 11,
        ylabelsize = 11,
        xticklabelsize = 9,
        yticklabelsize = 9,
        xticksize = 9,
        yticksize = 9,
        topspinevisible = false,
        rightspinevisible = false,
        spinewidth = 1.7,
    ),
    Legend = (
        titlesize = 10,
        titlehalign = :left,
        labelsize = 9,
        framevisible = false,
        rowgap = 1,
        colgap = 9,
    )

)