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
        titlesize = 13,
        xlabelsize = 12,
        ylabelsize = 12,
        xticklabelsize = 11,
        yticklabelsize = 11,
        xticksize = 9,
        yticksize = 9,
        topspinevisible = false,
        rightspinevisible = false,
        spinewidth = 1.7,
    ),
    Legend = (
        titlesize = 12,
        titlehalign = :left,
        labelsize = 12,
        framevisible = false,
        rowgap = 1,
        colgap = 9,
    )

)