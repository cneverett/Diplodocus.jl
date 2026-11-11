const inch = 96
const pt = 4/3
const cm = inch / 2.54 

#= 
Latex double column width is 236pt which is 3.29 inches or 8.36cm

Deault figure size will be designed for single column with a width of 3.25 inches or 8.255cm and a aspect ratio of 4:3, therefore height of 2.4375 inches or 6.19125cm

=#

function DiplodocusDark()
    # inspired from Makie's theme_black but with latex fonts and some modifications
    return Theme(
        #figure_padding = 10.0,
        backgroundcolor=(:black,0.0),
        textcolor=:white,
        linecolor=:white,
        size=(3.25inch,2.4375inch),
        fontsize = 9pt,
        palette = (color = Makie.ColorScheme(get(Makie.ColorSchemes.batlow10,range(0.9,0.1,length=7))),),
        colormap = Makie.ColorScheme(get(Makie.ColorSchemes.batlow,range(0.1,0.9,length=512))), # use a subset to avoid very dark or very light colors
        colormap_var = Makie.ColorScheme(get(Makie.ColorSchemes.devon,range(0.0,1.0,length=512))), # use a subset to avoid very dark or very light colors
        fonts = Attributes(
            :bold => Makie.texfont(:bold),
            :bolditalic => Makie.texfont(:bolditalic),
            :italic => Makie.texfont(:italic),
            :regular => Makie.texfont(:regular)
        ),
        Figure = (
            backgroundcolor = (:black,0.0),#:transparent,
            fontsize = 9pt
        ),
        Axis = (
            fontsize = 9pt,
            backgroundcolor=(:black,0.0),
            bottomspinecolor = :white,
            topspinecolor = :white,
            leftspinecolor = :white,
            rightspinecolor = :white,
            xgridcolor = :grey45,
            ygridcolor = :grey45,
            xtickcolor = :white,
            ytickcolor = :white,
            xtickwidth = 2.0,
            ytickwidth = 2.0,
            # modification from theme_black
            titlesize = 11pt,
            #xlabelsize = 10pt, # defults if :fontsize
            #ylabelsize = 10pt, # defults if :fontsize
            titlecolor= :white,
            ylabelcolor= :white,
            xlabelcolor= :white,
            spinewidth = 2.0,
            xgridwidth = 1.0,
            ygridwidth = 1.0,
            xminorgridcolor=:grey45,
            yminorgridcolor=:grey45,
            xminortickcolor=:white,
            yminortickcolor=:white,
            xminorgridwidth = 0.5,
            yminorgridwidth = 0.5,
            xticklabelcolor=:white,
            yticklabelcolor=:white,
            xticklabelsize = 8pt,
            yticklabelsize = 8pt,
        ),
        Legend = (
            fontsize = 9pt,
            labelsize = 8pt,
            framecolor = :white,
            backgroundcolor = (:black),
            framewidth = 2.0,
            rowgap = -5,
            padding = (4f0,4f0,2f0,2f0),
        ),
        Colorbar = (
            fontsize = 9pt,
            labelsize = 9pt,
            ticklabelsize = 8pt,
            tickcolor = :white,
            spinecolor = :white,
            topspinecolor = :white,
            bottomspinecolor = :white,
            leftspinecolor = :white,
            rightspinecolor = :white,
            # modifications from theme_black
            spinewidth= 2.0,
        ),
        ScatterLines = (
            linewidth = 2.0,
            joinstyle = :round,
            linecap = :round
        ),
        Lines = (
            linewidth = 2.0,
            joinstyle = :round,
            linecap = :round
        ),
        PolarAxis = (
            fontsize = 9pt,
            labelsize = 9pt,
            rticklabelsize = 8pt,
            thetaticklabelsize = 8pt,
            backgroundcolor = (:black,0.0),
            rticksvisible = true,
            rticksmirrored = true,
            thetaticksvisible = true,
            rtickcolor = :white,
            rtickwidth = 2.0,
            thetatickcolor = :white,
            thetatickwidth = 2.0,
            gridcolor = :grey45,
            gridwidth = 1.0,
            tickcolor = :white,
            tickwidth = 2.0,
            labelcolor = :white,
            spinewidth = 2.0,
            spinecolor = :white,
            clip = false,
        )
    )
end

function DiplodocusLight()
    # slight modifications from the defa1ult theme but with Latex fonts and some modifications
    return Theme(
        #figure_padding = 10.0,
        backgroundcolor=(:white,0.0),
        textcolor=:black,
        linecolor=:black,
        size=(3.25inch,2.4375inch),
        fontsize = 9pt,
        palette = (color = Makie.ColorScheme(get(Makie.ColorSchemes.batlow10,range(0.1,0.9,length=7))),),
        colormap = Makie.ColorScheme(get(Makie.ColorSchemes.batlow,range(0.9,0.1,length=512))),
        colormap_var = Makie.ColorScheme(get(Makie.ColorSchemes.devon,range(1.0,0.0,length=512))),
        fonts = Attributes(
            :bold => Makie.texfont(:bold),
            :bolditalic => Makie.texfont(:bolditalic),
            :italic => Makie.texfont(:italic),
            :regular => Makie.texfont(:regular)
        ),
        Figure = (
            backgroundcolor = (:white,0.0),
            fontsize = 9pt
        ),
        Axis = (
            fontsize = 9pt,
            backgroundcolor=(:white,0.0),
            titlesize = 11pt,
            xlabelsize = 9pt, # defults if :fontsize
            ylabelsize = 9pt, # defults if :fontsize
            titlecolor= :black,
            ylabelcolor= :black,
            xlabelcolor= :black,
            spinewidth = 2.0,
            xgridwidth = 1.0,
            ygridwidth = 1.0,
            xgridcolor = :grey45,
            ygridcolor = :grey45,
            xminorgridcolor=:grey45,
            yminorgridcolor=:grey45,
            xminorgridwidth = 0.5,
            yminorgridwidth = 0.5,
            xtickcolor=:black,
            ytickcolor=:black,
            xtickwidth = 2.0,
            ytickwidth = 2.0,
            xticklabelcolor=:black,
            yticklabelcolor=:black,
            xticklabelsize = 8pt,
            yticklabelsize = 8pt,
        ),
        Colorbar = (
            fontsize = 9pt,
            labelsize = 9pt,
            ticklabelsize = 8pt,
            backgroundcolor = (:white,0.0),
            spinewidth=1.5,
        ),
        Legend = (
            fontsize = 9pt,
            labelsize = 8pt,
            backgroundcolor = (:white),
            framecolor = :black,
            framewidth = 2.0,
            rowgap = -5,
            padding = (4f0,4f0,2f0,2f0),
        ),
        ScatterLines = (
            linewidth = 2.0,
            joinstyle = :round,
            linecap = :round
        ),
        Lines = (
            linewidth = 2.0,
            joinstyle = :round,
            linecap = :round
        ),
        PolarAxis = (
            fontsize = 9pt,
            labelsize = 9pt,
            rticklabelsize = 8pt,
            thetaticklabelsize = 8pt,
            backgroundcolor = (:white,0.0),
            rticksvisible = true,
            rticksmirrored = true,
            thetaticksvisible = true,
            rtickcolor = :black,
            rtickwidth = 2.0,
            thetatickcolor = :black,
            thetatickwidth = 2.0,
            gridcolor = :grey45,
            gridwidth = 1.0,
            tickcolor = :black,
            tickwidth = 2.0,
            labelcolor = :black,
            spinewidth = 2.0,
            spinecolor = :black,
            clip = false,
        )
    )
end