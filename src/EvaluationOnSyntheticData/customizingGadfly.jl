paperTheme = Gadfly.Theme( default_color="salmon",
                           point_size=1mm,  #
                           point_size_min=1mm, #
                           point_size_max=1mm,  #
                           #discrete_sizemap, 
                           #continuous_sizemap,
                           point_shapes=[Gadfly.square,],
                           line_width=1mm, 
                           line_style=[:solid,],
                           alphas=[1.0,], #
                           panel_fill="white",
                           panel_stroke="black",
                           #panel_opacity=0.0,
                           background_color="white",
                           plot_padding=[0mm,20mm,10mm,0mm],
                           grid_color="grey36",
                           grid_line_style=:dashdotdot,
                           grid_line_width=0.1mm,
                           minor_label_font="Arial",
                           minor_label_font_size=5mm,
                           minor_label_color="grey15",
                           major_label_font="sans-serif",
                           major_label_font_size=6mm,
                           major_label_color="black",
                           #point_label_font="Helvetica",
                           #point_label_font_size=8mm,
                           #point_label_color="white",
                           key_title_font="Arial",
                           key_title_font_size=4mm,
                           key_title_color="black",
                           key_label_font="Arial",
                           key_label_font_size=4mm,
                           key_label_color="black", 
                           #key_color_gradations,
                           bar_spacing=0mm, 
                           #boxplot_spacing=0mm,
                           #errorbar_cap_length,
                           highlight_width=0.2mm,
                           discrete_highlight_color=(x)->"black",
                           #continuous_highlight_color,
                           #lowlight_color,
                           #middle_color,
                           #middle_width,
                           guide_title_position=:center, 
                           colorkey_swatch_shape=:circle,
                           #key_swatch_shape,
                           #key_swatch_color="salmon",
                           #key_swatch_size,
                           key_position=:inside,
                           #bar_highlight, 
                           rug_size=0mm,
                           )

# using Colors, Measures

const halfpagewidth =  85mm;  
const fullpagewidth = 170mm; 
const maxplotheight = 225mm; 

struct custom_manual_color_key <: Gadfly.GuideElement
    title::AbstractString
    labels::Vector{String}
    pos::Vector
    colors::Vector{Colorant}
    shapes::Vector{Function}
    sizes::Vector{Measure}
    visible::Bool
end

# should the inputs be arrays?
function custom_manual_color_key(; title="", pos=[], labels=String[], color=Colorant[], shape=Function[], size=Measure[] )

    ncolors, nshapes, nsizes = length(color), length(shape), length(size)
    n = max(ncolors, nshapes, nsizes)
    ncolors == 1  &&  (color = repeat(color, n)) 
    nshapes == 1  &&  (shape = repeat(shape, n))
    nsizes  == 1  &&  (size  = repeat(size , n))
    # only repeat if length(arr) == 1, otherwise it would have n*length(arr) elements after repeat( ... ); 

    theme = paperTheme;
    
    CT, ST, SZT = eltype(color), eltype(shape), eltype(size)
    clrs = CT  <: Int ? theme.discrete_colormap(maximum(color))[color] : Gadfly.parse_colorant(color) 
    szs  = SZT <: Int ? theme.discrete_sizemap(maximum(size))[size] : size
    shps = ST  <: Int ? theme.point_shapes[shape] : shape
    # notes: 
    # theme.discrete_colormap(maximum(color))[color] <- interesting
    # Gadfly.parse_colorant is imported from Compose.jl
    # discrete_sizemap -> scales.default_discrete_sizes(n) -> range(0.45mm, 1.8mm, length=n)


    cataes   = [clrs, shps, szs]
    notempty = .!isempty.(cataes)
    if any(notempty)
        swatches = collect(Tuple, zip(cataes[notempty]...))
        !allunique(swatches) && error("Swatches should not be repeated in a manual key")
    end
    return custom_manual_color_key(title, labels, pos, clrs, shps, szs,  true)
end

function render(guide::custom_manual_color_key, theme::Gadfly.Theme, aes::Gadfly.Aesthetics)

    gpos = guide.pos

    (theme.key_position == :inside) && isempty(gpos) &&  (gpos = [0.7w, 0.25h])

    title_context, title_width = render_key_title2(guide.title, theme)

    ctxs = render_discrete_key(
                              guide.labels, 
                              title_context, 
                              title_width, 
                              theme, 
                              shapes=guide.shapes,
                              colors=guide.colors, 
                               sizes=guide.sizes)
    
    position, stackable = right_guide_position, true
    if !isempty(gpos)
        position, stackable = over_guide_position, false
        ctxs = [ compose(
                        context(), 
                        (context(gpos[1],gpos[2]), ctxs[1]) )]
    elseif theme.key_position == :left
        position = left_guide_position
    elseif theme.key_position == :top
        position = top_guide_position
    elseif theme.key_position == :bottom
        position = bottom_guide_position
    end

    return [PositionedGuide(ctxs, 0, position, stackable )]
end