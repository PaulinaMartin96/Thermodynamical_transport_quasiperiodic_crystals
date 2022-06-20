using LinearAlgebra
import Plots: plot, plot!

# Estructuras de elementos geométricos
mutable struct Segment{T <: Real}
    initial_point::Vector{T}
    final_point::Vector{T}
    length::T
end


mutable struct PolygonalLine{T <: Real}
    vertices::Array{Vector{T}}
    segments::Array{Segment{T}} # no ordenados
    length::T
end

mutable struct Polygon{T <: Real}
    vertices::Array{Vector{T}}
    segments::Array{Segment{T}} # no ordenados
    perimeter::T
    area::T
end

mutable struct Rhomboid{T <: Real}
    diagonal::Vector{T}
    vertices::Vector{Vector{T}}
    segments::Vector{Segment{T}}
    center::Vector{T}
end

mutable struct PolyhedricSurface{T <: Real}
    polygons::Array{Polygon{T}}
    surface::T
end

mutable struct Polyhedron{T <: Real}
    faces::Array{Polygon{T}}
    surface::PolyhedricSurface{T}
    volume::T
end

# Constructores de Segmento, LineaPoligonal y Poligono

segment(p1::Vector{T}, p2::Vector{T}) where T <: Real = Segment(p1, p2, LinearAlgebra.norm(p1 - p2))

function polygonalline(s::Array{Segment{T}}) where T <: Real
    vertices = unique(vcat([[segmento.initial_point, segmento.final_point] for segmento in s]...))
    line_length = sum([segmento.longitud for segmento in s])
    PolygonalLine(vertices, s, line_length)
end

function polygon_area(vertices::Vector{Vector{T}}) where T <: Real
    1/2*abs(sum(vcat([vertices[end][1]*vertices[1][2] - vertices[end][2]*vertices[1][1]], [vertices[k][1]*vertices[k+1][2] - vertices[k][2]*vertices[k+1][1] for k in 1:length(vertices)-1])))
end


function polygon(vertices::Vector{Vector{T}}) where T <: Real
    aristas = vcat(segment(vertices[1], vertices[end]), [segment(vertices[k], vertices[k+1]) for k in 1:length(vertices)-1])
    perimetro = sum([arista.length for arista in aristas])
    area = polygon_area([vertice for vertice in vertices])
    return Polygon(vertices, aristas, perimetro, area)
end

function rhomboid(diagonal::Vector{T}, center::Vector{T}, orientation::Val{:vertical}) where T <: Real
    a, b = [minimum(diagonal) / 2, 0.], [0., maximum(diagonal) / 2]
    vertices = [center .- a, center .- b, center .+ a, center .+ b]
    segments = vcat(segment(vertices[1], vertices[end]), [segment(vertices[k], vertices[k+1]) for k in 1:length(vertices)-1])
    Rhomboid(diagonal, vertices, segments, center)
end

rhomboid(diagonal, center, orientation::Symbol) = rhomboid(diagonal, center, Val(orientation))

# Extesión de métodos para graficar las estructuras previas
"""
La funcion ´plot` a partir de un segmento y genera una grafica
"""
function plot(s::Segment{T}; kw...) where T <: Real
    inicio = s.initial_point
    termino = s.final_point
    if length(inicio) == length(termino) == 2
        g = plot([inicio[1], termino[1]], [inicio[2], termino[2]]; kw...)
    else
        g = plot([inicio[1], termino[1]], [inicio[2], termino[2]], [inicio[3], termino[3]]; kw...)
    end

    return g
end

"""
La funcion ´plot!` a partir de un segmento y genera una grafica más elaborada 
"""
function plot!(s::Segment{T}; kw...) where T <: Real
    inicio = s.initial_point
    termino = s.final_point
    if length(inicio) == length(termino) == 2
        g = plot!([inicio[1], termino[1]], [inicio[2], termino[2]]; kw...)
    else
        g = plot!([inicio[1], termino[1]], [inicio[2], termino[2]], [inicio[3], termino[3]]; kw...)
    end

    return g
end

#### Ahora extendemos los métodos para LineaPoligonal

"""
La función ´plot` a partir de un objeto LineaPoligonal genera la grafica de esa línea.
"""
function plot(lp::PolygonalLine{T}; kw...) where T <: Real
    p = plot()
    for segment in lp.segments
        plot!(segment; kw...)
    end
    return p
end


"""
La función ´plot!` a partir de un objeto LineaPoligonal genera la grafica de las lineas poligonales que desee el usuario.
"""
function plot!(lp::PolygonalLine{T}; kw...) where T <: Real
    for segment in lp.segments
        plot!(segment; kw...)
    end
end


"""
La función ´plot` para este caso recibe un poligono y arroja la gráfica de dicho objeto.
"""
function plot(p::Union{Polygon{T}, Rhomboid{T}}; kw...) where T <: Real
    plot()
    for segment in p.segments
        plot!(segment; kw...)
    end
    plot!()
end


"""
La función ´plot` para este caso recibe un poligono y arroja la gráfica de dicho objeto sobre algun otro ya dibujado
"""
function plot!(p::Union{Polygon{T}, Rhomboid{T}}; kw...)  where T <: Real
    for segment in p.segments
        plot!(segment; kw...)
    end
    plot!()
end


# Function for intersecting geoemtric structures

function point_on_segment(p::Vector{T}, s::Segment{T}; epsilon::T = 1e-8, debug = false)  where T <: Real
    x1, y1 = s.initial_point
    x2, y2 = s.final_point
    x, y = p
    
    Δx1 = x - x1
    Δy1 = y - y1
    Δx2 = x2 - x
    Δy2 = y2 - y
    
    cross = Δx1 * Δy2 - Δy1 * Δx2
    debug && @show cross
    x_min = min(x1, x2)
    x_max = max(x1, x2)
    y_min = min(y1, y2)
    y_max = max(y1, y2)
    (x_min <= x <= x_max) & (y_min <= y <= y_max) & (abs(cross) < epsilon) ? true : false
end

point_on_segment(s::Segment{T}, p::Vector{T}; kwargs...) where T <: Real = point_on_segment(p, s; kwargs...)

function intersection(p::Vector{T}, s::Segment{T})  where T <: Real
    (min(s.initial_point[1], s.final_point[1]) <= p[1] <= max(s.initial_point[1], s.final_point[1])) && (min(s.initial_point[2], s.final_point[2]) <= p[2] <= max(s.initial_point[2], s.final_point[2]))
end

intersection(s::Segment{T}, p::Vector{T}) where T <: Real = intersection(p, s)

function cross_product(a::Vector{T}, b::Vector{T}, c::Vector{T})  where T <: Real
    vec = [a, b, c]
    p = [vcat(vec[k], [0.]) for k in 1:length(vec)]
    cross(p[1] - p[3], p[2] - p[3]) 
end

function intersection(s1::Segment{T}, s2::Segment{T})  where T <: Real
    p1 = s1.initial_point
    p2 = s1.final_point
    p3 = s2.initial_point
    p4 = s2.final_point
    
    d1 = cross_product(p1,p4,p3)[3]
    d2 = cross_product(p2,p4,p3)[3]
    d3 = cross_product(p3,p2,p1)[3]
    d4 = cross_product(p4,p2,p1)[3]
   
    if ((d1 > 0 && d2 < 0) || (d1 < 0 && d2 > 0)) && ((d3 > 0 && d4 < 0) || (d3 < 0 && d4 > 0))
         intersect = true
    elseif d1 == 0 && intersection(p1, segment(p3, p4))
        intersect = true
    elseif d2 == 0 && intersection(p2, segment(p3, p4))
        intersect =  true
    elseif d3 == 0 && intersection(p3, segment(p1, p2))
        intersect = true
    elseif d4 == 0 && intersection(p4, segment(p1, p2))
        intersect =  true
    else
        intersect = false # solo este es false
    end
   return intersect
end

function point_on_polygon_frontier(point::Vector{T}, polygon::Union{Polygon{T}, Rhomboid{T}}) where T <: Real
    on_polygon_frontier = [point_on_segment(segment, point) for segment in polygon.segments]
    any(on_polygon_frontier) ? true : false
end

point_on_polygon_frontier(polygon::Union{Polygon{T}, Rhomboid{T}}, point::Vector{T}) where T <: Real = point_on_polygon_frontier(point, polygon)

function point_inside_polygon(point::Vector{T}, polygon::Union{Polygon{T}, Rhomboid{T}}; x_inf::Real = 1e13, y_inf = 0.4, include_frontier = true)  where T <: Real
    on_polygon_frontier = [point_on_segment(segment, point) for segment in polygon.segments]
    if (any(on_polygon_frontier) && include_frontier)
        return true
    elseif (any(on_polygon_frontier) && include_frontier == false)
        return false
    else
        point_inf = [x_inf, rand()]
        semiray = segment(point, point_inf)
        intersections = [intersection(semiray, edge) for edge in polygon.segments]
        isodd(sum(intersections)) ? true : false
    end
end
    
#function point_inside_polygon(point::Vector{T}, polygon::Union{Polygon{T}, Rhomboid{T}}; x_inf::Real = 1e13, y_inf = 0.4, include_frontier = true)  where T <: Real
#    on_polygon_frontier = [point_on_segment(segment, point) for segment in polygon.segments]
#    any(on_polygon_frontier) && return true
#    point_inf = [x_inf, rand()]
#    semiray = segment(point, point_inf)
#    intersections = [intersection(semiray, edge) for edge in polygon.segments]
#    isodd(sum(intersections)) ? true : false
#end

point_inside_polygon(polygon::Union{Polygon{T}, Rhomboid{T}}, point::Vector{T}; kwargs...) where T <: Real = point_inside_polygon(point, polygon; kwargs...)