using LinearAlgebra
import Plots: plot, plot!
import Base: zero, ==

# Estructuras de elementos geométricos
mutable struct Segment{T <: Real}
    initial_point::Vector{T}
    final_point::Vector{T}
    length::T
end


mutable struct PolygonalLine{T <: Real}
    vertices::Vector{Vector{T}}
    segments::  Vector{Segment{T}} # no ordenados
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

mutable struct Mesh_Frontier{T <: Real}
    shape::Polygon{T}
    segments::Vector{Segment{T}}
    vertices::Vector{Vector{T}}
    right_side::PolygonalLine{T}
    left_side::PolygonalLine{T}
    upper_side::PolygonalLine{T}
    lower_side::PolygonalLine{T}
    length::T
end

mutable struct Mesh{T1 <: Real, T2 <: Real}
    #init_point::Vector{T}
    cells_per_dims::Vector{T1} # in cartesian
    cell_size::Vector{T2}
    cells::Union{Vector{Polygon{T2}}, Vector{Rhomboid{T2}}}
    cells_centers::Vector{Vector{T2}}
    frontier::Mesh_Frontier{T2} #frontier
    include_frontier::Bool # include polygons on frontier
    reflective_segments::Union{Vector{Segment}, Nothing}
end

# Para que el código sea mas rapido, cambiar Mesh{T1 <: Real, T2 <: Real} a Mesh{T <: Real}


## Constructores de Segmento, LineaPoligonal y Poligono

segment(p1::Vector{T}, p2::Vector{T}) where T <: Real = Segment(p1, p2, LinearAlgebra.norm(p1 - p2))

zero(s::Segment{T}) where T <: Real = Segment(zero(s.initial_point), zero(s.final_point), zero(s.length))

function ==(s1::Segment, s2::Segment)
    s1p1 = s1.initial_point
    s1p2 = s1.final_point
    s2p1 = s2.initial_point
    s2p2 = s2.final_point
    
    s1p1, s1p2 = sort([s1p1, s1p2]; by = x -> x[1])
    s2p1, s2p2 = sort([s2p1, s2p2]; by = x -> x[1])
    
    (s1p1 == s2p1) & (s1p2 == s2p2)
end

function polygonalline(s::Vector{Segment{T}}) where T <: Real
    vertices = unique(vcat([[segmento.initial_point, segmento.final_point] for segmento in s]...))
    line_length = sum([segmento.length for segmento in s])
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

function polygon(segments::Vector{Segment{T}}) where T <: Real
    vertices = vcat([[segment.initial_point, segment.final_point] for segment in segments]...)
    perimeter = sum(segment.length for segment in segments)
    area = polygon_area([vertice for vertice in vertices])
    return Polygon(vertices, segments, perimeter, area)
end

function rhomboid(diagonal::Vector{T}, center::Vector{T}, orientation::Val{:vertical}) where T <: Real
    a, b = [minimum(diagonal) / 2, 0.], [0., maximum(diagonal) / 2]
    vertices = [center .- a, center .- b, center .+ a, center .+ b]
    segments = vcat(segment(vertices[1], vertices[end]), [segment(vertices[k], vertices[k+1]) for k in 1:length(vertices)-1])
    Rhomboid(diagonal, vertices, segments, center)
end

rhomboid(diagonal::Vector{T}, center::Vector{T}, orientation::Symbol) where T <: Real  = rhomboid(diagonal, center, Val(orientation))

function mesh_frontier(mesh_shape::Polygon{T}) where T <: Real
    (; vertices, segments, perimeter, area) = mesh_shape
    segment_types = identify_segment_type(mesh_shape)
    segment_types_names = [:right_side, :left_side, :upper_side, :lower_side]
    which_segments = Vector{Vector{Segment{T}}}(undef, length(segment_types_names))
    for (i, s_type) in enumerate(segment_types_names)
        idx = findall(x -> x == s_type, segment_types)
        length(idx) == 0 ? (which_segments[i] = [zero(segments[1])]) : (which_segments[i] = segments[idx])
    end
    polygonal_lines = polygonalline.(which_segments)
    Mesh_Frontier(mesh_shape, segments, vertices, polygonal_lines[1], polygonal_lines[2], polygonal_lines[3], polygonal_lines[4], perimeter)
end

Mesh(cells_per_dims, cell_size, cells, cells_centers, frontier, include_frontier) = Mesh(
    cells_per_dims, cell_size, cells, cells_centers, frontier, include_frontier, nothing
)

## Methods extension for plotting previous structures
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
    plot!()
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


function plot(mesh::Mesh; plot_cell_center::Bool = false, color_cell::Symbol = :black, color_mesh_shape::Symbol = :green, color_cell_center::Symbol = :blue, kwargs...)
    (;cells_per_dims, cell_size, cells, cells_centers, frontier, include_frontier) = mesh
    plot(ratio = :equal, legend = false)
    for cell in cells
        plot!(cell; color = color_cell, kwargs...)
    end
    if plot_cell_center 
        scatter!(Tuple.(cells_centers); color = color_cell_center, kwargs...)
    end
    plot!(frontier.shape; color = color_mesh_shape, kwargs...)
end


function plot!(mesh::Mesh; plot_cell_center::Bool = false, color_cell::Symbol = :black, color_mesh_shape::Symbol = :green, color_cell_center::Symbol = :blue, kwargs...)
    (;cells_per_dims, cell_size, cells, cells_centers, frontier, include_frontier) = mesh
    plot!(ratio = :equal, legend = false)
    for cell in cells
        plot!(cell; color = color_cell, kwargs...)
    end
    if plot_cell_center 
        scatter!(Tuple.(cells_centers); color = color_cell_center, kwargs...)
    end
    plot!(frontier.shape; color = color_mesh_shape, kwargs...)
end


# Function for intersecting geoemtric structures
function cross_product(a::Vector{T}, b::Vector{T}, c::Vector{T})  where T <: Real
    vec = [a, b, c]
    p = [vcat(vec[k], [0.]) for k in 1:length(vec)]
    cross(p[1] - p[3], p[2] - p[3]) 
end

function intersection(p::Vector{T}, s::Segment{T})  where T <: Real
    (min(s.initial_point[1], s.final_point[1]) <= p[1] <= max(s.initial_point[1], s.final_point[1])) && (min(s.initial_point[2], s.final_point[2]) <= p[2] <= max(s.initial_point[2], s.final_point[2]))
end

intersection(s::Segment{T}, p::Vector{T}) where T <: Real = intersection(p, s)

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

vector_segment_intersection(vec::Vector{T}, vec_initial_point::Vector{T}, seg::Segment{T}, dim::Val{:x}) where T <: Real = (seg.initial_point[1] - vec_initial_point[1]) / (vec[1] - seg.final_point[1] + seg.initial_point[1])
vector_segment_intersection(vec::Vector{T}, vec_initial_point::Vector{T}, seg::Segment{T}, dim::Val{:y}) where T <: Real = (seg.initial_point[2] - vec_initial_point[2]) / (vec[2] - seg.final_point[2] + seg.initial_point[2])
vector_segment_intersection(vec::Vector{T}, vec_initial_point::Vector{T}, seg::Segment{T}, dim::Symbol) where T <: Real = vector_segment_intersection(vec, vec_initial_point, seg, Val(dim))

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

function polygon_on_frontier(polygon::Rhomboid{T}, frontier::Union{Polygon{T}, Rhomboid{T}}) where T <: Real
    center_on_frontier = point_on_polygon_frontier(polygon.center, frontier)
end

function find_idx_center_polygons_on_frontier(points::Vector{Vector{T}}, frontier::Union{Polygon{T}, Rhomboid{T}}, polygon_type::Val{:rhomboid}) where T <: Real
    # points are the centers, polygon is the mesh shape
    idx_points_on_polygon_frontier = findall([point_on_polygon_frontier(point, frontier) for point in points])
    points_on_polygon_frontier = points[idx_points_on_polygon_frontier]
    return (idx_points_on_polygon_frontier, points_on_polygon_frontier)
end

find_idx_center_polygons_on_frontier(points::Vector{Vector{T}}, frontier::Union{Polygon{T}, Rhomboid{T}}, polygon_type::Symbol) where T <: Real = find_idx_center_polygons_on_frontier(points, frontier, Val(polygon_type))

function find_polygons_on_frontier(mesh::Mesh{T2, T1}) where {T1<:Real, T2<:Real}
    (; cells_per_dims, cell_size, cells, cells_centers, frontier, include_frontier) = mesh
    idx, center_points = find_idx_center_polygons_on_frontier(cells_centers, frontier.shape, :rhomboid)
    cells_on_frontier = cells[idx]
    centers = [c.center for c in cells_on_frontier]
    min_y = minimum(x -> x[2], centers)
    max_y = maximum(x -> x[2], centers)
    cells_frontier_dict = Dict()
    cells_frontier_dict[:lower] = [cell for cell in cells_on_frontier if cell.center[2] == min_y]
    cells_frontier_dict[:upper] = [cell for cell in cells_on_frontier if cell.center[2] == max_y]
    return idx, center_points, cells_frontier_dict
end

function get_frontier(mesh::Mesh{T1,T2}; get_array = false) where {T1<:Real,T2<:Real}
    (; cells_per_dims, cell_size, cells, cells_centers, frontier, include_frontier) = mesh
    centers_min_x = minimum(x -> x[1], cells_centers)
    centers_max_x = maximum(x -> x[1], cells_centers)
    centers_min_y = minimum(x -> x[2], cells_centers)
    centers_max_y = maximum(x -> x[2], cells_centers)
    cells_frontier_dict = Dict()
    cells_frontier_dict[:lower] = [cell for cell in cells if cell.center[2] == centers_min_y]
    cells_frontier_dict[:upper] = [cell for cell in cells if cell.center[2] == centers_max_y]
    cells_frontier_dict[:left] = [cell for cell in cells if cell.center[1] == centers_min_x]
    cells_frontier_dict[:right] = [cell for cell in cells if cell.center[1] == centers_max_x]
    if get_array
        return [cell for cell_arr in values(cells_frontier_dict) for cell in cell_arr]
    else
        return cells_frontier_dict
    end
end

## Mesh generation auxiliary functions

# Definir metodos apra estas funciones donde solo se tenga un tipo de Real
find_cell(position::Vector{T1}, mesh::Mesh{T2, T1}) where {T1 <: Real, T2 <: Real} = findall([point_inside_polygon(polygon, position) for polygon in mesh.cells])
find_cell(position::Vector{T1}, mesh::Mesh{T2, T1}) where {T1 <: Real, T2 <: Real} = findall([point_inside_polygon(polygon, position) for polygon in mesh.cells])
find_cell(mesh::Mesh{T2, T1}, position::Vector{T1}) where {T1 <: Real, T2 <: Real} = find_cell(position, mesh)
find_cell(mesh::Mesh{T2, T1}, position::Vector{T1}) where {T1 <: Real, T2 <: Real} = find_cell(position, mesh)
find_cell(position::Vector{T}, cells::Vector{Rhomboid{T}}) where T <: Real = findall([point_inside_polygon(polygon, position) for polygon in cells])[1]
find_cell(position::Vector{T}, cells::Vector{Polygon{T}}) where T <: Real = findall([point_inside_polygon(polygon, position) for polygon in cells])[1]
find_cell(cells::Vector{Polygon{T}}, position::Vector{T}) where T <: Real = find_cell(position, cells)
find_cell(cells::Vector{Rhomboid{T}}, position::Vector{T}) where T <: Real = find_cell(position, cells)

""" 
    Identifica el numero de celda en un mesh dado. 
"""
function find_cell(cell::Rhomboid{T1}, mesh::Mesh{T2, T1}) where {T1 <: Real, T2 <: Real} #lo que identifica a cada rombo es su centro
    index = 0
    (; diagonal, vertices, segments, center) = cell
    for (i, cell_center) in enumerate(mesh.cells_centers)
        if cell_center == center
            index = i
        end
    end
    return index
end

function generate_rhomboid_mesh(rhomboids_per_dims::Vector{T1}, rhomboid_diagonal::Vector{T2}, mesh_shape::Polygon{T2}; init_point::Vector{T2} = [0., 0.], include_frontier::Bool = true, min_val::Vector{Int} = [0, 0,]) where {T1 <: Real, T2 <: Real}# side represents minor diagonal
    # init_point == [0., 0.] ? init_point = [-rhomboid_diagonal[1] * (rhomboids_per_dims[1] - 0.5), 0.] : init_point = init_point
    L = rhomboids_per_dims .* rhomboid_diagonal
    a = [rhomboid_diagonal[1] * 0.5, 0.]
    b = [0, rhomboid_diagonal[2] * 0.5]
    v1 = a .+ b
    v2 = a .- b
    
    rhomboid_center_points = [init_point .+ (m * v1) .+ (n .* v2) for m in min_val[1]:maximum(rhomboids_per_dims) + 1 for n in min_val[2]:maximum(rhomboids_per_dims) + 1]
    idx_center_points_inside = findall([point_inside_polygon(point, mesh_shape; include_frontier = include_frontier) for point in rhomboid_center_points])
    rhomboid_center_points_inside = rhomboid_center_points[idx_center_points_inside]
    rhomboids = [rhomboid(rhomboid_diagonal, center, :vertical) for center in rhomboid_center_points_inside];
    frontier = mesh_frontier(mesh_shape)
    mesh = Mesh(rhomboids_per_dims, rhomboid_diagonal, rhomboids, rhomboid_center_points_inside, frontier, include_frontier)
    mesh_segments = [s for cell in mesh.cells for s in cell.segments]
    reflective_segments = [s for s in mesh_segments if s in mesh_shape.segments]
    mesh.reflective_segments = reflective_segments
    return mesh
end

function generate_rhomboid_mesh(rhomboids_per_dims::Vector{T}, rhomboid_diagonal::Vector{T}, mesh_shape::Polygon{T}; init_point::Vector{T} = [0., 0.], include_frontier::Bool = true, min_val::Vector{Int} = [0, 0,]) where T <: Real# side represents minor diagonal
    # init_point == [0., 0.] ? init_point = [-rhomboid_diagonal[1] * (rhomboids_per_dims[1] - 0.5), 0.] : init_point = init_point
    L = rhomboids_per_dims .* rhomboid_diagonal
    a = [rhomboid_diagonal[1] * 0.5, 0.]
    b = [0, rhomboid_diagonal[2] * 0.5]
    v1 = a .+ b
    v2 = a .- b
    
    rhomboid_center_points = [init_point .+ (m * v1) .+ (n .* v2) for m in min_val[1]:maximum(rhomboids_per_dims) + 1 for n in min_val[2]:maximum(rhomboids_per_dims) + 1]
    idx_center_points_inside = findall([point_inside_polygon(point, mesh_shape; include_frontier = include_frontier) for point in rhomboid_center_points])
    rhomboid_center_points_inside = rhomboid_center_points[idx_center_points_inside]
    rhomboids = [rhomboid(rhomboid_diagonal, center, :vertical) for center in rhomboid_center_points_inside];
    frontier = mesh_frontier(mesh_shape)
    mesh = Mesh(rhomboids_per_dims, rhomboid_diagonal, rhomboids, rhomboid_center_points_inside, frontier, include_frontier)
    mesh_segments = [s for cell in mesh.cells for s in cell.segments]
    reflective_segments = [s for s in mesh_segments if s in mesh_shape.segments]
    mesh.reflective_segments = reflective_segments
    return mesh
end


function generate_rectangular_stripe(cells_per_dims::Vector{T1}, cell_size::Vector{T2}, init_point::Vector{T2}) where {T1 <: Real, T2 <: Real}#cell dims repressents the dimensions of an unitary cell 
    # init_point == [0., 0.] ? init_point = [-cell_size[1] * (cells_per_dims[1] - 0.5), 0.] : init_point = init_point
    L = cells_per_dims .* cell_size
    p1 = init_point .+ [0., - cell_size[2]]
    p2 = init_point .+ [0., cell_size[2]]
    p4 = init_point .+ [cell_size[1] + L[1], - cell_size[2]]
    p3 = init_point .+ [cell_size[1] + L[1], cell_size[2]]

    rectangular_stripe = polygon([p1, p2, p3, p4])     
end

function generate_rectangular_stripe(cells_per_dims::Vector{T}, cell_size::Vector{T}, init_point::Vector{T}) where T <: Real #cell dims repressents the dimensions of an unitary cell 
    # init_point == [0., 0.] ? init_point = [-cell_size[1] * (cells_per_dims[1] - 0.5), 0.] : init_point = init_point
    L = cells_per_dims .* cell_size
    p1 = init_point .+ [0., - cell_size[2]]
    p2 = init_point .+ [0., cell_size[2]]
    p4 = init_point .+ [cell_size[1] + L[1], - cell_size[2]]
    p3 = init_point .+ [cell_size[1] + L[1], cell_size[2]]

    rectangular_stripe = polygon([p1, p2, p3, p4])     
end

# This only work for rhomboid mesh and Lorentz gas like systemsm with two rows of rhomboidal cells
function generate_frontier_lorentz_gas(
    rhomboids_per_dims::Vector{T1},
    rhomboid_diagonal::Vector{T2};
    init_point::Vector{T2} = [0., 0.],
    include_frontier::Bool = false,
    min_val::Vector{Int} = [0, 0,]
) where {T1 <: Real, T2 <: Real}# side represents minor diagonal
# init_point == [0., 0.] ? init_point = [-rhomboid_diagonal[1] * (rhomboids_per_dims[1] - 0.5), 0.] : init_point = init_point
rectangular_stripe = generate_rectangular_stripe(rhomboids_per_dims, rhomboid_diagonal, init_point)
mesh = generate_rhomboid_mesh(rhomboids_per_dims, rhomboid_diagonal, rectangular_stripe; init_point = init_point, include_frontier = include_frontier, min_val = min_val)
cells = mesh.cells
left_side = [cells[1].segments[1], cells[1].segments[2], cells[2].segments[1], cells[2].segments[2]]
right_side = [cells[end].segments[3], cells[end].segments[4], cells[end-1].segments[3], cells[end- 1].segments[4]]
upper_side = [segment(cells[2].segments[1].final_point, cells[end].segments[4].final_point)]
lower_side = [segment(cells[1].segments[2].final_point, cells[end-1].segments[3].initial_point)]
segments = vcat(left_side, right_side, upper_side, lower_side)
return polygon(segments)
end

function identify_segment_type(frontier::Polygon{T}) where T <: Real
    (; vertices, segments, perimeter, area) = frontier
    min_x = findmin(x -> x[1], vertices)[1]
    max_x = findmax(x -> x[1], vertices)[1]
    min_y = findmin(x -> x[2], vertices)[1]
    max_y = findmax(x -> x[2], vertices)[1]
    segment_types_dict = Dict(:min_x => :left_side, :max_x => :right_side, :min_y => :lower_side, :max_y => :upper_side)
    segment_types = Array{Symbol}(undef, length(segments))
    for(i, segment) in enumerate(segments)
        if ((segment.initial_point[1] == min_x) || (segment.final_point[1] == min_x))
            segment_types[i] = segment_types_dict[:min_x]
        elseif ((segment.initial_point[1] == max_x) || (segment.final_point[1] == max_x))
            segment_types[i] = segment_types_dict[:max_x]
        elseif ((segment.initial_point[2] == min_y) || (segment.final_point[2] == min_y))
            segment_types[i] = segment_types_dict[:min_y]
        elseif ((segment.initial_point[2] == max_y) || (segment.final_point[2] == max_y))
            segment_types[i] = segment_types_dict[:max_y]
        else
        
        end     
    end
    return segment_types
end


## Other functions

is_horizontal(seg::Segment{T}) where T <: Real = (seg.initial_point[2] == seg.final_point[2])
is_vertical(seg::Segment{T}) where T <: Real = (seg.initial_point[1] == seg.final_point[1])

function is_vertical(line_seg::PolygonalLine{T}) where T <: Real
    (; vertices, segments, length) = line_seg
    length(unique(x -> x[1], vertices)) < length(vertices) ? true : false
end


function translate_mesh(mesh::Mesh, trans_vec::Vector)
    cells_per_dims = mesh.cells_per_dims
    cell_size = mesh.cell_size
    include_frontier = mesh.include_frontier
    frontier = generate_frontier_lorentz_gas(cells_per_dims, cell_size; init_point = trans_vec)
    mesh_trans = generate_rhomboid_mesh(cells_per_dims, cell_size, frontier; include_frontier = include_frontier, init_point = trans_vec)
    
    return mesh_trans
end

function get_extended_mesh(mesh::Mesh)
    m, n = mesh.cells_per_dims
    d1, d2 = mesh.cell_size
    trans_vec_arr = [[-m * d1, 0.], [m * d1, 0.], [0., -2 * d2], [0., 2 * d2]] # only works for 2 rows
    
    [(trans_vec, translate_mesh(mesh, trans_vec)) for trans_vec in trans_vec_arr]
end

function get_trans_vectors(diagonal)
    a = [diagonal[1] * 0.5, 0.]
    b = [0., diagonal[2] * 0.5]
    v1 = a + b
    v2 = a - b
    return v1, v2
end

function find_neighbors(cell::Rhomboid{T2}, mesh::Mesh{T1, T2}) where {T1 <: Real, T2 <: Real}
    (; diagonal, vertices, segments, center) = cell
    v1, v2 = get_trans_vectors(diagonal)
    neighbors_centers = [
        center + (m * v1) + (n * v2)
        for m in -1:1, n in -1:1
        if !((m == 0) & (n == 0))
    ]
    indices = Vector{Int}(undef, length(neighbors_centers))
    for (idx, center) in enumerate(neighbors_centers)
        cell_idx = find_cell(center, mesh)
        length(cell_idx) == 0 ? (indices[idx] = -1) : (indices[idx] = cell_idx[1])
    end
    
    on_extended_mesh_indices = findall(x -> x == -1, indices)
    on_extended_mesh_centers = neighbors_centers[on_extended_mesh_indices]
    
    extended_mesh = get_extended_mesh(mesh)
    for idx in on_extended_mesh_indices
        center = neighbors_centers[idx]
        for (trans_vec, ex_mesh) in extended_mesh
            extended_idx = find_cell(center, ex_mesh)
            if length(extended_idx) != 0
                extended_center = ex_mesh.cells_centers[extended_idx[1]]
                inv_trans_center = extended_center - trans_vec
                original_idx = find_cell(inv_trans_center, mesh)
                if length(original_idx) != 0
                    indices[idx] = original_idx[1]
                    neighbors_centers[idx] = inv_trans_center
                end
            end
        end
    end
    
    indices, neighbors_centers
end
