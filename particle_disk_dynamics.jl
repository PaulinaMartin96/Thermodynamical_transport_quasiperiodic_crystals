import Base.zero
using LinearAlgebra

## Structures
mutable struct Hard_Disk{T <: Real}
    mass::T
    radius::T
    position::Vector{T}
    #velocity::Vector{T}
    angular_velocity::Vector{T}
    cell::Int
end

mutable struct Particle{T <: Real}
    mass::T
    position::Vector{T}
    velocity::Vector{T}
    normal_velocity::Vector{T}
    tangential_velocity::Vector{T}
    cell::Int
end


## Hard_disk constructors
 # We used two different functions instead of one function with Union in the arguments because a vector of Union of two types makes the code slower
 # This could be solved by defining an Abstract type called Polygon that includes both Polygon and Rhomboid
function hard_disk(mass::T1, radius::T1, position::Vector{T1}, angular_velocity::Vector{T1}, mesh::Mesh{T2, T1}) where {T1 <: Real, T2 <: Real}
    n_cell = find_cell(position, mesh)[1]
    Hard_Disk(mass, radius, position, angular_velocity, n_cell)
end

function hard_disk(mass::T, radius::T, position::Vector{T}, angular_velocity::Vector{T}, mesh::Mesh{T}) where T <: Real
    n_cell = find_cell(position, mesh)[1]
    Hard_Disk(mass, radius, position, angular_velocity, n_cell)
end

function zero(d::Hard_Disk{T}) where T <: Real
    x = d.mass
    y = d.position
    Hard_Disk(zero(x), zero(x), zero(y), zero(y), 1)
end


## Particle constructors
#Agregar emtodos con T1 y T2 comutados para los argumentos. Esto es para Mesh

function particle(mass::T1, position::Vector{T1}, velocity::Vector{T1}, mesh::Mesh{T2, T1}) where {T1 <: Real, T2 <: Real}
    normal_vec = find_normal_vector(velocity) # In two dimensions, the nullspace() retunrs 
    tangential_vec = find_tangential_vector(velocity)
    normal_velocity = projection(velocity, normal_vec) .* normal_vec
    tangential_velocity = projection(velocity, tangential_vec) .* tangential_vec
    n_cell = find_cell(position, mesh)[1]
    Particle(mass, position, velocity, normal_velocity, tangential_velocity, n_cell)
end

function particle(mass::T, position::Vector{T}, velocity::Vector{T}, mesh::Mesh{T}) where T <: Real
    normal_vec = find_normal_vector(velocity) # In two dimensions, the nullspace() retunrs 
    tangetial_vec = find_tangential_vector(velocity)
    normal_velocity = projection(velocity, normal_vec) .* normal_vec
    tangential_velocity = projection(velocity, tangential_vec) .* tangential_vec
    n_cell = find_cell(position, mesh)[1]
    Particle(mass, position, velocitiy, normal_velocity, tangetial_velocity, n_cell)
end

function particle(mass::T1, position::Vector{T1}, normal_velocity::Vector{T1}, tangential_velocity::Vector{T1}, mesh::Mesh{T2, T1}) where {T1 <: Real, T2 <: Real}
    velocity = normal_velocity .+ tangential_velocity
    n_cell = find_cell(position, mesh)[1]
    Particle(mass, position, velocitiy, normal_velocity, tangetial_velocity, n_cell)
end

function particle(mass::T, position::Vector{T}, normal_velocity::Vector{T}, tangential_velocity::Vector{T}, mesh::Mesh{T}) where T <: Real
    velocity = normal_velocity .+ tangential_velocity
    n_cell = find_cell(position, mesh)[1]
    Particle(mass, position, velocitiy, normal_velocity, tangetial_velocity, n_cell)
end

function zero(p::Particle{T}) where T <: Real
    x = p.mass
    y = p.position
    Particle(zero(x), zero(x), zero(y), zero(y), 1)
end

function generate_hard_disk_vector(mass::T1, radius::T1, positions::Vector{Vector{T1}}, angular_velocities::Vector{Vector{T1}}, mesh::Mesh{T2, T1}) where {T1 <: Real, T2 <: Real}
    n = length(positions)
    hard_disks_vec = Vector{Hard_Disk{T1}}(undef, n)
    for i in 1:n
        hard_disks_vec[i] = hard_disk(mass, radius, positions[i], angular_velocities[i], mesh)
    end
    return hard_disks_vec
end

## Initialize a vector of hard disks or particles

function generate_hard_disk_vector(mass::T, radius::T, positions::Vector{Vector{T}}, angular_velocities::Vector{Vector{T}}, mesh::Mesh{T}) where T <: Real
    n = length(positions)
    hard_disks_vec = Vector{Hard_Disk{T}}(undef, n)
    for i in 1:n
        hard_disks_vec[i] = hard_disk(mass, radius, positions[i], angular_velocities[i], mesh)
    end
    return hard_disks_vec
end

function generate_particle_vector(mass::T1, positions::Vector{Vector{T1}}, velocities::Vector{Vector{T1}}, mesh::Mesh{T2, T1}) where {T1 <: Real, T2 <: Real}
    n = length(positions)
    particles_vec = Vector{Particle{T1}}(undef, n)
    for i in 1:n
        particles_vec[i] = particle(mass, positions[i], velocities[i], mesh)
    end
    return particles_vec
end

function generate_particle_vector(mass::T, positions::Vector{Vector{T}}, velocities::Vector{Vector{T}}, mesh::Mesh{T}) where T <: Real
    n = length(positions)
    particles_vec = Vector{Particle{T}}(undef, n)
    for i in 1:n
        particles_vec[i] = particle(mass, positions[i], velocities[i], mesh)
    end
    return particles_vec
end


## Dynamics auxiliary functions

MRU(initial_position::Vector{T}, velocity::Vector{T}, time::T, dim::Val{:x}) where T <: Real = (velocity[1] * time) + initial_position[1]
MRU(initial_position::Vector{T}, velocity::Vector{T}, time::T, dim::Val{:y}) where T <: Real = (velocity[2] * time) + initial_position[2]
MRU(initial_position::Vector{T}, velocity::Vector{T}, time::T, s::Symbol) where T <: Real = MRU(time, velocity, initial_position, Val(s))

# Ecuaci??n param??trica de la circunferencia
x_disk(center::Vector{T}, time::T, radius::T) where T <: Real = center[1] + 2*radius*cos(time)
y_disk(center::Vector{T}, time::T, radius::T) where T <: Real = center[2] + 2*radius*sin(time)

# Ecuaciones de movimeinto de la part??cula (MRU)
x_particle(position::Vector{T}, velocity::Vector{T}, time::T) where T <: Real = velocity[1] * time + position[1]
y_particle(position::Vector{T}, velocity::Vector{T}, time::T) where T <: Real = velocity[2] * time + position[2]


# Ecuaciones de movimiento
#x_collision(center::Vector{T}, velocity::Vector{T}, position::Vector{T}, t::T) where T <: Real = center[1] + 2*r*cos(t) - velocity[1]*t - position[1] # x denota el vector de posiciones iniciales
#y_collision(center::Vector{T}, velocity::Vector{T}, position::Vector{T}, t::T) where T <: Real= center[2] + 2*r*sin(t) - velocity[2]*t - position[2]

function projection(a::Vector{T}, b::Vector{T}) where T <: Real # projection of vector a over vector b
    proj = LinearAlgebra.dot(a, b) ./ LinearAlgebra.norm(b)
end

function find_normal_vector(v::Vector{T}) where T <: Real
    if length(v) == 1  
        error("There is no normal_vector for a one dimensional vector.")
    elseif length(v) == 2
        normal_vec = vec(LinearAlgebra.nullspace(v'))
    else
        normal_vec = LinearAlgebra.nullspace(v)
    end
end

find_tangential_vector(v::Vector{T}) where T <: Real =  v ./ LinearAlgebra.norm(v)

normal_velocity(velocity::Vector{T}) where T <: Real = projection(velocity, find_normal_vector(velocity)) .* find_normal_vector(velocity)

tangential_velocity(velocity::Vector{T}) where T <: Real = projection(velocity, find_tangential_vector(velocity)) .* find_tangential_vector(velocity)

## 

??(??, mass::T, radius::T; coord::Symbol = :z) where T <: Real = ??(mass, radius, :z) / (mass * radius^2)
??(mass::T, radius::T, coord::Val{:z}) where T <: Real = 0.5 * mass * (radius^2)
??(mass::T, radius::T, coord::Val{:x}) where T <: Real = 0.25 * mass * (radius^2)
??(mass::T, radius::T, coord::Val{:y}) where T <: Real = ??(mass, radius, Val{:x})
??(mass::T, radius::T, coord::Symbol) where T <: Real = ??(mass, radius, Val(coord)) 