abstract type FieldProperty end
struct ChargeDensity <: FieldProperty end
struct CurrentDensity{D} <: FieldProperty end

mutable struct Particles{N}
    positions::Array{Float,2}
    velocities::Array{Float,2}
    accelerations::Array{Float,2}
    mass::Float
    charge::Float
    dim::Dim{N}

    function Particles{N}(positions::Array{<:Any,2}, mass, charge) where {N}
        dims, count = size(positions)
        @checkeq(dims, N)
        velocities = zeros(N, count)
        accelerations = zeros(N, count)
        new(positions, velocities, accelerations, mass, charge, Dim{N}())
    end
end

count(particles::Particles) = size(particles.positions)[2]
mass(particles::Particles) = particles.mass
charge(particles::Particles) = particles.charge

function advance_positions!(particles::Particles, Δt)
    particles.positions .+= particles.velocities .* Δt
end

function advance_velocities!(particles::Particles, Δt)
    particles.velocities .+= particles.accelerations .* Δt
end

function update_accelerations!(
    particles::Particles{N},
    fields::Fields{N},
) where {N}
    E_interp = VectorFieldInterpolator(fields, ElectricField)
    B_interp = VectorFieldInterpolator(fields, MagneticField)

    q_over_m = charge(particles) / mass(particles)
    for p = 1:count(particles)
        r = particles.positions[:, p]
        v = particles.velocities[:, p]
        E = interpolate(E_interp, r)
        B = interpolate(B_interp, r)
        particles.accelerations[:, p] .=
            lorentz_force(q_over_m, v, E, B, Dim{N}())
    end
end

function lorentz_force(q_over_m, v, E, B, dim::Dim{N}) where {N}
    q_over_m .* (0 .* E[1:N] .+ v_cross_B(v, B, dim) .+ [1.0, 0])
end

v_cross_B(v, B, ::Dim{2}) = [v[2] * B[3], -v[1] * B[3]]

v_cross_B(v, B, ::Dim{3}) = cross(v, B)

function update_j!(fields::Fields{N}, particles::Particles{N}) where {N}
    for n = 1:N
        sample_current_density!(fields.j[n], fields.grid, particles, Dim{n}())
    end
end

function sample_charge_density!(
    field::Array{Float,N},
    grid::Grid{N},
    particles::Particles{N},
) where {N}
    sample_property!(field, grid, particles, fill(LowerEdges, N), ChargeDensity)
end

function sample_current_density!(
    field::Array{Float,N},
    grid::Grid{N},
    particles::Particles{N},
    ::Dim{D},
) where {N,D}
    coordinate_references = fill(LowerEdges, N)
    coordinate_references[D] = Centers
    sample_property!(
        field,
        grid,
        particles,
        coordinate_references,
        CurrentDensity{D},
    )
end

function sample_property!(
    field::Array{Float,N},
    grid::Grid{N},
    particles::Particles{N},
    coordinate_references,
    property_type::Type{<:FieldProperty},
) where {N}
    @dbgasserteq(size(field), size(grid))

    field .= 0

    shape = size(field)
    weights = zeros(2, N)

    for p = 1:count(particles)
        position = particles.positions[:, p]
        indices = find_indices_at(grid, position, coordinate_references)
        if any((indices .< 1) .| (indices .>= shape))
            continue
        end

        weights .= 0
        for n = 1:N
            coords = coordinates(grid, coordinate_references[n])
            weights[2, n] =
                (position[n] - coords[n][indices[n]]) * cell_extents⁻¹(grid)[n]
            weights[1, n] = 1 - weights[2, n]
        end

        distribute_property!(
            field,
            compute_property(particles, grid, property_type, p),
            indices,
            weights,
        )
    end
end

function compute_property(
    particles::Particles{N},
    grid::Grid{N},
    ::Type{ChargeDensity},
    ::Integer,
) where {N}
    charge(particles) * cell_volume⁻¹(grid)
end

function compute_property(
    particles::Particles{N},
    grid::Grid{N},
    ::Type{CurrentDensity{D}},
    index::Integer,
) where {N,D}
    charge(particles) * particles.velocities[D, index] * cell_volume⁻¹(grid)
end

function distribute_property!(field::Array{Float,2}, property, indices, weights)
    for i = 1:2, j = 1:2
        field[indices[1]+i-1, indices[2]+j-1] +=
            property * weights[i, 1] * weights[j, 2]
    end
end

function distribute_property!(field::Array{Float,3}, property, indices, weights)
    for i = 1:2, j = 1:2, k = 1:2
        field[indices[1]+i-1, indices[2]+j-1, indices[3]+k-1] +=
            property * weights[i, 1] * weights[j, 2] * weights[k, 3]
    end
end
