abstract type BoundaryCondition end
struct ZeroGradient <: BoundaryCondition end
struct ZeroLaplacian <: BoundaryCondition end
struct Free <: BoundaryCondition end

abstract type PoissonSolver end

struct Jacobi <: PoissonSolver
    max_deviation::Float
    max_iterations::UInt64
    Jacobi(; max_deviation = 1e-6, max_iterations = 100000) =
        @check(max_deviation > 0.0) && new(max_deviation, max_iterations)
end

struct GaussSeidel <: PoissonSolver
    relaxation_factor::Float
    max_deviation::Float
    max_iterations::UInt64
    GaussSeidel(;
        relaxation_factor = 1.0,
        max_deviation = 1e-6,
        max_iterations = 100000,
    ) =
        @check(relaxation_factor >= 1.0 && relaxation_factor < 2.0) &&
        @check(max_deviation > 0.0) &&
        new(relaxation_factor, max_deviation, max_iterations)
end

struct Analytical <: PoissonSolver end

function solve_poisson!(
    field::Array{<:Any,N},
    source::Array{<:Any,N},
    grid::Grid{N},
    boundary_condition::Type{<:BoundaryCondition},
    solver::Jacobi,
) where {N}
    old_field = field
    new_field = similar(field)

    total_source = sum(abs.(source))
    if total_source == 0
        field .= 0
        return
    end

    set_poisson_boundaries!(old_field, source, grid, boundary_condition)

    n_iter = 0
    deviation_scale = 1 / compute_solution_scale(total_source, grid, Dim{N}())
    while true
        max_difference = update_poisson!(
            new_field,
            old_field,
            source,
            grid,
            boundary_condition,
            solver,
        )
        n_iter += 1

        if max_difference * deviation_scale < solver.max_deviation
            field[:] = new_field
            break
        elseif n_iter >= solver.max_iterations
            println(stderr, "Warning: Jacobi Poisson solver did not converge")
            field[:] = new_field
            break
        end

        tmp = old_field
        old_field = new_field
        new_field = tmp
    end
end

function solve_poisson!(
    field::Array{<:Any,N},
    source::Array{<:Any,N},
    grid::Grid{N},
    boundary_condition::Type{<:BoundaryCondition},
    solver::GaussSeidel,
) where {N}
    total_source = sum(abs.(source))
    if total_source == 0
        field .= 0
        return
    end

    set_poisson_boundaries!(field, source, grid, boundary_condition)

    n_iter = 0
    deviation_scale = 1 / compute_solution_scale(total_source, grid, Dim{N}())
    while true
        max_difference =
            update_poisson!(field, source, grid, boundary_condition, solver)
        n_iter += 1
        if max_difference * deviation_scale < solver.max_deviation
            break
        elseif n_iter >= solver.max_iterations
            println(
                stderr,
                "Warning: Gauss Seidel Poisson solver did not converge",
            )
            break
        end
    end
end

function solve_poisson!(
    field::Array{<:Any,3},
    source::Array{<:Any,3},
    grid::Grid{3},
    boundary_condition::Type{Free},
    solver::Type{Analytical},
)
    source_indices = findall(!iszero, source)
    if length(source_indices) == 0
        field .= 0
        return
    end
    shape = size(grid)
    coords = lower_edges(grid)
    grid_cell_volume = cell_volume(grid)
    for k = 1:shape[3]
        for j = 1:shape[2]
            for i = 1:shape[1]
                field[i, j, k] = compute_free_potential(
                    SVector(coords[1][i], coords[2][j], coords[3][k]),
                    source,
                    grid_cell_volume,
                    coords,
                    source_indices,
                )
            end
        end
    end
end

function solve_poisson!(
    field::Array{<:Any,2},
    source::Array{<:Any,2},
    grid::Grid{2},
    boundary_condition::Type{Free},
    solver::Type{Analytical},
)
    source_indices = findall(!iszero, source)
    shape = size(grid)
    coords = lower_edges(grid)
    grid_cell_volume = cell_volume(grid)
    for j = 1:shape[2]
        for i = 1:shape[1]
            field[i, j] = compute_free_potential(
                SVector(coords[1][i], coords[2][j]),
                source,
                grid_cell_volume,
                coords,
                source_indices,
            )
        end
    end
end

function compute_solution_scale(total_source, grid::Grid{3}, ::Dim{3})
    total_source * cell_volume(grid) * length(extents(grid)) /
    sum(extents(grid))
end

function compute_solution_scale(total_source, grid::Grid{2}, ::Dim{2})
    total_source *
    cell_volume(grid) *
    log(ℯ - 1 + sum(extents(grid)) / length(extents(grid)))
end

function update_poisson!(
    dest::Array{<:Any,N},
    field::Array{<:Any,N},
    source::Array{<:Any,N},
    grid::Grid{N},
    boundary_condition::Type{<:BoundaryCondition},
    solver::Jacobi,
) where {N}
    @dbgasserteq(size(dest), size(grid))
    @dbgasserteq(size(source), size(grid))
    @dbgassert (dest !== field) && (dest !== source)

    scale = 0.5 / sum(cell_extents⁻²(grid))
    coefs = cell_extents⁻²(grid) .* scale

    dest[interior(grid)...] .=
        poisson_jacobi_update_field(field, grid, coefs) .-
        source[interior(grid)...] .* scale

    update_poisson_boundaries!(dest, field, boundary_condition)

    max_difference = maximum(abs.(dest - field))
    max_difference
end

function poisson_jacobi_update_field(
    field::Array{<:Any,3},
    grid::Grid{3},
    coefs,
)
    shiftedsum(field, grid, UpDown, DIM1) .* coefs[1] .+
    shiftedsum(field, grid, UpDown, DIM2) .* coefs[2] .+
    shiftedsum(field, grid, UpDown, DIM3) .* coefs[3]
end

function poisson_jacobi_update_field(
    field::Array{<:Any,2},
    grid::Grid{2},
    coefs,
)
    shiftedsum(field, grid, UpDown, DIM1) .* coefs[1] .+
    shiftedsum(field, grid, UpDown, DIM2) .* coefs[2]
end

function update_poisson!(
    field::Array{<:Any,N},
    source::Array{<:Any,N},
    grid::Grid{N},
    boundary_condition::Type{<:BoundaryCondition},
    solver::GaussSeidel,
) where {N}
    @dbgasserteq(size(field), size(grid))
    @dbgasserteq(size(source), size(grid))
    @dbgassert field !== source

    source_scale = 0.5 * solver.relaxation_factor / sum(cell_extents⁻²(grid))
    coefs = cell_extents⁻²(grid) .* source_scale
    SOR_scale = 1.0 - solver.relaxation_factor

    max_difference = poisson_gauss_seidel_update_field!(
        field,
        source,
        coefs,
        source_scale,
        SOR_scale,
    )
    update_poisson_boundaries!(field, boundary_condition)
    max_difference
end

function poisson_gauss_seidel_update_field!(
    field::Array{<:Any,3},
    source::Array{<:Any,3},
    coefs,
    source_scale,
    SOR_scale,
)
    shape = size(field)
    max_difference = 0.0
    for k = 2:shape[3]-1
        for j = 2:shape[2]-1
            for i = 2:shape[1]-1
                new_value =
                    (field[i-1, j, k] + field[i+1, j, k]) * coefs[1] +
                    (field[i, j-1, k] + field[i, j+1, k]) * coefs[2] +
                    (field[i, j, k-1] + field[i, j, k+1]) * coefs[3] -
                    source[i, j, k] * source_scale + field[i, j, k] * SOR_scale

                max_difference =
                    max(max_difference, abs(new_value - field[i, j, k]))

                field[i, j, k] = new_value
            end
        end
    end
    max_difference
end

function poisson_gauss_seidel_update_field!(
    field::Array{<:Any,2},
    source::Array{<:Any,2},
    coefs,
    source_scale,
    SOR_scale,
)
    shape = size(field)
    max_difference = 0.0
    for j = 2:shape[2]-1
        for i = 2:shape[1]-1
            new_value =
                (field[i-1, j] + field[i+1, j]) * coefs[1] +
                (field[i, j-1] + field[i, j+1]) * coefs[2] -
                source[i, j] * source_scale + field[i, j] * SOR_scale

            max_difference = max(max_difference, abs(new_value - field[i, j]))

            field[i, j] = new_value
        end
    end
    max_difference
end

function set_poisson_boundaries!(
    field::Array{<:Any,3},
    source::Array{<:Any,3},
    grid::Grid{3},
    ::Type{Free},
)
    source_indices = findall(!iszero, source)
    shape = size(grid)
    coords = lower_edges(grid)
    grid_cell_volume = cell_volume(grid)
    for k in [1, shape[3]]
        for j = 1:shape[2]
            for i = 1:shape[1]
                field[i, j, k] = compute_free_potential(
                    SVector(coords[1][i], coords[2][j], coords[3][k]),
                    source,
                    grid_cell_volume,
                    coords,
                    source_indices,
                )
            end
        end
    end
    for k = 2:shape[3]-1
        for j in [1, shape[2]]
            for i = 1:shape[1]
                field[i, j, k] = compute_free_potential(
                    SVector(coords[1][i], coords[2][j], coords[3][k]),
                    source,
                    grid_cell_volume,
                    coords,
                    source_indices,
                )
            end
        end
    end
    for k = 2:shape[3]-1
        for j = 2:shape[2]-1
            for i in [1, shape[1]]
                field[i, j, k] = compute_free_potential(
                    SVector(coords[1][i], coords[2][j], coords[3][k]),
                    source,
                    grid_cell_volume,
                    coords,
                    source_indices,
                )
            end
        end
    end
end

function set_poisson_boundaries!(
    field::Array{<:Any,2},
    source::Array{<:Any,2},
    grid::Grid{2},
    ::Type{Free},
)
    source_indices = findall(!iszero, source)
    shape = size(grid)
    coords = lower_edges(grid)
    grid_cell_area = cell_volume(grid)
    for j in [1, shape[2]]
        for i = 1:shape[1]
            field[i, j] = compute_free_potential(
                SVector(coords[1][i], coords[2][j]),
                source,
                grid_cell_area,
                coords,
                source_indices,
            )
        end
    end
    for j = 2:shape[2]-1
        for i in [1, shape[1]]
            field[i, j] = compute_free_potential(
                SVector(coords[1][i], coords[2][j]),
                source,
                grid_cell_area,
                coords,
                source_indices,
            )
        end
    end
end

function compute_free_potential(
    potential_position::SVector{N,<:Any},
    source::Array{<:Any,N},
    grid_cell_volume,
    coordinates::SVector{N,Vector{T}},
    source_indices::Vector{CartesianIndex{N}},
) where {N,T<:Any}
    potential = 0.0
    for source_index in source_indices
        distance² = 0.0
        for n = 1:N
            distance² +=
                (potential_position[n] - coordinates[n][source_index[n]])^2
        end
        potential += compute_free_potential_contribution(
            grid_cell_volume * source[source_index],
            sqrt(distance²),
            Dim{N}(),
        )
    end
    potential
end

compute_free_potential_contribution(charge, distance, ::Dim{3}) =
    -charge / (4π * distance)

compute_free_potential_contribution(charge, distance, ::Dim{2}) =
    charge * log(distance) / 2π

function update_poisson_boundaries!(::Array{<:Any,N}, ::Type{Free}) where {N} end

function update_poisson_boundaries!(
    dest::Array{<:Any,3},
    field::Array{<:Any,3},
    ::Type{Free},
)
    dest[1:end, 1:end, 1] .= field[1:end, 1:end, 1]
    dest[1:end, 1:end, end] .= field[1:end, 1:end, end]
    dest[1:end, 1, 2:end] .= field[1:end, 1, 2:end]
    dest[1:end, end, 2:end] .= field[1:end, end, 2:end]
    dest[1, 2:end, 2:end] .= field[1, 2:end, 2:end]
    dest[end, 2:end, 2:end] .= field[end, 2:end, 2:end]
end

function update_poisson_boundaries!(
    dest::Array{<:Any,2},
    field::Array{<:Any,2},
    ::Type{Free},
)
    dest[1:end, 1] .= field[1:end, 1]
    dest[1:end, end] .= field[1:end, end]
    dest[1, 2:end] .= field[1, 2:end]
    dest[end, 2:end] .= field[end, 2:end]
end

set_poisson_boundaries!(
    field::Array{<:Any,N},
    ::Array{<:Any,N},
    ::Grid{N},
    boundary_condition::Type{ZeroGradient},
) where {N} = update_poisson_boundaries!(field, boundary_condition)

set_poisson_boundaries!(
    field::Array{<:Any,N},
    ::Array{<:Any,N},
    ::Grid{N},
    boundary_condition::Type{ZeroLaplacian},
) where {N} = update_poisson_boundaries!(field, boundary_condition)

function update_poisson_boundaries!(
    field::Array{<:Any,N},
    boundary_condition::Type{<:BoundaryCondition},
) where {N}
    update_poisson_boundaries!(field, field, boundary_condition)
end

function update_poisson_boundaries!(
    dest::Array{<:Any,3},
    field::Array{<:Any,3},
    ::Type{ZeroGradient},
)
    dest[2:end-1, 2:end-1, 1] .=
        ghost_cell.(field[2:end-1, 2:end-1, 2], ZeroGradient)
    dest[2:end-1, 2:end-1, end] .=
        ghost_cell.(field[2:end-1, 2:end-1, end-1], ZeroGradient)

    dest[2:end-1, 1, 2:end-1] .=
        ghost_cell.(field[2:end-1, 2, 2:end-1], ZeroGradient)
    dest[2:end-1, end, 2:end-1] .=
        ghost_cell.(field[2:end-1, end-1, 2:end-1], ZeroGradient)

    dest[1, 2:end-1, 2:end-1] .=
        ghost_cell.(field[2, 2:end-1, 2:end-1], ZeroGradient)
    dest[end, 2:end-1, 2:end-1] .=
        ghost_cell.(field[end-1, 2:end-1, 2:end-1], ZeroGradient)
end

function update_poisson_boundaries!(
    dest::Array{<:Any,2},
    field::Array{<:Any,2},
    ::Type{ZeroGradient},
)
    dest[2:end-1, 1] .= ghost_cell.(field[2:end-1, 2], ZeroGradient)
    dest[2:end-1, end] .= ghost_cell.(field[2:end-1, end-1], ZeroGradient)

    dest[1, 2:end-1] .= ghost_cell.(field[2, 2:end-1], ZeroGradient)
    dest[end, 2:end-1] .= ghost_cell.(field[end-1, 2:end-1], ZeroGradient)
end

function update_poisson_boundaries!(
    dest::Array{<:Any,3},
    field::Array{<:Any,3},
    ::Type{ZeroLaplacian},
)
    dest[2:end-1, 2:end-1, 1] .=
        ghost_cell.(
            field[2:end-1, 2:end-1, 2],
            field[2:end-1, 2:end-1, 3],
            ZeroLaplacian,
        )
    dest[2:end-1, 2:end-1, end] .=
        ghost_cell.(
            field[2:end-1, 2:end-1, end-1],
            field[2:end-1, 2:end-1, end-2],
            ZeroLaplacian,
        )

    dest[2:end-1, 1, 2:end-1] .=
        ghost_cell.(
            field[2:end-1, 2, 2:end-1],
            field[2:end-1, 3, 2:end-1],
            ZeroLaplacian,
        )
    dest[2:end-1, end, 2:end-1] .=
        ghost_cell.(
            field[2:end-1, end-1, 2:end-1],
            field[2:end-1, end-2, 2:end-1],
            ZeroLaplacian,
        )

    dest[1, 2:end-1, 2:end-1] .=
        ghost_cell.(
            field[2, 2:end-1, 2:end-1],
            field[3, 2:end-1, 2:end-1],
            ZeroLaplacian,
        )
    dest[end, 2:end-1, 2:end-1] .=
        ghost_cell.(
            field[end-1, 2:end-1, 2:end-1],
            field[end-2, 2:end-1, 2:end-1],
            ZeroLaplacian,
        )
end

function update_poisson_boundaries!(
    dest::Array{<:Any,2},
    field::Array{<:Any,2},
    ::Type{ZeroLaplacian},
)
    dest[2:end-1, 1] .=
        ghost_cell.(field[2:end-1, 2], field[2:end-1, 3], ZeroLaplacian)
    dest[2:end-1, end] .=
        ghost_cell.(field[2:end-1, end-1], field[2:end-1, end-2], ZeroLaplacian)

    dest[1, 2:end-1] .=
        ghost_cell.(field[2, 2:end-1], field[3, 2:end-1], ZeroLaplacian)
    dest[end, 2:end-1] .=
        ghost_cell.(field[end-1, 2:end-1], field[end-2, 2:end-1], ZeroLaplacian)
end

function ghost_cell(u₂, ::Type{ZeroGradient})
    u₂
end

function ghost_cell(u₂, u₃, ::Type{ZeroGradient})
    (4 * u₂ - u₃) / 3
end

function ghost_cell(u₂, u₃, ::Type{ZeroLaplacian})
    2 * u₂ - u₃
end

function ghost_cell(u₂, u₃, u₄, ::Type{ZeroLaplacian})
    0.5 * (5 * u₂ - 4 * u₃ + u₄)
end
