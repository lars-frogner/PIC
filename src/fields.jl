VectorField{T,N} = SVector{3,Array{T,N}}
VectorFieldF{N} = VectorField{Float,N}

new_vectorfield(shape) = SVector(ntuple(_ -> zeros(shape), 3))

mutable struct Fields{N}
    grid::Grid{N}
    E::VectorFieldF{N} # Defined in the middle of cell edges at half time steps
    B::VectorFieldF{N} # Defined in the middle of cell faces at whole time steps
    curl_E::VectorFieldF{N} # Defined in the same place as B at whole time steps
    curl_B::VectorFieldF{N} # Defined in the same place as E at half time steps
    j::VectorFieldF{N} # Defined in the same place as E at whole time steps

    function Fields(grid::Grid{N}) where {N}
        E = new_vectorfield(size(grid))
        B = new_vectorfield(size(grid))
        curl_E = new_vectorfield(size(grid))
        curl_B = new_vectorfield(size(grid))
        j = new_vectorfield(size(grid))
        new{N}(grid, E, B, curl_E, curl_B, j)
    end
end

function advance_E!(fields::Fields{N}, Δt) where {N}
    for n = 1:N
        fields.E[n] .+=
            (Δt * SPEED_OF_LIGHT²) .* fields.curl_B[n] .-
            (Δt * VACUUM_PERMITTIVITY⁻¹) .* fields.j[n]
    end
end

function advance_B!(fields::Fields{N}, Δt) where {N}
    for n = 1:N
        fields.B[n] .-= Δt .* fields.curl_E[n]
    end
end

function update_curl_E!(fields::Fields)
    curl!(fields.curl_E, fields.E, fields.grid, Up)
end

function update_curl_B!(fields::Fields)
    curl!(fields.curl_B, fields.B, fields.grid, Down)
end

abstract type FieldType end
struct ElectricField <: FieldType end
struct MagneticField <: FieldType end

struct FieldInterpolator{N}
    interpolator::Interpolator
    grid::Grid{N}
    coordinate_references::NTuple{N,Type{<:ReferenceCoords}}

    function FieldInterpolator(
        fields::Fields{N},
        ::Type{MagneticField},
        ::Dim{D},
    ) where {N,D}
        interpolator = Interpolator(fields.B[D])
        coordinate_references = ntuple(n -> (n == D) ? LowerEdges : Centers, N)
        new{N}(interpolator, fields.grid, coordinate_references)
    end

    function FieldInterpolator(
        fields::Fields{N},
        ::Type{ElectricField},
        ::Dim{D},
    ) where {N,D}
        interpolator = Interpolator(fields.E[D])
        coordinate_references = ntuple(n -> (n == D) ? Centers : LowerEdges, N)
        new{N}(interpolator, fields.grid, coordinate_references)
    end
end

function interpolate(field_interpolator::FieldInterpolator, position)
    continuous_indices = find_indices_at(
        field_interpolator.grid,
        position,
        field_interpolator.coordinate_references,
        index_type = Continuous,
    )
    interpolate(field_interpolator.interpolator, continuous_indices)
end

struct VectorFieldInterpolator{N}
    interpolators::NTuple{3,FieldInterpolator{N}}

    function VectorFieldInterpolator(
        fields::Fields{N},
        field_type::Type{<:FieldType},
    ) where {N}
        interpolators =
            ntuple(n -> FieldInterpolator(fields, field_type, Dim{n}()), 3)
        new{N}(interpolators)
    end
end

function interpolate(
    vectorfield_interpolator::VectorFieldInterpolator,
    position,
)
    [
        interpolate(vectorfield_interpolator.interpolators[n], position)
        for n = 1:3
    ]
end

function curl!(
    dest::VectorField{<:Any,3},
    field::VectorField{<:Any,3},
    grid::Grid{3},
    direction::Type{<:SliceShift},
)
    @dbgasserteq(size(dest[1]), size(grid))
    @dbgassert dest !== field
    for n = 1:3
        dest[n][interior(grid)...] .= curlcomp(field, grid, direction, Dim{n}())
    end
end

function curl!(
    dest::VectorField{<:Any,2},
    field::VectorField{<:Any,2},
    grid::Grid{2},
    direction::Type{<:SliceShift},
)
    @dbgasserteq(size(dest[1]), size(grid))
    @dbgassert dest !== field
    dest[3][interior(grid)...] .= curlcomp(field, grid, direction, DIM3)
end

function curlcomp(
    field::VectorField{<:Any,3},
    grid::Grid{3},
    direction::Type{<:SliceShift},
    ::Dim{1},
)
    shifteddiff(field[3], grid, direction, DIM2) .* cell_extents⁻¹(grid)[2] .-
    shifteddiff(field[2], grid, direction, DIM3) .* cell_extents⁻¹(grid)[3]
end

function curlcomp(
    field::VectorField{<:Any,3},
    grid::Grid{3},
    direction::Type{<:SliceShift},
    ::Dim{2},
)
    shifteddiff(field[1], grid, direction, DIM3) .* cell_extents⁻¹(grid)[3] .-
    shifteddiff(field[3], grid, direction, DIM1) .* cell_extents⁻¹(grid)[1]
end

function curlcomp(
    field::VectorField{<:Any,N},
    grid::Grid{N},
    direction::Type{<:SliceShift},
    ::Dim{3},
) where {N}
    shifteddiff(field[2], grid, direction, DIM1) .* cell_extents⁻¹(grid)[1] .-
    shifteddiff(field[1], grid, direction, DIM2) .* cell_extents⁻¹(grid)[2]
end

function grad!(
    dest::VectorField{<:Any,N},
    field::Array{<:Any,N},
    grid::Grid{N},
    direction::Type{<:SliceShift},
) where {N}
    @dbgasserteq(size(dest[1]), size(grid))
    for n = 1:N
        @dbgassert dest[n] !== field
        dest[n][interior(grid)...] .=
            shifteddiff(field, grid, direction, Dim{n}()) .*
            cell_extents⁻¹(grid)[n]
    end
end

function grad_padded!(
    dest::VectorField{<:Any,N},
    padded_field::Array{<:Any,N},
    padded_grid::Grid{N},
    direction::Type{<:SliceShift},
) where {N}
    @dbgasserteq(size(dest[1]) .+ 2, size(padded_grid))
    for n = 1:N
        dest[n] .=
            shifteddiff(padded_field, padded_grid, direction, Dim{n}()) .*
            cell_extents⁻¹(padded_grid)[n]
    end
end

function div!(
    dest::Array{<:Any,3},
    field::VectorField{<:Any,3},
    grid::Grid{3},
    direction::Type{<:SliceShift},
)
    @dbgasserteq(size(dest), size(grid))
    @dbgassert (dest !== field[1]) && (dest !== field[2]) && (dest !== field[3])
    dest[interior(grid)...] .=
        shifteddiff(field[1], grid, direction, DIM1) .*
        cell_extents⁻¹(grid)[1] .+
        shifteddiff(field[2], grid, direction, DIM2) .*
        cell_extents⁻¹(grid)[2] .+
        shifteddiff(field[3], grid, direction, DIM3) .* cell_extents⁻¹(grid)[3]
end

function div!(
    dest::Array{<:Any,2},
    field::VectorField{<:Any,2},
    grid::Grid{2},
    direction::Type{<:SliceShift},
)
    @dbgasserteq(size(dest), size(grid))
    @dbgassert (dest !== field[1]) && (dest !== field[2])
    dest[interior(grid)...] .=
        shifteddiff(field[1], grid, direction, DIM1) .*
        cell_extents⁻¹(grid)[1] .+
        shifteddiff(field[2], grid, direction, DIM2) .* cell_extents⁻¹(grid)[2]
end

function div(
    field::VectorField{<:Any,N},
    grid::Grid{N},
    direction::Type{<:SliceShift},
) where {N}
    div_field =
        shifteddiff(field[1], grid, direction, DIM1) .* cell_extents⁻¹(grid)[1]
    for n = 2:N
        div_field .+=
            shifteddiff(field[n], grid, direction, Dim{n}()) .*
            cell_extents⁻¹(grid)[n]
    end
    div_field
end
