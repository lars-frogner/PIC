struct Dim{N}
    Dim{1}() = new()
    Dim{2}() = new()
    Dim{3}() = new()
end

const DIM1 = Dim{1}()
const DIM2 = Dim{2}()
const DIM3 = Dim{3}()

otherdims(::Dim{1}) = [2, 3]
otherdims(::Dim{2}) = [3, 1]
otherdims(::Dim{3}) = [1, 2]

abstract type ReferenceCoords end
struct LowerEdges <: ReferenceCoords end
struct Centers <: ReferenceCoords end

abstract type SliceShift end
struct Up <: SliceShift end
struct Down <: SliceShift end
struct UpDown <: SliceShift end

abstract type IndexType end
struct Discrete <: IndexType end
struct Continuous <: IndexType end

struct Grid{N}
    shape::SVectorU{N}
    extents::SVectorF{N}
    origins::SVectorF{N}
    cell_extents::SVectorF{N}
    cell_extents⁻¹::SVectorF{N}
    cell_extents⁻²::SVectorF{N}
    centers::SVector{N,Vector{Float}}
    lower_edges::SVector{N,Vector{Float}}
    interior::SVector{N,UnitRange}
    interior_up::SVector{N,SVector{N,UnitRange}}
    interior_dn::SVector{N,SVector{N,UnitRange}}

    function Grid(
        shape::Union{NTuple{N},SVectorU{N}},
        extents::Union{NTuple{N},SVectorF{N}},
        origins::Union{NTuple{N},SVectorF{N}},
    ) where {N}
        @check all(shape .> 2)
        @check all(extents .> 0.0)
        cell_extents = extents ./ shape
        cell_extents⁻¹ = 1 ./ cell_extents
        cell_extents⁻² = cell_extents⁻¹ .* cell_extents⁻¹
        centers = ntuple(
            n -> range(
                origins[n] + 0.5 * cell_extents[n],
                origins[n] + extents[n] - 0.5 * cell_extents[n],
                length = shape[n],
            ),
            N,
        )
        lower_edges = ntuple(
            n -> range(
                origins[n],
                origins[n] + extents[n] - cell_extents[n],
                length = shape[n],
            ),
            N,
        )
        interior = ntuple(n -> 2:shape[n]-1, N)
        interior_up = ntuple(
            n -> SVector(ntuple(
                m -> ((m == n) ? (3:shape[m]) : (2:shape[m]-1)),
                N,
            )),
            N,
        )
        interior_dn = ntuple(
            n -> SVector(ntuple(
                m -> ((m == n) ? (1:shape[m]-2) : (2:shape[m]-1)),
                N,
            )),
            N,
        )
        new{N}(
            shape,
            extents,
            origins,
            cell_extents,
            cell_extents⁻¹,
            cell_extents⁻²,
            centers,
            lower_edges,
            interior,
            interior_up,
            interior_dn,
        )
    end

    function Grid(
        shape::Union{NTuple{N},SVectorU{N}},
        extents::Union{NTuple{N},SVectorF{N}},
    ) where {N}
        Grid(shape, extents, ntuple(_ -> 0.0, N))
    end
end

Base.size(grid::Grid) = grid.shape.data

extents(grid::Grid) = grid.extents

origins(grid::Grid) = grid.origins

interior(grid::Grid) = grid.interior

interior_dn(grid::Grid) = grid.interior_dn

interior_up(grid::Grid) = grid.interior_up

lower_edges(grid::Grid) = grid.lower_edges

centers(grid::Grid) = grid.centers

cell_extents(grid::Grid) = grid.cell_extents

cell_extents⁻¹(grid::Grid) = grid.cell_extents⁻¹

cell_extents⁻²(grid::Grid) = grid.cell_extents⁻²

cell_volume(grid::Grid) = prod(cell_extents(grid))

cell_volume⁻¹(grid::Grid) = 1 / cell_volume(grid)

coordinates(grid::Grid, ::Type{LowerEdges}) = lower_edges(grid)

coordinates(grid::Grid, ::Type{Centers}) = centers(grid)

function padded(grid::Grid; padding = 1)
    new_shape = size(grid) .+ 2 .* padding
    new_extents = extents(grid) .+ 2 .* padding .* cell_extents(grid)
    new_origins = origins(grid) .- padding .* cell_extents(grid)
    Grid(new_shape, new_extents, new_origins)
end

function find_indices_at(
    grid::Grid{N},
    position,
    reference::Type{<:ReferenceCoords};
    index_type::Type{<:IndexType} = Discrete,
) where {N}
    [
        find_index_at(
            grid,
            position[n],
            Dim{n}(),
            reference,
            index_type = index_type,
        ) for n = 1:N
    ]
end

function find_indices_at(
    grid::Grid{N},
    position,
    references;
    index_type::Type{<:IndexType} = Discrete,
) where {N}
    [
        find_index_at(
            grid,
            position[n],
            Dim{n}(),
            references[n],
            index_type = index_type,
        ) for n = 1:N
    ]
end

function find_index_at(
    grid::Grid,
    coord,
    ::Dim{D},
    ::Type{LowerEdges};
    index_type::Type{<:IndexType} = Discrete,
) where {D}
    compute_index_unform(
        size(grid)[D],
        lower_edges(grid)[D][1],
        extents(grid)[D],
        coord,
        index_type,
    )
end

function find_index_at(
    grid::Grid,
    coord,
    ::Dim{D},
    ::Type{Centers};
    index_type::Type{<:IndexType} = Discrete,
) where {D}
    compute_index_unform(
        size(grid)[D],
        centers(grid)[D][1],
        extents(grid)[D],
        coord,
        index_type,
    )
end

function compute_index_unform(
    length,
    lower_edge,
    extent,
    coord,
    ::Type{Discrete},
)
    floor(
        Int,
        compute_index_unform(length, lower_edge, extent, coord, Continuous),
    )
end

function compute_index_unform(
    length,
    lower_edge,
    extent,
    coord,
    ::Type{Continuous},
)
    1 + (coord - lower_edge) * (length - 1) / extent
end

function shifteddiff(
    field::Array{<:Any,N},
    grid::Grid{N},
    ::Type{Up},
    ::Dim{D},
) where {N,D}
    @dbgasserteq(size(field), size(grid))
    field[interior_up(grid)[D]...] .- field[interior(grid)...]
end

function shifteddiff(
    field::Array{<:Any,N},
    grid::Grid{N},
    ::Type{Down},
    ::Dim{D},
) where {N,D}
    @dbgasserteq(size(field), size(grid))
    field[interior(grid)...] .- field[interior_dn(grid)[D]...]
end

function shiftedsum(
    field::Array{<:Any,N},
    grid::Grid{N},
    ::Type{UpDown},
    ::Dim{D},
) where {N,D}
    @dbgasserteq(size(field), size(grid))
    field[interior_up(grid)[D]...] .+ field[interior_dn(grid)[D]...]
end
