module PIC
using StaticArrays
using Interpolations
using Plots
import Base

gr()

include("error.jl")
include("debug.jl")
include("num.jl")
include("constants.jl")
include("interp.jl")
include("grid.jl")
include("fields.jl")
include("particles.jl")
include("poisson.jl")
include("timestepping.jl")
include("plotting.jl")

include("test.jl")

#
# positions = zeros(2, 1)
# positions[:, 1] = [0.0, 0.0]
# particles = Particles{2}(positions, 1.0)
# particles.velocities[:] = [0.2, 0.07]
# grid = Grid((30, 30), (10.0, 10.0))
# charge_density = zeros(size(grid))
# potential = zeros(size(grid))
#
# n_steps = 100
# n_skips = 1
# anim = @animate for _ = 1:n_steps
#     charge_density .= 0
#     #sample_charge_density!(charge_density, grid, particles)
#     sample_current_density!(charge_density, grid, particles, DIM1)
#     positions .+= particles.velocities
#     @show sum(charge_density) * cell_volume(grid)
#     # potential .= 0
#     # solve_poisson!(
#     #     potential,
#     #     charge_density,
#     #     grid,
#     #     Free,
#     #     Jacobi(max_deviation = 1e-6),
#     # )
#     scalarheatmap(charge_density)#, clims = (0, 8))#(-0.3, 0.3))
# end every n_skips
# gif(anim, "test.gif", fps = 15)

#
# grid = Grid((30, 30, 3), (10.0, 10.0, 0.1))
# charge_density = zeros(size(grid))
# positions = []
# pt = [60, 60, 2]
# charge_density[pt...] = 1.0

# potential = zeros(size(grid))
# solve_poisson!(
#     potential,
#     charge_density,
#     grid,
#     Free,
#     Jacobi(max_deviation = 1e-6),
# )
# p1 = scalarheatmap(potential, clims = (-1e-3, 0))
# savefig(p1, "jac.png")
#
# p6 = plot(lower_edges(grid)[1], potential[1:end, pt[2]], ylims = (-1e-3, 0))
# savefig(p6, "jac_1d.png")


# p2 = scalarheatmap(potential[:, :, 2], clims = (-1e-3, 0))
# savefig(p2, "gs.png")
#
# p5 = plot(
#     lower_edges(grid)[1],
#     potential[1:end, pt[2:end]...],
#     ylims = (-1e-3, 0),
# )
# savefig(p5, "gs_1d.png")
#
# potential = zeros(size(grid))
# source_indices = findall(!iszero, charge_density)
# for k = 1:size(grid)[3]
#     for j = 1:size(grid)[2]
#         for i = 1:size(grid)[1]
#             potential[i, j, k] = compute_free_potential(
#                 SVector(
#                     lower_edges(grid)[1][i],
#                     lower_edges(grid)[2][j],
#                     lower_edges(grid)[3][k],
#                 ),
#                 charge_density,
#                 cell_volume(grid),
#                 lower_edges(grid),
#                 source_indices,
#             )
#         end
#     end
# end
# p2 = scalarheatmap(potential[:, :, 2], clims = (-1e-3, 0))
# savefig(p2, "correct.png")
#
# p3 = plot(
#     lower_edges(grid)[1],
#     potential[1:end, pt[2:end]...],
#     ylims = (-1e-3, 0),
# )
# savefig(p3, "correct_1d.png")
#
# analytical(charge, ::Dim{3}) =
#     -charge ./ (4π .* abs.(lower_edges(grid)[1] .- lower_edges(grid)[1][pt[1]]))
#
# analytical(charge, ::Dim{2}) =
#     charge .* log.(abs.(lower_edges(grid)[1] .- lower_edges(grid)[1][pt[1]])) ./
#     2π
#
# p4 = plot(
#     lower_edges(grid)[1],
#     analytical(
#         charge_density[pt...] * cell_volume(grid),
#         Dim{length(size(grid))}(),
#     ),
#     ylims = (-1e-3, 0),
# )
# savefig(p4, "analytical_1d.png")

# potential = zeros(size(grid)[1])
# analytical1 = analytical(charge_density[pt...] * cell_extents(grid)[1])
# potential[1] = analytical1[1]
# potential[end] = analytical1[end]
# for n = 1:10000
#     for i = 2:length(potential)-1
#         potential[i] =
#             0.5 * (
#                 potential[i-1] + potential[i+1] -
#                 cell_extents(grid)[1]^2 * charge_density[i, pt[2]]
#             )
#     end
# end
# p7 = plot(lower_edges(grid)[1], potential)
# savefig(p7, "gs1_1d.png")
#
# p8 = plot(lower_edges(grid)[1], analytical1)
# savefig(p8, "analytical1_1d.png")

# n_steps = 20
# n_skips = 1
# set_poisson_boundaries!(potential, charge_density, grid, Free)
# anim = @animate for _ = 1:n_steps
#     @show update_poisson!(
#         potential,
#         charge_density,
#         grid,
#         Free,
#         GaussSeidel(relaxation_factor = 1.5),
#     ) / (prod(extents(grid)) * sum(abs.(charge_density)))
#     scalarheatmap(potential)
# end every n_skips
# gif(anim, "test.gif", fps = 15)

end
