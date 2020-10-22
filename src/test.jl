using Test

function test_grid()
    shape = (4, 4, 4)
    grid = Grid(shape, (1.0, 1.0, 1.0))
    all(size(grid) .== shape)
end

function test_fields()
    grid = Grid((4, 5, 6), (1.0, 1.0, 1.0))
    fields = Fields(grid)
    update_curl_E!(fields)
    true
end

function test_poisson()
    grid = Grid((30, 30, 3), (10.0, 10.0, 0.1))
    charge_density = zeros(size(grid))
    charge_density[6, 6, 2] = 1.0
    gs_potential = zeros(size(grid))
    solve_poisson!(
        gs_potential,
        charge_density,
        grid,
        Free,
        GaussSeidel(max_deviation = 1e-6, relaxation_factor = 1.6),
    )
    correct_potential = zeros(size(grid))
    solve_poisson!(correct_potential, charge_density, grid, Free, Analytical)
    maximum(abs.(gs_potential - correct_potential)[isfinite.(
        correct_potential,
    )]) / compute_solution_scale(sum(abs.(charge_density)), grid, DIM3) < 0.06
end

function test_particles()
    positions = zeros(2, 1)
    positions[:, 1] = [0.0, 0.0]
    grid = Grid((30, 30), (10.0, 10.0))
    charge_density = zeros(size(grid))
    potential = zeros(size(grid))

    n_steps = 20
    total_charge = zeros(n_steps)

    for i = 1:n_steps
        particles = Particles{2}(positions, 1.0, 1.0)
        charge_density .= 0
        sample_charge_density!(charge_density, grid, particles)
        total_charge[i] = sum(charge_density) * cell_volume(grid)
        positions .+= [0.1, 0.07]
    end

    maximum(abs.(total_charge .- 1.0)) < 1e-12
end

function test_interpolation()
    grid = Grid((4, 5, 6), (1.0, 1.0, 1.0))
    fields = Fields(grid)
    E_interp = FieldInterpolator(fields, ElectricField, DIM1)
    interpolate(E_interp, [0.3, 0.81, 0.02])
    true
end

function test_timestepping()
    grid = Grid((50, 50), (1.0, 1.0))
    positions = zeros(2, 1)
    positions[:, 1] = [0.2, 0.5]
    #positions[:, 2] = [0.5, 0.5]
    particles = Particles{2}(positions, 1.0, 1.0)
    fields = initialize_fields(particles, grid)
    savefig(
        scalarheatmap(sqrt.(
            fields.E[1] .^ 2 .+ fields.E[2] .^ 2 .+ fields.E[3] .^ 2,
        )),
        "Enorm.png",
    )
    n_steps = 100
    n_skips = 1
    anim = @animate for _ = 1:n_steps
        vectorfieldquiver(centers(grid)[1], centers(grid)[2], fields.E)
        advance!(particles, fields, 1e-2)
        @show(maximum(fields.E[1]))
    end every n_skips
    gif(anim, "E.gif", fps = 15)

    true
end

export test_timestepping

function runtests()

    @testset "grid.jl" begin
        @test test_grid()
    end

    @testset "fields.jl" begin
        @test test_fields()
    end

    @testset "poisson.jl" begin
        @test test_poisson()
    end

    @testset "particles.jl" begin
        @test test_particles()
    end

    @testset "interp.jl" begin
        @test test_interpolation()
    end

    @testset "timestepping.jl" begin
        @test test_timestepping()
    end
end
export runtests
