function initialize_fields(particles::Particles{N}, grid::Grid{N}) where {N}
    padded_grid = padded(grid)

    potential_source = zeros(size(padded_grid))
    sample_charge_density!(potential_source, padded_grid, particles)
    potential_source .*= VACUUM_PERMITTIVITY⁻¹

    potential = zeros(size(padded_grid))
    solve_poisson!(
        potential,
        potential_source,
        padded_grid,
        Free,
        GaussSeidel(),
    )
    potential .*= -1

    fields = Fields(grid)
    grad_padded!(fields.E, potential, padded_grid, Up)
    fields
end

function advance!(particles::Particles{N}, fields::Fields{N}, Δt) where {N}
    update_accelerations!(particles, fields)
    advance_positions!(particles, 0.5 * Δt)
    advance_velocities!(particles, 0.5 * Δt)
    update_j!(fields, particles)
    update_curl_E!(fields)
    advance_B!(fields, 0.5 * Δt)
    update_curl_B!(fields)
    advance_E!(fields, Δt)
    advance_B!(fields, 0.5 * Δt)
    advance_positions!(particles, 0.5 * Δt)
    advance_velocities!(particles, 0.5 * Δt)
end
