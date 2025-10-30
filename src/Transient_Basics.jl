"""
    get_internal_state_vec_pbc_profile(lattice, elongation_vector, k₋, k₊, mobile, jammed, paused)

Initialize the internal state vector for the lattice assuming all particles are active (non-paused).

- Sets internal state of every occupied site (ribosome) to the active rate `k₊`.
- Returns the `internal_state_vec` array where each occupied lattice site has value `k₊`, else zero.
"""
function get_internal_state_vec_pbc_profile(lattice, elongation_vector, k₋, k₊, mobile, jammed, paused)
    internal_state_vec = zeros(length(lattice)) # Create vector for internal states

    for i in eachindex(lattice)
        if lattice[i] == 1
            internal_state_vec[i] = k₊  # All particles are active, no pausing
        end
    end

    return internal_state_vec
end


"""
    elongation_process_pbc_profile(
        elongation_vector, internal_state_vec, posR, lattice, rates, current, k₊, mobile, jammed, l_ribosome, num_clust
    )

Perform one elongation step in the Gillespie simulation with periodic boundary conditions (PBC),
assuming all particles are active (no pausing).

- Selects a particle to move based on elongation rates weighted by `elongation_vector`.
- Updates the lattice by moving the particle forward by one site.
- Updates position vector `posR` that tracks particle positions.
- Updates `current` count for the site from which the particle moved (tracking flow).
- Updates internal state vector to follow particle movement.
- Updates elongation propensity vector accordingly by calling helper functions to adjust possible moves.

Returns the position of the moving particle.
"""
function elongation_process_pbc_profile(
    elongation_vector, internal_state_vec, posR, lattice, 
    rates, current, k₊, mobile, jammed, l_ribosome, num_clust
    )

    # Choose moving particle weighted by elongation rates
    moving_particle = sample(posR, Weights(elongation_vector))

    if 1 + l_ribosome <= moving_particle <= length(lattice) - l_ribosome
        # Bulk hopping inside lattice

        # Update lattice sites: vacate current, occupy next
        lattice[moving_particle] = 0
        lattice[moving_particle + 1] = 1

        # Update position vector accordingly
        posR[moving_particle + 1] = posR[moving_particle] + 1
        posR[moving_particle] = 0

        # Increment current flow count for site
        current[moving_particle] += 1

        # Transfer internal state to new position
        internal_state_vec[moving_particle + 1] = internal_state_vec[moving_particle]
        internal_state_vec[moving_particle] = 0

        # Update elongation propensity vector
        elongation_vector[moving_particle] = 0
        change_elong_vector_bulk_pbc(elongation_vector, moving_particle, lattice, rates, internal_state_vec, k₊, jammed, mobile, l_ribosome, num_clust)

    elseif moving_particle <= l_ribosome
        # Initiation step: particle hops from start

        lattice[moving_particle] = 0
        lattice[moving_particle + 1] = 1

        posR[moving_particle + 1] = posR[moving_particle] + 1
        posR[moving_particle] = 0

        current[moving_particle] += 1

        internal_state_vec[moving_particle + 1] = internal_state_vec[moving_particle]
        internal_state_vec[moving_particle] = 0

        elongation_vector[moving_particle] = 0
        change_elong_vector_init_pbc(elongation_vector, moving_particle, lattice, rates, internal_state_vec, k₊, mobile, jammed, l_ribosome, num_clust)

    elseif length(lattice) - l_ribosome + 1 <= moving_particle <= length(lattice)
        # Termination step with periodic boundary condition

        change_elong_vector_term_pbc(elongation_vector, moving_particle, lattice, rates, internal_state_vec, k₊, mobile, jammed, l_ribosome, num_clust)

        if moving_particle == length(lattice)
            current[moving_particle] += 1

            # Particle hops from last site back to first (p.b.c.)
            lattice[moving_particle] = 0
            lattice[1] = 1

            posR[moving_particle] = 0
            posR[1] = 1

            internal_state_vec[1] = internal_state_vec[moving_particle]
            internal_state_vec[moving_particle] = 0

        else
            # Particle moves forward within termination region (not last site)

            lattice[moving_particle] = 0
            lattice[moving_particle + 1] = 1

            posR[moving_particle + 1] = posR[moving_particle] + 1
            posR[moving_particle] = 0

            current[moving_particle] += 1

            internal_state_vec[moving_particle + 1] = internal_state_vec[moving_particle]
            internal_state_vec[moving_particle] = 0

        end

        elongation_vector[moving_particle] = 0 # Clear old site propensity
    end

    return moving_particle
end
