module GillespieTransientSimulations

using Random
using Distributions

export Gillespie_pbc_transient, Gillespie_pbc_transient_text_file, get_transient_profiles

"""
Run a Gillespie simulation with periodic boundary conditions (PBC) on a lattice of length `L`
with ribosome density `ρ`, hopping rate `γ`, over a total time `run_time`, measuring every `delta_t`.
Internal states switch with rates `k₋` (paused → active) and `k₊` (active → paused).
Ribosome length is `l_ribosome`.
Returns time-resolved matrices of current and density and the measurement times.
"""
function Gillespie_pbc_transient(L, ρ, γ, run_time, delta_t, k₋, k₊, l_ribosome; kymo=false)
    
    number_of_measurements = Int(run_time / delta_t)
    density_matrix = Matrix{Float64}(undef, number_of_measurements+1, L)
    current_matrix = Matrix{Float64}(undef, number_of_measurements+1, L)

    # Initialize counts of mobile, jammed, and paused particles as Float arrays of length 1 (mutable)
    mobile = [0.0]
    jammed = [0.0]
    paused = [0.0]

    # Initialize vectors for weighted statistics over time (not used in return, but for internal tracking)
    ρ_weighted = [0.0]
    mobile_weighted = [0.0]
    jammed_weighted = [0.0]
    paused_weighted = [0.0]

    cluster_num_bucket = initiate_cluster_buckets(ρ, L, l_ribosome) # assumed external function
    paused_dist_bucket = initiate_cluster_buckets(ρ, L, l_ribosome) # assumed external function

    current = zeros(L)  # Tracks hops (currents) per lattice site
    J = [0.0]          # Total current counter at the last site

    # Initialize lattice and rates
    lattice = create_lattice_l_stuck(L, ρ, l_ribosome) # assumed external function
    rates = create_rates(L, γ)                         # assumed external function
    posR = positionRibosomes(lattice)                  # get ribosome positions vector
    elongation_vector = get_elongation_vector_pbc(lattice, rates, mobile, jammed, l_ribosome)
    internal_state_vec = get_internal_state_vec_pbc_profile(lattice, elongation_vector, k₋, k₊, mobile, jammed, paused)

    # Calculate initial counts of mobile, paused, and jammed particles
    m_start = sum(elongation_vector .== γ)
    p_start = sum(internal_state_vec .== k₋)
    j_start = sum(lattice) - m_start - p_start
    mobile, paused, jammed = [m_start], [p_start], [j_start]

    density = zeros(L)
    time_vec = [0.0]
    collecting_time = delta_t

    density_matrix[1, :] = lattice 
    current_matrix[1, :] = current

    # Initialize cluster count
    if l_ribosome == 1
        num_clust = [length(cluster(lattice))]  # cluster is assumed external
    else
        num_clust = [1.0]
    end

    # Avoid exact equality in rates to prevent division issues
    if k₋ == k₊
       k₊ *= 0.99999
    end

    t = 0.0  # simulation time
    J[1] = 0
    i = 2    # measurement index

    while t <= run_time
        elong_sum = sum(elongation_vector)
        state_sum = sum(internal_state_vec)
        total_sum = elong_sum + state_sum

        RV = rand()
        time_added = rand(Exponential(1/total_sum))
        t += time_added

        # Weight populations by time step (not returned)
        mobile_weighted[1] += (mobile[1]/L) * time_added
        paused_weighted[1] += (paused[1]/L) * time_added
        jammed_weighted[1] += (jammed[1]/L) * time_added
        
        density += lattice .* time_added

        # Evolve system by elongation or internal state change based on random choice
        if RV <= elong_sum / total_sum
            elongation_process_pbc_profile(elongation_vector, internal_state_vec, posR, lattice, rates, current, k₊, mobile, jammed, l_ribosome, num_clust)
        else
            change_internal_state_pbc(internal_state_vec, posR, lattice, rates, elongation_vector, k₋, k₊, mobile, jammed, paused, l_ribosome)
        end

        # Record averages at fixed delta_t intervals
        if t >= collecting_time
            push!(time_vec, t)
            time_past = time_vec[end] - time_vec[end-1]

            density_matrix[i, :] = density ./ time_past
            current_matrix[i, :] = current ./ time_past

            current .= 0.0
            density = lattice .* time_added
            t = collecting_time
            collecting_time += delta_t
            i += 1
        end
    end

    return current_matrix, density_matrix, time_vec
end


"""
Run Gillespie simulation as above, but write transient density profiles to a text file at each delta_t.
Returns the number of measurements written.
"""
function Gillespie_pbc_transient_text_file(
    filename, L, ρ, γ, run_time, delta_t, k₋, k₊, l_ribosome; kymo=false)

    number_of_measurements = Int(run_time / delta_t)
    density_matrix = Matrix{Float64}(undef, number_of_measurements+1, L)
    current_matrix = Matrix{Float64}(undef, number_of_measurements+1, L)

    mobile = [0.0]
    jammed = [0.0]
    paused = [0.0]

    ρ_weighted = [0.0]
    mobile_weighted = [0.0]
    jammed_weighted = [0.0]
    paused_weighted = [0.0]

    cluster_num_bucket = initiate_cluster_buckets(ρ, L, l_ribosome)
    paused_dist_bucket = initiate_cluster_buckets(ρ, L, l_ribosome)

    current = zeros(L)
    J = [0.0]

    lattice = create_lattice_l_stuck(L, ρ, l_ribosome)
    rates = create_rates(L, γ)
    posR = positionRibosomes(lattice)
    elongation_vector = get_elongation_vector_pbc(lattice, rates, mobile, jammed, l_ribosome)
    internal_state_vec = get_internal_state_vec_pbc_profile(lattice, elongation_vector, k₋, k₊, mobile, jammed, paused)

    m_start = sum(elongation_vector .== γ)
    p_start = sum(internal_state_vec .== k₋)
    j_start = sum(lattice) - m_start - p_start
    mobile, paused, jammed = [m_start], [p_start], [j_start]

    density = zeros(L)
    time_vec = [0.0]
    collecting_time = delta_t

    density_matrix[1, :] = lattice
    current_matrix[1, :] = current

    if l_ribosome == 1
        num_clust = [length(cluster(lattice))]
    else
        num_clust = [1.0]
    end

    if k₋ == k₊
       k₊ *= 0.99999
    end

    t = 0.0
    J[1] = 0
    i = 1

    open(filename, "w") do file
        # Write initial lattice state
        println(file, "$t: $(join(lattice, ','))")
        while t <= run_time
            elong_sum = sum(elongation_vector)
            state_sum = sum(internal_state_vec)
            total_sum = elong_sum + state_sum
            RV = rand()
            time_added = rand(Exponential(1/total_sum))
            t += time_added

            mobile_weighted[1] += (mobile[1]/L) * time_added
            paused_weighted[1] += (paused[1]/L) * time_added
            jammed_weighted[1] += (jammed[1]/L) * time_added

            density += lattice .* time_added

            if RV <= elong_sum / total_sum
                elongation_process_pbc_profile(elongation_vector, internal_state_vec, posR, lattice, rates, current, k₊, mobile, jammed, l_ribosome, num_clust)
            else
                change_internal_state_pbc(internal_state_vec, posR, lattice, rates, elongation_vector, k₋, k₊, mobile, jammed, paused, l_ribosome)
            end

            # Write densities to file at each delta_t interval
            if t >= i * delta_t
                tmp_den = density ./ t
                println(file, "$t: $(join(tmp_den, ','))")
                i += 1
            end
        end
    end

    return i  # number of measurements written
end


"""
Run multiple Gillespie transient simulations and average the resulting
current and density profiles over `number_of_runs`.
"""
function get_transient_profiles(
    number_of_runs, L, ρ, γ, run_time, delta_t, k₋, k₊, l_ribosome; kymo=false)

    J_matrix = zeros(Int(run_time/delta_t), L)
    ρ_matrix = zeros(Int(run_time/delta_t), L)
    time_vector = zeros(Int(run_time/delta_t))

    for i in 1:number_of_runs
        tmpJ, tmpρ, tmp_t = Gillespie_pbc_transient(L, ρ, γ, run_time, delta_t, k₋, k₊, l_ribosome; kymo=false)
        J_matrix += tmpJ
        ρ_matrix += tmpρ
        time_vector += tmp_t
    end

    J_matrix ./= number_of_runs
    ρ_matrix ./= number_of_runs
    time_vector ./= number_of_runs
    
    return J_matrix, ρ_matrix, time_vector
end


end # module
