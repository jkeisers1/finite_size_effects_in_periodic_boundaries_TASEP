module GillespieClusterAnalysis

using StatsBase

export sort_clusters, initiate_cluster_buckets, initiate_clustertip_bucket, find_all_cluster_lengths,
       cluster, back_counting, cluster_tip, count_backwards_pbc, 
       sort_clusters_length_and_dist_pbc, paused_particle_distribution_single_cluster_pbc

"""
Update cluster buckets with the counts of hole clusters (`h_clusters`) and particle clusters (`p_clusters`),
weighted by the elapsed `time_added`. Returns `true` if the sum of paused clusters times
their size does not match the total number of particles (error check).
"""

function sort_clusters(h_clusters, h_cluster_buckets, p_clusters, p_cluster_buckets, time_added, lattice)
    test = false
    h_clusters_classes = sort(unique(h_clusters))
    p_clusters_classes = sort(unique(p_clusters))

    # Count clusters of each size multiplied by elapsed time (weights)
    number_of_h_clusters = [count(x -> x == item, h_clusters) * time_added for item in h_clusters_classes]
    number_of_p_clusters = [count(x -> x == item, p_clusters) * time_added for item in p_clusters_classes]

    # Error check: total paused particles should equal total particles on lattice
    if sum(number_of_p_clusters .* p_clusters_classes) != sum(lattice)
        test = true
    end

    # Accumulate counts into buckets for holes and particles
    for cluster in eachindex(h_clusters_classes)
        h_cluster_buckets[h_clusters_classes[cluster]] += number_of_h_clusters[cluster]
    end

    for cluster in eachindex(p_clusters_classes)
        p_cluster_buckets[p_clusters_classes[cluster]] += number_of_p_clusters[cluster]
    end
    
    return test
end

"""
Initialize a dictionary to keep counts of paused clusters sized from -1 up to max expected cluster size.
"""
function initiate_cluster_buckets(ρ, L, l_ribosome)
    p_cluster_buckets = Dict{Float64, Float64}()
    max_cluster_size = ρ * L / l_ribosome
    for i in -1:1:max_cluster_size
        p_cluster_buckets[i] = 0.0
    end
    return p_cluster_buckets
end

"""
Initialize a dictionary bucket for cluster tips from -1 up to expected maximum number of tips.
"""
function initiate_clustertip_bucket(ρ, L)
    c_tip_bucket = Dict{Float64, Float64}()
    max_tip_count = ρ * L
    for i in -1:max_tip_count
        c_tip_bucket[i] = 0.0
    end
    return c_tip_bucket
end

"""
Find cluster lengths in a binary vector (1 = particle, 0 = hole),
updating `cluster_lengths` and `cluster_numbers` weighted by `time_added` and normalized by `num_clust`.
Returns the number of clusters found.
"""
function find_all_cluster_lengths(vector, cluster_lengths, cluster_numbers, time_added, num_clust)
    current_length = 0
    current_number = 0

    first_cluster = true
    first_cluster_length = 0

    for digit in vector
        if digit == 1
            current_length += 1
        elseif current_length > 0
            cluster_lengths[current_length] += time_added / num_clust[1]
            current_number += 1
            if first_cluster
                first_cluster_length += current_length
                first_cluster = false
            end
            current_length = 0
        end
    end

    if current_length > 0
        cluster_lengths[current_length] += time_added / num_clust[1]
        current_number += 1
    end

    # Account for periodic boundary: if cluster connects end and start, merge lengths
    if vector[end] == 1 && vector[1] == 1
        if current_length > 0
            cluster_lengths[current_length + first_cluster_length] += time_added / num_clust[1]
            cluster_lengths[first_cluster_length] -= time_added / num_clust[1]
            cluster_lengths[current_length] -= time_added / num_clust[1]
            current_number -= 1
        end
    end

    cluster_numbers[current_number] += time_added
    if current_number != num_clust[1]
        println("Warning: cluster count mismatch: $current_number vs old: $(num_clust[1])")
    end

    return current_number
end

"""
Find hole and particle clusters on the lattice.
Returns the vector of hole cluster sizes.
"""
function cluster(lattice)
    holes = findall(x -> x == 0, lattice)
    particles = findall(x -> x == 1, lattice)

    p_clusters = [holes[i+1] - holes[i] - 1 for i in 1:length(holes)-1]
    h_clusters = [particles[i+1] - particles[i] - 1 for i in 1:length(particles)-1]

    # Periodic boundary conditions cluster extension
    pbc_p_cluster = length(lattice) - holes[end] + holes[1] - 1
    pbc_h_cluster = length(lattice) - particles[end] + particles[1] - 1

    push!(p_clusters, pbc_p_cluster)
    push!(h_clusters, pbc_h_cluster)

    # Filter clusters smaller than 1
    filter!(x -> x >= 1, p_clusters)
    filter!(x -> x >= 1, h_clusters)

    # Avoid empty cluster lists (return 0 instead)
    if isempty(p_clusters)
        p_clusters = [0.0]
    end
    if isempty(h_clusters)
        h_clusters = [0.0]
    end

    return p_clusters
end

"""
Helper function to count backwards along the lattice from the leading particle's position,
tracking cluster length under periodic boundary conditions.
"""
function back_counting(pos_leading_part, distance_count, pbc_counter, k₊, lattice, internal_state_vec)
    if pos_leading_part != nothing && pos_leading_part - 1 - distance_count[1] <= 0
        if lattice[end - pbc_counter[1]] == 1 && internal_state_vec[end - pbc_counter[1]] == k₊
            distance_count[1] += 1
            pbc_counter[1] += 1
        end
    elseif pos_leading_part != nothing && lattice[pos_leading_part - 1 - distance_count[1]] == 1 && internal_state_vec[pos_leading_part - 1 - distance_count[1]] == k₊
        distance_count[1] += 1
    end
end

"""
Calculate the cluster tip size, counting backwards from the last paused particle, 
considering periodic boundaries.
"""
function cluster_tip(k₊, lattice, k₋, internal_state_vec, paused)
    distance_count = [0]
    pbc_counter = [0]
    last_hole = findlast(x -> x == 0, lattice)
    first_hole = findfirst(x -> x == 0, lattice)

    if sum(lattice[1:first_hole]) + sum(lattice[last_hole:end]) == length(lattice) - (last_hole - first_hole + 1) && lattice[1] == 1 && lattice[end] == 1
        pos_leading_part = findlast(x -> x == k₋, internal_state_vec[1:findlast(x -> x == 0, lattice)])
    else
        pos_leading_part = findlast(x -> x == k₋, internal_state_vec)
    end

    if pos_leading_part == length(lattice)
        pbc_counter[1] += 1
    end

    if paused[1] == 0
        distance_count[1] = -1
    else
        while back_counting(pos_leading_part, distance_count, pbc_counter, k₊, lattice, internal_state_vec) != nothing
            # Loop until no more backward counting possible
        end
    end

    return distance_count[1]
end

"""
Counts the length of the paused particle cluster going backwards from a start position
considering periodic boundaries.
"""
function count_backwards_pbc(lattice, internal_state_vec, pos_paused_particles, paused_dist_cluster, start_pos, k₊, time_added)
    count = 0
    stop = false

    while !stop
        if start_pos == pos_paused_particles[1] && sum(lattice[1:start_pos]) == start_pos
            pbc_count = 0
            pbc_stop = false
            while !pbc_stop
                if start_pos - count != 1 && lattice[start_pos - 1 - count] == 1 && internal_state_vec[start_pos - 1 - count] == k₊
                    count += 1
                elseif lattice[end - pbc_count] == 1 && internal_state_vec[end - pbc_count] == k₊
                    pbc_count += 1
                else
                    paused_dist_cluster[count + pbc_count] += time_added
                    pbc_stop = true
                    stop = true
                end
            end
        elseif lattice[start_pos - 1 - count] == 1 && internal_state_vec[start_pos - 1 - count] == k₊
            count += 1
        else
            paused_dist_cluster[count] += time_added
            stop = true
        end
    end

    return count
end

"""
Sort clusters by length and paused particle distribution,
tracking paused particle cluster lengths weighted by elapsed time.
"""
function sort_clusters_length_and_dist_pbc(
    lattice, internal_state_vec, paused_dist_cluster, 
    k₊, time_added, k₋)

    pos_paused_particles = findall(x -> x == k₋, internal_state_vec)
    paused_particle_distribution_single_cluster_pbc(lattice, internal_state_vec, paused_dist_cluster, pos_paused_particles, k₊, time_added)
end

"""
For each paused particle in the lattice, count backwards the cluster length and
accumulate cluster length distribution in `paused_dist_cluster`.
"""
function paused_particle_distribution_single_cluster_pbc(
    lattice, internal_state_vec, paused_dist_cluster, 
    pos_paused_particles, k₊, time_added)

    if !isempty(pos_paused_particles)
        for start_pos in reverse(pos_paused_particles)
            count_backwards_pbc(lattice, internal_state_vec, pos_paused_particles, paused_dist_cluster, start_pos, k₊, time_added)
        end
    else
        # When no paused particles, increment bucket for -1 (special case)
        paused_dist_cluster[-1] += time_added
    end
end

end # module
