"""
    find_max_stdv(Matrix)

Finds the maximum standard deviation value inside a nested matrix structure.

- `Matrix` is expected to be a 2D array where each element contains pairs `(stdC, stdt)`.
- Returns a tuple `(maximal_stdv, index)` where `maximal_stdv` is the largest standard deviation found
  and `index` is a 2-element vector indicating its position in `Matrix`.
"""
function find_max_stdv(Matrix)
    maximal_stdv = 0
    index = [0,0]
    for i in axes(Matrix,1)
        for j in axes(Matrix,2)
            max_std = maximum([max(stdC,stdt) for (stdC,stdt) in Matrix[i,j]])
            if max_std > maximal_stdv
                maximal_stdv = max_std
                index[1] = i
                index[2] = j
            end
        end
    end

    return maximal_stdv, index
end


"""
    check_mobile(elongation_vector, mobile)

Checks consistency between the elongation vector and the count of mobile particles.

- Returns `true` if the number of nonzero elongation rates equals `mobile[1]`.
- Otherwise, returns `false`.
"""
function check_mobile(elongation_vector, mobile)
    check_counter = 0.0

    for item in elongation_vector
        if item != 0
            check_counter += 1
        end
    end
    return check_counter == mobile[1]
end


"""
    check_jammed(lattice, jammed, internal_state_vec, k₊)

Checks if the count of jammed particles is consistent with the lattice configuration.

- A particle is considered jammed if it is active (`internal_state_vec[i] == k₊`) and the site ahead is occupied.
- Also considers periodic boundary condition between the last and first site.
- Returns `true` if jammed count matches `jammed[1]`, else `false`.
"""
function check_jammed(lattice, jammed, internal_state_vec, k₊)
    check_counter = 0.0
    for i in eachindex(lattice[1:end-1])
        if lattice[i] == 1 && internal_state_vec[i] == k₊ && lattice[i+1] == 1
            check_counter += 1
        end
    end

    # Check PBC condition between last and first site
    if lattice[1] == 1 && lattice[end] == 1 && internal_state_vec[end] == k₊
        check_counter += 1
    end

    return check_counter == jammed[1]
end


"""
    check_cluster(lattice, num_clust, l_ribosome)

Counts the number of clusters of ribosomes on the lattice.

- A cluster is a pair of occupied sites separated by `l_ribosome` sites.
- Checks bulk and periodic boundary conditions.
- Returns the count of clusters and verifies it equals `num_clust[1]`.
"""
function check_cluster(lattice, num_clust, l_ribosome)
    count = 0

    for index in eachindex(lattice[1:end-l_ribosome])
        if lattice[index] == 1 && lattice[index+l_ribosome] == 1
            count += 1
        end
    end

    for index in eachindex(lattice[end-l_ribosome+1:end])
        if lattice[end+1-index] == 1 && lattice[l_ribosome+1-index] == 1
            count += 1
        end
    end

    return count
end


"""
    check_clust(internal_state_vec, lattice, l_ribosome, mobile)

Verifies the consistency of cluster counts including mobile and paused particles near boundaries.

- Counts special cases near the lattice edges to accommodate PBC.
- Returns total count of mobile particles plus detected cluster edges.
"""
function check_clust(internal_state_vec, lattice, l_ribosome, mobile)
    count1 = 0
    count2 = 0
    count3 = 0
    pos_paused = findall(x->x == k₋, internal_state_vec)
    for pos in pos_paused
        distance_to_end = length(lattice) - pos

        overlap = l_ribosome == 1 ? 0 : l_ribosome - distance_to_end - 1

        if pos == length(lattice) && lattice[l_ribosome] == 0
            count1 += 1
        elseif l_ribosome > 1 && pos > length(lattice) - l_ribosome && lattice[1+overlap] == 0
            count2 += 1
        elseif 1 <= pos <= length(lattice) - l_ribosome && lattice[pos+l_ribosome] == 0
            count3 += 1
        end
    end

    return mobile[1] + count1 + count2 + count3
end


"""
    plot_kymo(lattice_matr, internal_state_matr, k₊, k₋, t_vec)

Generates a kymograph plot of ribosome lattice occupancy over time.

- Uses colors to represent mobile (lightgreen), jammed (lightblue), and paused (black) ribosomes.
- `lattice_matr` and `internal_state_matr` are matrices with lattice states and internal states at each timepoint.
- `k₊` and `k₋` represent active and paused state rates.
- `t_vec` is the vector of time points corresponding to each row in the matrices.

Returns a plot object.
"""
function plot_kymo(lattice_matr, internal_state_matr, k₊, k₋, t_vec)
    if k₊ == k₋
        k₊ *= 0.999999
    end
    plo = scatter(legend=false, xlabel="lattice site", ylabel="t") # initiate kymo plot
    i = 1
    for time in t_vec
        lattice = lattice_matr[i, 1:end]
        internal_state_vec = internal_state_matr[i, 1:end]
        t = time
        for index in eachindex(lattice[1:end-1])
            # mobile (green square)
            if lattice[index] == 1 && lattice[index+1] == 0 && internal_state_vec[index] == k₊
                scatter!([index], [t], markercolor="lightgreen", markersize=2.5, markershape=:rect, markerstrokewidth=0)
            end

            # jammed (blue square)
            if lattice[index] == 1 && lattice[index+1] == 1 && internal_state_vec[index] == k₊
                scatter!([index], [t], markercolor="lightblue", markersize=2.5, markershape=:rect, markerstrokewidth=0)
            end

            # paused (black square)
            if lattice[index] == 1 && internal_state_vec[index] == k₋
                scatter!([index], [t], markercolor="black", markersize=2.5, markershape=:rect, markerstrokewidth=0)
            end
        end

        # periodic boundary conditions for last site
        if lattice[end] == 1 && lattice[1] == 1 && internal_state_vec[end] == k₊
            scatter!([length(lattice)], [t], markercolor="lightgrey", markersize=2.5, markershape=:rect, markerstrokewidth=0)
        end

        if lattice[end] == 1 && lattice[1] == 0 && internal_state_vec[end] == k₊
            scatter!([length(lattice)], [t], markercolor="lightgreen", markersize=2.5, markershape=:rect, markerstrokewidth=0)
        end

        if lattice[end] == 1 && internal_state_vec[end] == k₋
            scatter!([length(lattice)], [t], markercolor="black", markersize=2.5, markershape=:rect, markerstrokewidth=0)
        end

        i += 1
    end
    
    return plo
end
