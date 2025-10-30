
function J_ρ_pbc_transient(number_of_runs, L, ρ, γ , k₋, k₊, run_time, delta_t, l_ribosome)

    worker_typ = Tuple{Matrix{Float64},Matrix{Float64}, Vector{Float64}}
    futures = Array{Future}(undef, nworkers())
    density = zeros(Float64, Int(run_time/delta_t)+1, L)
    current = zeros(Float64, Int(run_time/delta_t)+1, L)
    time = zeros(Float64, Int(run_time/delta_t)+1)

    for i in 1:number_of_runs

        for (i,id) in enumerate(workers())
            futures[i] = @spawnat id Gillespie_pbc_transient(L, ρ, γ, run_time, delta_t, k₋, k₊, l_ribosome; kymo=false)
        end   

        #w1::worker_typ, w2::worker_typ,w3::worker_typ, w4::worker_typ, w5::worker_typ= fetch.(futures)
        w1, w2 = fetch.(futures)

        w1::worker_typ, w2::worker_typ = fetch.(futures)

        # current::Vector{Float64} = vcat(w1[1], w2[1], w3[1], w4[1])
        # density::Matrix{Float64}= vcat(w1[2], w2[2], w3[2], w4[2])
        # clusters::Vector{Dict{Float64, Float64}}= vcat(w1[3], w2[3], w3[3], w4[3])
        tmp = 1
        current += w1[tmp] .+ w2[tmp] #.+ w3[tmp] .+ w4[tmp] .+ w5[tmp]# .+ w6[tmp] .+w7[tmp] .+ w8[tmp] .+w9[tmp] .+ w10[tmp] 
        tmp = 2
        density += w1[tmp] .+ w2[tmp] #.+ w3[tmp] .+ w4[tmp] .+ w5[tmp]# .+ w6[tmp] .+w7[tmp] .+ w8[tmp] .+w9[tmp] .+ w10[tmp] 
        tmp = 3
        time += w1[tmp] .+ w2[tmp] #.+ w3[tmp] .+ w4[tmp] .+ w5[tmp]# .+ w6[tmp] .+w7[tmp] .+ w8[tmp] .+w9[tmp] .+ w10[tmp] 
    end

    current = current ./ (number_of_runs * nworkers())
    density = density ./ (number_of_runs * nworkers())
    time = time ./ (number_of_runs * nworkers())
    return current, density, time
end

function J_ρ_pbc_transient_vs_ρ_list(number_of_runs, L, ρ_list, γ , k₋, k₊, run_time, delta_t, l_ribosome)
    current_vec = []
    density_vec = []
    time_vec = []

    for ρ in ρ_list
        current, density, time = J_ρ_pbc_transient(number_of_runs, L, ρ, γ , k₋, k₊, run_time, delta_t, l_ribosome)
        push!(current_vec, current)
        push!(density_vec, density)
        push!(time_vec, time)
    end
    
    return current_vec, density_vec, time_vec
end

function j_pbc(L, ρ_list, γ, k₋,k₊, number_of_measurements, steps_to_steady_state, number_of_steps, l_ribosome)
    
    # vectors with preset length equal to batch size / part of the rate list they are working on
    current = Matrix{Float64}(undef, length(ρ_list), 3) 
    density = Matrix{Float64}(undef, length(ρ_list) , 4) # gives mobile, jammed, paused, total density
    clusters = Vector{Dict{String, Dict{Float64, Float64}}}(undef, length(ρ_list))
    stdev = Vector{Tuple{Float64, Float64}}(undef, length(ρ_list))
    time = Matrix{Float64}(undef, length(ρ_list),3)
    i = 1 # counter for adding to preset length above
    
    for item in eachindex(ρ_list) # goes through rate of partial rate_list

        ρ = ρ_list[item] # in LD set the initiation rate
        tmpJ, tmpD, tmpC,tmpS, tmpT = Gillespie_pbc(L, ρ, γ, k₋,k₊, number_of_measurements, steps_to_steady_state, number_of_steps, l_ribosome)
        
        current[i,1:end] = tmpJ
        density[i,1:end] = tmpD
        clusters[i] = tmpC
        stdev[i] = tmpS
        time[i, 1:end] = tmpT
        i += 1
    end

    return current, density, clusters, stdev, time
end
 
function J_k₊_chunks(L, ρ, init, term, elong, k₋, k₊_list, run_time, starting_t, delta_t)
    
    # vectors with preset length equal to batch size / part of the rate list they are working on
    current = Vector{Float64}(undef, length(k₊_list)) 
    density = Matrix{Float64}(undef, length(k₊_list), 4) # gives mobile, jammed, paused, total density
    clusters = Matrix{Float64}(undef, length(k₊_list), 2)

    i = 1 # counter for adding to preset length above
    for item in eachindex(k₊_list) # goes through rate of partial rate_list

       k₊= k₊_list[item] # in LD set the initiation rate
        tmpJ, tmpD, tmpC = Gillespie_pbc(L, ρ, init, term, elong, k₋,k₊, run_time, starting_t, delta_t)
        
        current[item] = tmpJ
        density[item,1:end] = tmpD
        clusters[item,1:end] = tmpC
    
        i += 1
    end

    return current, density, clusters
end

function J_k₋_chunks(L, ρ, init, term, elong, k₋_list ,f, run_time, starting_t, delta_t)
    
    # vectors with preset length equal to batch size / part of the rate list they are working on
    current = Vector{Float64}(undef, length(k₋_list)) 
    density = Matrix{Float64}(undef, length(k₋_list), 4) # gives mobile, jammed, paused, total density
    clusters = Matrix{Float64}(undef, length(k₋_list), 2)

    i = 1 # counter for adding to preset length above
    for item in eachindex(k₋_list) # goes through rate of partial rate_list

        k₋::Float64 = k₋_list[item] # in LD set the initiation rate
        tmpJ, tmpD, tmpC = Gillespie_pbc(L, ρ, init, term, elong, k₋, k₊, run_time, starting_t, delta_t)
        
        current[i] = tmpJ
        density[i,1:end] = tmpD
        clusters[i,1:end] = tmpC
    
        i += 1
    end

    return current, density, clusters
end

function J_k₋_pbc(L, ρ, init, term, elong, k₋_list, k₊, run_time, starting_t, delta_t)

    worker_typ = Tuple{Vector{Float64},Matrix{Float64}, Matrix{Float64}}
    futures = Array{Future}(undef, nworkers())
    
    for (i,id) in enumerate(workers())
        batch_size = Int64(length(k₋_list)/(nworkers()))
        batch = 1:batch_size
        futures[i] = @spawnat id J_k₋_chunks(L, ρ, init, term, elong, k₋_list[batch .+ (i-1)*batch_size],k₊, run_time, starting_t, delta_t)
    end   
    
    w1::worker_typ, w2::worker_typ, w3::worker_typ, w4::worker_typ = fetch.(futures)

    current::Vector{Float64} = vcat(w1[1], w2[1], w3[1], w4[1])
    density::Matrix{Float64}= vcat(w1[2], w2[2], w3[2], w4[2])
    clusters::Matrix{Float64}= vcat(w1[3], w2[3], w3[3], w4[3])

    return current, density, clusters
end

function J_ρ_pbc(L, ρ_list, γ , k₋, k₊, number_of_measurements, steps_to_steady_state, number_of_steps, l_ribosome)

    worker_typ = Tuple{ Matrix{Float64}, Matrix{Float64}, Vector{Dict{String, Dict{Float64, Float64}}},Vector{Tuple{Float64, Float64}}, Matrix{Float64}}
    futures = Array{Future}(undef, nworkers())
    
    for (i,id) in enumerate(workers())
        batch_size = Int64(length(ρ_list)/(nworkers()))
        batch = 1:batch_size
        futures[i] = @spawnat id j_pbc(L, ρ_list[batch .+ (i-1)*batch_size],γ , k₋,k₊, number_of_measurements, steps_to_steady_state, number_of_steps,l_ribosome)
    end   
    
    #w1::worker_typ, w2::worker_typ, w3::worker_typ, w4::worker_typ,w5::worker_typ, w6::worker_typ, w7::worker_typ, w8::worker_typ,w9::worker_typ, w10::worker_typ = fetch.(futures)
    w1::worker_typ, w2::worker_typ, w3::worker_typ, w4::worker_typ,w5::worker_typ = fetch.(futures)
    
    tmp = 1
    current::Matrix{Float64} = vcat(w1[tmp], w2[tmp], w3[tmp], w4[tmp], w5[tmp], w6[tmp], w7[tmp], w8[tmp],w9[tmp], w10[tmp])
    tmp = 2
    density::Matrix{Float64}= vcat(w1[tmp], w2[tmp], w3[tmp], w4[tmp],w5[tmp], w6[tmp], w7[tmp], w8[tmp],w9[tmp], w10[tmp])
    tmp = 3
    clusters::Vector{Dict{String, Dict{Float64,Float64}}}= vcat(w1[tmp], w2[tmp], w3[tmp], w4[tmp],w5[tmp], w6[tmp], w7[tmp], w8[tmp],w9[tmp], w10[tmp])
    tmp = 4
    stdev::Vector{Tuple{Float64, Float64}} = vcat(w1[tmp], w2[tmp], w3[tmp], w4[tmp],w5[tmp], w6[tmp], w7[tmp], w8[tmp],w9[tmp], w10[tmp])
    tmp = 5
    time::Matrix{Float64} = vcat(w1[tmp], w2[tmp], w3[tmp], w4[tmp],w5[tmp], w6[tmp], w7[tmp], w8[tmp],w9[tmp], w10[tmp])

    return current, density, clusters, stdev, time
end

function J_k₊_pbc(L, ρ, init, term, elong, k₋, k₊_list, run_time, starting_t, delta_t)

    worker_typ = Tuple{Vector{Float64},Matrix{Float64}, Matrix{Float64}}
    futures = Array{Future}(undef, nworkers())
    
    for (i,id) in enumerate(workers())
        batch_size = Int64(length(k₊_list)/(nworkers()))
        batch = 1:batch_size
        futures[i] = @spawnat id J_k₊_chunks(L, ρ, init, term, elong, k₋, k₊_list[batch .+ (i-1)*batch_size], run_time, starting_t, delta_t)
    end   
    
    w1::worker_typ, w2::worker_typ, w3::worker_typ, w4::worker_typ,w5::worker_typ, w6::worker_typ, w7::worker_typ, w8::worker_typ,w9::worker_typ, w10::worker_typ = fetch.(futures)
    tmp = 1
    current::Vector{Float64} = vcat(w1[tmp], w2[tmp], w3[tmp], w4[tmp],w5[tmp], w6[tmp], w7[tmp], w8[tmp],w9[tmp], w10[tmp])
    tmp = 2
    density::Matrix{Float64}= vcat(w1[tmp], w2[tmp], w3[tmp], w4[tmp],w5[tmp], w6[tmp], w7[tmp], w8[tmp],w9[tmp], w10[tmp])
    tmp = 3
    clusters::Matrix{Float64}= vcat(w1[tmp], w2[tmp], w3[tmp], w4[tmp],w5[tmp], w6[tmp], w7[tmp], w8[tmp],w9[tmp], w10[tmp])

    return current, density, clusters
end

function J_ρ_pbc_k₋_k₊(L, ρ_list, γ , k₋_list, k₊_list, l_ribosome, number_of_measurements, steps_to_steady_state, number_of_steps)

    #J_vec::Vector{Vector{Float64}} = []
    #ρ_vec::Vector{Matrix{Float64}} = []
    #c_vec = []
    
    data_J = Matrix{Matrix{Float64}}(undef, length(k₋_list), length(k₊_list))
    data_ρ = Matrix{Matrix{Float64}}(undef, length(k₋_list), length(k₊_list))
    data_c = Matrix{Vector{Dict{String,Dict{Float64,Float64}}}}(undef, length(k₋_list), length(k₊_list))
    data_t = Matrix{Matrix{Float64}}(undef, length(k₋_list), length(k₊_list))
    data_s = Matrix{Vector{Tuple{Float64, Float64}}}(undef, length(k₋_list), length(k₊_list))
    
    for i in eachindex(k₋_list)
        k₋ = k₋_list[i]
        for j in eachindex(k₊_list)
            k₊ = k₊_list[j]

            if k₋ == k₊
                k₊ *= 0.99999999
            end

            J,ρ,c,s,t = J_ρ_pbc(L, ρ_list, γ , k₋, k₊, number_of_measurements, steps_to_steady_state, number_of_steps, l_ribosome)
            
            data_J[i,j] = J
            data_ρ[i,j] = ρ
            data_c[i,j] = c
            data_s[i,j] = s
            data_t[i,j] = t
            #push!(c_vec, c)
        end
    end

    return data_J, data_ρ, data_c, data_s, data_t
end



function J_ρ_profile_pbc(ρ_list, L, init, term, elong, k₋,k₊, run_time, starting_t, delta_t)

    J_vec = []
    a_vec = []
    p_vec = []

    i = 0
    for item in ρ_list
        i += 1
        if  i > length(ρ_list)
            break
        else
            ρ = ρ_list[i]
            r = remotecall(Gillespie_pbc, 2, L, ρ, init, term, elong, k₋,k₊, run_time, starting_t, delta_t)
        end

        i += 1
        if  i > length(ρ_list)
            break
        else
            ρ = ρ_list[i]
            s = remotecall(Gillespie_pbc, 3, L, ρ, init, term, elong, k₋,k₊, run_time, starting_t, delta_t)
        end
        
        i += 1 
        if  i > length(ρ_list)
            break
        else
            ρ = ρ_list[i]
            t = remotecall(Gillespie_pbc, 4, L, ρ, init, term, elong, k₋,k₊, run_time, starting_t, delta_t)
        end
    
        i +=1
        if  i > length(ρ_list)
            break
        else
            ρ = ρ_list[i]
            u = remotecall(Gillespie_pbc, 5, L, ρ, init, term, elong, k₋,k₊, run_time, starting_t, delta_t)
        end
       
        (Jr, ar, pr), (Js, as, ps), (Jt, at, pt), (Ju, au, pu) = fetch(r), fetch(s), fetch(t), fetch(u)
        
        push!(J_vec, Jr)
        push!(a_vec, ar)
        push!(p_vec, pr)

        push!(J_vec, Js)
        push!(a_vec, as)
        push!(p_vec, ps)

        push!(J_vec, Jt)
        push!(a_vec, at)
        push!(p_vec, pt)

        push!(J_vec, Ju)
        push!(a_vec, au)
        push!(p_vec, pu)
    end
    return hcat(ρ_list, J_vec), hcat(a_vec, p_vec)
end

# r = remotecall(Gillespie_pbc, 4, L, ρ, init, term, elong, k₋,k₊, run_time, starting_t, delta_t)
# u = remotecall(Gillespie_pbc, 3, L, ρ, init, term, elong, k₋,k₊, run_time, starting_t, delta_t)
# (x1, x2, x3),(x4,x5,x6) = fetch(r), fetch(u)

function k₋vs_k₊(k₋_list, k₊_list, L, init, term, elong, run_time, starting_t, delta_t)
    # element_type = Tuple(Vector{Array{Float64}})
    master_matrix = Matrix(undef,length(k₋_list), length(k₊_list)) # preset matrix for overall measurement stuff
    i= 1 # i is the row counter and j the column one
    for k₋ in k₋_list # start by setting the unpausing rate
        j = 1
        for k₊ in k₊_list # cycle through entering paused state rate
            
            ρ_max = Float64(round(ρ_max_fct(k₋,k₊, γ), digits=2)) # find rho max
            
            if k₋ == k₊
               k₊*= 0.9999
            end      
            ρ_list = collect(range(ρ_max/10, length = 20, stop=1-(ρ_max/10))) # ten measurement from the rho max in both directions
            tmpJ::Vector{Float64}, tmpD::Matrix{Float64} = J_ρ_pbc(L, ρ_list, init, term, elong, k₋, k₊, run_time, starting_t, delta_t)
            
            master_matrix[i,j] = [tmpJ, tmpD]
            j += 1
        end
        i += 1
    end
    return master_matrix
end

function J_vs_k₊(k₊_list,ρ_list, k₋, L, init, term, elong, run_time, starting_t, delta_t)

    master_vector = []
    
    for rate in k₊_list
    
        k₊=  rate # set the unpausing rate
        k₋ = 1 / k₋
        if k₊== k₋
           k₊*= 0.99999
        end
        tmpJ, tmpD = J_ρ_pbc(L, ρ_list, init, term, elong, k₋,k₊, run_time, starting_t, delta_t)
        
        push!(master_vector, (tmpJ, tmpD)) 
        
    end
    master_matrix = hcat(k₊_list, master_vector) 
    return master_matrix
end