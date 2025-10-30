using Blink: Window, body!
using Interact: @manipulate
using Plots
using JLD
using Colors


data_path = "C:\\Users\\johan\\OneDrive\\Desktop\\Work\\Montpellier\\Code\\Data\\"
#plotlyjs()
# include("analyitcal_functions.jl")
begin
  # @load "data_k₋2_10_k₊005_1.jld"
  # @load "data_k₋2_10_k₊005_1.jld"
  # @load "data_k₋005_1_k₊2_10.jld"
  # @load "data_k₋2_10_k₊2_10.jld"
  # @load "matrix_k₋01_1_k₊01_1.jld"

  # @load "data_16_01.jld"
  # data_k₋2_10_k₊01_1 = data_k₋2_10_k₊005_1[1:end,2:2:end]
  # data_k₋01_1_k₊2_10 = data_k₋005_1_k₊2_10[2:2:end,1:end]
  # data_k₋01_1_k₊01_1 = matrix_k₋01_1_k₊01_1
  # data_k₋2_10_k₊2_10
  # data = vcat(hcat(data_k₋01_1_k₊01_1,  data_k₋01_1_k₊2_10), hcat(data_k₋2_10_k₊01_1, data_k₋2_10_k₊2_10))

  # k₋_list = vcat(collect(0.1:0.1:1), collect(2:1:10))
  # k₊_list = vcat(collect(0.1:0.1:1), collect(2:1:10))
  # colors = distinguishable_colors(10)[5:10]
  

  # jammed_test(ρ,k₋, k₊) = (ρ * (elong * mobile_test(ρ,k₋,k₊) + k₋*ρ*k₊/(k₊ + k₋)- elong*mobile_test(ρ,k₋,k₊)))/k₊
  # paused_test(ρ, k₋, k₊) = ρ * k₊ / (k₊ + k₋)
  # mobile_test(ρ,k₋,k₊) = k₋*ρ*(1-ρ) *Φ1(k₋, k₊)/(k₊+elong*ρ)

  
end

function sort_and_norm_dict(dict)
  filtered_dict = Dict(key => value for (key, value) in dict if value !=0)
  sorted_dict = sort(filtered_dict)
  total_sum = sum(values(sorted_dict))
  norm_dic = Dict(k => v / total_sum for (k,v) in sorted_dict)
  # mean_dic = sum(Dict(k => v / total_sum for (k,v) in sorted_dict))
  return norm_dic
end


function plot_J_vs_ρ(J_data,ρ_data)
  w = Window()
  ui = @manipulate for ρ in ρ_list
    index1 = findfirst(x-> x==ρ, ρ_list)
    plot(ylabel = "current", xlabel = "density", xtickfont = 16, ytickfont=16, guidefont = 16, xlims  =(0,1))
    scatter!(ρ_list , data[index1, index2][1], label ="numerics", markersize = 3)
  end
end

function plot_number_of_clusters(ρ_list, k₋_list, k₊_list, L, data)
  w = Window()
  ui = @manipulate for k₋ in k₋_list, k₊ in k₊_list, ρ in ρ_list

      index1 = findfirst(x-> x==k₊, k₊_list)
      index2 = findfirst(x-> x==k₋, k₋_list)
      index3 = findfirst(x-> x == ρ, ρ_list)

      dict1 = data[index2, index1][index3]["NumberOfClusters"]
      norm_dict1 = sort_and_norm_dict(dict1)
      mean_dict1 = round(sum(k * v for (k,v) in norm_dict1), digits= 2)

      max_particles = floor(ρ * L) 

      if max_particles in keys(norm_dict1)
        max_key = maximum(keys(norm_dict1)) 
        max_length = round(norm_dict1[max_key], digits = 4)
      else 
        max_length = 0
      end

      p1 = bar(
        norm_dict1, 
        xlabel = "Number Of Clusters", ylabel = "Distribution",
        label = "mean:$mean_dict1" , title = "$max_length"
        )
      
      dict2 = data[index2, index1][index3]["ClusterLengthDistribution"]
      norm_dict2 = sort_and_norm_dict(dict2)
      
      if 1 in keys(norm_dict2)  
        one_cluster = round(norm_dict2[1], digits = 4)
      else
        one_cluster = 0
      end

      mean_dict2 = round(sum(k * v for (k,v) in norm_dict2), digits = 2)

      p2 = bar(
        dict2, 
        xlabel = "Cluster Length", ylabel = "Distribution",
        label = "mean:$mean_dict2", title = "$one_cluster" 
        )

      plot(p1, p2)
  end 
  body!(w, ui)
end

function plot_current_whoosh_cluster(data, ρ_list, L; error = true)
  w = Window()
  γ = 1
  ui = @manipulate for k₋ in k₋_list, k₊ in k₊_list
      index2 = findfirst(x-> x==k₋, k₋_list)
      index1 = findfirst(x-> x==k₊, k₊_list)

      xaxis_inp = ρ_list[1:end]

      p1 = Plots.scatter(
          xaxis_inp, data[index2, index1][1:end,1], label =false, markercolor = colors[1], 
          ylabel = "J", xlabel="ρ", title = "total current", xtickfont = 10, ytickfont=10, guidefont = 10, xlims  =(0,1),
          legendfont = 7
          )
          #Plots.plot!(xaxis_inp, greulich.(xaxis_inp, k₋, k₊,γ),label = "Greulich", linewidth = 6, linecolor = colors[3],alpha = 0.7)
          Plots.plot!(xaxis_inp, J_whoosh.(γ, xaxis_inp, L, k₋, k₊), label = "pTASEP+TASEP", linewidth = 6, linecolor = colors[4],alpha = 0.7)
          Plots.plot!(xaxis_inp, J_transient.(xaxis_inp,k₊,k₋,L), label = "Single Cluster", linewidth = 6,  linecolor = colors[5],alpha = 0.7)
          #Plots.plot!(xaxis_inp,  Wang.(xaxis_inp, k₋, k₊, γ), label = "pTASEP", linewidth = 6, linecolor = colors[2], alpha = 0.7)

      p2 = Plots.scatter(
          xaxis_inp, data[index2, index1][1:end,2],label =false,  markercolor = colors[1],
          ylabel = "JW", xlabel="ρ", title = "Pw × Jw ", xtickfont = 10, ytickfont=10, guidefont = 10, xlims  =(0,1)
          )
          plot!(xaxis_inp, γ .* xaxis_inp .* (1 .- xaxis_inp) .* pw_fct.(xaxis_inp, L, k₋, k₊), linewidth = 6, label = "pTASEP+TASEP",  linecolor = colors[4])
          plot!(xaxis_inp, J_transient_w.(xaxis_inp,k₊,k₋,L), linewidth = 6,  linecolor = colors[5], label = "Single Cluster")

      p3 = Plots.scatter(
          xaxis_inp, data[index2, index1][1:end,3],label =false, markercolor = colors[1],
          ylabel = "JC", xlabel="ρ", title = "(1-Pw) × Jc", xtickfont = 10, ytickfont=10, guidefont = 10, xlims  =(0,1)
          )
          plot!(xaxis_inp, Wang.(xaxis_inp, k₋, k₊, γ) .* (1 .- pw_fct.(xaxis_inp, L, k₋, k₊)), linewidth = 6, label = "pTASEP+TASEP",  linecolor = colors[4])   
          plot!(xaxis_inp, J_transient_c.(xaxis_inp,k₊,k₋,L), linewidth = 6, label = "Single Cluster", linecolor = colors[5])

      if error == false
        p4 = bar(
                xaxis_inp, data[index2,index1][1:end,3]./ data[index2,index1][1:end,1], label = "(1-Pw) × Jc", alpha = 0.7, color = colors[5]
            )
            bar!(
            xaxis_inp, data[index2,index1][1:end,2] ./ data[index2,index1][1:end,1]  , label = "Pw × Jw", color = colors[end], alpha = 0.7,
            ylabel = "Distribution", xlabel="ρ", title = "Current Contribution", xtickfont = 10, ytickfont=10, guidefont = 10, xlims  =(0,1)
            )
      else
          numerical = data[index2, index1][1:end,1]
          ana_W = Wang.(xaxis_inp, k₋, k₊, γ)
          ana_WW = J_whoosh.(γ, xaxis_inp, L, k₋, k₊)
          ana_G = greulich.(xaxis_inp, k₋, k₊,γ)
          ana_Our = J_transient.(xaxis_inp,k₊,k₋,L)

          sum_rel_err_W = round(sum(relative_error(numerical, ana_W)), digits = 2)
          RMSD_W = round(sqrt(mean_squared_error(numerical, ana_W)), digits = 4)

          sum_rel_err_WW = round(sum(relative_error(numerical, ana_WW)), digits = 2)
          RMSD_WW = round(sqrt(mean_squared_error(numerical, ana_WW)), digits = 4)

          sum_rel_err_G = round(sum(relative_error(numerical, ana_G)), digits = 2)
          RMSD_G = round(sqrt(mean_squared_error(numerical, ana_G)), digits = 4)

          sum_rel_err_Our = round(sum(relative_error(numerical, ana_Our)), digits = 2)
          RMSD_Our = round(sqrt(mean_squared_error(numerical, ana_Our)), digits = 4)

          p4 = Plots.scatter(
            xaxis_inp, relative_error(numerical, ana_W), 
            label="Sum:$sum_rel_err_W",#\nRMSD:$RMSD_W",
            xlabel = "ρ", ylabel = "relative error", color = colors[2],
            markershape =:square,
            title = "Error"
            )
            Plots.scatter!(
              xaxis_inp, relative_error(numerical, ana_G), 
              label="Sum:$sum_rel_err_G",#\nRMSD:$RMSD_G",
              xlabel = "ρ", ylabel = "relative error", color = colors[3],
              markershape = :hexagon
            )
            Plots.scatter!(
              xaxis_inp, relative_error(numerical, ana_WW), 
              label="Sum:$sum_rel_err_WW",#\nRMSD:$RMSD_WW",
              xlabel = "ρ", ylabel = "relative error", color = colors[4],
              markershape = :diamond
            )
            Plots.scatter!(
              xaxis_inp, relative_error(numerical, ana_Our), 
              label="Sum:$sum_rel_err_Our",#\nRMSD:$RMSD_Our",
              color = colors[5], xlabel = "ρ", ylabel = "relative error", legendfont = 7, legendalpha = 0.7
            )
      end
          #plot!(ρ_list, Cluster_Wang_fraction.(γ, ρ_list, L, k₋, k₊), linewidth = 6, label = false, linecolor = colors[5])
          #plot!(ρ_list, Woosh_fraction.(γ, ρ_list, L, k₋, k₊),  linewidth = 6, label =false, linecolor = colors[end])
      # plot!(t,  test2_0.(t,k₋,k₊, γ), label = "Mobile Particles")
      # plot!(t, greulich.(t, k₋, k₊), label = "Greulich")
      #scatter!(t,  mobile_test.(t,k₋,k₊), label = "New Analytical Solution")
      #scatter!(t, test2_0.(t, k₋, k₊))

      Plots.plot(p1,p2,p3,p4, layout= grid(2,2))
  end
  body!(w, ui)
end

function plot_current(data,k₋_list, k₊_list)
  w = Window()
  t = 0.05:0.02:0.95
  ui = @manipulate for k₋ in k₋_list, k₊ in k₊_list
    index1 = findfirst(x-> x==k₋, k₋_list)
    index2 = findfirst(x-> x==k₊, k₊_list)
    plot(ylabel = "J", xlabel="ρ", title = "k₋=$k₋,k₊=$k₊", xtickfont = 16, ytickfont=16, guidefont = 16, xlims  =(0,1))
    #plot!(t, J_ana.(γ,t,k₋, k₊), label = "Wang + ν + ψ ", linecolor= colors[1])
    scatter!(data[index1, index2][2][1:end,4], data[index1, index2][1], label ="numerics", markersize = 3)
    τ,f = 1/k₋, k₊
    plot!(t, Wang.(t, k₋, k₊, γ), label = "Wang")
    # plot!(t,  test2_0.(t,k₋,k₊, γ), label = "Mobile Particles")
    # plot!(t, greulich.(t, k₋, k₊), label = "Greulich")
    #scatter!(t,  mobile_test.(t,k₋,k₊), label = "New Analytical Solution")
    #scatter!(t, test2_0.(t, k₋, k₊))
  end
  body!(w, ui)
end

function plot_current_obc(data, γ,k₋_list, k₊_list)
  w = Window()
  # k₊_list =[0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 2.0, 4.0, 6.0, 8.0, 10.0] #[0.1, 0.2, 0.4, 0.6, 0.8]#,[ 1.0, 2, 4, 6, 8, 10]
  # k₋_list =[0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 2.0, 4.0, 6.0, 8.0, 10.0]
  
  ui = @manipulate for k₋ in k₋_list, k₊ in k₊_list
    index1 = findfirst(x-> x==k₋, k₋_list)
    index2 = findfirst(x-> x==k₊, k₊_list)
    plot(ylabel = "J", xlabel="ρ", title = "k₋=$k₋,k₊=$k₊", xtickfont = 16, ytickfont=16, guidefont = 16, xlims  =(0,1))
    #plot!(t, J_ana.(γ,t,k₋, k₊), label = "Wang + ν + ψ ", linecolor= colors[1])
    scatter!(data[index1, index2][2][1:8,4], data[index1, index2][1][1:8,2], label ="LD", markersize = 3, color = "blue")
    scatter!(data[index1, index2][2][9:end,4], data[index1, index2][1][9:end,2], label ="HD", markersize = 3, color = "red")
    # scatter!(data_pbc[index1, index2][2][1:end,4], data_pbc[index1, index2][1][1:end], label ="PBC", markersize = 3, color = "green")
    τ,f = 1/k₋, k₊
    # ρ_step_HD = ρ_HD(β_crit_fct(k₋, k₊, elong), 1):0.001:1
    # ρ_step_LD = 0:0.01:ρ_LD(α_crit_fct(k₋,k₊, elong), k₋, k₊)
    β_step = 0:0.0001:β_crit_fct(k₋, k₊, γ)
    α_step = 0:0.0001:α_crit_fct(k₋, k₊ , γ)

    # plot!(t, Wang.(t, k₋, k₊, elong), label = "Wang")
    # plot!(ρ_step_HD, Wang.(ρ_step_HD, k₋, k₊, γ), label = false, color = "red")
    # plot!(ρ_step_LD, Wang.(ρ_step_LD, k₋, k₊, γ), label = false, color = "blue")
    
    plot!(ρ_HD.(β_step,1), Wang_HD.(β_step, k₋, k₊, 1), label = false, color = "red")
    plot!(ρ_LD.(α_step, k₋, k₊), Wang_LD.(α_step, k₋, k₊, 1), label = false, color = "blue")

    # plot!(ρ_step, Wang.(ρ_step,k₋, k₊, γ))
    #plot!(t,  test2_0.(t,k₋,k₊, γ), label = "Mobile Particles")
    #plot!(t, greulich.(t, k₋, k₊), label = "Greulich")
    #scatter!(t,  mobile_test.(t,k₋,k₊), label = "New Analytical Solution")
    #scatter!(t, test2_0.(t, k₋, k₊))
  end
  body!(w, ui)
end

function plot_particle_classes(data, k₋_list, k₊_list)
  w = Window()
  # k₋_list = vcat(collect(0.1:0.1:1), collect(2:1:10))
  # k₊_list = vcat(collect(0.1:0.1:1), collect(2:1:10))
  t = 0:0.02:0.96
  ui = @manipulate for k₋ in k₋_list, k₊ in k₊_list
    index1 = findfirst(x-> x==k₋, k₋_list)
    index2 = findfirst(x-> x==k₊, k₊_list)
    τ,f = 1/k₋, k₊
    p2 = plot(
      ylabel = "particle density", 
      xlabel="ρ", 
      legend=:topleft, 
      xlims=(0,1), 
      title = "k₋ = $k₋, k₊ = $k₊", 
      xtickfont = 16, 
      ytickfont=16, 
      guidefont = 16,
      ylims = (-0.04,1),
      yticks = (0:0.2:1)
      )#, ylims = (0,1))

    scatter!(data[index1, index2][2][1:end,4], data[index1, index2][2][1:end,3], label ="jammed", markercolor = colors[1], markershape=:circle)
    plot!(t, ν.(t, k₋, k₊, γ) , linecolor= colors[1], label = false, linestyle = :dashdotdot, linewidth= 3)
    # plot!(t,ν.(γ,t,k₋, k₊) , linecolor= colors[1], label = false, linestyle = :dashdotdot, linewidth= 3)

    #scatter!(data[index1, index2][2][1:end,4], data[index1, index2][2][1:end,4] .- data[index1, index2][2][1:end,2] .- data[index1, index2][2][1:end,3])

    scatter!(data[index1, index2][2][1:end,4], data[index1, index2][2][1:end,2], label ="paused", markercolor = colors[3], markershape=:cross, markerstrokewidth = 2)
    plot!(t, ψ.(t,k₋, k₊), label = false, linecolor = colors[3],linestyle = :dot, linewidth= 3)

    scatter!(data[index1, index2][2][1:end,4], data[index1, index2][2][1:end,1], label ="mobile", markercolor = colors[2], markershape=:star5)
    plot!(t, Wang.(t, k₋, k₊, γ)./γ, label=false, linecolor = colors[2], linestyle = :dash, linewidth= 3)
    # plot!(t, test2_0.(t, k₋, k₊, γ), label=false, linecolor = colors[2])
    scatter!(data[index1, index2][2][1:end,4], data[index1, index2][2][1:end,1] .+ data[index1, index2][2][1:end,2] .+ data[index1, index2][2][1:end,3], label = false, markershape=:dtriangle )
    plot!(t, Wang.(t, k₋, k₊, γ)./γ .+ ψ.(t,k₋, k₊) .+ ν.(t, k₋, k₊, γ), label="Sum", linecolor = colors[4])
    # plot!(t, test2_0.(t, k₋, k₊, γ) .+ ψ.(t,k₋, k₊) .+ jammed.(t, k₋, k₊, γ), label="Sum", linecolor = colors[4])
    # plot!(data[index1, index2][2][1:end,4], data[index1, index2][1], label ="current")
    # plot!(data[index1, index2][2][1:end,4], data[index1, index2][1], label ="current")
  end
  
  body!(w, ui)
end

function plot_J_vs_k₊(data)
  t = k₊_list[1]:0.01:k₊_list[end]
  w = Window()
  ρ_list = collect(range(0.04, stop = 0.96, length=19))

  ui = @manipulate for k₋ in k₋_list, ρ in ρ_list
    index1 = findfirst(x-> x==k₋, k₋_list)
    index2 = findfirst(x-> x==ρ, ρ_list)
    # index3 = findfist(x-> x== k₊, k₊_list)
    nice_ρ = round(ρ, digits=3)
    τ = 1/k₋
    p2 = plot(
      ylabel = "J", 
      xlabel="k₊", 
      legend=:topleft, 
      xlims=(k₊_list[1], k₊_list[end] *1.1), 
      title = "k₋ = $k₋, ρ = $nice_ρ", 
      xtickfont = 16, 
      ytickfont=16, 
      guidefont = 16,
      # ylims = (-0.04,1),
      # yticks = (0:0.2:1)
      )#, ylims = (0,1))
    # J_k₊ = [for item in ]
    J_k₊ = [data[index1, item][1][index2] for item in axes(data, 2)]
    scatter!(k₊_list, J_k₊, label = false, markercolor = colors[1])

    plot!(
      t, 
      Wang.(data[index1, index2][2][1:end,4][index2], t, τ), 
      label = false, linecolor = colors[1])
  end
  body!(w, ui)
end

function plot_interact_kymo(data, ρ_list, k₊_list)
  w = Window()
   
  ui = @manipulate for k₊ in k₊_list, ρ in ρ_list
    index1 = findfirst(x-> x==k₊, k₊_list)
    index2 = findfirst(x-> x==ρ, ρ_list)
    J = round(data[index2, index1][1], digits=6)
    # index3 = findfist(x-> x== k₊, k₊_list)

    plot_kymo(data[index2, index1][4][1:end,1:2:end]', data[index2, index1][5][1:end,1:2:end]', k₊, 1/10^-4, data[index2, index1][3][1:2:end])
    plot!(title = "J = $J", )
    
    # p2 = plot(
    #   ylabel = "J", 
    #   xlabel="k₊", 
    #   legend=:topleft, 
    #   xlims=(k₊_list[1], k₊_list[end] *1.1), 
    #   title = "k₋ = $k₋, ρ = $nice_ρ", 
    #   xtickfont = 16, 
    #   ytickfont=16, 
    #   guidefont = 16,
    #   # ylims = (-0.04,1),
    #   # yticks = (0:0.2:1)
    #   )#, ylims = (0,1))
    # J_k₊ = [for item in ]
  
  end
  body!(w, ui)
end

function plot_transient_density(data, time_vec, ρ)
  w = Window()
  time_vec_rounded = round.(time_vec, digits = 6)
  ui = @manipulate for time in time_vec_rounded
    index1 = findfirst(x-> x==time, time_vec_rounded)
    
    scatter(data[index1, 1:end], ylabel= "lattice site", xlabel = "density", title = "$time", label = false, ylims = (-0.,1.05), yticks=0:0.2:1)
    hline!([ρ], label = false, linewidth = 5)
  end
  body!(w, ui)
end

function plot_transient_current_profile(data, time_vec, ρ)
  w = Window()
  J = ρ * (1-ρ)
  time_vec_rounded = round.(time_vec, digits = 6)
  ui = @manipulate for time in time_vec_rounded
    index1 = findfirst(x-> x==time, time_vec_rounded)
    
    scatter(data[index1,1:end], ylabel= "Current", xlabel = "lattice", title = "$time", label = false, ylims = (-0.05,0.32))
    hline!([J], label = false, linewidth = 5)
  end
  body!(w, ui)
end

