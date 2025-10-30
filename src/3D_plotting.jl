#using Plots
using PlotlyJS
using Colors
using Blink: Window, body!
using JLD
Pkg.add("PlotlyJS")

include("analyitcal_functions.jl")


@load "pbc_plotting/data/c_data_250.jld"
@load "pbc_plotting/data/J_data_250.jld"
@load "pbc_plotting/data/ρ_data_250.jld"
@load "pbc_plotting/data/t_data_250.jld"

# load L250

#@load "pbc_plotting\\data2nd_run\\c_data_250.jld"
#@load "pbc_plotting\\data2nd_run\\J_data_250.jld"
#@load "pbc_plotting\\data2nd_run\\ρ_data_250.jld"
#@load "pbc_plotting\\data2nd_run\\t_data_250.jld"

c250, J250, ρ250, t250 = c, J[1:end-1, 1:end-1], ρ, t


#load L500
@load "pbc_plotting/data/c_data_500.jld"
@load "pbc_plotting/data/J_data_500.jld"
@load "pbc_plotting/data/ρ_data_500.jld"
@load "pbc_plotting/data/t_data_500.jld"

c500, J500, ρ500, t500 = c, J[1:end-1, 1:end-1], ρ, t

#load L750
@load "pbc_plotting/data/c_data_750.jld"
@load "pbc_plotting/data/J_data_750.jld"
@load "pbc_plotting/data/ρ_data_750.jld"
@load "pbc_plotting/data/t_data_750.jld"

c750, J750, ρ750, t750 = c, J[1:end-1, 1:end-1], ρ, t

function create_error_matrix(data,ρ_list, k₊_list, k₋_list, L)
   
    error_matrix_W = zeros(length(k₊_list), length(k₋_list))
    error_matrix_WW = zeros(length(k₊_list), length(k₋_list))
    error_matrix_G = zeros(length(k₊_list), length(k₋_list))
    error_matrix_L = zeros(length(k₊_list), length(k₋_list))
    xaxis_inp = ρ_list[1:end] 
    i, j = length(k₋_list), 1

    for k₋ in k₋_list
        j = 1
        index2 = findfirst(x-> x==k₋, k₋_list)
        
        for k₊ in k₊_list
            index1 = findfirst(x-> x==k₊, k₊_list)

            numerical = data[index2, index1][1:end,1]
            ana_W = Wang.(xaxis_inp, k₋, k₊, γ)
            ana_WW = J_whoosh.(γ, xaxis_inp, L, k₋, k₊)
            ana_G = greulich.(xaxis_inp, k₋, k₊,γ)
            ana_L = J_transient.(xaxis_inp,k₊,k₋,L)

            #rho_max = density_at_max(numerical, xaxis_inp)
            #ana_W = Wang(rho_max, k₋, k₊, γ)
            #ana_WW = J_whoosh(γ, rho_max, L, k₋, k₊)
            #ana_G = greulich(rho_max, k₋, k₊,γ)
            #ana_L = J_transient(rho_max,k₊,k₋,L)

            #sum_rel_err_W = round(sum(relative_error(numerical, ana_W)), digits = 2)
            #sum_rel_err_WW = round(sum(relative_error(numerical, ana_WW)), digits = 2)
            #sum_rel_err_G = round(sum(relative_error(numerical, ana_G)), digits = 2)
            #sum_rel_err_L = round(sum(relative_error(numerical, ana_L)), digits = 2)
            
            sum_rel_err_W = round(relative_max_error(numerical, ana_W), digits = 2)
            sum_rel_err_WW = round(relative_max_error(numerical, ana_WW), digits = 2)
            sum_rel_err_G = round(relative_max_error(numerical, ana_G), digits = 2)
            sum_rel_err_L = round(relative_max_error(numerical, ana_L), digits = 2)

            #sum_rel_err_W = round(absolute_bounded_peak_error(numerical, ana_W), digits = 2)
            #sum_rel_err_WW = round(absolute_bounded_peak_error(numerical, ana_WW), digits = 2)
            #sum_rel_err_G = round(absolute_bounded_peak_error(numerical, ana_G), digits = 2)
            #sum_rel_err_L = round(absolute_bounded_peak_error(numerical, ana_L), digits = 2)
            
            #sum_rel_err_W = round(positional_peak_error(numerical, ana_W,xaxis_inp), digits = 2)
            #sum_rel_err_WW = round(positional_peak_error(numerical, ana_WW,xaxis_inp), digits = 2)
            #sum_rel_err_G = round(positional_peak_error(numerical, ana_G,xaxis_inp), digits = 2)
            #sum_rel_err_L = round(positional_peak_error(numerical, ana_L,xaxis_inp), digits = 2)
          
            error_matrix_W[i,j] = sum_rel_err_W
            error_matrix_WW[i,j] = sum_rel_err_WW
            error_matrix_G[i,j] = sum_rel_err_G
            error_matrix_L[i,j] = sum_rel_err_L
            j += 1
        end
        i -= 1
    end
    return error_matrix_W, error_matrix_WW, error_matrix_G, error_matrix_L
end

k₋, k₊ = 10^-5, 10^-4#10^-2, 10^-3
index2 = findfirst(x-> x==k₋, k₋_list)
index1 = findfirst(x-> x==k₊, k₊_list)

A = [1 2 ; 3 5]
A[2,2]

newJ250 = flip_matrix(J250)

xaxis_inp = ρ_list[1:end]
numerical = J250[index2, index1][1:end,1]

Plots.scatter(xaxis_inp, numerical)

ana_W = Wang.(xaxis_inp, k₋, k₊, γ)
ana_WW = J_whoosh.(γ, xaxis_inp, L, k₋, k₊)
ana_G = greulich.(xaxis_inp, k₋, k₊,γ)
ana_L = J_transient.(xaxis_inp,k₊,k₋,L)

plot!(xaxis_inp, ana_WW)
plot!(xaxis_inp, ana_G)
plot!(xaxis_inp, ana_L)


relative_max_error(numerical, ana_L)
            
sum_rel_err_W = round(relative_max_error(numerical, ana_W), digits = 2)

sum_rel_err_WW = round(relative_max_error(numerical, ana_WW), digits = 2)

sum_rel_err_G = round(relative_max_error(numerical, ana_G), digits = 2)
sum_rel_err_L = round(relative_max_error(numerical, ana_L), digits = 2)


function flip_matrix(matrix, k₊_list, k₋_list)

    new_matrix = Matrix{Any}(undef, length(k₊_list), length(k₋_list))
    j = 0
    for i in 1:length(k₊_list)
        new_matrix[i,1:end] = matrix[end-j,1:end]
        j += 1
    end
    return new_matrix
end

function set_layouts(name_of_plot)
    layoutW = Layout(

        title = attr(
            text = "Error:"*name_of_plot,
            x= 0.5,
            y = 1.4,
            xanchor="center",
            yanchor = "top",
            font_size = 30
        ),
        xaxis = attr(
            title = "k₊ (log)",
            tickfont_size = 18,
            tick0 = -6,
            dtick = 1,
        ),

        yaxis = attr(
            title = "k₋ (log)",
            tickfont_size = 18,
            tick0 = -6,
            dtick = 1
        ),

        font = attr(
            family = "Cambria Math",
            #size = 22,
        )
    )
    layoutL = Layout(

        title = attr(
            text = "Error:"*name_of_plot,
            x= 0.5,
            y = 1.4,
            xanchor="center",
            yanchor = "top",
            font_size = 30
        ),
        xaxis = attr(
            title = "k₊ (log)",
            tickfont_size = 18,
            tick0 = -6,
            dtick = 1,
        ),

        yaxis = attr(
            titlefont_size = 22,
            tickfont_size = 18,
            tick0 = -6,
            dtick = 1
        ),

        font = attr(
            family = "Cambria Math",
            #size = 22,
        )
    )
    layoutWW = Layout(

        title = attr(
            text = "Error:"*name_of_plot,
            x= 0.5,
            y = 1.4,
            xanchor="center",
            yanchor = "top",
            font_size = 30
        ),
        xaxis = attr(
            tickfont_size = 18,
            tick0 = -6,
            dtick = 1,
        ),

        yaxis = attr(
            tickfont_size = 18,
            tick0 = -6,
            dtick = 1
        ),

        font = attr(
            family = "Cambria Math",
            #size = 22,
        )
    )

    layoutG = Layout(

        title = attr(
            text = "Error:"*name_of_plot,
            x= 0.5,
            y = 1.4,
            xanchor="center",
            yanchor = "top",
            font_size = 30
        ),
        xaxis = attr(
            tickfont_size = 18,
            tick0 = -6,
            dtick = 1,
        ),

        yaxis = attr(
            title = "k₋ (log)",
            tickfont_size = 18,
            tick0 = -6,
            dtick = 1
        ),

        font = attr(
            family = "Cambria Math",
            #size = 22,
        )
    )


    if name_of_plot == "Wang"
        return layoutW
    elseif name_of_plot == "Lorenzo"
        return layoutL
    elseif name_of_plot == "Whoosh and Wang"
        return layoutWW
    else
        return layoutG
    end


end

function plot_error_contour(err_matrix, name_of_plot,k₊_list, k₋_list)
    
    matrix = flip_matrix(err_matrix, k₊_list, k₋_list)
    
    x_axis = collect(range(start=-6, stop=0, length=length(k₋_list)))
    y_axis = collect(range(start=-6, stop=0, length=length(k₊_list)))

    
    trace1 =  PlotlyJS.contour(
        z = matrix,
        x = x_axis,
        y = y_axis,
        zauto = true,
        #zmin = 0,
        #zmax = 1,

        contours= attr(
            coloring = :heatmap,
            contours_start = 0.0,
            contours_size = 0.1,
            contours_end = 1.0,
            showlabels = true,
            labelfont = attr(
                size = 12,
                color = "black"
            )
        )
    )
    layout = Layout(

        title = attr(
            text = "Error:"*name_of_plot,
            x= 0.5,
            y = 1.4,
            xanchor="center",
            yanchor = "top",
            font_size = 30
        ),
        xaxis = attr(
            title = "k₊ (log)",
            tickfont_size = 18,
            titlefont_size = 20,
            tick0 = -6,
            dtick = 1,
        ),

        yaxis = attr(
            title = "k₋ (log)",
            tickfont_size = 18,
            titlefont_size = 20,
            tick0 = -6,
            dtick = 1
        ),

        font = attr(
            family = "Cambria Math",
            #size = 22,
        )
    )
    #=
    if name_of_plot == "Wang"
        layoutW = set_layouts(name_of_plot)
        p = PlotlyJS.plot(trace1, layoutW)
    elseif name_of_plot == "Lorenzo"
        layoutL = set_layouts(name_of_plot)
        p = PlotlyJS.plot(trace1, layoutL)
    elseif name_of_plot == "Whoosh and Wang"
        layoutWW = set_layouts(name_of_plot)
        p = PlotlyJS.plot(trace1, layoutWW)
    else
        layoutG = set_layouts(name_of_plot)
        p = PlotlyJS.plot(trace1, layoutG)
    end
    =#
    p = PlotlyJS.plot(trace1, layout)
    return p  
end

using PlotlyJS

function plot_error_isolines(err_matrices, labels, k₊_list, k₋_list; threshold=0.5)
    x_axis = collect(range(-6, 0; length=length(k₋_list)))
    y_axis = collect(range(-6, 0; length=length(k₊_list)))

    traces = AbstractTrace[]
    colors = ["red", "blue", "green", "purple", "orange", "black"]

    for (i, err_matrix) in enumerate(err_matrices)
        M = flip_matrix(err_matrix, k₊_list, k₋_list)

        push!(traces, PlotlyJS.contour(
            x = x_axis,
            y = y_axis,
            z = M,

            # turn off auto‐contours
            autocontour = false,

            # force exactly one contour at threshold
            ncontours = 1,
            showscale = false,

            contours = attr(
                contours_start     = threshold,   # only level = threshold
                contours_end       = threshold,
                size      = 1,           # irrelevant when start==end
                coloring  = "lines",
                showlabels= true
            ),

            line = attr(
                width = 3,
                color = colors[i]
            ),

            name = labels[i]
        ))
    end

    layout = Layout(
        title = "Error isolines at $(threshold*100)% error",
        xaxis = attr(title="log k₋"),
        yaxis = attr(title="log k₊"),
        legend = attr(x=1.02, y=1)
    )

    return PlotlyJS.plot(traces, layout)
end

function plot_single_contour(err_matrix, k₊_list, k₋_list; threshold=0.5)
    # Axes
    x_axis = collect(range(start=-6, stop=0, length=length(k₋_list)))
    y_axis = collect(range(start=-6, stop=0, length=length(k₊_list)))


    # Flip if needed so rows ↦ y, cols ↦ x
    M = flip_matrix(err_matrix, k₊_list, k₋_list)

    # Quick check
    @info "Data range: [$(minimum(M)), $(maximum(M))], threshold = $threshold"

    # One‐level contour trace
    trace = PlotlyJS.contour(
        x = x_axis,
        y = y_axis,
        z = M,
        showscale   = false,    # no colorbar
        autocontour = false,    # disable auto levels
        contours = attr(
            cvalues    = [threshold],  # <-- exactly this one level
            coloring  = "lines",      # draw only the line
            showlabels= true
        ),
        line = attr(
            width = 3,
            color = "red"
        )
    )

    layout = Layout(
        title = "Isoline at $(threshold*100)% error",
        xaxis = attr(title="log k₊"),
        yaxis = attr(title="log k₋")
    )

    return PlotlyJS.plot(trace, layout)
end

ρ_list = collect(range(0.01, stop= 0.99, length = 50))
k₊_list = [10^-6, 10^-5, 10^-4, 10^-3, 10^-2, 10^-1, 10^0]#, 10]
k₋_list = [10^-6, 10^-5, 10^-4, 10^-3, 10^-2, 10^-1, 10^0]#, 10]
γ = 1


# ask for user input
println("What system size do you want to look at?type(250,500,750): ")
#L = parse(Int64,readline())

newJ250 = flip_matrix(J250,k₊_list, k₋_list)

L = 250
errW250, errWW, errG, errL = create_error_matrix(J250,ρ_list, k₊_list, k₋_list, L)
errW250

traceW = plot_error_contour(errW250, "Wang",k₊_list, k₋_list)
plot_single_contour(errW250, k₊_list, k₋_list)


L = 500
errW500, errWW, errG, errL = create_error_matrix(J500,ρ_list, k₊_list, k₋_list, L)
traceW = plot_error_contour(errW500, "Wang",k₊_list, k₋_list)
plot_single_contour(errW500, k₊_list, k₋_list)

L = 750
errW750, errWW, errG, errL = create_error_matrix(J750,ρ_list, k₊_list, k₋_list, L)
traceW = plot_error_contour(errW750, "Wang",k₊_list, k₋_list)
plot_single_contour(errW750, k₊_list, k₋_list)

errW

plot_error_isolines(
    [errW250], 
    ["L = 250", "L = 500", "L = 750"], 
    k₊_list, 
    k₋_list
    )


PlotlyJS.savefig(traceWW, "Whoosh_and_Wang_error_contour_L_$L.svg")

errG
errL
errWW
errW

L = 250
if L == 250
    errW, errWW, errG, errL = create_error_matrix(J250,ρ_list, k₊_list, k₋_list, L )

    traceL = plot_error_contour(errL, "Lorenzo",k₊_list, k₋_list)
    traceW = plot_error_contour(errW, "Wang",k₊_list, k₋_list)
    traceWW = plot_error_contour(errWW, "Whoosh and Wang",k₊_list, k₋_list)
    traceG = plot_error_contour(errG, "Greulich",k₊_list, k₋_list)


    p = [
        traceWW traceL
        traceW traceG
    ]

    w = Window()
    body!(w,p)
elseif L == 500
    errW, errWW, errG, errL = create_error_matrix(J500,ρ_list, k₊_list, k₋_list, L )

    traceL = plot_error_contour(errL, "Lorenzo",k₊_list, k₋_list)
    traceW = plot_error_contour(errW, "Wang",k₊_list, k₋_list)
    traceWW = plot_error_contour(errWW, "Whoosh and Wang",k₊_list, k₋_list)
    traceG = plot_error_contour(errG, "Greulich",k₊_list, k₋_list)


        p = [
        traceWW traceL
        traceW traceG
    ]

    w = Window()
    body!(w,p)

elseif L == 750
    errW, errWW, errG, errL = create_error_matrix(J750,ρ_list, k₊_list, k₋_list, L )

    traceL = plot_error_contour(errL, "Lorenzo",k₊_list, k₋_list)
    traceW = plot_error_contour(errW, "Wang",k₊_list, k₋_list)
    traceWW = plot_error_contour(errWW, "Whoosh and Wang",k₊_list, k₋_list)
    traceG = plot_error_contour(errG, "Greulich",k₊_list, k₋_list)

            p = [
        traceWW traceL
        traceW traceG
    ]

    w = Window()
    body!(w,p)
end


errW, errWW, errG, errL = create_error_matrix(J250,ρ_list, k₊_list, k₋_list, L )

    traceL = plot_error_contour(errL, "Lorenzo",k₊_list, k₋_list)
    traceW = plot_error_contour(errW, "Wang",k₊_list, k₋_list)
    traceWW = plot_error_contour(errWW, "Whoosh and Wang",k₊_list, k₋_list)
    traceG = plot_error_contour(errG, "Greulich",k₊_list, k₋_list)


    p = [
        traceWW traceL
        traceW traceG
    ]

traceL = plot_error_contour(errL, "Lorenzo",k₊_list, k₋_list)

traceW = plot_error_contour(errW, "Wang",k₊_list, k₋_list)
traceWW = plot_error_contour(errWW, "Whoosh and Wang",k₊_list, k₋_list)
traceG = plot_error_contour(errG, "Greulich",k₊_list, k₋_list)

PlotlyJS.savefig(traceW, "pT_error_L_$L.svg")


