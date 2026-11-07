# # title: Making nice plots
# author: alexander stein


# Plot the growth size trajectory for the two-type model: S and R
# Input: time_points, populatoin size of sensitive cells and resistance cells
#       and finally also the name of the plot
function plot_trajectory(time_points, sensitive_pop, resistant_pop, name)
    f = Figure(size = (600,600), fontsize = 24)
    ax = Axis(f[1,1], xlabel="Time", ylabel="Number of cells")
    lines!(ax, time_points, sensitive_pop, color=:red, linewidth=5, label="S")
    lines!(ax, time_points, resistant_pop, color=:blue, linewidth=5, label="R")
    lines!(ax, time_points, sensitive_pop.+resistant_pop, color=:black, linewidth=5, label="N")
    axislegend(position = :rc)
    save("./PhenoModelling/figures/"*name, f)
end

# Plot heatmap of an arbitray matrix
# Input: time_points, populatoin size of sensitive cells and resistance cells
#       and finally also the name of the plot
function plot_heatmap(matrix, xaxis, yaxis, name)
    f = Figure(size = (600,600), fontsize = 24)
    ax = Axis(f[1,1][1,1], xlabel = "x", ylabel = "y")
    hm = heatmap!(ax, xaxis, yaxis, matrix)
    Colorbar(f[1,1][1,2], hm)
    save("./PhenoModelling/figures/"*name, f)
end

# Plot time to progression including theoretical predictions
function plot_heatmap_theory(matrix, wSR_list, wRS_list, wRS_extinction_list, wRS_stabilisation_list, name)
    f = Figure(size = (600,600), fontsize = 24)
    ax = Axis(f[1,1][1,1], xlabel = "w_SR", ylabel = "w_RS")
    hm = heatmap!(ax, wSR_list, wRS_list, matrix)
    Colorbar(f[1,1][1,2], hm)
    lines!(ax, wSR_list, wRS_extinction_list, color=:red, linewidth=5)
    lines!(ax, wSR_list, wRS_stabilisation_list, color=:red, linewidth=5)
    xlims!(ax, minimum(wSR_list), maximum(wSR_list))
    ylims!(ax, minimum(wRS_list), maximum(wRS_list))
    save("./PhenoModelling/figures/"*name, f)
end



function plot_heatmap_categorial(mat, xaxis, yaxis, threshold_values, name)

    # build category matrix
    cat = zeros(Int, size(mat))
    for thr in threshold_values
        cat .+= mat .> thr
    end

    f = Figure()
    ax = Axis(f[1, 1], xlabel = "x", ylabel = "y")

    hm = heatmap!(
        ax, xaxis, yaxis, cat;
        colormap = cgrad(:viridis, 4; categorical=true),
        colorrange = (0, maximum(cat))
    )

    Colorbar(f[1, 2], hm, label = "Category")
    save("./PhenoModelling/figures/"*name, f)
end



