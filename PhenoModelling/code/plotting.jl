# # title: Making nice plots
# author: alexander stein


# Plot the growth size trajectory
# Input: time_points, populatoin size of sensitive cells and resistance cells
#       and finally also the name of the plot
function plot_trajectory(time_points, sensitive_pop, resistant_pop, name)
    f = Figure(size = (600,650))
    ax = Axis(f[1, 1])
    lines!(ax, time_points, sensitive_pop, label="S")
    lines!(ax, time_points, resistant_pop, label="R")
    lines!(ax, time_points, sensitive_pop.+resistant_pop, label="N")
    save("./PhenoModelling/figures/"*name, f)
end


function plot_heatmap(matrix, xaxis, yaxis, name)
    f = Figure(size = (600,800))
    ax, hm = heatmap(f[1,1][1,1], xaxis, yaxis, matrix)
    Colorbar(f[1,1][1,2], hm)
    save("./PhenoModelling/figures/"*name, f)
end

