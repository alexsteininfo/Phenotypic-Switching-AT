# title: Solving coupled ODEs that simulate phenotypic switching during adaptive therapy
# author: alexander stein

#######################
### Import Packages ###
#######################

using DifferentialEquations
using CairoMakie

### Include functions that implement ODE systems
include("./ODEsystems.jl")

### Include functions that implement treatment schedules
include("./treatmentschedules.jl")

### Include functions to make plots 
include("./plotting.jl")

#########################
### Model  parameters ###
#########################

### Baseline parameters
rS = 1.0        # Growth rate of sensitive cells
rR = 1.0        # Growth rate of resistant cells
K = 4.0e11      # Carrying capacity
wSR = 0.5       # Switching rate from sensitive to resistant
wRS = 0.3       # Switching rate from resistant to sensitive

### Treatment implemention

# Additive death rate function based on drug dose
#   = 0.0 if we are interested in multiplicate death rate
#   = d_max * m for linear scaling
#   = ... for more Michaelis–Menten kinetics
d_max = 2.0
dm(m) = d_max * m

# Multiplicative growth rate reduction function based on drug dose
#   = 1.0 if we are only interested in additive death rate
#   = ... for linear scaling
#   = ... for more Michaelis–Menten kinetics
gamma(m) = 1.0

### Treatment dependent growth and switching rates
rS_f(m) = 1.0       # currently not used
rR_f(m) = 1.0       # currently not used
wSR_f(m) = 2.0      # currently not used
wRS_f(m) = 2.0      # currently not used


#######################################
### Run simulations of growth phase ###
#######################################

### Tumour growth
#   Starting from a single sensitive cell
#   growing to size = size_detection

# Initial conditions
S0_GR = 1
R0_GR = 0
u0_GR = [S0_GR, R0_GR]

# Simulation parameters 
size_detection = 1e11
tmin = 0
tmax = Inf
param_simulation_growth = (size_detection, tmin, tmax)

# Run single simulation with baseline parameters
param_GR = (rS, rR, K, wSR, wRS, dm(0.0), gamma(0.0))   # m = 0.0 during growth
time_points_GR, sensitive_pop_GR, resistant_pop_GR = tumourgrowth(SR!, param_simulation_growth, u0_GR, param_GR)

# Define variable parameter values
wSR_list = range(0,2,51)
wRS_list = range(0,2,51)

# Initialize lists to save relevant simulation information
T_detection_list = Matrix{Float64}(undef, length(wSR_list), length(wRS_list))
S_detection_list = Matrix{Float64}(undef, length(wSR_list), length(wRS_list))
R_detection_list = Matrix{Float64}(undef, length(wSR_list), length(wRS_list))

# Run simulations over the parameter space
for (i,wSR) in enumerate(wSR_list)
    for (j,wRS) in enumerate(wRS_list)
        # Run simulations with corresponding parameterisation
        local param_local = (rS, rR, K, wSR, wRS, dm(0.0), gamma(0.0))  # m = 0.0 during growth
        local T, S, R = tumourgrowth(SR!, param_simulation_growth, u0_GR, param_local)
        
        # Save relevant simulation data
        global T_detection_list[i,j] = T[end]
        global S_detection_list[i,j] = S[end]
        global R_detection_list[i,j] = R[end]
    end
end

### Plotting
plot_trajectory(time_points_GR, sensitive_pop_GR, resistant_pop_GR, "trajectory_GR.png")
plot_heatmap(R_detection_list/size_detection, wSR_list, wRS_list, "prop_R_detection.png")

##############################
### Maximum tolerated dose ###
##############################

### Maximum tolerated dose
#   starting with initial conditions given at detection
#   apply treatment at maximum tolerated dose
#   simulation stops at progression or maximum treatment time

# Simulation parameters 
size_progression = 1.2*size_detection
tmin = 0
tmax = 600.0
param_simulation_MTD = (size_progression, tmin, tmax)

# Run a single simulation with baseline parameters
u0 = [sensitive_pop_GR[end], resistant_pop_GR[end]]
param_MTD = (rS, rR, K, wSR, wRS, dm(1.0), gamma(1.0))   # m = 1.0 during MTD
time_points_MTD, sensitive_pop_MTD, resistant_pop_MTD = maxtoldose(SR!, param_simulation_MTD, u0, param_MTD)

# Initialize lists to save relevant simulation information
drugusage_MTD = Matrix{Float64}(undef, length(wSR_list), length(wRS_list))
ttp_MTD =       Matrix{Float64}(undef, length(wSR_list), length(wRS_list))
S_MTD_list =    Matrix{Float64}(undef, length(wSR_list), length(wRS_list))
R_MTD_list =    Matrix{Float64}(undef, length(wSR_list), length(wRS_list))
FS_MTD_list =   Matrix{Float64}(undef, length(wSR_list), length(wRS_list))

# Run simulations over the parameter space
for (i,wSR) in enumerate(wSR_list)
    for (j,wRS) in enumerate(wRS_list)
        # Define initial conditions
        local u0_detect = [S_detection_list[i,j], R_detection_list[i,j]]
        #local u0_detect = [0.5*size_detection, 0.5*size_detection]
        
        # Run simulations with corresponding parameterisation
        local param_local = (rS, rR, K, wSR, wRS, dm(1.0), gamma(1.0)) # m = 1 during MTD
        local T, S, R = maxtoldose(SR!, param_simulation_MTD, u0_detect, param_local)
        
        # Save relevant simulation data
        global drugusage_MTD[i,j] = 1.0 * T[end]
        global ttp_MTD[i,j] =    T[end]
        global S_MTD_list[i,j] = S[end]
        global R_MTD_list[i,j] = R[end]
        global FS_MTD_list[i,j] = S[end] + R[end]
    end
end


### Plotting
# Trajectory of a single realisation
plot_trajectory(time_points_MTD, sensitive_pop_MTD, resistant_pop_MTD, "trajectory_MTD.png")

# Time to progression over parameter space
ttp_MTD_NaN = copy(ttp_MTD)
ttp_MTD_NaN[ttp_MTD_NaN .== tmax] .= NaN
plot_heatmap(ttp_MTD_NaN, wSR_list, wRS_list, "ttp_MTD.png")

# Cumulative drug usage over parameter space
drugusage_MTD_NaN = copy(drugusage_MTD)
drugusage_MTD_NaN[drugusage_MTD_NaN .== tmax] .= NaN
plot_heatmap(drugusage_MTD, wSR_list, wRS_list, "drugusage_MTD.png")

# Whether the population goes extinct, stabilizes or progresses
plot_heatmap_categorial(FS_MTD_list, wSR_list, wRS_list, [0.5*size_detection, 1.0*size_detection, 1.2*size_detection], "finalsizes_MTD.png")

########################################
### Adaptive therapy (dose skipping) ###
########################################

### Adaptive therapy (dose skipping)
#   starting with initial conditions given at detection
#   apply treatment in cycles according to the Zhang schedule
#   simulation stops at progression or maximum treatment time

# Simulation parameters 
#size_progression = 1.2*size_detection  # keep same as for MTD
#tmin = 0.0                             # same as for MTD
#tmax = 200.0                           # same as for MTD
lt = 0.5*size_detection   # lower threshold when treatment is turned off
ut = 1.0*size_detection   # upper threshold when treatment is turned on
param_simulation_AT = (size_progression, tmin, tmax, lt, ut)

# Run a single simulation with baseline parameters
u0 = [sensitive_pop_GR[end], resistant_pop_GR[end]]
#u0 = [0.99,0.01]*size_detection
param = (rS, rR, K, wSR, wRS, dm(1.0), gamma(1.0))   # m = 1.0 initially (but will be altered)
time_points_AT, sensitive_pop_AT, resistant_pop_AT = at_doseskipping(SR!, param_simulation_AT, u0, param)

# Initialize lists to save relevant simulation information
drugusage_AT = Matrix{Float64}(undef, length(wSR_list), length(wRS_list))
ttp_AT =       Matrix{Float64}(undef, length(wSR_list), length(wRS_list))
S_AT_list =    Matrix{Float64}(undef, length(wSR_list), length(wRS_list))
R_AT_list =    Matrix{Float64}(undef, length(wSR_list), length(wRS_list))

# Run simulations over the parameter space
for (i,wSR) in enumerate(wSR_list)
    for (j,wRS) in enumerate(wRS_list)
        # Define initial conditions
        local u0_detect = [S_detection_list[i,j], R_detection_list[i,j]]
        #local u0_detect = [0.5*size_detection, 0.5*size_detection]
        
        # Run simulations with corresponding parameterisation
        local param_local = (rS, rR, K, wSR, wRS, dm(1.0), gamma(1.0)) # m = 1.0 initially (but will be altered)
        local T, S, R = at_doseskipping(SR!, param_simulation_AT, u0_detect, param_local)
        
        # Save relevant simulation data
        global drugusage_AT[i,j] = 1.0 * T[end]
        global ttp_AT[i,j] =    T[end]
        global S_AT_list[i,j] = S[end]
        global R_AT_list[i,j] = R[end]
    end
end

### Plotting
plot_trajectory(time_points_AT, sensitive_pop_AT, resistant_pop_AT, "trajectory_AT.png")

ttp_AT_NaN = copy(ttp_AT)
ttp_AT_NaN[ttp_AT_NaN .== tmax] .= NaN
plot_heatmap(ttp_AT_NaN, wSR_list, wRS_list, "ttp_AT.png")

drugusage_AT_NaN = copy(drugusage_AT)
drugusage_AT_NaN[drugusage_AT_NaN .== tmax] .= NaN
plot_heatmap(drugusage_AT_NaN, wSR_list, wRS_list, "drugusage_AT.png")

#############################
### Higher level analysis ###
#############################

ttp_diff = ttp_AT - ttp_MTD
#ttp_diff[abs.(ttp_diff) .< 1e-2] .= 0.0

ttp_diff_NaN = copy(ttp_diff)
ttp_diff_NaN[ttp_diff_NaN .== tmax] .= NaN

#ttp_diff_NaN = ttp_diff_NaN[2:end,1:end]  # Avoid the case with no resistant cells
plot_heatmap(ttp_diff_NaN, wSR_list, wRS_list, "ttp_difference.png")

# Theoretical extinction threshold for wRS
wRS_extinction(wSR) = ( (rS-d_max)*rR - (wSR*rR) ) / (rS-d_max)
wRS_extinction_list = [wRS_extinction(wSR) for wSR in wSR_list]

# Theoretical stabilisation threshold for wRS
Nstar = 1.2*size_detection
wRS_stabilisation(wSR) = (K - Nstar) * rR*(1/K + wSR /(d_max*K - K*rS + Nstar*rS));
wRS_stabilisation_list = [wRS_stabilisation(wSR) for wSR in wSR_list]

# Plot difference in time to progression
plot_heatmap_theory(ttp_diff_NaN, wSR_list, wRS_list, wRS_extinction_list, wRS_stabilisation_list, "ttp_difference.png")

# Plot difference in drug usage
drugusage_diff = drugusage_AT - drugusage_MTD
drugusage_diff_NaN = copy(drugusage_diff)
drugusage_diff_NaN[drugusage_diff_NaN .== tmax] .= NaN
plot_heatmap_theory(drugusage_diff_NaN, wSR_list, wRS_list, wRS_extinction_list, wRS_stabilisation_list, "drugusage_difference.png")



