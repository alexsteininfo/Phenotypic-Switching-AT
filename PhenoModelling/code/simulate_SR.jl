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
K = 1.5e11      # Carrying capacity
wSR = 0.0       # Switching rate from sensitive to resistant
wRS = 0.0       # Switching rate from resistant to sensitive

### Treatment implemention

# Additive death rate function based on drug dose
#   = 0.0 if we are interested in multiplicate death rate
#   = m_alpha * m for linear scaling
#   = ... for more Michaelis–Menten kinetics
dm(m) = m_alpha * m

# Multiplicative growth rate reduction function based on drug dose
#   = 1.0 if we are only interested in additive death rate
#   = ... for linear scaling
#   = ... for more Michaelis–Menten kinetics
gamma(m) = 1.0

### Treatment dependent growth and switching rates
rS_f(m) = 1.0       # currently not used
rR_f(m) = 1.0       # currently not used
wSR_f(m) = 1.0      # currently not used
wRS_f(m) = 1.0      # currently not used


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
wSR_list = range(0,1,11)
wRS_list = range(0,1,11)

# Initialize lists to save relevant simulation information
T_detection_list = Matrix{Float64}(undef, length(wSR_list), length(wRS_list))
S_detection_list = Matrix{Float64}(undef, length(wSR_list), length(wRS_list))
R_detection_list = Matrix{Float64}(undef, length(wSR_list), length(wRS_list))

# Run simulations over the parameter space
for (i,wSR) in enumerate(wSR_list)
    for (j,wRS) in enumerate(wRS_list)
        # Run simulations with corresponding parameterisation
        local param_local = (rS, rR, K, wSR, wRS, dm(m), gamma(m))
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
#   simulation stops at progression, virutal extinction or maximum treatment time

# Simulation parameters 
size_progression = 1.2*size_detection
tmin = 0
tmax = 200.0
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

# Run simulations over the parameter space
for (i,wSR) in enumerate(wSR_list)
    for (j,wRS) in enumerate(wRS_list)
        # Define initial conditions
        local u0_detect = [S_detection_list[i,j], R_detection_list[i,j]]
        #local u0_detect = [0.5*size_detection, 0.5*size_detection]
        
        # Run simulations with corresponding parameterisation
        local param_local = (rS, rR, K, wSR, wRS, dm(1.0), gamma(1.0)) # m = 1 by definition
        local T, S, R = maxtoldose(SR!, param_simulation_MTD, u0_detect, param_local)
        
        # Save relevant simulation data
        global drugusage_MTD[i,j] = 1.0 * T[end]
        global ttp_MTD[i,j] =    T[end]
        global S_MTD_list[i,j] = S[end]
        global R_MTD_list[i,j] = R[end]
    end
end


### Plotting
plot_trajectory(time_points_MTD, sensitive_pop_MTD, resistant_pop_MTD, "trajectory_MTD.png")
plot_heatmap(ttp_MTD, wSR_list, wRS_list, "ttp_MTD.png")
plot_heatmap(drugusage_MTD, wSR_list, wRS_list, "drugusage_MTD.png")

########################################
### Adaptive therapy (dose skipping) ###
########################################

### Adaptive therapy (dose skipping)
#   starting with initial conditions given at detection
#   apply treatment in cycles according to the Zhang schedule
#   simulation stops at progression, 'virutal' extinction or maximum treatment time

# Simulation parameters 
size_progression = 1.2*size_detection
tmin = 0.0
tmax = 200.0
lt = 0.5*size_progression   # lower threshold when treatment is turned off
ut = 1.0*size_progression   # upper threshold when treatment is turned on
param_simulation_AT = (size_progression, tmin, tmax, lt, ut)

# Run a single simulation with baseline parameters
u0 = [sensitive_pop_GR[end], resistant_pop_GR[end]]
param = (rS, rR, K, wSR, wRS, dm(1.0), gamma(1.0))   # m = 1.0 initially (but will be altered)
time_points_AT, sensitive_pop_AT, resistant_pop_AT = maxtoldose(SR!, param_simulation_AT, u0, param)

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
        local param_local = (rS, rR, K, wSR, wRS, dm(1.0), gamma(1.0)) # m = 1 by definition
        local T, S, R = maxtoldose(SR!, param_simulation_AT, u0_detect, param_local)
        
        # Save relevant simulation data
        global drugusage_AT[i,j] = 1.0 * T[end]
        global ttp_AT[i,j] =    T[end]
        global S_AT_list[i,j] = S[end]
        global R_AT_list[i,j] = R[end]
    end
end

### Plotting
plot_trajectory(time_points_AT, sensitive_pop_AT, resistant_pop_AT, "trajectory_AT.png")
plot_heatmap(ttp_AT, wSR_list, wRS_list, "ttp_AT.png")
plot_heatmap(drugusage_AT, wSR_list, wRS_list, "drugusage_AT.png")


