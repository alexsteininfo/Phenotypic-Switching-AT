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

### Treatment implemention
# Additive death rate function based on drug dose
dm(m) = m
# Multiplicative growth rate reduction function based on drug dose
gamma(m) = 1.0

#################################
### Model baseline parameters ###
#################################

r1 = 1.0        # Growth rate of sensitive cells
r2 = 1.0        # Growth rate of resistant cells
K = 1e12        # Carrying capacity
wSR = 0.1       # Switching rate from sensitive to resistant
wRS = 0.5       # Switching rate from resistant to sensitive
m = 0.0         # Drug dose (0 to infinity)
param = (r1, r2, K, m, wSR, wRS, dm(m), gamma(m))


#######################################
### Run simulations of growth phase ###
#######################################

### Tumour growth
#   Starting from a single sensitive cell
#   growing to detction_size

# Initial conditions
S0 = 1
R0 = 0
u0 = [S0, R0]

# Simulation parameters 
size_detection = 1e10
tmin = 0
tmax = Inf
param_simulation_growth = (size_detection, tmin, tmax)

# Run single simulation with baseline parameters
param = (r1, r2, K, m, wSR, wRS, dm(m), gamma(m))
time_points_1, sensitive_pop_1, resistant_pop_1 = tumourgrowth(SR!, param_simulation_growth, u0, param)

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
        local param_local = (r1, r2, K, wSR, wRS, dm(m), gamma(m))
        local T, S, R = tumourgrowth(SR!, param_simulation_growth, u0, param_local)
        
        # Save relevant simulation data
        global T_detection_list[i,j] = T[end]
        global S_detection_list[i,j] = S[end]
        global R_detection_list[i,j] = R[end]
    end
end


##########################################
### Run simulations of treatment phase ###
##########################################


### Maximum tolerated dose
#   starting with initial conditions given at detection
#   simulation stops at progression or 'virutal' extinction

# Simulation parameters 
size_progression = 1.2*size_detection
tmin = 0
tmax = 200
param_simulation_MTD = (size_progression, tmin, tmax)

# Run a single simulation with baseline parameters
u0 = [sensitive_pop_1[end], resistant_pop_1[end]]
param = (r1, r2, K, wSR, wRS, dm(1.0), gamma(1.0))   # m = 1.0
time_points_2, sensitive_pop_2, resistant_pop_2 = maxtoldose(SR!, param_simulation_MTD, u0, param)

# Initialize lists to save relevant simulation information
cum_drugusage = Matrix{Float64}(undef, length(wSR_list), length(wRS_list))
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
        local param_local = (r1, r2, K, wSR, wRS, dm(1.0), gamma(1.0))
        local T, S, R = maxtoldose(SR!, param_simulation_MTD, u0_detect, param_local)
        
        # Save relevant simulation data
        global cum_drugusage[i,j] = m * T[end]
        global ttp_MTD[i,j] =    T[end]
        global S_MTD_list[i,j] = S[end]
        global R_MTD_list[i,j] = R[end]
    end
end

### Adaptive therapy (dose skipping)
#   starting with initial conditions given at detection
#   simulation stops at progression or 'virutal' extinction








