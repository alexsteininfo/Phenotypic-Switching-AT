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

#########################################
### Parameters and initial conditions ###
#########################################

r1 = 1.0        # Growth rate of sensitive cells
r2 = 1.2        # Growth rate of resistant cells
K = 1e11        # Carrying capacity
wSR = 0.1       # Switching rate from sensitive to resistant
wRS = 0.5       # Switching rate from resistant to sensitive
m = 0.0         # Drug dose (0 to infinity)
param = (r1, r2, K, m, wSR, wRS, dm(m), gamma(m))

# Initial conditions
S0 = 1
R0 = 0
u0 = [S0, R0]

#######################
### Run simulations ###
#######################

detection_size = 1e9
tmin = 0
tmax = 200

time_points, sensitive_pop, resistant_pop = tumourgrowth(SR!, detection_size, tmin, tmax, u0, param)














