# title: implementation of treatment schedules
# author: alexander stein

#####################
### Tumour growth ###
#####################

function tumourgrowth(ODEsystem!, detection_size, tmin, tmax, u0, param)
    
    # Define stopping size criterion
    cond_up(u,t,integrator) = u[1] + u[2] - detection_size
    affect_up!(integrator) = terminate!(integrator)
    cb_detect = ContinuousCallback(cond_up, affect_up!)

    # Time span for the simulation
    tspan = (tmin, tmax)
    # Define the ODE problem
    prob = ODEProblem(ODEsystem!, u0, tspan, param)
    # Solve the ODE problem
    sol = solve(prob, Tsit5(), callback = cb_detect)

    # Extract time points and populations
    time_points = sol.t
    sensitive_pop = sol[1, :]
    resistant_pop = sol[2, :]  
    total_pop = sensitive_pop .+ resistant_pop

    return time_points, sensitive_pop, resistant_pop
end

############################## 
### Maximum tolerated dose ###
##############################

function maxtoldose()
    # Implement simulation stopping criteria based on population size

    # Upper bound threshold
    threshold_up = 0.8 * K
    cond_up(u,t,integrator) = u[1] + u[2] - threshold_up
    affect_up!(integrator) = terminate!(integrator)
    cb_up = ContinuousCallback(cond_up, affect_up!)

    # Lower bound threshold
    threshold_low = 0.4 * K
    cond_low(u,t,integrator) = u[1] + u[2] - threshold_low
    affect_low!(integrator) = terminate!(integrator)
    cb_low = ContinuousCallback(cond_low, affect_low!)

    cbs = CallbackSet(cb_up, cb_low)

    # Extract time points and populations
    time_points = sol.t
    sensitive_pop = sol[1, :]
    resistant_pop = sol[2, :]  
    total_pop = sensitive_pop .+ resistant_pop

    return time_points, sensitive_pop, resistant_pop

end

########################################
### Adaptive therapy (dose skipping) ###
########################################

function at_dskipping()
    return 0
end

##########################################
### Adaptive therapy (dose modulation) ###
##########################################

function at_dmodulation()
    return 0
end




