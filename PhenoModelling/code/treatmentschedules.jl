# title: Implementation of treatment schedules
# author: alexander stein

#####################
### Tumour growth ###
#####################

function tumourgrowth(ODEsystem!, param_simulation, u0, param_model)
    # Unpack simulation parameters
    size_detection, tmin, tmax = param_simulation
    
    # Define stopping size criterion
    cond_up(u,t,integrator) = sum(u) - size_detection
    affect_up!(integrator) = terminate!(integrator)
    cb_detect = ContinuousCallback(cond_up, affect_up!)

    # Time span for the simulation
    tspan = (tmin, tmax)
    # Define the ODE problem
    prob = ODEProblem(ODEsystem!, u0, tspan, param_model)
    # Solve the ODE problem
    sol = solve(prob, Tsit5(), callback = cb_detect)

    # Extract time points and populations
    if(length(u0)==2)
        return sol.t, sol[1, :], sol[2, :] 
    else
        return sol.t, sol[1, :], sol[2, :], sol[3, :]
    end
end

################# 
### Treatment ###
#################

### Maximum tolerated dose
#   terminates when 
#       (i)     t = tmax --> maximum treatment time "tmax-tmin" is reached
#       (ii)    size = 1 --> virtual extinction
#       (iii)   size = size_progression --> progression

function maxtoldose(ODEsystem!, param_simulation, u0, param_model)
    # Unpack simulation parameters
    size_progression, tmin, tmax = param_simulation
    
    # Progression threshold
    cond_prog(u,t,integrator) = sum(u) - size_progression   # activated when function crosses 0
    affect_prog!(integrator) = terminate!(integrator)
    cb_prog = ContinuousCallback(cond_prog, affect_prog!)

    # Extinction threshold
    cond_extinct(u,t,integrator) = sum(u) - 1
    affect_extinct!(integrator) = terminate!(integrator)
    cb_extinct = ContinuousCallback(cond_extinct, affect_extinct!)
    
    cbs = CallbackSet(cb_prog, cb_extinct)

    # Time span for the simulation
    tspan = (tmin, tmax)
    # Define the ODE problem
    prob = ODEProblem(ODEsystem!, u0, tspan, param_model)
    # Solve the ODE problem
    sol = solve(prob, Tsit5(), callback = cbs)

    # Extract time points and populations
    if(length(u0)==2)
        return sol.t, sol[1, :], sol[2, :] 
    else
        return sol.t, sol[1, :], sol[2, :], sol[3, :]
    end
end

### Adaptive therapy (dose skipping)
#   starting with treatment m = 1
#   changes treatment when
#       size = lower_threshold: m --> 0
#       size = upper_threshold: m --> 1
#   terminates when 
#       (i)     t = tmax --> maximum treatment time "tmax-tmin" is reached
#       (ii)    size = 1 --> virtual extinction
#       (iii)   size = size_progression --> progression

function at_dskipping(ODEsystem!, param_simulation, u0, param_model)
    # Unpack simulation parameters
    size_progression, tmin, tmax, lt, ut = param_simulation

    # Progression threshold
    cond_prog(u,t,integrator) = sum(u) - size_progression   # activated when function crosses 0
    affect_prog!(integrator) = terminate!(integrator)
    cb_prog = ContinuousCallback(cond_prog, affect_prog!)

    # Extinction threshold
    cond_extinct(u,t,integrator) = sum(u) - 1
    affect_extinct!(integrator) = terminate!(integrator)
    cb_extinct = ContinuousCallback(cond_extinct, affect_extinct!)
    
    # Upper bound threshold
    cond_upper(u,t,integrator) = sum(u) - ut
    affect_upper!(integrator) = begin
        println("Changing parameter at t=$(integrator.t)")
        integrator.param[6] = dm(1.0)         # modify treatment parameter, dm
        integrator.param[7] = gamma(1.0)      # modify treatment parameter, gamma
    end
    cb_upper = ContinuousCallback(cond_upper,affect_upper!)

    # Lower bound threshold
    cond_lower(u,t,integrator) = sum(u) - lt
    affect_lower!(integrator) = begin
        println("Changing parameter at t=$(integrator.t)")
        integrator.param[6] = dm(0.0)   # modify treatment parameter, dm
        integrator.param[7] = gamma(0.0)   # modify treatment parameter, gamma
    end
    cb_lower = ContinuousCallback(cond_lower,affect_lower!)

    cbs = CallbackSet(cb_prog, cb_extinct, cb_upper, cb_lower)

    # Time span for the simulation
    tspan = (tmin, tmax)
    # Define the ODE problem
    prob = ODEProblem(ODEsystem!, u0, tspan, param_model)
    # Solve the ODE problem
    sol = solve(prob, Tsit5(), callback = cbs)

    # Extract time points and populations
    if(length(u0)==2)
        return sol.t, sol[1, :], sol[2, :] 
    else
        return sol.t, sol[1, :], sol[2, :], sol[3, :]
    end
end

##########################################
### Adaptive therapy (dose modulation) ###
##########################################

function at_dmodulation()
    return 0
end




