# title: Implementation of treatment schedules
# author: alexander stein

#####################
### Tumour growth ###
#####################

# The population starts growing from a single sensitive cell
#   until it reaches the detection size
#   no treatment is applied (m=0)

function tumourgrowth(ODEsystem!, param_simulation, u0, param_model; solver=Tsit5(), kwargs...)
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
    sol = solve(prob, solver, callback = cb_detect)

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
#       (ii)   size = size_progression --> progression

function maxtoldose(ODEsystem!, param_simulation, u0, param_model; solver=Tsit5(), kwargs...)
    # Unpack simulation parameters
    size_progression, tmin, tmax = param_simulation
    
    # Progression threshold
    cond_prog(u,t,integrator) = sum(u) - size_progression   # activated when function crosses 0
    affect_prog!(integrator) = terminate!(integrator)
    cb_prog = ContinuousCallback(cond_prog, affect_prog!)

    cbs = CallbackSet(cb_prog)

    # Time span for the simulation
    tspan = (tmin, tmax)
    # Define the ODE problem
    prob = ODEProblem(ODEsystem!, u0, tspan, param_model)
    # Solve the ODE problem
    sol = solve(prob, solver, callback = cbs)

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
#       (ii)   size = size_progression --> progression

function at_doseskipping(ODEsystem!, param_simulation, u0, param_model; solver=Tsit5(), kwargs...)
    # Unpack simulation parameters
    size_progression, tmin, tmax, lt, ut = param_simulation

    # Progression threshold
    cond_prog(u,t,integrator) = sum(u) - size_progression   # activated when function crosses 0
    affect_prog!(integrator) = terminate!(integrator)
    cb_prog = ContinuousCallback(cond_prog, affect_prog!)
    
    # Upper bound threshold
    cond_upper(u,t,integrator) = sum(u) - ut
    affect_upper!(integrator) = begin
        p = integrator.p
        integrator.p = (p[1], p[2], p[3], p[4], p[5], dm(1.0), gamma(1.0))
    end
    cb_upper = ContinuousCallback(cond_upper,affect_upper!)

    # Lower bound threshold
    cond_lower(u,t,integrator) = sum(u) - lt
    affect_lower!(integrator) = begin
        p = integrator.p
        integrator.p = (p[1], p[2], p[3], p[4], p[5], dm(0.0), gamma(0.0))
    end
    cb_lower = ContinuousCallback(cond_lower,affect_lower!)

    cbs = CallbackSet(cb_prog, cb_upper, cb_lower)

    # Time span for the simulation
    tspan = (tmin, tmax)
    # Define the ODE problem
    prob = ODEProblem(ODEsystem!, u0, tspan, param_model)
    # Solve the ODE problem
    sol = solve(prob, solver, callback = cbs)

    # Extract time points and populations
    if(length(u0)==2)
        return sol.t, sol[1, :], sol[2, :] 
    else
        return sol.t, sol[1, :], sol[2, :], sol[3, :]
    end
end

function at_doseskipping2(ODEsystem!, param_simulation, u0, param_model)
    # Unpack simulation parameters
    size_progression, tmin, tmax, lt, ut = param_simulation

    # Variables for treatment implementation
    t = tmin
    u = copy(u0)

    # Vectors for time and population sizes
    ts = Float64[]
    us = Vector{typeof(u)}()

    next_p = Ref(p)
    reason = Ref(:none)
    size_progression, tmin, tmax, lt, ut = param_simulation
    while t < tmax
            

    end

    return sol.t, sol[1, :], sol[2, :], t_cycles, n_cycles
end


### Another code version suggested by chatGPT
#=
function at_doseskipping(ODEsystem!, param_simulation, u0, p0; solver=Tsit5(), kwargs...)
    size_progression, tmin, tmax, lt, ut = param_simulation

    t = tmin
    u = copy(u0)
    p = p0

    ts = Float64[]
    us = Vector{typeof(u)}()

    next_p = Ref(p)
    reason = Ref(:none)

    while t < tmax
        reason[] = :none
        next_p[] = p

        cond_prog(u,t,_)  = sum(u) - size_progression
        affect_prog!(int) = (reason[] = :prog; terminate!(int))

        cond_upper(u,t,_)  = sum(u) - ut
        affect_upper!(int) = begin
            p = int.p
            next_p[] = (p[1], p[2], p[3], p[4], p[5], dm(1.0), gamma(1.0))
            reason[] = :upper
            terminate!(int)
        end

        cond_lower(u,t,_)  = sum(u) - lt
        affect_lower!(int) = begin
            p = int.p
            next_p[] = (p[1], p[2], p[3], p[4], p[5], dm(0.0), gamma(0.0))
            reason[] = :lower
            terminate!(int)
        end

        cbs = CallbackSet(
            ContinuousCallback(cond_prog,  affect_prog!),
            ContinuousCallback(cond_upper, affect_upper!),
            ContinuousCallback(cond_lower, affect_lower!)
        )

        sol = solve(ODEProblem(ODEsystem!, u, (t, tmax), p), solver; callback=cbs, kwargs...)

        # append segment (skip first point if continuing)
        if !isempty(ts)
            append!(ts, sol.t[2:end])
            append!(us, sol.u[2:end])
        else
            append!(ts, sol.t)
            append!(us, sol.u)
        end

        t = sol.t[end]
        u = sol.u[end]
        p = next_p[]

        reason[] == :prog && break
        t = nextfloat(t)  # move past the event time to avoid immediately retriggering
    end

    U = reduce(hcat, us)  # states × time
    return ts, eachrow(U)...
end
=#

##########################################
### Adaptive therapy (dose modulation) ###
##########################################

function at_dosemodulation()
    return 0
end




