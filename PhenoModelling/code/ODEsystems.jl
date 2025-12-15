# title: Implementation of ODE systems
# author: alexander stein


### 2-type Lotka-Volterra model with switching and treatment
#   one sensitive (S) and one resistant population (R)
function SR!(du, u, param, t)
    # Unpack parameters and variables
    r1, r2, K, wSR, wRS, dm, gamma = param
    S, R = u

    # Define the ODEs
    dS = gamma * r1 * S * (1 - (S + R) / K) - dm * S - wSR * S + wRS * R
    dR =         r2 * R * (1 - (S + R) / K)          + wSR * S - wRS * R

    du[1] = dS
    du[2] = dR
end


### 3-type Lotka-Volterra model with switching and treatment
#   two sensitive (S1, S2) and one resistant population (R)
#   where only S2 and R have phenotypic switching
function S1S2R!(du, u, param, t)
    # Unpack parameters and variables
    rS1, rS2, rR, K, wSR, wRS, dm, gamma = param
    S1, S2, R = u

    # Define the ODEs
    dS1 = gamma * r1 * S1 * (1 - (S1 + S2 + R) / K) - dm * S - wSR * S + wRS * R
    dR =         r2 * R * (1 - (S + R) / K)          + wSR * S - wRS * R

    du[1] = dS
    du[2] = dR
    du[3] = dR
end


### 3-type Lotka-Volterra model with switching and treatment
#   one sensitive (S) and two resistant population (R1,R2)
#   where R2 is an irreversible resistant state
function SR1R2!(du, u, param, t)
    # Unpack parameters and variables
    r1, r2, K, wSR, wRS, dm, gamma = param
    S, R = u

    # Define the ODEs
    dS = gamma * r1 * S * (1 - (S + R) / K) - dm * S - wSR * S + wRS * R
    dR =         r2 * R * (1 - (S + R) / K)          + wSR * S - wRS * R

    du[1] = dS
    du[2] = dR
    du[3] = dR
end

