# Simplest setup ####

# Parameters
gravity_exponent <- -2
birth_rate <- 2
pest_area <- 1
human_area <- 1

param <- list(gravity_exponent=gravity_exponent, birth_rate=birth_rate, pest_area=pest_area, human_area=human_area)

# Environment
N <- 5
area <- const_area(N, 10)
pos <- unif_pos(N)
dist <- euclidean_dist(pos)
human <- const_human (N, 1000)
connectivity <- compute_connectivity (dist, human, gravity_exponent)
environment <- list(N=N, area=area, pos=pos, dist=dist, human=human, connectivity=connectivity)

# State
initial_state <- const_initial_state (N, tree_density=0.1, pest_density=1, area=area)
state <- array (dim=c(dim(initial_state), 1), dimnames=dimnames(initial_state))
state[,,1] <- as.matrix(initial_state)

# Simulation
sim <- list(param=param, environment=environment, state=state)

# Deterministic simulation ####
det_results <- step_engine(sim, eab_contagion_det, added_time=100)

# dynamic_graph(det_results, display="trees")
# dynamic_graph(det_results, display="pests")

state_ts_plot(det_results)
state_ts_plot(det_results, normalize=F)

# Stochastic simulation ####
sto_results <- step_engine(sim, eab_contagion_sto, added_time=100)


# dynamic_graph(sto_results, display="trees")
# dynamic_graph(sto_results, display="pests")

state_ts_plot(sto_results)
state_ts_plot(sto_results, normalize=F)
