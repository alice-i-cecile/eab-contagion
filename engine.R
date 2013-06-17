# Engines ####

# A simple engine for discrete models
step_engine <- function(sim, update_func, added_time=50)
{
  # Expand the state array to accomadate new data
  initial_time <- dim(sim$state)[3]
  
  empty_state <- array(NA, c(dim(sim$state)[1], dim(sim$state)[2], initial_time + added_time), dimnames(sim$state))
  empty_state[,,1:initial_time] <- sim$state
  
  sim$state <- empty_state 
  
  # Run the simulation
  for (time in (initial_time + 1):(initial_time+added_time))
  {
    sim$state[,,time] <- update_func(sim)
  }
  
  return (sim)
  
}

# Deterministic simulation ####

# Deterministic modelling of emerald ash borer spread as a contagion
eab_contagion_det <- function(sim)
{
  # Resume simulation from the last time step
  last_time <- max(which(!is.na(sim$state), arr.ind=T)[,3])
  last_state <- sim$state[,,last_time]
  
  # Extract number of patches
  N <- sim$environment$N

  
  # Calculate transport
  transport_det <- function(i)
  {
    # Check if human transport is in effect
    if (sum(sim$environment$human)==0)
    {
      return(rep(0, times=N))
    }
    
    # Extract relevant parameters
    human_area <- sim$param$human_area
    patch_area <- sim$environment$area[i]
    human_i <- sim$environment$human[i]
    
    pests_i <- last_state[i, "pests"]
    trees_i <- last_state[i, "pests"]
    
    connectivity_i <- sim$environment$connectivity[i, ]
    
    
    # Compute chance for a single human to find a single pest
    p <- min(human_area/patch_area, 1)
    
    # Compute emigration rate
    transport_i <-  sum(dneedles(n=pests_i, p=p, K=human_i)*(0:pests_i))
    
    # Split emigrants by destination according to the connectivity matrix
    transport_ij <- transport_i * connectivity_i
    
    return(transport_ij)
  }
  
  transport <- t(sapply(1:N, transport_det))
  
  # Move pests
  migration_det <- function(i)
  {
    
    pests_i <- last_state[i, "pests"]
    immigration <- sum(transport[,i])
    emigration <- sum(transport[i,]) 
    
    migrated_pest_i <- pests_i + immigration - emigration
    
    return(migrated_pest_i)
    
  }
  
  migrated_pests <- sapply(1:N, migration_det)
  
  # Compute local infection rate
  infection_det <- function(i)
  {
    
    # Extract relevant parameters
    patch_area <- sim$environment$area[i]
    pest_area <- sim$param$pest_area
    trees_i <- last_state[i, "trees"]
    migrated_pests_i <- migrated_pests [i]
    
    # Compute chance for a single pest to find a single tree
    p <- min(pest_area/patch_area, 1)
    
    # Compute the number of successful infections
    infection_i <- sum(dneedles(n=migrated_pests_i, p=p, K=trees_i)*(0:migrated_pests_i))
    
    return(infection_i)
  }
  
  infection <- sapply(1:N, infection_det)
  
  # Kill old generation of pests and give birth to new pests
  pest_det <- function(i)
  {
    # Extract relevant parameters
    birth_rate <- sim$param$birth_rate
    
    births <- birth_rate*infection[i]    
    pest_i <- births
    
    return (pest_i)
  }
  
  pests <- sapply (1:N, pest_det)
  
  # Kill infected trees
  tree_det <- function(i)
  {
    last_state[i, "trees"] - infection[i]
  }

  trees <- sapply(1:N, tree_det)
  
  # Generate new state
  new_state <- as.matrix(data.frame (trees=trees, pests=pests))
  
  return (new_state)
}

# Stochastic simulation ####

# Stochastic modelling of emerald ash borer spread as a contagion
eab_contagion_sto <- function(sim)
{
  # Resume simulation from the last time step
  last_time <- max(which(!is.na(sim$state), arr.ind=T)[,3])
  last_state <- sim$state[,,last_time]
  
  # Extract number of patches
  N <- sim$environment$N
  
  # Calculate transport
  transport_sto <- function(i)
  {
    # Check if human transport is in effect
    if (sum(sim$environment$human)==0)
    {
      return(rep(0, times=N))
    }
      
    # Extract relevant parameters
    human_area <- sim$param$human_area
    patch_area <- sim$environment$area[i]
    human_i <- sim$environment$human[i]
    
    pests_i <- last_state[i, "pests"]
    trees_i <- last_state[i, "pests"]
    
    connectivity_i <- sim$environment$connectivity[i, ]
    
    # Compute chance for a single human to find a single pest
    p <- min(human_area/patch_area, 1)
    
    # Compute emigration rate
    transport_i <- rneedles(n=pests_i, p=p, K=human_i)
    
    # Split emigrants by destination according to the connectivity matrix
    transport_ij <- rmultinom(1, transport_i, connectivity_i)
    
    return(transport_ij)
  }
  
  transport <- t(sapply(1:N, transport_sto))
  
  # Move pests
  migration_sto <- function(i)
  {
    
    pests_i <- last_state[i, "pests"]
    immigration <- sum(transport[,i])
    emigration <- sum(transport[i,]) 
    
    migrated_pest_i <- pests_i + immigration - emigration
    
    return(migrated_pest_i)
    
  }
  
  migrated_pests <- sapply(1:N, migration_sto)
  
  # Compute local infection rate
  infection_sto <- function(i)
  {
    
    # Extract relevant parameters
    patch_area <- sim$environment$area[i]
    pest_area <- sim$param$pest_area
    trees_i <- last_state[i, "trees"]
    migrated_pests_i <- migrated_pests [i]
    
    # Compute chance for a single pest to find a single tree
    p <- min(pest_area/patch_area, 1)
    
    # Compute the number of successful infections
    infection_i <- rneedles(n=migrated_pests_i, p=p, K=trees_i)
    
    return(infection_i)
  }
  
  infection <- sapply(1:N, infection_sto)
  
  # Kill old generation of pests and give birth to new pests
  pest_sto <- function(i)
  {
    # Extract relevant parameters
    birth_rate <- sim$param$birth_rate
    
    births <- sum(rpois(n=infection[i], lambda=birth_rate))
    pest_i <- births
    
    return (pest_i)
  }
  
  pests <- sapply (1:N, pest_sto)
  
  # Kill infected trees
  tree_sto <- function(i)
  {
    last_state[i, "trees"] - infection[i]
  }
  
  trees <- sapply(1:N, tree_sto)
  
  # Generate new state
  new_state <- as.matrix(data.frame (trees=trees, pests=pests))
  
  return (new_state)
}