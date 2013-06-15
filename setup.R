# Initial state (pests and trees) ####

const_initial_state <- function (N, tree_density, pest_density, area) 
{
  trees <- tree_density*area
  pests <- pest_density*area
  patch <- 1:N
  
  state <- data.frame (trees=trees, pests=pests)
  rownames (state) <- patch
  
  return (state)
}

# Area ####
const_area <- function(N, area)
{
  return(rep.int(area, N))
}

rlnorm_area <- function(N, meanlog, sdlog)
{
  return(rlnorm(N, meanlog, sdlog))
}

# Position ####
unif_pos <- function (N, bounding_box=data.frame(x=c(0,1), y=c(0, 1)))
{
  x <- runif(n=N, min=bounding_box$x[1], max=bounding_box$x[2])
y <- runif(n=N, min=bounding_box$y[1], max=bounding_box$y[2])
pos <- data.frame (x=x, y=y)
return (pos)
}

norm_pos <- function (N, sd_x=1, sd_y=1)
{
  x <- rnorm(n=N, mean=0, sd=sd_x)
  y <- rnorm(n=N, mean=0, sd=sd_y)
  pos <- data.frame (x=x, y=y)
  return (pos)
}

# Distance ####
const_dist <- function(N, dist)
{
  return(rep.int(human, N))
}

euclidean_dist <- function(pos)
{
  return (as.matrix(dist(pos,  method = "euclidean", diag=TRUE, upper=TRUE)))
}

manhattan_dist <- function(pos)
{
  return (as.matrix(dist(pos,  method = "manhattan", diag=TRUE, upper=TRUE)))
}

# Human use ####
const_human <- function(N, human=1)
{
  return(rep.int(human, N))
}

prop_area_human <- function(area, human_0)
{
  return (area*human_0)
}                      

# Connectivity ####
compute_connectivity <- function (dist, human, gravity_exponent=-2)
{
  pair_connectivity <- function(dist_cell, human_target)
  {
    connectivity <- dist_cell^gravity_exponent*human_target
    return (connectivity)
  }
  
  N <- length(human)
  
  connectivity_matrix <- matrix(0, N, N)
  
  # rows (i) are origin, columns (j) are target
  for (i in 1:N)
  {
    for (j in 1:N)
    {
      # Connectivity from patches to self is 0
      if (i != j)
      {
        connectivity_matrix[i, j] <- pair_connectivity(dist[i, j], human[j])
      }
    }
  }
  
  # Normalize to ensure total outgoing fraction is 1 for each patch
  connectivity_matrix <- connectivity_matrix/rowSums(connectivity_matrix)
  
  
  return(connectivity_matrix)
}
