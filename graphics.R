# Libraries ####
library(ggplot2)
#library(igraph)
library(animation)
library(reshape2)

# Data cleaning ####
normalize_state <- function(sim, by_area = TRUE)
{
  # Normalize state by area
  if (by_area)
  {
    sim$state <- sweep(sim$state, MARGIN=1, STATS=sim$environment$area, FUN="/")
  }
  
  # Ensure the values for the state variables range between 0 and 1
  sim$state <- sweep(sim$state, MARGIN=2, STATS=apply(sim$state, MARGIN=2, FUN=max), FUN="/")
  
  return (sim)
}

# Network graphs ####

# Make a static graph showing the state of each patch
static_graph <- function(time=1, sim, display="trees", normalize=TRUE, by_area=TRUE)
{
  # Clean up the data by normalizing it first?
  if (normalize)
  {
    sim <- normalize_state (sim, by_area)
  }
  
  # Extract relevant data
  x <- sim$environment$pos$x
  y <- sim$environment$pos$y
  area <- sim$environment$area
  trees <- sim$state[,,time][,"trees"]
  pests <- sim$state[,,time][,"pests"]
  
  df <- data.frame(x, y, area, trees, pests)
  
  # Make graphs 
  base_graph <- ggplot(df, aes(x=x, y=y)) + theme_bw() + theme (legend.position="none")
  #base_graph <- ggplot(df, aes(x=x, y=y, size=sqrt(area)*100)) + theme_bw() + theme (legend.position="none")
  
  if (display=="trees")
  {
    final_graph <- base_graph + geom_point (size=10, shape=20, aes(colour=trees)) +  scale_colour_gradient (low="white", high="forestgreen", limits = c(0,1))
  }
  else if (display=="pests")
  {
    final_graph <- base_graph + geom_point (size=10, shape=20, aes(colour=pests)) +  scale_colour_gradient (low="white", high="red", limits = c(0,1))
  }
  else
  {
    final_graph <- base_graph + geom_point(size=10, shape=20)
  }
  return (final_graph)
}

# Make a stack of static graph showing the state of each patch
dynamic_graph <- function (sim, display="trees",  normalize=TRUE, by_area=TRUE)
{
  # Clean up the data by normalizing it first?
  if (normalize)
  {
    sim <- normalize_state (sim, by_area)
  }
  
  # Make graphs
  time <- dim(sim$state)[3]
  graphs <- lapply (1:time, static_graph, sim=sim, display=display, normalize=FALSE)
  
  return(graphs)
}

# Turn a stack of graphs into an animation
anim_graph <- function (graphs, ...)
{
  return(saveHTML(print(graphs), ...))
}

# Time series graphs ####

state_ts_plot <- function(sim, normalize=TRUE, by_area=TRUE)
{
  # Clean up the data by normalizing it first?
  if (normalize)
  {
    sim <- normalize_state (sim, by_area)
  }
  
  df <- melt(sim$state)
  names(df) <- c("patch", "var", "time", "value")
  
  ggplot(df, aes(x=time, y=value, colour=var)) + facet_grid(patch~var) + geom_line() +theme_bw() + theme (legend.position="none")
  
}
