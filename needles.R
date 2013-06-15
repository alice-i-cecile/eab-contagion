# Libraries ####
library(plyr)

# Searching for rare objects ####

# We're looking for rare objects
# And have n chances
# Chance to find an object is p

# We start with K objects to be found

# If we've already found k objects, the chance to find an object is:
# 1-(1-p)^(K-k)
# We can't find objects if they've all been taken
found <- function(k, p, K)
{
  if (k < K){
    return (1-(1-p)^(K-k))
  }
  else
  {
    return(0)
  }  
}

# For each trial, we either find an object or we don't
find_paths <- function(n)
{
  paths <- as.matrix(expand.grid(rep(list(c(0, 1)), n)))
  return (paths)
}

# The probability of each branch is dependent on the number of events before it
find_k_branches <- function (paths)
{
  # Compute the number of objects found for each branch
  # Finding an object increases the count by 1
  t(sapply (1:nrow(paths), 
    function(i)
    {
      sapply(1:ncol(paths), function(j, x)
      {
        ifelse(j==1, 0, sum(x[1:(j-1)]))
      }
      , x=paths[i,])
    }
  ))
}

# We can find the probability of each path by multiplying the probability of each branch
find_prob_paths <- function (paths, k_branches, p, K)
{
  prob_branch <- function (path_found, k_branch)
  {
    # If the event was a success, the probability is given by found()
    if (path_found)
    {
      return (found(k_branch, p, K))
    }
    # Otherwise, the probability is given by !found()
    else
    {
      return (1 - found(k_branch, p, K))
    }
  }
  
  # Create a matrix to store the probabilities
  prob_matrix <- paths
  
  # Find the probability of each branch
  prob_matrix[,] <- mapply (FUN=prob_branch, path_found=paths, k_branch=k_branches)
  
  # Find the probabilities of each path
  prob_paths <- apply (prob_matrix, 1, prod)
  
  return (prob_paths)
}

# Now we need to group the paths by the number of objects found
# Summing across the group is the probability mass for that outcome 
find_k_paths <- function (paths)
{
  rowSums(paths)
}

k_probability <- function(prob_paths, k_paths)
{
  df <- data.frame(prob=prob_paths, k=k_paths)
  k_probs <- ddply(df, "k", summarise, total_prob=sum(prob))
  return (k_probs$total_prob)
}

# Putting it all together to make a probability mass function
dneedles2 <- function(n, p, K)
{
  paths <- find_paths(n)
  k_branches <- find_k_branches(paths)
  prob_paths <- find_prob_paths (paths, k_branches, p, K)
  
  k_paths <- find_k_paths (paths)
  pmf <- k_probability (prob_paths, k_paths)
  return (pmf)
}

# And then a RNG function
rneedles2 <- function(n, p, K)
{
  pmf <- dneedles2(n, p, K)
  random <- sample(0:n, 1, prob=pmf)
  return(random)
}

# A faster approach ####

# The pmf approach is extremely slow! Finding the pmf takes O(2^n) time :(
# A tree traversal approach should be doable in O(n) time to generate random samples from the distribution
# Bootstrapping from this rng can be used to approximate the correct pmf if needed in O(n^2) time

# Let's just loop through
rneedles <- function(n, p, K)
{
  # Start with 0 successes
  k <- 0
  
  # Make a decision at each branch
  for (i in 1:n)
  {
    grab <- rbinom(1, 1, found(k, p, K))
    
    if (grab){k = k + 1}
  }
  
  return (k)
}

dneedles <- function(n, p, K, samples=100*(n+1))
{
    results <- replicate(samples, rneedles(n, p, K))
    counts <- count(results)
    counts$freq <- counts$freq / sum(counts$freq)
    
    # Add 0 probability for unseen outcomes
    all_outcomes <- 0:n
    found_outcomes <- counts$x
    missing_outcomes <- setdiff(all_outcomes, found_outcomes)
    
    if (length(missing_outcomes) > 0)
    {
      missing_rows <- data.frame(x=missing_outcomes, freq=0)
      counts <- rbind (missing_rows, counts)
      counts <- counts[order(counts$x), ]
    }
    
    
    pmf <- counts$freq
    
    return(pmf)
}

# Comparing the approaches ####

# Test parameters
n <- 3
p <- 0.1
K <- 2

# Timing 
system.time(dneedles(n, p, K))
system.time(dneedles2(n, p, K))

# Comparing outcomes
pmf1 <- dneedles(n, p, K)
pmf2 <- dneedles2(n, p, K)

plot(pmf1)
lines(pmf2)

plot(pmf1, pmf2)
abline(a=0, b=1)