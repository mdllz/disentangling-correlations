# load packages and data generating function -----------------------------------

library(mvtnorm)
library(rtdists)

# function to simulate DDM data
make_data <- function(parameters, n_trials){
  # First for the easy condition
  v_easy <- parameters['v']
  a_easy <- parameters['a']
  t0_easy <- parameters['t0']
  data_easy <- rdiffusion(n = n_trials,
                          v = v_easy,
                          a = a_easy,
                          t0 = t0_easy,
                          stop_on_error = TRUE)
  data_easy$difficulty <- 'easy'
  
  # Then for the hard condition, first reconstructing parameters.
  v_hard <- v_easy + parameters['v_diff']
  a_hard <- a_easy + parameters['a_diff']
  t0_hard <- t0_easy + parameters['t0_diff']
  data_hard <- rdiffusion(n = n_trials,
                          v = v_hard,
                          a = a_hard,
                          t0 = t0_hard,
                          stop_on_error = TRUE)
  data_hard$difficulty <- 'hard'
  
  data_sim <- rbind(data_easy, data_hard)
  
  return(data_sim)
}

# simulate data with 1000 trials per condition without correlations between parameters -----------------------------

n_subj <- 1000 # Number of subjects
n_trials <- 1000 # Number of trials per condition

# we will use the difference parameterization used in simulation 5 in Grange & Schuch (2023)
# e.g. a_hard = a_easy + a_d

# group effects or not (TRUE means no difference in the group means between conditions)
centered_group_effects <- FALSE # needs to be run twice, once with TRUE, once with FALSE

set.seed(5464)

# Define the group-level means for our main and difference parameters
if(centered_group_effects){ 
  # no difference in the population means between conditions
  #v #v_d #a  #a_d #t0 #t0_d
  group_means <- c(2.8, 0, 1.25, 0, .35, 0) # easy condition from table 6
} else { 
  # Specify group-level differences (table 6; large effect)
  #v #v_d #a  #a_d #t0 #t0_d
  group_means <- c(2.8, -0.8, 1.25, 0.32, .35, 0.08) 
}

# This will be the median standard deviation for the effects (50%)
max_effect_sds <- c(.5, .1, .05)
names(max_effect_sds) <-c("v", "a", "t0")

# multipliers to get the varying extents of individual differences
sd_multipliers <- seq(0, 2,length.out = 11)

# count the number of times that resampling is needed
attempt <- 1
count_subj_resamples <- NA
multiple_attempts <- numeric()

for(sd_multiplier in sd_multipliers){
  # We generate data without any correlation between the parameters
  sigma <- diag(c(0.8^2, 
                  (sd_multiplier*max_effect_sds['v'])^2, 
                  0.3^2, 
                  (sd_multiplier*max_effect_sds['a'])^2,
                  0.05^2, 
                  (sd_multiplier*max_effect_sds['t0'])^2))
  
  # Check whether our parameters are valid, (i.e. a > 0 and t0>0)
  # There are better ways to do this, with log transformations
  # But we wanted to keep the simulations as close to the original paper as possible.
  # Thus instead we resample until we no longer have any of those values
  data_bad <- TRUE
  while(data_bad){
    # simulate individual parameters
    subj_params <- mvtnorm::rmvnorm(n_subj,mean=group_means,sigma=sigma)
    colnames(subj_params) <- c("v", "v_diff", "a", "a_diff", "t0", "t0_diff")
    # Make a short output that can be checked, also reconstructing v_hard, a_hard and t0_hard
    for_print <- as.data.frame(subj_params)
    for_print$v_hard <- for_print$v + for_print$v_diff
    for_print$a_hard <- for_print$a + for_print$a_diff
    for_print$t0_hard <- for_print$t0 + for_print$t0_diff
    if(all(for_print$a > 0 & for_print$t0 > 0 & for_print$a_hard > 0 & for_print$t0_hard > 0)) {
      data_bad <- FALSE
    } else{
      attempt <- attempt + 1
      multiple_attempts <- c(multiple_attempts, sd_multiplier) # store for which multiplier resampling occurred
      count_subj_resamples <- c(count_subj_resamples, sum(!(for_print$a > 0 & for_print$t0 > 0 & for_print$a_hard > 0 & for_print$t0_hard > 0)))
      print(paste0("attempt: ", attempt))
    }
  }
  print(apply(for_print, 2, quantile, probs = c(.025, .5, .975))) # Print 95% credible interval to check if they're reasonable
  
  # store the data as a list, each artificial participant as a list element
  data <- vector("list", length = n_subj) 
  for(i in 1:n_subj){
    tmp <- make_data(subj_params[i,], n_trials)
    tmp$subjects <- i
    data[[i]] <- tmp
  }
  save(data, subj_params, file = paste0("simulated-data-1000/SDmult_", sd_multiplier, "_centr_", centered_group_effects, "_.RData"))
}

## simulate data with 10000 trials per condition without correlation ----------

n_subj <- 1000 # Number of subjects
n_trials <- 10000 # Number of trials per condition

# we will use the difference parameterization used in simulation 5 in Grange & Schuch (2023)
# e.g. a_hard = a_easy + a_d

# group effects or not (TRUE means no difference in the group means between conditions)
centered_group_effects <- TRUE # needs to be run twice, once with TRUE, once with FALSE

set.seed(5464)

# Define the group-level means for our main and difference parameters
if(centered_group_effects){ 
  # no difference in the population means between conditions
  #v #v_d #a  #a_d #t0 #t0_d
  group_means <- c(2.8, 0, 1.25, 0, .35, 0) # easy condition from table 6
} else { 
  # Specify group-level differences (table 6; large effect)
  #v #v_d #a  #a_d #t0 #t0_d
  group_means <- c(2.8, -0.8, 1.25, 0.32, .35, 0.08) 
}

# This will be the median standard deviation for the effects (50%)
max_effect_sds <- c(.5, .1, .05)
names(max_effect_sds) <-c("v", "a", "t0")

# multipliers to get the varying extents of individual differences
sd_multipliers <- seq(0, 2,length.out = 11)

# count the number of times that resampling is needed
attempt <- 1
count_subj_resamples <- NA
multiple_attempts <- numeric()

for(sd_multiplier in sd_multipliers){
  # We generate data without any correlation between the parameters
  sigma <- diag(c(0.8^2, 
                  (sd_multiplier*max_effect_sds['v'])^2, 
                  0.3^2, 
                  (sd_multiplier*max_effect_sds['a'])^2,
                  0.05^2, 
                  (sd_multiplier*max_effect_sds['t0'])^2))
  
  # Check whether our parameters are valid, (i.e. a > 0 and t0>0)
  # There are better ways to do this, with log transformations
  # But we wanted to keep the simulations as close to the original paper as possible.
  # Thus instead we resample until we no longer have any of those values
  data_bad <- TRUE
  while(data_bad){
    # simulate individual parameters
    subj_params <- mvtnorm::rmvnorm(n_subj,mean=group_means,sigma=sigma)
    colnames(subj_params) <- c("v", "v_diff", "a", "a_diff", "t0", "t0_diff")
    # Make a short output that can be checked, also reconstructing v_hard, a_hard and t0_hard
    for_print <- as.data.frame(subj_params)
    for_print$v_hard <- for_print$v + for_print$v_diff
    for_print$a_hard <- for_print$a + for_print$a_diff
    for_print$t0_hard <- for_print$t0 + for_print$t0_diff
    if(all(for_print$a > 0 & for_print$t0 > 0 & for_print$a_hard > 0 & for_print$t0_hard > 0)) {
      data_bad <- FALSE
    } else{
      attempt <- attempt + 1
      multiple_attempts <- c(multiple_attempts, sd_multiplier) # store for which multiplier resampling occurred
      count_subj_resamples <- c(count_subj_resamples, sum(!(for_print$a > 0 & for_print$t0 > 0 & for_print$a_hard > 0 & for_print$t0_hard > 0)))
      print(paste0("attempt: ", attempt))
    }
  }
  print(apply(for_print, 2, quantile, probs = c(.025, .5, .975))) # Print 95% credible interval to check if they're reasonable
  
  # store the data as a list, each artificial participant as a list element
  data <- vector("list", length = n_subj) 
  for(i in 1:n_subj){
    tmp <- make_data(subj_params[i,], n_trials)
    tmp$subjects <- i
    data[[i]] <- tmp
  }
  save(data, subj_params, file = paste0("simulated-data-10000/SDmult_", sd_multiplier, "_centr_", centered_group_effects, "_.RData"))
}



## simulate data with 200 trials per condition without correlation ----------

n_subj <- 1000 # Number of subjects
n_trials <- 200 # Number of trials per condition

# we will use the difference parameterization used in simulation 5 in Grange & Schuch (2023)
# e.g. a_hard = a_easy + a_d

# group effects or not (TRUE means no difference in the group means between conditions)
centered_group_effects <- TRUE # needs to be run twice, once with TRUE, once with FALSE

set.seed(5464)

# Define the group-level means for our main and difference parameters
if(centered_group_effects){ 
  # no difference in the population means between conditions
  #v #v_d #a  #a_d #t0 #t0_d
  group_means <- c(2.8, 0, 1.25, 0, .35, 0) # easy condition from table 6
} else { 
  # Specify group-level differences (table 6; large effect)
  #v #v_d #a  #a_d #t0 #t0_d
  group_means <- c(2.8, -0.8, 1.25, 0.32, .35, 0.08) 
}

# This will be the median standard deviation for the effects (50%)
max_effect_sds <- c(.5, .1, .05)
names(max_effect_sds) <-c("v", "a", "t0")

# multipliers to get the varying extents of individual differences
sd_multipliers <- seq(0, 2,length.out = 11)

# count the number of times that resampling is needed
attempt <- 1
count_subj_resamples <- NA
multiple_attempts <- numeric()

for(sd_multiplier in sd_multipliers){
  # We generate data without any correlation between the parameters
  sigma <- diag(c(0.8^2, 
                  (sd_multiplier*max_effect_sds['v'])^2, 
                  0.3^2, 
                  (sd_multiplier*max_effect_sds['a'])^2,
                  0.05^2, 
                  (sd_multiplier*max_effect_sds['t0'])^2))
  
  # Check whether our parameters are valid, (i.e. a > 0 and t0>0)
  # There are better ways to do this, with log transformations
  # But we wanted to keep the simulations as close to the original paper as possible.
  # Thus instead we resample until we no longer have any of those values
  data_bad <- TRUE
  while(data_bad){
    # simulate individual parameters
    subj_params <- mvtnorm::rmvnorm(n_subj,mean=group_means,sigma=sigma)
    colnames(subj_params) <- c("v", "v_diff", "a", "a_diff", "t0", "t0_diff")
    # Make a short output that can be checked, also reconstructing v_hard, a_hard and t0_hard
    for_print <- as.data.frame(subj_params)
    for_print$v_hard <- for_print$v + for_print$v_diff
    for_print$a_hard <- for_print$a + for_print$a_diff
    for_print$t0_hard <- for_print$t0 + for_print$t0_diff
    if(all(for_print$a > 0 & for_print$t0 > 0 & for_print$a_hard > 0 & for_print$t0_hard > 0)) {
      data_bad <- FALSE
    } else{
      attempt <- attempt + 1
      multiple_attempts <- c(multiple_attempts, sd_multiplier) # store for which multiplier resampling occurred
      count_subj_resamples <- c(count_subj_resamples, sum(!(for_print$a > 0 & for_print$t0 > 0 & for_print$a_hard > 0 & for_print$t0_hard > 0)))
      print(paste0("attempt: ", attempt))
    }
  }
  print(apply(for_print, 2, quantile, probs = c(.025, .5, .975))) # Print 95% credible interval to check if they're reasonable
  
  # store the data as a list, each artificial participant as a list element
  data <- vector("list", length = n_subj) 
  for(i in 1:n_subj){
    tmp <- make_data(subj_params[i,], n_trials)
    tmp$subjects <- i
    data[[i]] <- tmp
  }
  save(data, subj_params, file = paste0("simulated-data-200/SDmult_", sd_multiplier, "_centr_", centered_group_effects, "_.RData"))
}

## simulate data with 100 trials per condition without correlation -----------------------------------

n_subj <- 1000 # Number of subjects
n_trials <- 100 # Number of trials per condition

# we will use the difference parameterization used in simulation 5 in Grange & Schuch (2023)
# e.g. a_hard = a_easy + a_d

# group effects or not (TRUE means no difference in the group means between conditions)
centered_group_effects <- TRUE # needs to be run twice, once with TRUE, once with FALSE

set.seed(5464)

# Define the group-level means for our main and difference parameters
if(centered_group_effects){ 
  # no difference in the population means between conditions
  #v #v_d #a  #a_d #t0 #t0_d
  group_means <- c(2.8, 0, 1.25, 0, .35, 0) # easy condition from table 6
} else { 
  # Specify group-level differences (table 6; large effect)
  #v #v_d #a  #a_d #t0 #t0_d
  group_means <- c(2.8, -0.8, 1.25, 0.32, .35, 0.08) 
}

# This will be the median standard deviation for the effects (50%)
max_effect_sds <- c(.5, .1, .05)
names(max_effect_sds) <-c("v", "a", "t0")

# multipliers to get the varying extents of individual differences
sd_multipliers <- seq(0, 2,length.out = 11)

# count the number of times that resampling is needed
attempt <- 1
count_subj_resamples <- NA
multiple_attempts <- numeric()

for(sd_multiplier in sd_multipliers){
  # We generate data without any correlation between the parameters
  sigma <- diag(c(0.8^2, 
                  (sd_multiplier*max_effect_sds['v'])^2, 
                  0.3^2, 
                  (sd_multiplier*max_effect_sds['a'])^2,
                  0.05^2, 
                  (sd_multiplier*max_effect_sds['t0'])^2))
  
  # Check whether our parameters are valid, (i.e. a > 0 and t0>0)
  # There are better ways to do this, with log transformations
  # But we wanted to keep the simulations as close to the original paper as possible.
  # Thus instead we resample until we no longer have any of those values
  data_bad <- TRUE
  while(data_bad){
    # simulate individual parameters
    subj_params <- mvtnorm::rmvnorm(n_subj,mean=group_means,sigma=sigma)
    colnames(subj_params) <- c("v", "v_diff", "a", "a_diff", "t0", "t0_diff")
    # Make a short output that can be checked, also reconstructing v_hard, a_hard and t0_hard
    for_print <- as.data.frame(subj_params)
    for_print$v_hard <- for_print$v + for_print$v_diff
    for_print$a_hard <- for_print$a + for_print$a_diff
    for_print$t0_hard <- for_print$t0 + for_print$t0_diff
    if(all(for_print$a > 0 & for_print$t0 > 0 & for_print$a_hard > 0 & for_print$t0_hard > 0)) {
      data_bad <- FALSE
    } else{
      attempt <- attempt + 1
      multiple_attempts <- c(multiple_attempts, sd_multiplier) # store for which multiplier resampling occurred
      count_subj_resamples <- c(count_subj_resamples, sum(!(for_print$a > 0 & for_print$t0 > 0 & for_print$a_hard > 0 & for_print$t0_hard > 0)))
      print(paste0("attempt: ", attempt))
    }
  }
  print(apply(for_print, 2, quantile, probs = c(.025, .5, .975))) # Print 95% credible interval to check if they're reasonable
  
  # store the data as a list, each artificial participant as a list element
  data <- vector("list", length = n_subj) 
  for(i in 1:n_subj){
    tmp <- make_data(subj_params[i,], n_trials)
    tmp$subjects <- i
    data[[i]] <- tmp
  }
  save(data, subj_params, file = paste0("simulated-data-fewer-trials/SDmult_", sd_multiplier, "_centr_", centered_group_effects, "_.RData"))
}

