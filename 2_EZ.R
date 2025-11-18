library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

### ez diffusion fitting function: Code from Grange & Schuch (2023) https://osf.io/2h4jt/
fit_ez <- function(Pc, VRT, MRT, s = 1, n_trials){
  s2 = s^2
  # The default value for the scaling parameter s equals .1
  
  if (Pc == 0)
    cat("Oops, Pc == 0!\n")
  if (Pc == 0.5)
    cat("Oops, Pc == .5!\n")
  if (Pc == 1){
    Pc = 1 - (1 / (2 * n_trials))
  }
  # If Pc equals 0, .5, or 1, the method will not work, and
  # an edge-correction is required.
  
  L = qlogis(Pc)
  # The function "qlogis" calculates the logit.
  x = L*(L*Pc^2 - L*Pc + Pc - 0.5)/VRT
  v = sign(Pc-0.5)*s*x^(1/4)
  # This gives drift rate.
  
  a = s2*qlogis(Pc)/v
  # This gives boundary separation.
  
  y   = -v*a/s2
  MDT = (a/(2*v))*(1-exp(y))/(1+exp(y))
  Ter = MRT-MDT
  # This gives nondecision time.
  
  return(list(v, a, Ter))
}


# Fit the EZ diffusion model to data without correlation 1000 trials  ---------------------------------------------------
data_names <- list.files("simulated-data-1000")
ez_results <- list()
true_params <- list()

for(n in 1:length(data_names)){
  ## load dataset --------------------
  load(paste0("simulated-data-1000/", data_names[n]))
  n_trials <- length(data)
  
  # create summary statistics -------------------
  dat_summary <- bind_rows(data) %>% 
    group_by(subjects, difficulty) %>% 
    summarize(accuracy = mean(response == "upper")*100,
              mean_rt = mean(rt),
              var_rt = var(rt))
  
  # fit EZ diffusion model ---------------
  # correction method
  ez <- dat_summary %>% 
    mutate(v = 0, a = 0, t0 = 0)
  
  for(i in 1:nrow(ez)){
    
    # apply edge correction to accuracy
    if(ez$accuracy[i] == 0){
      ez$accuracy[i] <- (1 / (2 * n_trials)) * 100
    }
    if(ez$accuracy[i] == 50){
      ez$accuracy[i] <- (0.5 + 1 / (2 * n_trials)) * 100
    }
    if(ez$accuracy[i] == 100){
      ez$accuracy[i] <- (1 - 1 / (2 * n_trials)) * 100
    }
    
    fit <- fit_ez(Pc = ez$accuracy[i] / 100, 
                  VRT = ez$var_rt[i], 
                  MRT = ez$mean_rt[i], 
                  s = 1, 
                  n_trials = n_trials)
    
    ez$v[i] <- fit[[1]]
    ez$a[i] <- fit[[2]]
    ez$t0[i] <- fit[[3]]
  }
  ez_results[[n]] <- ez
  #names(ez_results[[n]]) <- data_names[n]
  true_params[[n]] <- subj_params
}

## compute correlations -------------------------------------------------------

# store correlations in a list
ez_corr <- list()
cor_t0diff_adiff <- list()
res <- list()

for(n in 1:length(ez_results)){
  
res[[n]] <- ez_results[[n]] %>% 
    pivot_wider(id_cols = subjects, 
                names_from = difficulty, 
                values_from = c(v, a, t0)) %>% 
    mutate(a_difference = a_hard - a_easy, 
           v_difference = v_hard - v_easy, 
           t0_difference = t0_hard - t0_easy) %>% 
  ungroup(subjects) %>% 
    select(a_easy, a_hard, v_easy, v_hard, t0_easy, t0_hard, 
           a_difference, v_difference, t0_difference)
  
  ez_corr[[n]] <- round(cor(res[[n]]),3)
  cor_t0diff_adiff[[n]] <- cor.test(res[[n]]$t0_difference,res[[n]]$a_difference)
  
}

names(cor_t0diff_adiff) <- data_names

## save results 
#save(true_params, res, ez_results, ez_corr, cor_t0diff_adiff, file = "ez-results/ez_output_1000.RData")

# Fit the EZ diffusion model to data without true correlation and 10000 trials ----------------------------

data_names <- list.files("simulated-data-10000")
ez_results <- list()
true_params <- list()

for(n in 1:length(data_names)){
  ## load dataset --------------------
  load(paste0("simulated-data-10000/", data_names[n]))
  n_trials <- length(data)
  
  # create summary statistics -------------------
  dat_summary <- bind_rows(data) %>% 
    group_by(subjects, difficulty) %>% 
    summarize(accuracy = mean(response == "upper")*100,
              mean_rt = mean(rt),
              var_rt = var(rt))
  
  # fit EZ diffusion model ---------------
  # correction method
  ez <- dat_summary %>% 
    mutate(v = 0, a = 0, t0 = 0)
  
  for(i in 1:nrow(ez)){
    
    # apply edge correction to accuracy
    if(ez$accuracy[i] == 0){
      ez$accuracy[i] <- (1 / (2 * n_trials)) * 100
    }
    if(ez$accuracy[i] == 50){
      ez$accuracy[i] <- (0.5 + 1 / (2 * n_trials)) * 100
    }
    if(ez$accuracy[i] == 100){
      ez$accuracy[i] <- (1 - 1 / (2 * n_trials)) * 100
    }
    
    fit <- fit_ez(Pc = ez$accuracy[i] / 100, 
                  VRT = ez$var_rt[i], 
                  MRT = ez$mean_rt[i], 
                  s = 1, 
                  n_trials = n_trials)
    
    ez$v[i] <- fit[[1]]
    ez$a[i] <- fit[[2]]
    ez$t0[i] <- fit[[3]]
  }
  ez_results[[n]] <- ez
  #names(ez_results[[n]]) <- data_names[n]
  true_params[[n]] <- subj_params
}

## compute correlations -------------------------------------------------------

# store correlations in a list
ez_corr <- list()
cor_t0diff_adiff <- list()
res <- list()

for(n in 1:length(ez_results)){
  
  res[[n]] <- ez_results[[n]] %>% 
    pivot_wider(id_cols = subjects, 
                names_from = difficulty, 
                values_from = c(v, a, t0)) %>% 
    mutate(a_difference = a_hard - a_easy, 
           v_difference = v_hard - v_easy, 
           t0_difference = t0_hard - t0_easy) %>% 
    ungroup(subjects) %>% 
    select(a_easy, a_hard, v_easy, v_hard, t0_easy, t0_hard, 
           a_difference, v_difference, t0_difference)
  
  ez_corr[[n]] <- round(cor(res[[n]]),3)
  cor_t0diff_adiff[[n]] <- cor.test(res[[n]]$t0_difference,res[[n]]$a_difference)
  
}

names(cor_t0diff_adiff) <- data_names

## save results ----------------------------------------------------------------
#save(true_params, res, ez_results, ez_corr, cor_t0diff_adiff, file = "ez-results/ez_output_10000.RData")

# Fit the EZ diffusion model to data without true correlation and 200 trials ----------------------------

data_names <- list.files("simulated-data-200")
ez_results <- list()
true_params <- list()

for(n in 1:length(data_names)){
  ## load dataset --------------------
  load(paste0("simulated-data-200/", data_names[n]))
  n_trials <- length(data)
  
  # create summary statistics -------------------
  dat_summary <- bind_rows(data) %>% 
    group_by(subjects, difficulty) %>% 
    summarize(accuracy = mean(response == "upper")*100,
              mean_rt = mean(rt),
              var_rt = var(rt))
  
  # fit EZ diffusion model ---------------
  # correction method
  ez <- dat_summary %>% 
    mutate(v = 0, a = 0, t0 = 0)
  
  for(i in 1:nrow(ez)){
    
    # apply edge correction to accuracy
    if(ez$accuracy[i] == 0){
      ez$accuracy[i] <- (1 / (2 * n_trials)) * 100
    }
    if(ez$accuracy[i] == 50){
      ez$accuracy[i] <- (0.5 + 1 / (2 * n_trials)) * 100
    }
    if(ez$accuracy[i] == 100){
      ez$accuracy[i] <- (1 - 1 / (2 * n_trials)) * 100
    }
    
    fit <- fit_ez(Pc = ez$accuracy[i] / 100, 
                  VRT = ez$var_rt[i], 
                  MRT = ez$mean_rt[i], 
                  s = 1, 
                  n_trials = n_trials)
    
    ez$v[i] <- fit[[1]]
    ez$a[i] <- fit[[2]]
    ez$t0[i] <- fit[[3]]
  }
  ez_results[[n]] <- ez
  #names(ez_results[[n]]) <- data_names[n]
  true_params[[n]] <- subj_params
}

## compute correlations -------------------------------------------------------

# store correlations in a list
ez_corr <- list()
cor_t0diff_adiff <- list()
res <- list()

for(n in 1:length(ez_results)){
  
  res[[n]] <- ez_results[[n]] %>% 
    pivot_wider(id_cols = subjects, 
                names_from = difficulty, 
                values_from = c(v, a, t0)) %>% 
    mutate(a_difference = a_hard - a_easy, 
           v_difference = v_hard - v_easy, 
           t0_difference = t0_hard - t0_easy) %>% 
    ungroup(subjects) %>% 
    select(a_easy, a_hard, v_easy, v_hard, t0_easy, t0_hard, 
           a_difference, v_difference, t0_difference)
  
  ez_corr[[n]] <- round(cor(res[[n]]),3)
  cor_t0diff_adiff[[n]] <- cor.test(res[[n]]$t0_difference,res[[n]]$a_difference)
  
}

names(cor_t0diff_adiff) <- data_names

## save results ----------------------------------------------------------------
#save(true_params, res, ez_results, ez_corr, cor_t0diff_adiff, file = "ez-results/ez_output_200.RData")


# Fit the EZ diffusion model to data without true correlation and 100 trials ----------------------------

data_names <- list.files("simulated-data-100")
ez_results <- list()
true_params <- list()

for(n in 1:length(data_names)){
  ## load dataset --------------------
  load(paste0("simulated-data-100/", data_names[n]))
  n_trials <- length(data)
  
  # create summary statistics -------------------
  dat_summary <- bind_rows(data) %>% 
    group_by(subjects, difficulty) %>% 
    summarize(accuracy = mean(response == "upper")*100,
              mean_rt = mean(rt),
              var_rt = var(rt))
  
  # fit EZ diffusion model ---------------
  # correction method
  ez <- dat_summary %>% 
    mutate(v = 0, a = 0, t0 = 0)
  
  for(i in 1:nrow(ez)){
    
    # apply edge correction to accuracy
    if(ez$accuracy[i] == 0){
      ez$accuracy[i] <- (1 / (2 * n_trials)) * 100
    }
    if(ez$accuracy[i] == 50){
      ez$accuracy[i] <- (0.5 + 1 / (2 * n_trials)) * 100
    }
    if(ez$accuracy[i] == 100){
      ez$accuracy[i] <- (1 - 1 / (2 * n_trials)) * 100
    }
    
    fit <- fit_ez(Pc = ez$accuracy[i] / 100, 
                  VRT = ez$var_rt[i], 
                  MRT = ez$mean_rt[i], 
                  s = 1, 
                  n_trials = n_trials)
    
    ez$v[i] <- fit[[1]]
    ez$a[i] <- fit[[2]]
    ez$t0[i] <- fit[[3]]
  }
  ez_results[[n]] <- ez
  #names(ez_results[[n]]) <- data_names[n]
  true_params[[n]] <- subj_params
}

## compute correlations -------------------------------------------------------

# store correlations in a list
ez_corr <- list()
cor_t0diff_adiff <- list()
res <- list()

for(n in 1:length(ez_results)){
  
  res[[n]] <- ez_results[[n]] %>% 
    pivot_wider(id_cols = subjects, 
                names_from = difficulty, 
                values_from = c(v, a, t0)) %>% 
    mutate(a_difference = a_hard - a_easy, 
           v_difference = v_hard - v_easy, 
           t0_difference = t0_hard - t0_easy) %>% 
    ungroup(subjects) %>% 
    select(a_easy, a_hard, v_easy, v_hard, t0_easy, t0_hard, 
           a_difference, v_difference, t0_difference)
  
  ez_corr[[n]] <- round(cor(res[[n]]),3)
  cor_t0diff_adiff[[n]] <- cor.test(res[[n]]$t0_difference,res[[n]]$a_difference)
  
}

names(cor_t0diff_adiff) <- data_names

## save results ----------------------------------------------------------------
#save(true_params, res, ez_results, ez_corr, cor_t0diff_adiff, file = "ez-results/ez_output_100.RData")

