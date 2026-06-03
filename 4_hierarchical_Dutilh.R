
# load packages and read data ----------------------------------------------
library(dplyr)
library(tidyr)
library(EMC2) # fitted using EMC2 v2.1.0: https://github.com/ampl-psych/EMC2/releases/tag/v2.1.0

# Function for EZ diffusion: code taken from Grange & Schuch (2023) https://osf.io/2h4jt/: 
### ez diffusion fit code
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

# Simple DDM without RT cleaning ------------------------------------------------------
full_data <- read.csv("data-DutilhEtAl-exp1/raw_factorial_data.csv")
set.seed(1234)

data <- full_data %>%
  filter(ease == 1,
         sa == 1,
         bias == "no") %>%
  select(id = pp,
         response = correct,
         rt = rt) %>%
  mutate(join_id = 1:length(rt))

# create the two conditions by randomly splitting the data into two

a_data <- data %>%
  group_by(id) %>%
  sample_n(length(rt) / 2, replace = FALSE)

b_data <- data %>%
  anti_join(a_data, by = "join_id")

a_data <- a_data %>%
  select(-join_id) %>%
  ungroup() %>%
  mutate(condition = "A")

b_data <- b_data %>%
  select(-join_id) %>%
  ungroup() %>%
  mutate(condition = "B")

# combine data
all_data <- rbind(a_data, b_data)


all_data <- bind_rows(a_data, b_data) %>%
  transmute(difficulty = as.factor(condition),
            R = as.factor(ifelse(response == 0, "lower", "upper")),
            subjects = as.factor(id),
            rt = rt) %>%
  as.data.frame()

Dmat <- matrix(c(-0.5,0.5), nrow = 2,dimnames=list(NULL,"_d"))
ddm_design <- design(formula = list(v ~ difficulty,
                                    t0 ~ difficulty,
                                    a ~ difficulty), dat = all_data, model = DDM,
                     contrasts =  list(difficulty = Dmat))

# use more informative priors
prior_input <- prior(ddm_design, mu_mean = c(2, 0, -0.9, 0, 1, 0), mu_sd = c(2,1,0.6,1,0.6, 1))
#plot_prior(prior_input, ddm_design, selection = "mu")

emc_object <- make_emc(all_data, ddm_design, prior_list = prior_input, n_chains = 4, compress = FALSE, rt_resolution = 1e-20)

fit_dutilh_exp1 <- fit(emc_object, cores_per_chain = 8, iter = 2000)

save(fit_dutilh_exp1, all_data, file = paste0("fit-objects/dutilh-exp1-hierarchical/dutilh_exp1.RData"))


# Simple DDM with RT cleaning ------------------------------------------------------
full_data <- read.csv("data-DutilhEtAl-exp1/raw_factorial_data.csv")
set.seed(1234)

# 10 repetitions
for(rep in 1:10){
  data <- full_data %>%
    filter(ease == 1,
           sa == 1,
           bias == "no") %>%
    select(id = pp,
           response = correct,
           rt = rt) %>%
    mutate(join_id = 1:length(rt))
  
  # create the two conditions by randomly splitting the data into two
  
  a_data <- data %>%
    group_by(id) %>%
    sample_n(length(rt) / 2, replace = FALSE)
  
  b_data <- data %>%
    anti_join(a_data, by = "join_id")
  
  a_data <- a_data %>%
    select(-join_id) %>%
    ungroup() %>%
    mutate(condition = "A")
  
  b_data <- b_data %>%
    select(-join_id) %>%
    ungroup() %>%
    mutate(condition = "B")
  
  # combine data
  all_data <- rbind(a_data, b_data)
  
  
  all_data <- bind_rows(a_data, b_data) %>%
    transmute(difficulty = as.factor(condition),
              R = as.factor(ifelse(response == 0, "lower", "upper")),
              subjects = as.factor(id),
              rt = rt) %>%
    # only keep RTs greater than 200 ms
    filter(rt > 0.2) %>%
    as.data.frame()
  
  Dmat <- matrix(c(-0.5,0.5), nrow = 2,dimnames=list(NULL,"_d"))
  ddm_design <- design(formula = list(v ~ difficulty,
                                      t0 ~ difficulty,
                                      a ~ difficulty), dat = all_data, model = DDM,
                       contrasts =  list(difficulty = Dmat))
  
  # use more informative priors
  prior_input <- prior(ddm_design, mu_mean = c(2, 0, -0.9, 0, 1, 0), mu_sd = c(2,1,0.6,1,0.6, 1))
  #plot_prior(prior_input, ddm_design, selection = "mu")
  
  emc_object <- make_emc(all_data, ddm_design, prior_list = prior_input, n_chains = 4, compress = FALSE, rt_resolution = 1e-20)
  
  fit_dutilh_exp1 <- fit(emc_object, cores_per_chain = 7, iter = 2000)
  save(fit_dutilh_exp1, all_data, file = paste0("fit-objects/dutilh-exp1-hierarchical/dutilh_exp1_cleaned_nr_", rep, ".RData"))
  
}

## Full DDM ---------

files <- paste0("fit-objects/dutilh-exp1-hierarchical/dutilh_exp1_cleaned_nr_", 1:10, ".RData")

for(f in 1:length(files)){
  load(paste0("fit-objects/dutilh-exp1-hierarchical/", files[f]))
  Dmat <- matrix(c(-0.5,0.5), nrow = 2,dimnames=list(NULL,"_d"))
  full_ddm_design <- design(data = all_data, 
                            model=DDM, 
                            formula = list(v ~ difficulty,
                                           t0 ~ difficulty,
                                           a ~ difficulty,
                                           sv ~ difficulty,
                                           SZ ~ difficulty,
                                           st0 ~difficulty,
                                           Z ~ difficulty),
                            contrasts = list(difficulty = Dmat)
  )
  
  prior_input <- prior(full_ddm_design, mu_mean = c(2, 0, -0.9, 0, 1, 0, 0, 0, 0, 0, -3, 0,0,0), 
                       mu_sd = c(2, 1, 0.6, 1, 0.6, 1, 1, 1, 0.7, 0.5, 1.1, 1,0.2,0.4)
  )
  
  #plot_prior(prior_input, full_ddm_design)
  
  emc_object <- make_emc(all_data, full_ddm_design, prior_list = prior_input, n_chains = 4, compress = FALSE, rt_resolution = 1e-20)
  
  
  
  fit_dutilh_exp1_full <- fit(emc_object, cores_per_chain = 8, iter = 2000)
  save(fit_dutilh_exp1_full, all_data, file = paste0("fit-objects/dutilh-exp1-hierarchical/full_", files[f]))
}

# Sanity check: effect of RT exclusion on two-step approach --------------------

# EZ: two-step approach with RT cleaning ---------------------------------------
full_data <- read.csv("data-DutilhEtAl-exp1/raw_factorial_data.csv")
set.seed(1234)

data <- full_data %>%
  filter(ease == 1,
         sa == 1,
         bias == "no") %>%
  select(id = pp,
         response = correct,
         rt = rt) %>%
  mutate(join_id = 1:length(rt))

# create the two conditions by randomly splitting the data into two

a_data <- data %>%
  group_by(id) %>%
  sample_n(length(rt) / 2, replace = FALSE)

b_data <- data %>%
  anti_join(a_data, by = "join_id")

a_data <- a_data %>%
  select(-join_id) %>%
  ungroup() %>%
  mutate(condition = "A")

b_data <- b_data %>%
  select(-join_id) %>%
  ungroup() %>%
  mutate(condition = "B")

# combine data

all_data <- bind_rows(a_data, b_data) %>%
  transmute(difficulty = as.factor(condition),
            response = as.factor(ifelse(response == 0, "lower", "upper")),
            subjects = as.factor(id),
            rt = rt) %>%
  as.data.frame()

# without RT cleaning
dat_summary <- all_data %>% 
  group_by(subjects, difficulty) %>% 
  summarize(accuracy = mean(response == "upper")*100,
            mean_rt = mean(rt),
            var_rt = var(rt))

# with RT cleaning
dat_summary_cleaned <- all_data %>% 
  # only keep RTs greater than 200 ms
  filter(rt > 0.2) %>%
  group_by(subjects, difficulty) %>% 
  summarize(accuracy = mean(response == "upper")*100,
            mean_rt = mean(rt),
            var_rt = var(rt))


## get ez estimates for data with RT cleaning ------
ez_cleaned <- dat_summary_cleaned %>% 
  mutate(v = 0, a = 0, t0 = 0)

for(i in 1:nrow(ez_cleaned)){
  
  # apply edge correction to accuracy
  if(ez_cleaned$accuracy[i] == 0){
    ez_cleaned$accuracy[i] <- (1 / (2 * n_trials)) * 100
  }
  if(ez_cleaned$accuracy[i] == 50){
    ez_cleaned$accuracy[i] <- (0.5 + 1 / (2 * n_trials)) * 100
  }
  if(ez_cleaned$accuracy[i] == 100){
    ez_cleaned$accuracy[i] <- (1 - 1 / (2 * n_trials)) * 100
  }
  
  fit <- fit_ez(Pc = ez_cleaned$accuracy[i] / 100, 
                VRT = ez_cleaned$var_rt[i], 
                MRT = ez_cleaned$mean_rt[i], 
                s = 1, 
                n_trials = n_trials)
  
  ez_cleaned$v[i] <- fit[[1]]
  ez_cleaned$a[i] <- fit[[2]]
  ez_cleaned$t0[i] <- fit[[3]]
}

res <- ez_cleaned %>% 
  pivot_wider(id_cols = subjects, 
              names_from = difficulty, 
              values_from = c(v, a, t0)) %>% 
  mutate(a_difference = a_A - a_B, 
         v_difference = v_A - v_B, 
         t0_difference = t0_A - t0_B) %>% 
  ungroup(subjects) %>% 
  select(a_B, a_A, v_B, v_A, t0_B, t0_A, 
         a_difference, v_difference, t0_difference)

ez_cleaned_corr <- round(cor(res),3)
cor_t0diff_adiff_cleaned <- cor.test(res$t0_difference,res$a_difference)

## without RT cleaning ----------

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

res <- ez %>% 
  pivot_wider(id_cols = subjects, 
              names_from = difficulty, 
              values_from = c(v, a, t0)) %>% 
  mutate(a_difference = a_A - a_B, 
         v_difference = v_A - v_B, 
         t0_difference = t0_A - t0_B) %>% 
  ungroup(subjects) %>% 
  select(a_B, a_A, v_B, v_A, t0_B, t0_A, 
         a_difference, v_difference, t0_difference)

ez_corr <- round(cor(res),3)
cor_t0diff_adiff <- cor.test(res$t0_difference,res$a_difference)

cor_t0diff_adiff
cor_t0diff_adiff_cleaned

