
# load packages and read data ----------------------------------------------
library(dplyr)
library(tidyr)
library(EMC2) # fitted using EMC2 v2.1.0: https://github.com/ampl-psych/EMC2/releases/tag/v2.1.0

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
