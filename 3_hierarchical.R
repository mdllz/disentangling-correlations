# Hierarchical modeling with full variance-covariance matrix ----------------------------------------------
library(EMC2) # fitted using EMC2 v2.1.0: https://github.com/ampl-psych/EMC2/releases/tag/v2.1.0
library(dplyr)


# proportions of SDs
sds <- seq(0.2, 2, by = .2)

for(s in 1:length(sds)){
  load(paste0("simulated-data-1000/SDmult_", sds[s], "_centr_TRUE_.RData")) # centered at zero, so no differences between hard and easy at the group-level
  
  # first 250 subjects
  dat <- dplyr::bind_rows(data) %>% 
    filter(subjects %in% 1:250) %>% 
    transmute(rt = rt,
              R = as.factor(response),
              difficulty = as.factor(difficulty),
              subjects = as.factor(subjects)) %>% 
    as.data.frame()
  
  
  Dmat <- matrix(c(-0.5,0.5), nrow = 2,dimnames=list(NULL,"_d"))
  ddm_design <- design(formula = list(v ~ difficulty,
                                      t0 ~ difficulty,
                                      a ~ difficulty), dat = dat, model = DDM,
                       contrasts =  list(difficulty = Dmat))
  
  # use default priors
  prior_input <- prior(ddm_design)
  #plot_prior(prior_input, ddm_design, selection = "correlation")
  
  emc_object <- make_emc(dat, ddm_design, prior_list = prior_input, n_chains = 4, compress = FALSE, rt_resolution = 1e-20)
  fit_emc <- fit(emc_object, cores_per_chain = 7, iter = 1000)
  
  save(dat, fit_emc, file = paste0("fit-objects/fully-hierarchical/fit_SDmult_", sds[s], "_centr_TRUE_.RData"))
  
}

