library(ggplot2)
library(patchwork)
library(viridis)
library(stringr)
library(dplyr)
library(EMC2) # EMC2 v2.1.0: https://github.com/ampl-psych/EMC2/releases/tag/v2.1.0
library(rtdists)
library(RWiener)
library(bayesplot)

# Check if "plots" folder exists; if not, create it
if (!dir.exists("plots")) {
  dir.create("plots")
}

# DDM visualization ---------
# Code adapted from: Lüken, M., Heathcote, A., Haaf, J.M. et al. Parameter 
# identifiability in evidence-accumulation models: The effect of error rates 
# on the diffusion decision model and the linear ballistic accumulator. 
# Psychon Bull Rev 32, 1411–1424 (2025). https://doi.org/10.3758/s13423-024-02621-1
# https://osf.io/6tpj9

simWienerWalk <- function(start) {
  e <- start
  walk <- c(e)
  
  while (e < 1) {
    inc <- e + rnorm(1, 0.01, 0.2)
    if (inc > -0.5) {
      e <- inc
      walk <- append(walk, e)
    }
  }
  
  out <- walk[walk <= 1]
  out[length(out)] <- 1
  
  return(out)
}

set.seed(342)
yUpper <- simWienerWalk(0) * 0.5
yLower <- -1*simWienerWalk(0) * 0.5

ddmWalkUpper <- data.frame(x = 1:length(yUpper)/1.6, y = yUpper)
ddmWalkLower <- data.frame(x = 1:length(yLower)/1.6, y = yLower)

a <- 1.2
v <- 0.05
t0 <- 0.2
n <- 1000
TEXT_SIZE <- 4

densStart <- 0.1*max(length(yUpper), length(yLower))
densEnd <- 2*max(length(yUpper), length(yLower))

ddmDensUpper <- data.frame(x = seq(densStart, densEnd, length.out = n),
                           y = 0.4*dwiener(seq(t0, 2, length.out = n),
                                           alpha = a, tau = t0, beta = a/2, delta = v, resp = "upper"))
ddmDensLower <- data.frame(x = seq(densStart, densEnd, length.out = n),
                           y = -0.4*dwiener(seq(t0, 2, length.out = n),
                                            alpha = a, tau = t0, beta = a/2, delta = v, resp = "lower"))

ddmDensUpper$y <- (ddmDensUpper$y - min(ddmDensUpper$y))/(max(ddmDensUpper$y) - min(ddmDensUpper$y))
ddmDensLower$y <- (ddmDensLower$y - min(ddmDensUpper$y))/(max(ddmDensUpper$y) - min(ddmDensUpper$y)) 

ddmDensLowerMax <- min(ddmDensLower$y)

p_ddm <- ggplot() +
  geom_line(aes(x = x, y = y), ddmWalkUpper, alpha = 0.5) +
  geom_line(aes(x = x, y = y), ddmWalkLower, alpha = 0.3) +
  geom_hline(yintercept = c(1, -1)*0.5, size = 0.5) +
  annotate(geom = "segment", x = -0.1*densEnd, xend = -0.1*densEnd, y = -0.5, yend = 0.5, linetype = "longdash") +
  annotate(geom = "segment", x = 0, xend = 0, y = -0.5, yend = 0.5, linetype = "longdash") +
  geom_line(aes(x = x, y = y + 0.5), ddmDensUpper, size = 0.5) +
  geom_line(aes(x = x, y = y - 0.5), ddmDensLower, size = 0.5) +
  annotate(geom = "text", x = 0.825*densEnd, y = 1.2*0.5, label = "Right", size = TEXT_SIZE*(3/4)) +
  annotate(geom = "text", x = 0.825*densEnd, y = -1.2*0.5, label = "Left", size = TEXT_SIZE*(3/4)) +
  annotate(geom = "segment", x = 0, xend = 0.25*densEnd, y = 0, yend = 0.75*0.5,
           arrow = arrow(length = unit(7.5, "pt"))) +
  annotate(geom = "text", x = 0.11*densEnd, y = 0.6*0.5, label = "v", size = TEXT_SIZE) +
  annotate(geom = "errorbarh", xmin = -0.1*densEnd, xmax = 0, y = 1.2*0.5, height = 0.1*0.5) +
  annotate(geom = "text", x = -0.05*densEnd, y = 1.35*0.5, label = "t[0]", parse = T, size = TEXT_SIZE) +
  annotate(geom = "text", x = -0.03*densEnd, y = 0, label = "z", size = TEXT_SIZE) +
  annotate(geom = "point", x = 0, y = 0, size = 0.75) +
  annotate(geom = "segment", x = 0.9*densEnd, xend = 0.9*densEnd, y = 0.1*0.5, yend = 0.9*0.5,
           arrow = arrow(length = unit(7.5, "pt"))) +
  annotate(geom = "segment", x = 0.9*densEnd, xend = 0.9*densEnd, y = -0.1*0.5, yend = -0.9*0.5,
           arrow = arrow(length = unit(7.5, "pt"))) +
  annotate(geom = "text", x = 0.9*densEnd, y = 0.02*0.5, label = "a", size = TEXT_SIZE) +
  stat_function(fun = "dnorm", xlim = c(-0.05, 0.05)*densEnd, args = list(mean = 0, sd = 2),
                position = position_nudge(x = 0.8*densEnd, y = 0.95)) +
  theme_void() +
  theme(plot.margin = margin(10, 10, 10, 10), plot.title = element_text(hjust = 0.5))

#ggsave(p_ddm, width =  3.44, height = 2.9, file="plots/DDM.pdf")

# Trade-off visualization ------------------------------------------------------

# Function to generate positively correlated data
generate_correlated_data <- function(mean_x, mean_y, n, correlation = -0.85) {
  x <- rnorm(n, mean = mean_x, sd = 1)
  y <- correlation * x + rnorm(n, mean = mean_y, sd = sqrt(1 - correlation^2))
  return(data.frame(x = x, y = y))
}

# Generate sample data for each group with a positive correlation and spread them across the plot
set.seed(234) # for reproducibility
n <- 300
data1 <- rbind(
  generate_correlated_data(2.5, 12, n),
  generate_correlated_data(8, 14, n),
  generate_correlated_data(2, 8, n),
  generate_correlated_data(5, 16, n)
)
data1$group <- factor(rep(1:4, each = n))


# Compute group means for data1
data1_means <- data1 %>%
  group_by(group) %>%
  summarise(mean_x = mean(x), mean_y = mean(y))

# Plot data1 with group means
p1 <- ggplot(data1, aes(x = x, y = y, color = group)) +
  geom_point(size = 0.4, alpha = 0.3) + # Transparent points
  geom_point(data = data1_means, aes(x = mean_x, y = mean_y, fill = group),
             size = 3, alpha = 0.9, shape = 21, stroke = 1, color = "black") +
  scale_fill_viridis(discrete = TRUE, option = "viridis") +
  stat_ellipse() +
  stat_ellipse(aes(group = 1), level = 0.95, color = "darkred", lwd = 1.5, linetype = 11) +
  theme_minimal() +
  labs(x = expression(t0[d]),
       y =  expression(a[d]),
       color = "Group") +
  xlim(-2, 10) +
  ylim(3, 15) +
  theme_classic() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  scale_color_viridis(discrete = TRUE, option = "viridis") +
  scale_fill_viridis(discrete = TRUE, option = "viridis")

# same but overlapping
data2 <- rbind(
  generate_correlated_data(3.5, 12.3, n),
  generate_correlated_data(3.8, 12.1, n),
  generate_correlated_data(4.2, 12, n),
  generate_correlated_data(4.6, 11.9, n)
)
data2$group <- factor(rep(1:4, each = n))

# Compute group means for data2
data2_means <- data2 %>%
  group_by(group) %>%
  summarise(mean_x = mean(x), mean_y = mean(y))


p2 <- ggplot(data2, aes(x = x, y = y, color = group)) +
  geom_point(size = 0.4, alpha = 0.3) +  # Keep color legend for this layer (Individual)
  scale_fill_viridis(discrete = TRUE, option = "viridis") +
  stat_ellipse() +                          # Add ellipses
  geom_point(data = data2_means, aes(x = mean_x, y = mean_y, fill = group),
             size = 3, shape = 21, stroke = 1, color = "black", show.legend = FALSE) +  # Remove fill legend (Group)
  stat_ellipse(aes(group = 1), level = 0.95, color = "darkred", lwd = 1.5, linetype = 11) + # Add an overall ellipse
  theme_minimal() +                         # Minimal theme for better visualization
  labs(x = expression(t0[d]),
       y =  expression(a[d]),
       color = "Individual") +
  xlim(-2, 10) +
  ylim(3, 15) +
  theme_classic() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  scale_color_viridis(discrete = TRUE, option = "viridis") +
  guides(fill = "none")  # Remove fill legend (Group)




p <-  (p1 | p2)
p
ggsave(plot = p, filename = "ind_variation.pdf", path = "plots", device = cairo_pdf,
       width = 5.44, height = 2.49)


# EZ plot --------------------

compute_observed_cor <- function(res, true_params, group_diff){
  cor_observed <- numeric()
  for(n in 1:length(res)){
    var_total_adiff <- var(res[[n]]$a_difference)
    var_total_t0diff <- var(res[[n]]$t0_difference)
    var_between_adiff <- var(true_params[[n]][,"a_diff"])
    var_between_t0diff <-  var(true_params[[n]][,"t0_diff"])
    var_within_a_diff <- var_total_adiff-var_between_adiff
    var_within_t0_diff <- var_total_t0diff-var_between_t0diff
    
    eta_adiff <- (var_between_adiff/var_total_adiff)
    eta_t0diff <- (var_between_t0diff/var_total_t0diff)
    
    cor_between <- ifelse(var_between_adiff==0,0,cor(true_params[[n]][,"a_diff"], true_params[[n]][,"t0_diff"]))
    cor_within <- ifelse(group_diff[n]=="no" & any(group_diff=="yes"), cor(res[[2]]$a_difference, res[[2]]$t0_difference), cor(res[[1]]$a_difference, res[[1]]$t0_difference))
    cor_observed[n] <- (sqrt(eta_adiff)*sqrt(eta_t0diff)*cor_between) + (sqrt(1-eta_adiff)*sqrt(1-eta_t0diff)*cor_within)  
  }  
  return(cor_observed)
}

create_plotdat <- function(path,ntrials_per_condition){
  load(path)
  group_diff <- ifelse(grepl("_TRUE_", names(cor_t0diff_adiff)), "no", "yes")
  cor_observed <- compute_observed_cor(res, true_params, group_diff)
  sd_multipliers <- as.numeric(stringr::str_extract(names(cor_t0diff_adiff), "[0-9]+(\\.[0-9]+)?"))
  
  # get correlation coef
  estimates <- sapply(cor_t0diff_adiff, function(x) x$estimate)
  
  # get confidence interval
  confint <- sapply(cor_t0diff_adiff, function(x) x$conf.int)
  
  p_dat <- data.frame(sd_multipliers= as.numeric(sd_multipliers),
                      group_diff = group_diff,
                      estimates = as.numeric(estimates),
                      conf.lower = confint[1,],
                      conf.upper = confint[2,],
                      cor_reconstructed = cor_observed,
                      ntrials_per_condition=factor(ntrials_per_condition, 
                                                   levels = c(10000, 1000, 200, 100)))
  
  return(p_dat)
}

p_dat <- bind_rows(create_plotdat("ez-results/ez_output_1000.RData", ntrials_per_condition=1000),
                   create_plotdat("ez-results/ez_output_10000.RData", ntrials_per_condition=10000),
                   create_plotdat("ez-results/ez_output_200.RData", ntrials_per_condition=200),
                   create_plotdat("ez-results/ez_output_100.RData", ntrials_per_condition=100))
sd_multipliers <- seq(0, 2, by = 0.2)
dodge_width <- 0.06  
p <- p_dat %>% 
  mutate(sd_multipliers = sd_multipliers/2) %>% # for plotting / readability, 0 to 100% instead of 0 to 200%
  ggplot(aes(x = sd_multipliers, group = interaction(group_diff, ntrials_per_condition))) +
  geom_hline(yintercept = 0, colour = "lightgrey") +
  geom_line(aes(y = estimates, 
                linetype = group_diff,  
                colour = factor(ntrials_per_condition)), 
            position = position_dodge(width = dodge_width)) +  
  geom_point(aes(y = estimates), position = position_dodge(width = dodge_width)) +  
  geom_errorbar(aes(ymin = conf.lower, ymax = conf.upper), 
                alpha = 0.6, width = 0.03, 
                position = position_dodge(width = dodge_width)) +  
  geom_point(aes(y = cor_reconstructed), colour = "black", size = 3, shape = 4, position = position_dodge(width = dodge_width)) +
  scale_x_continuous(labels = scales::percent_format(scale = 100), n.breaks = 11) +
  scale_y_continuous(n.breaks = 6) +
  
  scale_color_manual(
    values = RColorBrewer::brewer.pal(length(levels(p_dat$ntrials_per_condition)), "Set1"),
    labels = setNames(
      formatC(as.numeric(levels(p_dat$ntrials_per_condition)), format = "d", big.mark = ","),
      levels(p_dat$ntrials_per_condition)  
    )
  ) +  
  
  scale_linetype_manual(
    values = c("no" = "solid", "yes" = "dashed"),  
    labels = c("no" = "No", "yes" = "Yes")  
  ) +  
  scale_shape_manual(values = c(16, 17)) +  
  theme_classic() +
  labs(
    x = expression(paste("Percentage of max Population SD ", (sigma[B]))),
    y =  expression( r[a[d] * t0[d]] ), 
    colour = "Number of trials per condition",  
    shape = "Population-level difference",  
    linetype = "Population-level difference"  
  ) +
  theme(panel.grid = element_blank())
ggsave(p, file = "plots/EZ.pdf", device = cairo_pdf,width =  7, height=3)


# Hierarchical correlation plot ------------
sd_multipliers <- seq(0.2, 2, by = 0.2)

fit_objects <- list.files("fit-objects/fully-hierarchical/")

corr_t0diff_adiff <- list()
within_cor_real <- list()
for(m in 1:length(fit_objects)){
  load(paste0("fit-objects/fully-hierarchical/", fit_objects[m]))
  corr_t0diff_adiff[[m]] <- summary(fit_emc, "correlation")$a_difficulty_d["t0_difficulty_d",]
  sample_stage <- merge_chains(fit_emc)$samples$stage == "sample"
  merged <- merge_chains(fit_emc)
  alpha <- merged$samples$alpha[,,sample_stage]
  
  within_real <- rep(NA, 250)
  t0_d_real <- rep(NA, 250)
  a_d_real <- rep(NA, 250)
  for(s in 1:250){
    a_easy <- exp(alpha["a",s,] - 0.5 * alpha["a_difficulty_d",s,])
    a_hard <- exp(alpha["a",s,] + 0.5 * alpha["a_difficulty_d",s,])
    a_d_real <- a_hard-a_easy
    t0_easy <- exp(alpha["t0",s,] - 0.5 * alpha["t0_difficulty_d",s,])
    t0_hard <- exp(alpha["t0",s,] + 0.5 * alpha["t0_difficulty_d",s,])
    t0_d_real <- t0_hard-t0_easy
    within_real[s] <- cor(t0_d_real, a_d_real)
  }
  within_cor_real[[m]] <- within_real
}

p_hierarchical_corr <- bind_rows(corr_t0diff_adiff) %>% 
  mutate(type = "group-level") %>% 
  select(-Rhat, -ESS) %>% 
  bind_rows(bind_rows(lapply(within_cor_real, quantile, probs = c(0.025, 0.5, 0.975))) %>% 
              mutate(type = "within-subject")) %>% 
  mutate(proportion_of_max_sd = rep(as.numeric(unlist(stringr::str_extract_all(fit_objects, 
                                                                               pattern = "(?<=_SDmult_)\\d+(\\.\\d+)?"))), 2)/2) %>% 
  ggplot(aes(x = proportion_of_max_sd, y = `50%`)) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = -0.743) +
  geom_point( aes(colour = type)) +
  geom_errorbar(aes(ymax = `97.5%`, ymin = `2.5%`, colour = type),alpha = 0.6, width = 0.03,) +
  scale_y_continuous(n.breaks =  6) +
  theme_classic() +
  theme(panel.grid = element_blank()) +
  scale_x_continuous(breaks = sd_multipliers/2, labels = scales::percent_format(scale = 100))+
  labs(x = expression(paste("Percentage of max Population SD ", (sigma[B]))),
       y =  expression( hat(rho)[a[d] * t0[d]] ), colour = "") +
  #ggtitle("Correlation as a function of effect variability", subtitle =
  #          expression("Variability (SD) at 100%: " * v[diff] * " = 0.5, " * a[diff] * " = 0.1, " * t0[diff] * " = 0.05")) +
  scale_color_manual(values = c("darkred", "darkgrey"), labels = c("Between-subject", "Within-subject"))
p_hierarchical_corr


ggsave(p_hierarchical_corr, file = "plots/hierarchical_standard_within.pdf", 
       device = cairo_pdf, width = 5.8, height = 2)

# Dutilh et al. plots --------------

fit_objects_full <- list.files("fit-objects/dutilh-exp1-hierarchical/", pattern = "full")
fit_objects <- list.files("fit-objects/dutilh-exp1-hierarchical/", pattern = "cleaned_nr")[1:10]

numbers_full <- as.numeric(sub(".*_nr_(\\d+)\\.RData", "\\1", fit_objects_full))
numbers <- as.numeric(sub(".*_nr_(\\d+)\\.RData", "\\1", fit_objects))

corr_t0diff_adiff_full <- list()
for(m in 1:length(fit_objects_full)){
  load(paste0("fit-objects/dutilh-exp1-hierarchical/", fit_objects_full[m]))
  corr_t0diff_adiff_full[[m]] <- summary(fit_dutilh_exp1_full, "correlation")$a_difficulty_d["t0_difficulty_d",]
}

corr_t0diff_adiff <- list()
for(m in 1:length(fit_objects)){
  load(paste0("fit-objects/dutilh-exp1-hierarchical/", fit_objects[m]))
  corr_t0diff_adiff[[m]] <- summary(fit_dutilh_exp1, "correlation")$a_difficulty_d["t0_difficulty_d",]
}



p_hierarchical_corr <- bind_rows(corr_t0diff_adiff) %>% 
  select(-Rhat, -ESS) %>% 
  mutate(Repetition = numbers,
         type = "simple") %>%
  bind_rows(bind_rows(corr_t0diff_adiff_full) %>% 
              select(-Rhat, -ESS) %>% 
              mutate(Repetition = numbers_full,
                     type = "full")) %>%
  ggplot(aes(x = Repetition, y = `50%`)) +
  geom_pointrange(aes(ymax = `97.5%`, ymin = `2.5%`, linetype=type), size = 0.2, position = position_dodge(width = 0.4)) +
  scale_x_continuous(n.breaks =  4) +
  theme_classic() +
  theme(panel.grid = element_blank()) +
  geom_hline(yintercept = 0) +
  labs(x = "",
       y =  expression( hat(rho)[Ba[d] * t0[d]] ), 
       shape = "") +
  guides(linetype = guide_legend(title = NULL)) +
  scale_linetype_manual(
    values = c("full" = "solid", "simple" = "longdash"),
    labels = c("simple" = "Simple DDM", "full" = "Full DDM")) +
  theme(axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

p_hierarchical_corr


ggsave(p_hierarchical_corr, width = 5.81, height = 1.54, device = cairo_pdf, 
       file = "plots/DutilhExp1.pdf")

# Example fit stone vs. full DDM: repetition 5 
par(mfrow= c(2,2)) # 6.60 x 6.02

load("fit-objects/dutilh-exp1-hierarchical/dutilh_exp1_cleaned_nr_5.RData")
dat_5 <- all_data %>% 
  mutate(difficulty = as.factor(ifelse(difficulty == "A", "easy", "hard")),
         R = as.factor(ifelse(R == "upper", "correct", "incorrect")))
pp <- predict(fit_dutilh_exp1, n_post = 100, n_cores = 5) %>% 
  mutate(difficulty = as.factor(ifelse(difficulty == "A", "easy", "hard")),
         R = as.factor(ifelse(R == "upper", "correct", "incorrect")))

plot_fit(pp=pp, data=,dat_5, main = "Simple DDM:")

load("fit-objects/dutilh-exp1-hierarchical/full_dutilh_exp1_cleaned_nr_5.RData")
dat_5 <- all_data %>% 
  mutate(difficulty = as.factor(ifelse(difficulty == "A", "easy", "hard")),
         R = as.factor(ifelse(R == "upper", "correct", "incorrect")))
pp <- predict(fit_dutilh_exp1_full, n_post = 100, n_cores = 5) %>% 
  mutate(difficulty = as.factor(ifelse(difficulty == "A", "easy", "hard")),
         R = as.factor(ifelse(R == "upper", "correct", "incorrect")))

plot_fit(pp=pp, data=,dat_5, main = "Full DDM:")
