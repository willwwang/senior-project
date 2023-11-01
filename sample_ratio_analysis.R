library(randomizr)
library(tidyverse)
rm(list = ls())
set.seed(42)

# Find the (log) average ratios for each compound in each runtype. Joins
# compound type data. Pivot the data for ggplot
average_ratios <- function(ratios) {
  averaged_ratios <- ratios %>%
    group_by(runtype) %>% 
    summarize(across(starts_with("ratio"), mean)) %>% 
    pivot_longer(
      cols = !runtype,
      names_prefix = "ratio_",
      names_to = "compound",
      values_to = "ratio"
    ) %>% 
    pivot_wider(
      names_from = "runtype",
      names_prefix = "ratio_",
      values_from = "ratio",
    ) %>% 
    mutate(
      pct_change = (ratio_Shutdown / ratio_CUNY - 1) * 100,
      across(starts_with("ratio"), log10, .names = "log{.col}")
    )
  return(averaged_ratios)
}

# Get the overall percent change in average ratios over all compounds
get_mean_percent_change <- function(averaged_ratios) {
  mean_change <- averaged_ratios %>% 
    filter(compound != "benzene") %>% 
    summarize(mean_change = mean(pct_change)) %>% 
    as.numeric()
  return(mean_change)
}

ratios <- read_csv("data/sample_ratios.csv")
# Filter out days with extremely high benzene concentrations (== biomass burning)
ratios <- ratios %>% 
  slice_min(order_by = norm_benzene, n = nrow(ratios) - 4)

# Get averaged ratios
cpd_classes <- read_csv("info/compound_classes.csv")
averaged_ratios <- ratios %>%
  average_ratios() %>% 
  left_join(cpd_classes, by = "compound")
write_csv(averaged_ratios, "data/compound_avg_ratios.csv")

# Plot pre- and post-shutdown ratios against each other
ggplot(data = averaged_ratios, 
       aes(x = logratio_CUNY, y = logratio_Shutdown, col = type)) +
  geom_point(alpha = 0.65, size = 3) +
  labs(# title = "Mean VOC / benzene peak area ratios, pre- and post-shutdown",
       x = bquote('Pre-shutdown log(' *PA[VOC] ~ '/' ~ PA[benzene]* ')'),
       y = bquote('Post-shutdown log(' *PA[VOC] ~ '/' ~ PA[benzene]* ')'),
       col = "Compound type") +
  geom_abline(intercept = 0, slope = 1, linetype = "c5") +
  geom_abline(intercept = log10(2), slope = 1, lty = "44") +
  geom_abline(intercept = -log10(2), slope = 1, lty = "44") +
  coord_cartesian(xlim = c(-3, 0.25), ylim = c(-3, 0.25)) +
  scale_color_viridis_d() +
  theme_bw() +
  theme(aspect.ratio = 1, legend.position = "bottom")

ggsave("plots/pre_post_ratios.png", width = 6, height = 6)

# Plot distribution of percent changes from pre- to post-shutdown
ggplot(data = averaged_ratios, aes(x = pct_change, fill = type)) +
  geom_bar(position = "stack") +
  labs(# title = "Changes in VOC/benzene peak area ratios after the shutdown",
       x = "Percent change", y = "Count") +
  scale_x_binned() +
  scale_fill_viridis_d() +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave("plots/percent_change_histogram.png", width = 6, height = 3)

averaged_ratios %>% 
  group_by(type) %>% 
  summarize(avg_change = mean(pct_change))

################################## INFERENCE ##################################
# Get average percent change
mean_change <- averaged_ratios %>% 
  get_mean_percent_change()

# Prepare for and execute randomization of ratio measurements
N_samples <- nrow(ratios)
N_shutdown <- nrow(ratios %>% filter(runtype == "Shutdown"))
N_iter <- 10000
N_cpds <- ncol(ratios %>% select(starts_with("ratio")))
sim_mean_changes <- numeric(N_iter)

for (i in 1:N_iter) {
  sim_mean_changes[i] <- ratios %>% 
    mutate(
      runtype = complete_ra(N = N_samples, m = N_shutdown,
                            conditions = c("CUNY", "Shutdown"))
    ) %>% 
    average_ratios() %>% 
    get_mean_percent_change()
}
mean(sim_mean_changes >= mean_change)

ratios %>% average_ratios() %>% arrange(pct_change)
