library(lubridate)
library(randomizr)
library(tidyverse)
rm(list = ls())
set.seed(42)

# Combine response factors and total sample volume to get concentrations
get_concs <- function(ratios, std_cpds) {
  cpds <- std_cpds$compound
  rfs <- set_names(std_cpds$RF, cpds)
  re <- str_c("area_(", str_c(cpds, collapse = "|"), ")$")
  concs <- ratios %>% 
    select(runtype, QBTX, start_datetime, end_datetime, hours, moles_air,
           matches(re)) %>%
    mutate(across(
      matches(re), 
      ~ . / rfs[[str_sub(cur_column(), 6)]] / moles_air, 
      .names = "conc_{str_sub({.col}, 6)}")
    ) %>% 
    select(!c(moles_air, starts_with("area")))
  return(concs)
}

# Get the average concentration of a compound in a time period
avg_conc_in_span <- function(conc, datetimes, start, end) {
  included <- which((datetimes >= start) & (datetimes < end))
  return(mean(conc[included], na.omit = TRUE))
}

# Average PTR concentrations during each sample and add them to the data
add_ptr_concs <- function(concs, ptr) {
  all_concs <- concs %>% 
    select(!conc_furfural) %>% 
    filter(runtype == "CUNY") %>% 
    mutate(
      across(starts_with("conc"),
             ~ map2_dbl(start_datetime, end_datetime,
                        ~ avg_conc_in_span(
                          ptr[[str_to_title(str_sub(cur_column(), 6))]],
                          ptr$start_datetime, .x, .y)),
             .names = "ptr_{.col}")
    ) %>% 
    pivot_longer(
      cols = contains("conc"),
      names_to = c("source", "compound"),
      names_pattern = "(.*)_(.*)",
      values_to = "concentration"
    ) %>% 
    pivot_wider(
      names_from = "source",
      values_from = "concentration"
    )
  return(all_concs)
}

# Plot concentrations over time for a given compound
plot_conc <- function(concs, cpd) {
  conc_plot <- ggplot(
    data = concs, aes(x = start_datetime, y = .data[[str_c("conc_", cpd)]])
  ) +
    geom_point() +
    geom_smooth(method = "lm", formula = "y ~ 1", se = FALSE, lty = 2,
                col = "black", size = 0.5) +
    facet_grid(
      cols = vars(runtype), scales = "free_x", space = "free_x",
      labeller = labeller(
        runtype = c("CUNY" = "Pre-shutdown", "Shutdown" = "Post-shutdown")
      )) +
    scale_x_datetime(date_breaks = "1 week", date_labels = "%b %d") +
    labs(# title = str_c(str_to_title(cpd), " concentrations, early 2020"),
         x = "Time", y = str_c("[", cpd, "] (ppb)")) +
    theme_bw()
  return(conc_plot)
}

# Returns difference in average concentrations in each period for each compoound
get_conc_mean_diffs <- function(concs) {
  avg_concs <- concs %>% 
    summarize(
      across(
        starts_with("conc"),
        ~ {mean(.[runtype == "Shutdown"]) - mean(.[runtype == "CUNY"])}
      ),
    )
  return(avg_concs)
}

# Read in data
std_cpds <- read_csv("data/GasStd_summary.csv")
ratios <- read_csv("data/sample_ratios.csv")

# Get sample moles of air using flow rate and time info
ratios <- ratios %>% 
  mutate(moles_air = flow_rate_sccm * 60 * hours / 10^6 * 40.89)

# Get concentrations from areas and mark 4 high-benzene samples
concs <- get_concs(ratios, std_cpds) %>% 
  mutate(
    high_benzene = conc_benzene > conc_benzene[order(conc_benzene, decreasing = TRUE)[5]]
  )
write_csv(concs, "data/concentrations.csv")

# Read in PTR data
ptr <- read_csv("data/VOC PTR data Oct2021 version_update.csv")
ptr <- ptr %>% 
  mutate(start_datetime = mdy_hm(Timestamp))

all_concs <- add_ptr_concs(concs, ptr)

# Plot the calculated sample concentrations against the PTR concentrations taken
# at the same time
ggplot(data = all_concs %>% filter(compound != "naphthalene"), aes(x = conc, y = ptr_conc)) +
  geom_point(alpha = 0.5, size = 2.5, stroke = 0) +
  geom_abline(intercept = 0, slope = 1, lty = 2) +
  facet_wrap(vars(compound)) +
  labs(# title = "Sample and PTR concentrations, pre-shutdown",
       x = "Sample concentration (ppb)", y = "PTR-MS concentration (ppb)") +
  coord_cartesian(xlim = c(0, 2.5), ylim = c(0, 2.5)) +
  theme_bw() +
  theme(aspect.ratio = 1) 

ggsave("plots/sample_ptr_concs.png", width = 7, height = 3.5)

# Plot concentration time series for benzene, toluene, and furfural
cpds <- set_names(c("benzene", "toluene", "furfural"))
conc_plots <- map(cpds, ~ plot_conc(concs, .))
map2(conc_plots, cpds, ~ ggsave(plot = .x,
                                filename = str_c("plots/", .y, "_concs.png"),
                                width = 6, height = 2))
concs %>% 
  filter(high_benzene == FALSE) %>% 
  group_by(runtype) %>% 
  summarize(across(starts_with("conc"), ~ mean(.)))

################################## INFERENCE ##################################

# Read in concs, filter, and get observed differences
concs <- read_csv("data/concentrations.csv")
concs <- concs %>% 
  filter(high_benzene == FALSE)
mean_conc_diffs <- get_conc_mean_diffs(concs)

# Prepare for and execute randomization of concentrations across periods
N_samples <- nrow(concs)
N_shutdown <- nrow(concs %>% filter(runtype == "Shutdown"))
N_iter <- 10000
sim_diffs <- matrix(0, nrow = N_iter, ncol = ncol(mean_conc_diffs),
                    dimnames = list(NULL, names(mean_conc_diffs)))
for (i in 1:N_iter) {
  sim_concs <- concs %>% 
    mutate(
      runtype = complete_ra(N = N_samples, m = N_shutdown,
                            conditions = c("CUNY", "Shutdown"))
    )
  sim_diffs[i,] <- get_conc_mean_diffs(sim_concs) %>% as.numeric()
}

# Get p-values. One-sided for benzene and toluene, two-sided for furfural
p_benzene <- mean(sim_diffs[,"conc_benzene"] <= mean_conc_diffs[["conc_benzene"]])
p_toluene <- mean(sim_diffs[,"conc_toluene"] <= mean_conc_diffs[["conc_toluene"]])
p_furfural <- mean(
  abs(sim_diffs[,"conc_furfural"]) >= abs(mean_conc_diffs[["conc_furfural"]])
)
                     
