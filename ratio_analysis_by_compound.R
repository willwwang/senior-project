library(tidyverse)
rm(list = ls())

# Plot boxplot of peak area ratios before and after the shutdown 
plot_cpd_ratios <- function(selected_ratios, cpd) {
  cpd_ratios <- selected_ratios %>% 
    filter(compound == cpd)
  ratio_plot <- ggplot(data = cpd_ratios, 
                       aes(x = runtype, y = ratio, col = runtype)) +
    geom_boxplot() +
    labs(title = str_c("Ratio of ", cpd, " to\nbenzene"),
         x = "Time", y = "Ratio of peak area to benzene") +
    guides(col = "none") +
    theme_bw()
  return(ratio_plot)
}

ratios <- read_csv("data/sample_ratios.csv")

cpds <- c("styrene", "butyl_acetate", "limonene", "decamethylcyclopentasiloxane",
          "pentane_2_2_4_trimethyl", "methyl_methacrylate")
names(cpds) <- c("styrene", "butyl acetate", "limonene", "siloxane D5",
                  "isooctane", "methylmethacrylate")

# Select the compounds of interest and pivot the data for ggplot
selected_ratios <- ratios %>%
  slice_min(norm_benzene, n = nrow(ratios) - 4) %>%
  select(runtype, contains(str_c("ratio_", cpds))) %>% 
  pivot_longer(
    cols = contains("ratio"),
    names_to = "compound",
    names_prefix = "ratio_",
    values_to = "ratio"
  ) %>% 
  mutate(compound = factor(compound, levels = cpds) %>% 
           fct_recode(!!!cpds),
         runtype = factor(runtype) %>%
           fct_recode("Pre-shutdown" = "CUNY", "Post-shutdown" = "Shutdown"))

# Plot boxplots of peak area ratios, before and after the shutdown
ggplot(data = selected_ratios, aes(x = compound, y = ratio, col = runtype)) +
  geom_point(alpha = 0.7, position = position_jitterdodge()) +
  coord_cartesian(ylim = c(0, 0.6)) +
  labs(#title = "Selected VOC / benzene peak area ratios, pre- and post-shutdown",
       x = "Compound", y = "Ratio of peak area to benzene", col = "Time") +
  theme_bw() +
  theme(legend.position = "bottom") +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

ggsave("plots/selected_compound_ratios.png", width = 8, height = 4)

# Create separate plots for butyl acetate and siloxane D5
cpds_interest <- cpds[2:3]
ratio_plots <- map(names(cpds_interest), ~ plot_cpd_ratios(selected_ratios, .))
map2(ratio_plots, cpds_interest,
     ~ ggsave(plot = .x, filename = str_c("plots/", .y, "_ratio.png"),
              width = 3, height = 4.5))
