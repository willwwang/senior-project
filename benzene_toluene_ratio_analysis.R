library(tidyverse)
rm(list = ls())

# Read in concentrations and mark samples with high benzene concentrations
concs <- read_csv("data/concentrations.csv")


# Plot benzene and toluene concentrations in each sample
ggplot(data = concs, aes(x = conc_toluene, y = conc_benzene, 
                          col = high_benzene)) +
  geom_point(alpha = 0.65, size = 3) +
  geom_abline(intercept = 0, slope = 1, lty = 2) +
  facet_wrap(vars(runtype)) +
  guides(col = "none") +
  labs(# title = "Sample benzene and toluene concentrations, pre- and post-shutdown",
       x = "[toluene] (ppb)", y = "[benzene] (ppb)") +
       # caption = "Samples in blue were removed from all analysis because their benzene:toluene ratios reflected biomass burning.") +
  coord_cartesian(xlim = c(0, 2.25), ylim = c(0, 2.25)) +
  theme_bw() +
  theme(aspect.ratio = 1)

ggsave("plots/benzene_toluene_ratios.png", width = 6.5, height = 3)

# Read in peak areas and mark samples with high benzene areas
ratios <- read_csv("data/sample_ratios.csv")
ratios <- ratios %>% 
  mutate(
    high_benzene = norm_benzene > norm_benzene[order(norm_benzene, decreasing = TRUE)[5]]
  )

# Plot normalized benzene and isooctane peak areas in each sample
ggplot(data = ratios, aes(x = norm_pentane_2_2_4_trimethyl, y = norm_benzene, 
                               col = high_benzene)) +
  geom_point(alpha = 0.65, size = 3) +
  facet_wrap(vars(runtype)) +
  guides(col = "none") +
  labs(title = "Sample peak areas for benzene and isooctane, pre- and post-shutdown",
       subtitle = "These graphs aren't 1-to-1!",
       x = "Isooctane peak areas", y = "Benzene peak areas",
       caption = "Peak areas have been normalized for sample time.") +
  theme_bw()

ggsave("plots/benzene_isooctane_ratios.png", width = 7, height = 4)
