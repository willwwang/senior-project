library(lubridate)
library(tidyverse)
rm(list = ls())

# Reads in gas standard areas, selects compounds, and adds date and tube
# information. 
get_areas <- function(runtype) {
  dat <- read_csv(str_c("data/", runtype, "_peakareas.csv"))
  filenames <- read_csv(str_c("info/", runtype, "_filenames.csv"))
  cpds <- c("benzene", "furfural", "toluene", "naphthalene")

  dat <- dat %>% 
    mutate(
      filename = filenames$filename,
      date = mdy(str_extract(filename, "0[0-9]+")),
      tube = str_extract(filename, "QBTX[^ _\\-+]*")
    ) %>% 
    filter(!str_detect(filename, "blank|Cleaning")) %>% 
    select(date, tube, filename, all_of(cpds)) %>% 
    arrange(date)
  return(dat)
}

std_areas <- get_areas("GasStd")
write_csv(std_areas, "data/GasStd_peakareainfo.csv")

# Pivot the area for ggplot
std_areas_pivoted <- std_areas %>% 
  pivot_longer(
    cols = benzene:naphthalene,
    names_to = "compound",
    values_to = "area"
  )

# Get response factors using gas cylinder concentrations
cpd_summary <- std_areas_pivoted %>% 
  group_by(compound) %>% 
  summarize(mean_area = mean(area)) %>% 
  mutate(
    conc = c(205, 216, 99.7, 201), # from standard cylinder components file
    nmole = conc * 40.89 * 150 / 10^6,
    RF = mean_area / nmole
  )
write_csv(cpd_summary, "data/GasStd_summary.csv")

# Plot time series of standard peak areas over analysis period
ggplot(data = std_areas_pivoted, aes(x = date, y = area)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, col = "black", size = 0.5, lty = 2) + 
  facet_wrap(vars(compound), nrow = 1) +
  labs(#title = "Gas standard peak areas, January - February 2021",
       x = "Date of standard run", 
       y = "Peak area") +
  theme_bw()

ggsave("plots/gas_standard_areas.png", width = 9, height = 4)
