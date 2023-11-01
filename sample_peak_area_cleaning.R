library(lubridate)
library(tidyverse)
rm(list = ls())

# Reads in peak areas and adds tube names and sample log info. Drop weird
# compounds.
get_info <- function(runtype) {
  filenames <- list.files("Z:/CUNY Samples/_HDFs")
  runtype_names <- filenames %>% 
    str_subset(runtype)
  
  dat <- read_csv(str_c("data/", runtype, "_peakareas.csv"))
  info <- read_csv(str_c("info/", runtype, " Sample Log.csv"))
  
  dat_info <- dat %>% 
    select(!c(isopropyl_acetate, pentyl_acetate)) %>% 
    rename_with(~ str_c("area_", .), !DOY) %>% 
    mutate(
      filename = all_of(runtype_names),
      QBTX = str_match(filename, "QBTX([^_\\-+]*)")[,2]
    ) %>% 
    left_join(info, by = "QBTX") %>% 
    filter(is.na(Notes) | str_detect(Notes, "blank", negate = TRUE))
  return(dat_info)
}

# Normalize peak areas and find ratios of each compound's peak area to benzene's
# peak area in each sample
get_ratio <- function(info) {
  ratio <- info %>% 
    mutate(
      across(starts_with("area"), ~ replace(., is.na(.), 0)),
      across(starts_with("area"), ~ . / hours, 
             .names = "norm_{str_sub({.col}, 6)}"),
      across(starts_with("area"), ~ . / area_benzene,
             .names = "ratio_{str_sub({.col}, 6)}")
    ) %>% 
    relocate(starts_with("norm"), .after = notes) %>% 
    relocate(starts_with("ratio"), .after = notes)
  return(ratio)
}

runtypes <- set_names(c("CUNY", "Shutdown"))
infos <- map(runtypes, get_info)

# Clean up datetime columns and select in correct order
infos[["CUNY"]] <- infos[["CUNY"]] %>% 
  mutate(
    start_datetime = mdy_hms(str_c(`Start Date`, `Start Time (EST)`, sep = " ")),
    end_datetime = mdy_hms(str_c(`End Date`, `End Time (EST)`, sep = " ")),
    hours = `Length (hours)` %>% as.numeric(),
    flow_rate_sccm = `Flow Rate (sccm)`,
    notes = str_c(
      "Line: ", Line, ". ",
      "Big Pump Flow Rate (slpm): ", `Flow Rate Big Pump (slpm)`, ". ", Notes
    )
  ) %>% 
  select(QBTX, start_datetime:notes, starts_with("area"), filename)

infos[["Shutdown"]] <- infos[["Shutdown"]] %>% 
  mutate(
    start_datetime = mdy_hms(str_c(Date, `Start time (LT)`, sep = " ")),
    hours = time_length(`End time (LT)` - `Start time (LT)`, unit = "hour") %>% 
      map_dbl(~ . %% 24),
    end_datetime = start_datetime + dhours(hours),
    flow_rate_sccm = str_sub(`Flow rate`, end = 3) %>% as.numeric() * 1000,
    notes = Notes
  ) %>% 
  select(QBTX, start_datetime, end_datetime, hours, flow_rate_sccm:notes,
         starts_with("area"), filename)

ratios <- map(infos, get_ratio) %>% 
  bind_rows(.id = "runtype") %>% 
  arrange(start_datetime)
write_csv(ratios, "data/sample_ratios.csv")
