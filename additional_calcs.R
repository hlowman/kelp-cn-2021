# Additional Calculations Script
# Heili Lowman
# March 26, 2021

# The following script will run any remaining calculations for figures printed in the manuscript that are not covered by the remaining scripts in this project.

#### Setup ####

# Load packages.
library(tidyverse)
library(lubridate)
library(naniar)

# Load datasets from "data_tidying.R".
load("data_tidy/kelp_cn_data_clean.rda")

# Calculate mean annual C:N values.
cn_annual <- cn_full %>%
  group_by(YEAR) %>%
  summarize(annual_cn = mean(cn, na.rm=TRUE)) %>%
  ungroup()

# Rough plot of overall values
cn_date <- cn_full %>%
  mutate(date = mdy(DATE))

ggplot(cn_date, aes(x = date, y = cn)) +
  geom_point() +
  geom_hline(yintercept = 20)

# Examining the CN = 20 cutoff
cn_20 <- cn_date %>%
  filter(cn > 20)
