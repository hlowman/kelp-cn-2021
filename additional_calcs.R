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

# Calculating "starting" C:N values for discussion of relative increases.
# Mean C:N for the first 12 months of the dataset
cn_12mo <- cn_full[1:36,] %>% # pull out first 12 months of data
  mutate(logcn = log10(cn)) # add log(C:N) column

mean_cn_start <- mean(cn_12mo$cn, na.rm = TRUE) # calculate mean of C:N
# 12.71
mean_logcn_start <- mean(cn_12mo$logcn, na.rm = TRUE) # calculate mean of log(C:N)
# 1.08
slope_logcn <- 0.0097
final_logcn <- mean_logcn_start + (slope_logcn*19)
# 1.26
final_cn <- 10^(final_logcn)
# 18.28
delta_cn <- (final_cn - mean_cn_start)/mean_cn_start
# 0.44 = 44%
yearly_delta_cn <- (final_cn - mean_cn_start)/19
# 0.29

# Mean %c for the first 12 months of the dataset
mean_c_start <- mean(cn_12mo$c, na.rm = TRUE) # calculate mean of %C
# 28.39
slope_c <- 0.27
final_c <- mean_c_start + (slope_c*19)
# 33.52
delta_c <- (final_c - mean_c_start)/mean_c_start
# 0.18 = 18%

# Mean %n for the first 12 months of the dataset
mean_n_start <- mean(cn_12mo$n, na.rm = TRUE) # calculate mean of %C
# 2.48
slope_n <- -0.024
final_n <- mean_n_start + (slope_n*19)
# 2.02
delta_n <- (final_n - mean_n_start)/mean_n_start
# -0.18 = 18%

# some additional protein-related calculations
mean_prot_start <- mean_n_start*5
# 12.38
final_prot <- final_n*5
# 10.10
delta_prot <- (final_prot - mean_prot_start)/mean_prot_start
# -0.18 = 18%

# End of script.
