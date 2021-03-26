# Seasonal Kendall Test Script
# Heili Lowman
# March 24, 2021

# The following script will run Seasonal Kendall tests to examine for a monotonic trend in C:N values.

# Load packages.
library(tidyverse)
library(lubridate)
library(naniar)
library(EnvStats)

# Load datasets from "data_tidying.R".
load("data_tidy/kelp_cn_data_clean.rda")

# This workflow has been adapted from Marcus Beck's tutorial located at https://fawda123.github.io/swmp_workshop_2016/training_modules/module4_kendall/kendall.pdf

# For this analysis, summarize across all three sites.
# cn_summ <- cn_full %>%
#   group_by(YEAR, MONTH) %>%
#   summarize(meanc = mean(c, na.rm=TRUE),
#             meann = mean(n, na.rm=TRUE),
#             meancn = mean(cn, na.rm=TRUE)) %>%
#   ungroup()
# Upon further reading of function documentation, "there may be more than one observation within the same year within the same season." - https://www.rdocumentation.org/packages/EnvStats/versions/2.3.1/topics/kendallSeasonalTrendTest

# Before anything else, examine the distribution of values

hist(cn_full$c)
hist(cn_full$n)
hist(cn_full$cn)

# Log transform C:N values to follow advice of Isles 2020

cn_log <- cn_full %>%
  mutate(logcn = log10(cn)) # add log(CN) column

hist(cn_log$logcn)

# Run Seasonal Kendall Tests (due to known seasonal variability)
# nonparametric test for detecting a monotonic trend

# (1) log(C:N)
# function requires pulling from a numerical-only dataframe
cn_data <- cn_log %>%
  select(YEAR, MONTH, logcn)

cn_sk <- kendallSeasonalTrendTest(logcn ~ MONTH + YEAR, data = cn_data)

#Null Hypothesis:                 All 12 values of tau = 0
#
#Alternative Hypothesis:          The seasonal taus are not all equal
#                                 (Chi-Square Heterogeneity Test)
#                                 At least one seasonal tau != 0
#                                 and all non-zero tau's have the
#                                 same sign (z Trend Test)

# View results
cn_sk
# Chi-Square (Het) = 9.7155, z (Trend) = 9.4675, df = 11

cn_sk$estimate
#           tau          slope      intercept 
#   0.255779511    0.009707598  -17.471535315 
# tau is the direction/magnitude of the trend (between -1 and 1)
# slope is the rate of change in the interval

cn_sk$p.value
# Chi-Square (Het)        z (Trend) 
# 5.561451e-01            2.865211e-21  
# Heterogeneity p value is above 0.05, so there were no opposing trends in any season.
# Trend p value falls well below 0.05, so there is a positive, monotonic trend in the C:N data.

cn_sk$seasonal.estimates

# Let's examine each season with a non-seasonal test
for(i in 1:12){
  
  new_data <- cn_data %>%
    filter(MONTH == i) %>%
    na.omit()
  
  ktt <- kendallTrendTest(logcn ~ YEAR, data = new_data)
  
  print(ktt)
}

# April and October are the only months with no significant trend in log(C:N) values.

# (2) C
# function requires pulling from a numerical-only dataframe
c_data <- cn_log %>%
  select(YEAR, MONTH, c)

c_sk <- kendallSeasonalTrendTest(c ~ MONTH + YEAR, data = c_data)

# View results
c_sk
# Chi-Square (Het) = 10.2149, z (Trend) = 9.0174, df = 11

c_sk$estimate
#           tau          slope     intercept 
#     0.2468106      0.2748313  -574.1016439 
# tau is the direction/magnitude of the trend (between -1 and 1)
# slope is the rate of change in the interval

c_sk$p.value
# Chi-Square (Het)        z (Trend) 
# 5.111686e-01            1.926636e-19
# Heterogeneity p value is above 0.05, so there were no opposing trends in any season.
# Trend p value falls well below 0.05, so there is a positive, monotonic trend in the C data.

c_sk$seasonal.estimates

# Let's examine each season with a non-seasonal test
for(i in 1:12){
  
  new_data <- c_data %>%
    filter(MONTH == i) %>%
    na.omit()
  
  ktt <- kendallTrendTest(c ~ YEAR, data = new_data)
  
  print(ktt)
}

# September is the only month with no significant trend in C values.

# (3) N
# function requires pulling from a numerical-only dataframe
n_data <- cn_log %>%
  select(YEAR, MONTH, n)

n_sk <- kendallSeasonalTrendTest(n ~ MONTH + YEAR, data = n_data)

# View results
n_sk
# Chi-Square (Het) = 20.650, z (Trend) = -5.239, df = 11

n_sk$estimate
#           tau          slope     intercept 
#   -0.13846400    -0.02388889   44.14444592 
# tau is the direction/magnitude of the trend (between -1 and 1)
# slope is the rate of change in the interval

n_sk$p.value
# Chi-Square (Het)        z (Trend) 
# 3.718061e-02            1.614224e-07
# Heterogeneity p value is below 0.05, so there ARE opposing trends in a season.
# Trend p value falls well below 0.05, so there is a negative, monotonic trend in the N data.

n_sk$seasonal.estimates
# Appears that opposing season is April - is strong upwelling counteracting the overall upwards trend?

# Let's examine each season with a non-seasonal test
for(i in 1:12){
  
  new_data <- n_data %>%
    filter(MONTH == i) %>%
    na.omit()
  
  ktt <- kendallTrendTest(n ~ YEAR, data = new_data)
  
  print(ktt)
}

# February, March, April, May, June, October, and December are months with no significant trend in N values.

# Sample results text from : https://open.alberta.ca/dataset/8698d43c-d504-412b-a686-39726710feb7/resource/1df46420-b2ad-40a0-90d9-605018a09eb1/download/larp-analysis-of-water-quality-conditions-final-report-march-2015.pdf
# "The seasonal Mann-Kendall test found a highly significant increasing trend in total nitrogen (pvalue<0.001), with a magnitude of 0.005 mg/L/year and a 95% confidence interval of 0.002-0.008 mg/L/year. The test for heterogeneity between seasons was not significant (p-value=0.5), indicating no opposing seasonal trends."

# End of script.
