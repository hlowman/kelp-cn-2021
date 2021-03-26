# Linear Mixed Effects Models Script
# Heili Lowman
# March 26, 2021

# The following script will run a series of linear mixed effects models to examine for relationships between oceanographic indices and C:N values. These are an update from previously run linear models with sampling site and date in as a random effect to account for repeated sampling design and non-independent samples.

# Model creation begins with fixed effects and random effects using a random intercept structure. Then, model selection follows the protocol outlined by Zuur et al. (2009, Chapter 5), beginning with a linear model, accounting for variance structure, optimizing the fixed structure, and validating best model fit using distribution of residuals and AIC values. All 24 of the below models follow the same format.

#### Setup ####

# Load packages.
library(tidyverse)
library(lubridate)
library(naniar)
library(nlme)
library(multcomp)
library(emmeans)
library(GGally)
library(groupdata2)

# Load datasets from "data_tidying.R".
load("data_tidy/kelp_cn_data_clean.rda")

# Join together and tidy data.

# Summarize SST by month at each site, since monthly resolution is as fine as the C:N data gets.
sst_monthly <- sst_full %>%
  group_by(year, month, site) %>%
  summarize(temp_C_m = mean(temp_C)) %>%
  ungroup()

# Add log(C:N) to dataset
cn_ed <- cn_full %>%
  mutate(logcn = log10(cn))

# Put together SST and CN datasets
cn_sst <- cn_ed %>%
  left_join(sst_monthly, by = c("YEAR" = "year", "MONTH" = "month", "SITE" = "site")) %>%
  group(n = 3, method = "greedy") %>% # adds grouping column
  mutate(sampling = as.numeric(`.groups`)) %>%
  mutate(sitef = factor(SITE, levels = c("ABUR", "AQUE", "MOHK"))) %>% # factor-ize the sampling sites
  ungroup()

# Put together indices and CN datasets
cn_oi <- cn_ed %>%
  left_join(oi_full, by = c("YEAR" = "year", "MONTH" = "month")) %>%
  group(n = 3, method = "greedy") %>% # adds grouping column
  mutate(sampling = as.numeric(`.groups`)) %>%
  mutate(sitef = factor(SITE, levels = c("ABUR", "AQUE", "MOHK"))) %>% # factor-ize the sampling sites
  ungroup()

#### Sea Surface Temperature ####

# LMEM #1: C:N vs. SST

# Examine data
plot(cn ~ temp_C_m, data = cn_sst)
plot(logcn ~ temp_C_m, data = cn_sst) # transformed data also looks better, aside from being statistically more robust
boxplot(logcn ~ sitef, data = cn_sst)
boxplot(logcn ~ sampling, data = cn_sst) # will be adding this in as a crossed-fixed effect since the SKT informed us that it DOES drive a monotonic upwards trend
hist(cn_sst$logcn)
hist(cn_sst$temp_C_m)

#### STEP 1: Create a linear regression and check residuals.
# Based on boxplots & pairplots : 
# Response variable - log(C:N)
# Explanatory (fixed) variables - Sea Surface Temperature*Time (aka "sampling")
a1 <- lm(logcn ~ temp_C_m*sampling, data = cn_sst) # Creates the initial linear model.

ra1 <- data.frame(rstandard(a1)) # Assigns standardized residuals to ra1.
# Because there's missing data, here's a workaround to compare residuals to data.
ra1 <- ra1 %>%
  rownames_to_column("record")
ra1_ed <- ra1 %>%
  mutate(rn = as.numeric(record))
cn_sst_ed <- cn_sst %>%
  mutate(RECORD = seq(1,666))
cn_sst_res <- cn_sst_ed %>%
  left_join(ra1_ed, by = c("RECORD" = "rn"))

ggplot(data = cn_sst_res, aes(x = temp_C_m, y = rstandard.a1.)) +
  geom_point() +
  labs(x = "Temperature", y = "Standardised residuals") # Plots said residuals.

ggplot(data = cn_sst_res, aes(x = sampling, y = rstandard.a1.)) +
  geom_point() +
  labs(x = "Sampling Date", y = "Standardised residuals")

#### STEP 2: Fit the lm() with GLS and compare to lme().
# After all that work above, it seems we cannot proceed with missing data in the dataset so...
cn_sst_rmna <- cn_sst %>%
  drop_na(logcn)

a2 <- gls(logcn ~ temp_C_m*sampling, data = cn_sst_rmna) # Effectively a linear regression - no other calls.
a3 <- lme(logcn ~ temp_C_m*sampling, random =~1 | sitef, data = cn_sst_rmna) # Creates the first LMEM with random term.
anova(a2, a3) # Compares the two models. a3 preferred with AIC value of -768.0379.

# STEP 3: Decide on a variance structure (aka random terms).
plot(a3, col=1) # Check the residuals before jumping right to applying a variance transformation.
qqnorm(a3) # This actually looks pretty good - so I'm going to skip the added variance structure. Keep in mind, log-transforming from the start is more powerful than adding in a variance component on the back end.

# STEP 4: Fit the lme().

# Using a3 <- lme(logcn ~ temp_C_m*sampling, random =~1 | sitef, data = cn_sst_rmna)

# STEP 5: Compare the lm() and lme().

# See Step 2.

# STEP 6: Everything ok? Check residuals.

# See Step 3.

# STEP 7/8: Step-wise Optimal Fixed Structure

a3_ml <- lme(logcn ~ temp_C_m*sampling, 
             random =~1 | sitef,
             method = "ML", 
             data = cn_sst_rmna) # Switch over to ML for fixed component editing portion.

a4 <- lme(logcn ~ temp_C_m, 
          random =~1 | sitef,
          method = "ML", 
          data = cn_sst_rmna) # Remove "sampling" as a fixed factor.

a5 <- lme(logcn ~ sampling, 
          random =~1 | sitef,
          method = "ML", 
          data = cn_sst_rmna) # Remove "temp_C_m" as a fixed factor.

anova(a3_ml, a4, a5) # Compare the three models. a3_ml preferred with AIC value of -821.3963, so remove nothing.

# STEP 9: Refit with REML

afinal <- lme(logcn ~ temp_C_m*sampling, 
              random =~1 | sitef,
              method = "REML", 
              data = cn_sst_rmna)

# Output of the model.
summary(afinal)

# Checking residuals.
plot(afinal, col=1) # No pattern.
qqnorm(afinal) # Looks pretty good.

# Final results.
anova(afinal)

# STEP 10: What does this mean in WORDS?

# I applied a linear mixed effect modelling approach because the data were collected at each site (3) on multiple occasions (200+). My model suggests there is a significant effect of sea surface temperature on log(C:N) values of kelp tissue samples as well as a significant effect of sampling date on values; while the interactive effect is not significant, I have left the structure as such, because we know sea surface temperature is changing through time. Random intercepts by site were added.

# Equation: log(C:N) = 0.371 + 0.048[temp] + 0.001[sampling] - 0.00003[temp:sampling] + random

### Bakun ####

# LMEM #2: C:N vs. Bakun

# Examine data
plot(logcn ~ bakun, data = cn_oi)
hist(cn_oi$bakun)

#### STEP 1: Create a linear regression and check residuals.
# Based on boxplots & pairplots : 
# Response variable - log(C:N)
# Explanatory (fixed) variables - Bakun Index*Time (aka "sampling")
b1 <- lm(logcn ~ bakun*sampling, data = cn_oi) # Creates the initial linear model.

rb1 <- data.frame(rstandard(b1)) # Standard residuals
# Because there's missing data, here's a workaround to compare residuals to data.
rb1 <- rb1 %>%
  rownames_to_column("record")
rb1_ed <- rb1 %>%
  mutate(rn = as.numeric(record))
cn_oi_ed <- cn_oi %>%
  mutate(RECORD = seq(1,666))
cn_oi_res <- cn_oi_ed %>%
  left_join(rb1_ed, by = c("RECORD" = "rn"))

ggplot(data = cn_oi_res, aes(x = bakun, y = rstandard.b1.)) +
  geom_point() +
  labs(x = "Bakun Index", y = "Standardised residuals") # Plots said residuals.

ggplot(data = cn_oi_res, aes(x = sampling, y = rstandard.b1.)) +
  geom_point() +
  labs(x = "Sampling Date", y = "Standardised residuals")

#### STEP 2: Fit the lm() with GLS and compare to lme().
# We cannot proceed with missing data in the data set so...
cn_oi_rmna <- cn_oi %>%
  drop_na(logcn)

b2 <- gls(logcn ~ bakun*sampling, data = cn_oi_rmna) # Linear regression.
b3 <- lme(logcn ~ bakun*sampling, random =~1 | sitef, data = cn_oi_rmna) # First LMEM with random term.
anova(b2, b3) # Compares the two models. b2 preferred with AIC value of -506.3812, but going to keep in the random term to account for repeated sampling.

# STEP 3: Decide on a variance structure (aka random terms).
plot(b3, col=1) # Check the residuals.
qqnorm(b3) # This actually looks pretty good - so I'm going to skip the added variance structure.

# STEP 4: Fit the lme().

# Using b3 <- lme(logcn ~ bakun*sampling, random =~1 | sitef, data = cn_oi_rmna)

# STEP 5: Compare the lm() and lme().

# See Step 2.

# STEP 6: Everything ok? Check residuals.

# See Step 3.

# STEP 7/8: Step-wise Optimal Fixed Structure

b3_ml <- lme(logcn ~ bakun*sampling, 
             random =~1 | sitef,
             method = "ML", 
             data = cn_oi_rmna) # Switch over to ML for fixed component editing portion.

b4 <- lme(logcn ~ bakun, 
          random =~1 | sitef,
          method = "ML", 
          data = cn_oi_rmna) # Remove "sampling" as a fixed factor.

b5 <- lme(logcn ~ sampling, 
          random =~1 | sitef,
          method = "ML", 
          data = cn_oi_rmna) # Remove "bakun" as a fixed factor.

anova(b3_ml, b4, b5) # Compare the three models. b3_ml preferred with AIC value of -573.0749, so remove nothing.

# STEP 9: Refit with REML

bfinal <- lme(logcn ~ bakun*sampling, 
              random =~1 | sitef,
              method = "REML", 
              data = cn_oi_rmna)

# Output of the model.
summary(bfinal)
# Checking residuals.
plot(bfinal, col=1) # No pattern.
qqnorm(bfinal) # Looks pretty good.
# Final results.
anova(bfinal)

# STEP 10: What does this mean in WORDS?

# My model suggests there is a significant effect of the Bakun index on log(C:N) values of kelp tissue samples as well as a significant effect of sampling date on values; while the interactive effect is not significant, I have left the structure as such, because we know the two are related and AIC values suggest we do. Random intercepts by site were added.

# Equation: log(C:N) = 1.139 - 0.0003[bakun] + 0.0009[sampling] + 3.8e-07[bakun:sampling] + random

#### BEUTI ####

# LMEM #3: C:N vs. BEUTI
# Examine data
plot(logcn ~ beuti, data = cn_oi)
hist(cn_oi$beuti)

#### STEP 1: Create a linear regression and check residuals.
# Based on boxplots & pairplots : 
# Response variable - log(C:N)
# Explanatory (fixed) variables - BEUTI Index*Time (aka "sampling")
c1 <- lm(logcn ~ beuti*sampling, data = cn_oi) # Creates the initial linear model.

rc1 <- data.frame(rstandard(c1)) # Standard residuals
# Because there's missing data, here's a workaround to compare residuals to data.
rc1 <- rc1 %>%
  rownames_to_column("record")
rc1_ed <- rc1 %>%
  mutate(rn = as.numeric(record))
cn_oi_ed <- cn_oi %>%
  mutate(RECORD = seq(1,666))
cn_oi_res <- cn_oi_ed %>%
  left_join(rc1_ed, by = c("RECORD" = "rn"))

ggplot(data = cn_oi_res, aes(x = beuti, y = rstandard.c1.)) +
  geom_point() +
  labs(x = "BEUTI Index", y = "Standardised residuals") # Plots said residuals.
# Might need to add a variance structure later.

ggplot(data = cn_oi_res, aes(x = sampling, y = rstandard.c1.)) +
  geom_point() +
  labs(x = "Sampling Date", y = "Standardised residuals")

#### STEP 2: Fit the lm() with GLS and compare to lme().
c2 <- gls(logcn ~ beuti*sampling, data = cn_oi_rmna) # Linear regression.
c3 <- lme(logcn ~ beuti*sampling, random =~1 | sitef, data = cn_oi_rmna) # First LMEM with random term.
anova(c2, c3) # Compares the two models. c2 preferred with AIC value of -670.7479, but going to keep in the random term to account for repeated sampling.

# STEP 3: Decide on a variance structure (aka random terms).
plot(c3, col=1) # Check the residuals.
qqnorm(c3) # This looks pretty good, but I'm going to try adding in a variance structure by BEUTI index in response to the plots above.

c4 <- lme(logcn ~ beuti*sampling, 
          random =~1 | sitef,  
          data = cn_oi_rmna,
          weights = varIdent(form =~1 | beuti))

# Would not converge, so moving on with c3.

# STEP 4: Fit the lme().

# Using c3 <- lme(logcn ~ beuti*sampling, random =~1 | sitef, data = cn_oi_rmna)

# STEP 5: Compare the lm() and lme().

# See Step 2.

# STEP 6: Everything ok? Check residuals.

# See Step 3.

# STEP 7/8: Step-wise Optimal Fixed Structure

c3_ml <- lme(logcn ~ beuti*sampling, 
             random =~1 | sitef,
             method = "ML", 
             data = cn_oi_rmna)

c4 <- lme(logcn ~ beuti, 
          random =~1 | sitef,
          method = "ML", 
          data = cn_oi_rmna) # Remove "sampling" as a fixed factor.

c5 <- lme(logcn ~ sampling, 
          random =~1 | sitef,
          method = "ML", 
          data = cn_oi_rmna) # Remove "beuti" as a fixed factor.

anova(c3_ml, c4, c5) # Compare the three models. c3_ml preferred with AIC value of -724.1865, so remove nothing.

# STEP 9: Refit with REML

cfinal <- lme(logcn ~ beuti*sampling, 
              random =~1 | sitef,
              method = "REML", 
              data = cn_oi_rmna)

# Output of the model.
summary(cfinal)
# Checking residuals.
plot(cfinal, col=1) # No strong pattern.
qqnorm(cfinal) # Looks pretty good.
# Final results.
anova(cfinal)

# STEP 10: What does this mean in WORDS?

# My model suggests there is a significant effect of the BEUTI index on log(C:N) values of kelp tissue samples as well as a significant effect of sampling date on values; while the interactive effect is not significant, I have left the structure as such, because we know the two are related and AIC values suggest we do. Random intercepts by site were also included.

# Equation: log(C:N) = 1.147 - 0.019[beuti] + 0.001[sampling] - 5.1e-05[beuti:sampling] + random

### CUTI ####

# LMEM #4: C:N vs. CUTI
# Examine data
plot(logcn ~ `34N_CUTI`, data = cn_oi)
hist(cn_oi$`34N_CUTI`)

#### STEP 1: Create a linear regression and check residuals.
# Based on boxplots & pairplots : 
# Response variable - log(C:N)
# Explanatory (fixed) variables - CUTI Index*Time (aka "sampling")
d1 <- lm(logcn ~ `34N_CUTI`*sampling, data = cn_oi) # Creates the initial linear model.

rd1 <- data.frame(rstandard(d1)) # Standard residuals
# Because there's missing data, here's a workaround to compare residuals to data.
rd1 <- rd1 %>%
  rownames_to_column("record")
rd1_ed <- rd1 %>%
  mutate(rn = as.numeric(record))
cn_oi_ed <- cn_oi %>%
  mutate(RECORD = seq(1,666))
cn_oi_res <- cn_oi_ed %>%
  left_join(rd1_ed, by = c("RECORD" = "rn"))

ggplot(data = cn_oi_res, aes(x = `34N_CUTI`, y = rstandard.d1.)) +
  geom_point() +
  labs(x = "CUTI Index", y = "Standardised residuals") # Plots said residuals.

ggplot(data = cn_oi_res, aes(x = sampling, y = rstandard.d1.)) +
  geom_point() +
  labs(x = "Sampling Date", y = "Standardised residuals")

#### STEP 2: Fit the lm() with GLS and compare to lme().
# Need to rename the column for clarity
cn_oi_rmna <- cn_oi_rmna %>%
  rename(cuti = `34N_CUTI`)
d2 <- gls(logcn ~ cuti*sampling, data = cn_oi_rmna) # Linear regression.
d3 <- lme(logcn ~ cuti*sampling, random =~1 | sitef, data = cn_oi_rmna) # First LMEM.
anova(d2, d3) # Compares the two models. d2 preferred with AIC value of -544.0483, but going to keep in the random term to account for repeated sampling.

# STEP 3: Decide on a variance structure (aka random terms).
plot(d3, col=1) # Check the residuals.
qqnorm(d3) # This looks pretty good.

# STEP 4: Fit the lme().

# Using d3 <- lme(logcn ~ cuti*sampling, random =~1 | sitef, data = cn_oi_rmna)

# STEP 5: Compare the lm() and lme().

# See Step 2.

# STEP 6: Everything ok? Check residuals.

# See Step 3.

# STEP 7/8: Step-wise Optimal Fixed Structure

d3_ml <- lme(logcn ~ cuti*sampling, 
             random =~1 | sitef,
             method = "ML", 
             data = cn_oi_rmna)

d4 <- lme(logcn ~ cuti, 
          random =~1 | sitef,
          method = "ML", 
          data = cn_oi_rmna) # Remove "sampling" as a fixed factor.

d5 <- lme(logcn ~ sampling, 
          random =~1 | sitef,
          method = "ML", 
          data = cn_oi_rmna) # Remove "cuti" as a fixed factor.

anova(d3_ml, d4, d5) # Compare the three models. d3_ml preferred with AIC value of -587.4412, so remove nothing.

# STEP 9: Refit with REML

dfinal <- lme(logcn ~ cuti*sampling, 
              random =~1 | sitef,
              method = "REML", 
              data = cn_oi_rmna)

# Output of the model.
summary(dfinal)
# Checking residuals.
plot(dfinal, col=1) # No strong pattern.
qqnorm(dfinal) # Looks pretty good.
# Final results.
anova(dfinal)

# STEP 10: What does this mean in WORDS?

# My model suggests there is a significant effect of the CUTI index on log(C:N) values of kelp tissue samples as well as a significant effect of sampling date on values; while the interactive effect is not significant, I have left the structure as such, because we know the two are related and AIC values suggest we do. Random intercepts by site were also included.

# Equation: log(C:N) = 1.164 - 0.133[cuti] + 0.001[sampling] - 8.6e-05[cuti:sampling] + random

### ENSO ####

# LMEM #5: C:N vs. ENSO

# LMEM #6: C:N vs. MJO

# LMEM #7: C:N vs. NPGO

# LMEM #8: C:N vs. PDO