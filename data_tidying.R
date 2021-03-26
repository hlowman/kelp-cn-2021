# Data Assembly Script
# Heili Lowman
# March 24, 2021

# The following script will assemble the disparate datasets into a single main dataset to be used in the accompanying analyses and figures.

# Load packages.
library(tidyverse)
library(lubridate)
library(openxlsx)
library(readxl)
library(naniar)

#### Raw Data ####

# C:N data
# Note - this project also includes the raw Excel file, but there was one date (9/12/19) entered incorrectly, so this has been fixed manually and the file saved as a csv for import below.
cn <- read_csv("data_raw/CHN_Field_all_years.csv")

# SST datasets (3 need to be knitted together)
sst1 <- read_csv("data_raw/LTER_sites_20000101_20190501.csv")
sst2 <- read_csv("data_raw/LTER_sites_20190502_20191001.csv")
sst3 <- read_csv("data_raw/LTER_sites_20191001_20210301.csv")

# Oceanographic indices (7 unique)

bakun <- read_csv("data_raw/Bakun_monthly_upwelling.csv")
beuti <- read_csv("data_raw/BEUTI_Upwelling.csv")
cuti <- read_csv("data_raw/CUTI_Upwelling.csv")
enso <- read_csv("data_raw/MEIV2_ENSO.csv")
mjo <- read_csv("data_raw/MJO_Index6_120W_Daily.csv")
npgo <- read_csv("data_raw/npgo.csv")
pdo <- read_csv("data_raw/PDO.csv")

#### Tidy Data ####

period <- c(".")

cn_full <- cn[-c(1:2),] %>% # remove first two rows of data
  replace_with_na_all(condition = ~.x %in% period) %>% # remove all periods
  mutate(WET_WT = as.numeric(WET_WT),
         DRY_WET = as.numeric(DRY_WET)) %>% # need to change char to num
  group_by(YEAR, MONTH, DATE, SITE, `_N`) %>% # group and...
  summarize(wet_wt = mean(WET_WT, na.rm=TRUE),
            dry_wt = mean(DRY_WT, na.rm=TRUE),
            dry_wet = mean(DRY_WET, na.rm=TRUE),
            anal_wt = mean(ANAL_WT, na.rm=TRUE),
            c = mean(C, na.rm=TRUE),
            h = mean(H, na.rm=TRUE),
            n = mean(N, na.rm=TRUE),
            cn = mean(CN_RATIO, na.rm=TRUE)) %>% # summarize across both replicates so that future analyses don't include pseudoreplication
  ungroup() %>%
  replace_with_na_all(condition = ~.x == -99999) # remove all -99999s
  
# Knit together temperature file.
sst2_2 <- sst2 %>%
  mutate(DATE = lubridate::ymd(as.character(date))) %>% # convert to date
  mutate(day = day(DATE)) %>% # pull out day
  select(year, month, day, temp_C, site) # match sst1 formatting

sst3_2 <- sst3 %>%
  mutate(DATE = lubridate::ymd(as.character(date))) %>% # convert to date
  mutate(month = month(DATE), day = day(DATE)) %>% # pull out day & month
  select(year, month, day, temp_C, site) # match sst1 formatting

sst_full <- sst1 %>%
  rbind(sst2_2) %>%
  rbind(sst3_2)

# Create compiled oceanographic indices file

beuti_ed <- beuti %>%
  rename(beuti = `34N`)

# some raw files missing the month column
m0nth <- beuti_ed$month

bakun_ed <- bakun %>%
  mutate(month = m0nth) %>% # add month to bakun dataset
  select(year, month, bakun) # and reorder

enso_ed <- enso %>%
  mutate(month = m0nth) %>% # add month to enso dataset
  select(year, month, enso) # and reorder

# need to average MJO
mjo_ed <- mjo %>%
  group_by(year, month) %>%
  summarize(mjo = mean(MJO_120W)) %>%
  ungroup()

# need to tidy NPGO file a bit
# as of 3/24/2021, the most recent data available for npgo is 7/2020
# see http://www.o3d.org/npgo/
year_add <- c(2020, 2020, 2020, 2020, 2020, 2021, 2021)
month_add <- c(8, 9, 10, 11, 12, 1, 2)  
npgo_add <- c(NA, NA, NA, NA, NA, NA, NA)  
addition <- data.frame(year_add, month_add, npgo_add) %>%
  rename(year = year_add, month = month_add, npgo = npgo_add)

npgo_ed <- npgo[-c(1:4),] %>% # removes jan-april 2002
  rbind(addition)

# Joining all the indices into a single dataset (wide format)
oi_full <- full_join(bakun_ed, beuti_ed) %>%
  full_join(cuti) %>%
  full_join(enso_ed) %>%
  full_join(mjo_ed) %>%
  full_join(npgo_ed) %>%
  full_join(pdo)
  
#### Export Data ####

save(cn_full, sst_full, oi_full, file = "data_tidy/kelp_cn_data_clean.rda")

# End of script.
