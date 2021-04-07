# Figures Script
# Heili Lowman
# March 30, 2021

# The following script will create Figures 1 & 2 for the manuscript.
# All figures in the following sections have been exported to my desktop to circumvent any file size issues that could arise when pushing files to GitHub.

#### Setup ####

# Load packages.
library(tidyverse)
library(lubridate)
library(patchwork)
library(calecopal)

# Load datasets from "data_tidying.R".
load("data_tidy/kelp_cn_data_clean.rda")

#### Figure 1 ####

# Need to format dates & calculate log(C:N)
cn_ed <- cn_full %>%
  mutate(Date = mdy(DATE)) %>%
  mutate(logCN = log10(cn))

# Panel A (C:N)
fig1a <- ggplot(cn_ed, aes(x = Date, y = cn, fill = SITE)) + 
  geom_point(shape = 21) +
  scale_fill_manual(values = c("black", "gray60", "white")) +
  scale_x_date(breaks = seq(as.Date("2005-01-01"), as.Date("2020-01-01"), by="5 years"), date_labels = "%Y") +
  annotate('text', x = as.Date("2020-01-01"), y = 45, size = 8, label = "A", family = 'Times New Roman', fontface = "bold", parse = TRUE) +
  labs(x = "Date",
       y = "C:N") +
  theme_bw() +
  theme(text=element_text(family="Times New Roman", size = 20)) +
  theme(legend.title = element_blank(),
        legend.background=element_rect(fill = alpha("white", 0.1)),
        legend.position = c(0.11, 0.85))

fig1a

# Panel B (log(C:N))
fig1b <- ggplot(cn_ed, aes(x = Date, y = logCN, fill = SITE)) + 
  geom_point(shape = 21) +
  scale_fill_manual(values = c("black", "gray60", "white")) +
  scale_x_date(breaks = seq(as.Date("2005-01-01"), as.Date("2020-01-01"), by="5 years"), date_labels = "%Y") +
  annotate('text', x = as.Date("2004-01-01"), y = 1.6, size = 8, label = "B", family = 'Times New Roman', fontface = "bold", parse = TRUE) +
  labs(x = "Date",
       y = "log(C:N)") +
  theme_bw() +
  theme(text=element_text(family="Times New Roman", size = 20)) +
  theme(legend.position = "none")

fig1b

# Panel C (C)
fig1c <- ggplot(cn_ed, aes(x = Date, y = c, fill = SITE)) + 
  geom_point(shape = 21) +
  scale_fill_manual(values = c("black", "gray60", "white")) +
  scale_x_date(breaks = seq(as.Date("2005-01-01"), as.Date("2020-01-01"), by="5 years"), date_labels = "%Y") +
  annotate('text', x = as.Date("2020-01-01"), y = 42.5, size = 8, label = "C", family = 'Times New Roman', fontface = "bold", parse = TRUE) +
  labs(x = "Date",
       y = "% C") +
  theme_bw() +
  theme(text=element_text(family="Times New Roman", size = 20)) +
  theme(legend.position = "none")

fig1c

# Panel D (N)
fig1d <- ggplot(cn_ed, aes(x = Date, y = n, fill = SITE)) + 
  geom_point(shape = 21) +
  scale_fill_manual(values = c("black", "gray60", "white")) +
  scale_x_date(breaks = seq(as.Date("2005-01-01"), as.Date("2020-01-01"), by="5 years"), date_labels = "%Y") +
  annotate('text', x = as.Date("2020-01-01"), y = 4, size = 8, label = "D", family = 'Times New Roman', fontface = "bold", parse = TRUE) +
  labs(x = "Date",
       y = "% N") +
  theme_bw() +
  theme(text=element_text(family="Times New Roman", size = 20)) +
  theme(legend.position = "none")

fig1d

# Compile full figure
figure_1_full <- fig1a + fig1b + fig1c + fig1d

figure_1_full

# ggsave(("Figure_1.png"),
#      path = "/Users/heililowman/Desktop/R_Figures/Kelp_CN",
#      width = 30,
#      height = 15,
#      units = "cm"
#    )

#### Figure 2 ####

# Adding boxplot to demonstrate seasonality in values.

cn_fac <- cn_ed %>%
  mutate(Month = factor(MONTH))

# Make color palette
lake_pal <- cal_palette(name = "lake", n = 12, type = "continuous")

figure_2 <- ggplot(cn_fac, aes(x = Month, y = cn)) + 
  geom_boxplot(aes(fill = Month), alpha = 0.9) +
  scale_fill_manual(values = lake_pal) +
  labs(x = "Month",
       y = "C:N") +
  theme_bw() +
  theme(text=element_text(family="Times New Roman", size = 20)) +
  theme(legend.position = "none")

figure_2

# ggsave(("Figure_2.png"),
#        path = "/Users/heililowman/Desktop/R_Figures/Kelp_CN",
#        width = 30,
#        height = 15,
#        units = "cm"
# )

#### Figure 3 ####

# Summarize SST by month at each site, since monthly resolution is as fine as the C:N data gets.
sst_monthly <- sst_full %>%
  group_by(year, month, site) %>%
  summarize(temp_C_m = mean(temp_C)) %>%
  ungroup()

# Put together SST and CN datasets
cn_sst <- cn_ed %>%
  left_join(sst_monthly, by = c("YEAR" = "year", "MONTH" = "month", "SITE" = "site"))

# Put together indices and CN datasets
cn_sst_oi <- cn_sst %>%
  left_join(oi_full, by = c("YEAR" = "year", "MONTH" = "month"))

# Panel A (SST)
fig3a <- ggplot(cn_sst_oi , aes(x = temp_C_m, y = logCN)) + 
  geom_point() +
  annotate('text', x = 21.5, y = 1.6, size = 8, label = "A", family = 'Times New Roman', fontface = "bold", parse = TRUE) +
  labs(x = "Sea Surface Temperature\n(ºCelsius)",
       y = "log(C:N)") +
  theme_bw() +
  theme(text=element_text(family="Times New Roman", size = 20))

fig3a

# Panel B (Bakun)
fig3b <- ggplot(cn_sst_oi , aes(x = bakun, y = logCN)) + 
  geom_point() +
  annotate('text', x = 390, y = 1.6, size = 8, label = "B", family = 'Times New Roman', fontface = "bold", parse = TRUE) +
  labs(x = "Bakun index",
       y = "log(C:N)") +
  theme_bw() +
  theme(text=element_text(family="Times New Roman", size = 20))

fig3b

# Panel C (BEUTI)
fig3c <- ggplot(cn_sst_oi , aes(x = beuti, y = logCN)) + 
  geom_point() +
  annotate('text', x = 14, y = 1.6, size = 8, label = "C", family = 'Times New Roman', fontface = "bold", parse = TRUE) +
  labs(x = "BEUTI",
       y = "log(C:N)") +
  theme_bw() +
  theme(text=element_text(family="Times New Roman", size = 20))

fig3c

# Panel D (CUTI)
fig3d <- ggplot(cn_sst_oi , aes(x = `34N_CUTI`, y = logCN)) + 
  geom_point() +
  annotate('text', x = 1.25, y = 1.6, size = 8, label = "D", family = 'Times New Roman', fontface = "bold", parse = TRUE) +
  labs(x = "CUTI",
       y = "log(C:N)") +
  theme_bw() +
  theme(text=element_text(family="Times New Roman", size = 20))

fig3d

# Panel E (ENSO)
fig3e <- ggplot(cn_sst_oi , aes(x = enso, y = logCN)) + 
  geom_point() +
  annotate('text', x = 1.75, y = 1.6, size = 8, label = "E", family = 'Times New Roman', fontface = "bold", parse = TRUE) +
  labs(x = "ENSO index",
       y = "log(C:N)") +
  theme_bw() +
  theme(text=element_text(family="Times New Roman", size = 20))

fig3e

# Panel F (MJO)
fig3f <- ggplot(cn_sst_oi , aes(x = mjo, y = logCN)) + 
  geom_point(color = "gray60") +
  annotate('text', x = 1.25, y = 1.6, size = 8, label = "F", family = 'Times New Roman', fontface = "bold", parse = TRUE) +
  labs(x = "MJO index",
       y = "log(C:N)") +
  theme_bw() +
  theme(text=element_text(family="Times New Roman", size = 20))

fig3f

# Panel G (NPGO)
fig3g <- ggplot(cn_sst_oi , aes(x = npgo, y = logCN)) + 
  geom_point() +
  annotate('text', x = 2, y = 1.6, size = 8, label = "G", family = 'Times New Roman', fontface = "bold", parse = TRUE) +
  labs(x = "NPGO index",
       y = "log(C:N)") +
  theme_bw() +
  theme(text=element_text(family="Times New Roman", size = 20))

fig3g

# Panel H (PDO)
fig3h <- ggplot(cn_sst_oi , aes(x = PDO, y = logCN)) + 
  geom_point(color = "gray60") +
  annotate('text', x = 1.6, y = 1.6, size = 8, label = "H", family = 'Times New Roman', fontface = "bold", parse = TRUE) +
  labs(x = "PDO index",
       y = "log(C:N)") +
  theme_bw() +
  theme(text=element_text(family="Times New Roman", size = 20))

fig3h

# Compile full figure
figure_3_full <- fig3a + fig3b + fig3c + fig3d +
  fig3e + fig3f + fig3g + fig3h +
  plot_layout(ncol = 2)

figure_3_full

# ggsave(("Figure_3.png"),
#      path = "/Users/heililowman/Desktop/R_Figures/Kelp_CN",
#      width = 30,
#      height = 30,
#      units = "cm"
#    )

# End of script.