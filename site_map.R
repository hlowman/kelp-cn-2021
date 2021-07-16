# Supplemental Figure Script
# Kyle Emery
# July 16, 2021

# The following script will create Supplemental Figure 1 for the manuscript.
# All figures in the following sections have been exported to my desktop to circumvent any file size issues that could arise when pushing files to GitHub.

#### Setup ####

# Load packages.
library(ggplot2)
library(sf)
library(scales)
library(USAboundaries)
library(tidyverse)
library(wesanderson)
library(ggrepel)
library(ggspatial)
library(ggthemes)
library(ggmap)

# Load datasets.
data <- read.csv("data_raw/coordinates.csv", header = TRUE, sep = ',')

data2 <- st_as_sf(data, coords = c("lon", "lat"), remove = F,
                  crs = 4326, agr = "constant")

data3 <- read.csv("data_raw/coordinates2.csv", header = TRUE, sep = ',')

data4 <- st_as_sf(data3, coords = c("lon", "lat"), remove = F,
                  crs = 4326, agr = "constant")

plot(data2$geometry)

# Create regional basemap.
lats <- c(34.2, 34.6)
lons <- c(-120.5, -119.5)
bb <- make_bbox(lon = lons, lat = lats, f = 0.05)

sb_basemap <- get_stamenmap(bb, 
                      maptype = c("terrain-background"), 
                      source = 'stamen')

ggmap(sb_basemap)

# Create full map.
fullmap <- ggmap(sb_basemap) + # base google maps tile
  geom_label_repel(data=data2,aes(label = Site), size=4,segment.colour = NA,point.size = NA, fontface="bold",fill="black", colour="white") +
  geom_text_repel(data=data4,aes(label = Site), size=6,segment.colour = NA,point.size = NA, fontface="bold", fill=NA, colour="black") +
  annotation_scale( line_width = 2, height = unit(0.75, "cm"),tick_height = 0.6, text_cex = 1, location = "bl", style="ticks") +
  annotation_north_arrow(height = unit(3, "cm"), width = unit(3, "cm"),location = "br", which_north = "true", style = north_arrow_nautical) +
  labs(x = "Longitude",
       y = "Latitude") +
  geom_text(x = -119.85, y = 34.12, label = "Santa Barbara Channel", color = "black", size = 6, fontface = "bold") +
  geom_text(x = -119.78, y = 34.005, label = "Santa Cruz Island", color = "black", size = 6, fontface = "bold") +
  geom_text(x = -120.1, y = 33.96, label = "Santa Rosa Island", color = "black", size = 6, fontface = "bold") +
  theme_bw() +
  coord_sf(crs = st_crs(4326)) +
  theme(plot.subtitle = element_text(vjust = 1), 
          plot.caption = element_text(vjust = 1), 
          panel.background = element_rect(fill = "white")) + 
  theme(axis.line = element_line(linetype = "solid"), 
        panel.grid.major = element_line(colour = NA), 
        panel.grid.minor = element_line(colour = NA),
        panel.border = element_rect(colour = "black", fill=NA, size=2)) + 
  theme(panel.background = element_rect(fill = "aliceblue")) +
  theme(axis.title = element_text(size = 25), 
        axis.text = element_text(size = 20, colour = "black"), 
        legend.text = element_text(size = 20), 
        legend.title = element_text(size = 20)) + 
  theme(legend.position = "bottom", legend.direction = "horizontal")

# Examine full map.
fullmap

# Export map.
# ggsave(("Supp_Figure_1.png"),
#      path = "/Users/heililowman/Desktop/R_Figures/Kelp_CN",
#      width = 35,
#      height = 15,
#      units = "cm"
#    )

# End of script.
