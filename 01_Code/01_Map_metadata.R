# Info -------------------------------------------------------------------------

#
# Map of sample locations and metaData tables
# Locus:
#     - MiFishU
#     - COI
#     - 16S
# 
# Samples:
#     - Godbout 2021 and 2020
#     - VRoy 2019
#     - BDA 2020
#
# CL
# 2023-04
#

# Library ----------------------------------------------------------------------
rm(list = ls())
gc()

library(tidyverse)

library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)

library(marmap)
library(sf)

# Creat new_dir for results  ---------------------------------------------------

new_dir <- "06_Map_metadata"
dir.create(file.path(here::here(), "02_Results", new_dir),
           showWarnings = FALSE)

# Data -------------------------------------------------------------------------

metaData <- read.csv(file.path(here::here(), "00_Data", "00_FileInfos",
                               "metaData_2019_to_2021.csv"),
                     header = T, row.names = 1)
metaData <- metaData %>%
  mutate(group = ifelse(Project %in% c("BDA", "PPO.Leim.VRoy.2019"), 
                        metaData$Site,
                        ifelse(Date %in% c("2020-10-27", "2020-10-28"),
                               "Fall2020", 
                               ifelse(Date %in% c("2021-08-24", "2021-08-25"),
                                      "Summer2021",
                                      "Spring2021"))))

all <- subset(metaData, Project %in% c("PPO.Godbout.2021", "PPO.KMcGregor.2020",
                                        "BDA",
                                        "PPO.Leim.VRoy.2019"))

fine.scale <- subset(metaData, Project %in% c("PPO.Godbout.2021", "PPO.KMcGregor.2020"))

med.scale <- subset(metaData, Project %in% c("BDA")) %>%
  mutate(group = ifelse(group == "CrB" & Latitude >= 48.6, "CrA", group))

large.scale <- subset(metaData, Project %in% c("PPO.Leim.VRoy.2019"))

colvec <- RColorBrewer::brewer.pal(6, "Set1")
sites <- colvec[as.factor(large.scale$Site)]

large.xlim <- range(large.scale$Longitude) + c(-0.4, 0.4)
large.ylim <- range(large.scale$Latitude) + c(-0.2, 0.2)

med.xlim <- range(med.scale$Longitude) + c(-0.2, 0)
med.ylim <- range(med.scale$Latitude) + c(-0.2, 0.2)

fine.xlim <- range(fine.scale$Longitude) + c(-0.01, 0.01)
fine.ylim <- range(fine.scale$Latitude) + c(-0.01, 0.01) 

isobath <- read_sf(file.path(here::here(), "00_Data", "isobathe_10_50.shp")) %>%
  st_transform(crs = 4326)

box = c(xmin = -67.705, ymin = 49.290, xmax = -67.692, ymax = 49.35)
isobath2 <- st_crop(isobath, box)


## Maps ------------------------------------------------------------------------

pannel.a <- getNOAA.bathy(min(-71), max(-60), min(46), max(51), 
                          res = 1, keep = TRUE)

plot.A <- autoplot.bathy(pannel.a, geom="raster") +
  scale_fill_gradient(limits = c(-500, 0), low="steelblue4", high="white") +
  labs(y = "Latitude", x = "Longitude", fill = "Depth (m)") +
  coord_cartesian(expand = 0) +
  theme_bw() +
  geom_rect(#data = data.frame(),
            aes(xmin = large.xlim[1], xmax = large.xlim[2], 
                ymin = large.ylim[1], ymax = large.ylim[2]),
            colour = "red", fill = "transparent") +
  geom_rect(#data = data.frame(),
            aes(xmin = med.xlim[1], xmax = med.xlim[2], 
                ymin = med.ylim[1], ymax = med.ylim[2]),
            colour = "red", fill = "transparent") +
  geom_rect(##=data = data.frame(),
            aes(xmin = fine.xlim[1]-0.05, xmax = fine.xlim[2]+0.05, 
                ymin = fine.ylim[1]-0.05, ymax = fine.ylim[2]+0.05),
            colour = "red", fill = "transparent") +
  ggtitle("Study area") 
plot.A


### broad scale ----------------------------------------------------------------

broad.scale <- getNOAA.bathy(min(large.xlim), max(large.xlim), min(large.ylim), max(large.ylim), 
                             res = 0.1, keep = TRUE)

long <- large.scale$Longitude
lati <- large.scale$Latitude

plot.B <- autoplot.bathy(broad.scale, geom=c("raster")) +
  scale_fill_gradient(limits = c(-500, 0), low="steelblue4", high="white") +
  labs(y = "Latitude", x = "Longitude", fill = "Depth (m)") +
  coord_cartesian(expand = 0) +
  theme_bw() +
  annotate(geom = "point", x = long, y = lati, shape = 4, 
           color = alpha("#a30000", 0.5), 
           alpha = 0.5, size = 1) +
  annotate(geom = "label", 
           x = large.scale[!duplicated(large.scale$Site),]$Longitude,
           y = large.scale[!duplicated(large.scale$Site),]$Latitude, 
           label = large.scale[!duplicated(large.scale$Site),]$Site,
           color = "#a30000", size = 3, 
           vjust = 1.25, hjust = 0,
           fill = NA, 
           label.size = NA
           ) +
  ggtitle("Broad scale") 
plot.B

### intermediate scale ---------------------------------------------------------

inter.scale <- getNOAA.bathy(min(med.xlim), max(med.xlim), min(med.ylim), max(med.ylim), 
                             res = 0.1, keep = TRUE)

med.scale$Site <- gsub("Cr", "", med.scale$Site)
long <- med.scale$Longitude
lati <- med.scale$Latitude
colvec <- RColorBrewer::brewer.pal(3, "Set1")
sites <- colvec[as.factor(med.scale$Site)]

plot.c <- autoplot.bathy(inter.scale, geom=c("raster")) +
  scale_fill_gradient(limits = c(-500, 0), low="steelblue4", high="white") +
  labs(y = "Latitude", x = "Longitude", fill = "Depth (m)") +
  coord_cartesian(expand = 0) +
  theme_bw() +
  annotate(geom = "point", x = med.scale$Longitude, y = med.scale$Latitude, shape = 4, 
                     color = colvec[as.factor(med.scale$Site)],
                     alpha = 0.5, size = 2) +
  annotate(geom = "label", 
           x = med.scale[!duplicated(med.scale$Site),]$Longitude,
           y = med.scale[!duplicated(med.scale$Site),]$Latitude, 
           label = med.scale[!duplicated(med.scale$Site),]$Site,
           color = colvec, size = 4, fontface = "bold",
           vjust = c(0, 1.4, 1.6), hjust = c(-0.6, -0.2, 0.1),
           fill = alpha(c("white"),0.1), label.size = NA) +
  ggtitle("Intermediate scale") 
plot.c



### fine scale ------------------------------------------------------------------

long <- fine.scale$Longitude
lati <- fine.scale$Latitude

fine.scale$group <- factor(fine.scale$group, levels = c("Spring2021", "Summer2021", "Fall2020"))
levels(fine.scale$group) <- c("Spring 2021", "Summer 2021", "Fall 2020")

df.season <- data.frame(group = factor(c("Spring 2021", "Summer 2021", "Fall 2020")),
                        x = -67.704, y = 49.305)

plot.fine <- ggplot() +
  geom_sf(data = isobath2, aes(colour = Contour), size = .1) +
  geom_sf_label(data = isobath2, 
               aes(label = Contour, colour = Contour), 
               size = 3, label.size = NA, 
               fill = NA, hjust = -1, vjust = 0.2) +
  coord_sf(expand = FALSE, ylim = c(49.29, 49.306)) +
  scale_y_continuous(
    labels = scales::number_format(accuracy = 0.01),
    breaks = pretty(fine.scale$Latitude, n = 2)[c(1,3)]) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01),
                     breaks = pretty(fine.scale$Longitude, n = 1)[2])
plot.d <- plot.fine +
  geom_point(data = fine.scale, aes(x = Longitude, y = Latitude),
             col = "#a30000", shape = 4) +
  labs(y = "Latitude", x = "Longitude") +
  facet_grid(~group)+
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "none",
        panel.background = element_rect(fill = alpha("lightsteelblue1", 0.5))) +
  ggtitle("Fine scale") +
  geom_text(aes(label = group, x = x, y = y), 
            data = df.season, vjust = 1, hjust = 0, size = 3)
plot.d



fig.all <- ggpubr::ggarrange(plot.A, plot.B, plot.c, plot.d,
                             common.legend = TRUE, legend = "right",
                             align = "hv", labels = LETTERS) +
  theme(plot.background = element_rect(fill = "white"),
        plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "cm"))
fig.all

ggsave(file.path(here::here(), "02_Results", new_dir, "Samples_Map.png"),
       fig.all, width = 6.5, height = 4, scale = 1.5)



