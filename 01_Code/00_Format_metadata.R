# Info -------------------------------------------------------------------------

#
# Formating metadata as a function of the different subprojects
#
# Christelle Leung
# 2022-10
#

# Library ----------------------------------------------------------------------
rm(list = ls())
gc()

library(tidyverse)
library(readxl)

# Files ------------------------------------------------------------------------

SeqInfo_Biodiv <- read_xlsx("00_Data/00_FileInfos/SeqInfo_Biodiversite2021.xlsx")
SeqInfo_Marin <- read_xlsx("00_Data/00_FileInfos/SeqInfo_Marin2021.xlsx")

ppo <- read_xlsx("00_Data/00_FileInfos/20211201_ADNe_PPO_compilation 2018-2021_Version SV2.xlsx",
                    sheet = 2)

bda.samples <- read_xlsx("00_Data/00_FileInfos/20220830_BDA_ADNe_V02.xlsx",
                         sheet = "Echantillons_ADNe")
bda.extract <- read_xlsx("00_Data/00_FileInfos/20220830_BDA_ADNe_V02.xlsx",
                         sheet = "Extraits_ADNe")

bda <- bda.samples %>% 
  left_join(bda.extract, by = c("Nom_reception_echantillon" = "Nom_echantillon_reception"))

## Subset samples --------------------------------------------------------------

### Godbout 2021 

df <- subset(SeqInfo_Biodiv, 
                         ID_subprojet == "PPO.Godbout.2021" & 
               Type_echantillon == "ECH") %>%
  left_join(., ppo, by = c("ID_labo" = "ID_ADNe"))
df.godbout2021 <- data.frame(ID_ADNe = df$ID_labo,
                             Site = df$Site_echantillonnage.y,
                             ID_reception = df$ID_echantillon_reception, 
                             Latitude = df$Latitude,
                             Longitude = df$Longitude,
                             Depth_CTD = df$`Profondeur CTD (m)`,
                             Tide = df$`Marée (m)`,
                             Depth_m = as.numeric(df$`Profondeur CTD (m)`) - as.numeric(df$`Marée (m)`),
                             Temperature = df$`Temperature (°C)`,
                             Salinity = df$`Salinite ‰`,
                             Date = as.Date(as.numeric(df$Date_collecte), origin = "1899-12-30"),
                             Project = df$ID_subprojet,
                             Locus = df$Loci)
rm(df)

### Others PPO

df <- subset(SeqInfo_Marin, 
             ID_subproject %in% c("PPO.Leim.VRoy.2019",
                                  "PPO.PO.DLevesque.2019",
                                  "PPO.Leim.MJRoux.2019",
                                  "PPO.Kildir.DLevesque.2020",
                                  "PPO.Leim.MJRoux.2020",
                                  "PPO.KMcGregor.2020") &
             Type_echantillon == "ECH") %>%
  left_join(., ppo, by = c("ID_labo" = "ID_ADNe"))

df.ppo.other <- data.frame(ID_ADNe = df$ID_labo,
                          Site = df$Site_echantillonnage.y,
                          ID_reception = df$ID_echantillon_reception, 
                          Latitude = df$Latitude.y,
                          Longitude = df$Longitude.y,
                          Depth_CTD = df$`Profondeur CTD (m)`,
                          Tide = df$`Marée (m)`,
                          Depth_m = as.numeric(df$`Profondeur CTD (m)`) - as.numeric(df$`Marée (m)`),
                          Temperature = df$`Temperature (°C)`,
                          Salinity = df$`Salinite ‰`,
                          Date = as.Date(as.numeric(df$Date_collecte), origin = "1899-12-30"),
                          Project = df$ID_subproject,
                          Locus = df$Loci)
rm(df)

### PPO all
df.ppo <- rbind(df.godbout2021, df.ppo.other)

### BDA
df <- subset(SeqInfo_Marin, 
             ID_subproject == "BDA" &
               Type_echantillon == "ECH") %>%
  left_join(., bda, by = c("ID_labo" = "Numero_unique_extrait"))

df.bda <- data.frame(ID_ADNe = df$ID_labo,
                     Site = substr(df$Station, 1, 3),
                     ID_reception = df$Nom_reception_echantillon, 
                     Latitude = df$Latitude,
                     Longitude = df$Longitude,
                     Depth_CTD = NA,
                     Tide = NA,
                     Depth_m = as.numeric(df$Profondeur_m_echantillon),
                     Temperature = NA,
                     Salinity = NA,
                     Date = as.Date(as.numeric(df$Date_collected), origin = "1899-12-30"),
                     Project = df$ID_subproject,
                     Locus = df$Loci)
rm(df)

## metadata for ALL ----

metaData.all <- rbind(df.ppo, df.bda)
metaData.all <- data.frame(metaData.all)
write.csv(metaData.all, "00_Data/00_FileInfos/metaData_2019_to_2021.csv")


## Map of all sample sites ----
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)
library(sf)

world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

lati <- as.numeric(metaData.all$Latitude)
long <- as.numeric(metaData.all$Longitude) %>% abs()*(-1)
project <- metaData.all$Subproject

colvec <- rainbow(7)

map.fig <- ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-67, -69.5), ylim = c(48.5, 49.5),
           expand = F) +
  theme_bw() +
  annotation_scale(location = "br", width_hint = 0.25, bar_cols = c("grey", "white"),
                   line_col = "grey30") +
  annotation_north_arrow(location = "tr", which_north = "true", 
                         pad_x = unit(0, "in"), pad_y = unit(0.1, "in"),
                         style = north_arrow_fancy_orienteering(
                           line_col = "grey30", 
                           fill = c("white", "grey"),
                         )) +
  annotate(geom = "point", x = long, y = lati, shape = 19, 
           color = colvec[as.factor(project)], size = 3, alpha = 0.3) +
  # annotate(geom = "label", x = long, y = lati, label = sites,
  #          color = "black", size = 1, 
  #          vjust = 1.3, hjust = 0.5,
  #          fill = alpha(c("white"),0.75), label.size = NA) +
  xlab("") + ylab("") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

map.fig

### ggsave("02_Results/Map_Godbout2021.pdf", map.fig, width = 4, height = 2, scale = 1.5)

