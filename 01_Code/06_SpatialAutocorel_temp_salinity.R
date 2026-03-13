# Info -------------------------------------------------------------------------
#
# - Spatial autocorrelation for Temperature and Salinity 
# 
# CL
# 2024-05
#
# Library ----------------------------------------------------------------------
rm(list = ls())
gc()

library(tidyverse)
library(seqinr)
library(Biostrings)
library(msa)
library(ape)
library(ggtree)

library(ggpubr)
library(vegan)


# Spatial autocorrelation for Temperature and Salinity  ------------------------
rm(list = ls())
gc()

new_dir <- "99_Result_Temp_Salinity"


### Function - 3D dist
distance3D <- function (df) {
  planiDist <- geosphere::distm(df[,1:2])
  depthDist <- dist(df[3], method = "euclidean") %>% as.matrix()
  dist3D <- sqrt(planiDist^2 + depthDist^2) %>% as.dist()
  return(dist3D)
}

### metadata GODBOUT

list.locus <- c("COI", "MiFishU", "16Schord")

metaData.all <- mat.or.vec(0,0)

for(i in 1:length(list.locus)){
  
  project.name = "Godbout"
  locus.name = list.locus[i]
  
  metaData <- read.csv(file.path(here::here(), "00_Data", "00_ESVs_tables_clean",
                                 paste(project.name, locus.name, "metadata_ok.csv", sep = "_")),
                       header = T, row.names = 1)
  metaData <- metaData %>%
    mutate(group = ifelse(Project %in% c("BDA", "PPO.Leim.VRoy.2019"),
                          metaData$Site,
                          ifelse(Date %in% c("2020-10-27", "2020-10-28"),
                                 "Fall2020",
                                 ifelse(Date %in% c("2021-08-24", "2021-08-25"),
                                        "Summer2021",
                                        "Spring2021")))) %>%
    select(ID_ADNe, Latitude, Longitude, Depth_m, Temperature, Salinity, 
           Locus,
           Date, group)
  
  metaData.all <- rbind(metaData.all, metaData) %>% distinct
}


### Tous les sites ensemble ----------------------------------------------------

rm("metaData")

library(mpmcorrelogram)

res.mantel.all <- mat.or.vec(0,0)
list.group <- unique(metaData.all$group)


for(i in 1:length(list.group)){
  
  metaData <- metaData.all %>%
    select(-Locus) %>% 
    filter(group == list.group[i]) %>%
  distinct
  
  env.df <- metaData[, c("Longitude", "Latitude", "Depth_m")]
  rownames(env.df) <- metaData$ID_ADNe
  env.df <- env.df[rowSums(is.na(env.df)) == 0,]
  dist.depth <- vegdist(env.df$Depth_m, method = "euclidean")
  dist.2d <- geosphere::distm(env.df[, c("Longitude", "Latitude")])
  dist.3d <- distance3D(env.df)
  
  
  dist.data <- vegdist(metaData %>% select(Salinity, Temperature), method = "euclidean")
 
  correl.depth <- mpmcorrelogram(xdis = dist.data, 
                                 geodis = dist.depth,
                                 zdis = dist.2d, 
                                 print = F)
  
  correl.2d <- mpmcorrelogram(xdis = dist.data, 
                              geodis = dist.2d,
                              zdis = dist.depth,
                              print = F)
  
  correl.3d <- mpmcorrelogram(xdis = dist.data, 
                              geodis = dist.3d,
                              print = F)
  
  df.correl <- data.frame(class.index = c(correl.depth$breaks[-1],
                                          correl.2d$breaks[-1],
                                          correl.3d$breaks[-1]),
                          Mantel.cor = c(correl.depth$rM,
                                         correl.2d$rM,
                                         correl.3d$rM),
                          Pr.Mantel = c(correl.depth$pvalues,
                                        correl.2d$pvalues,
                                        correl.3d$pvalues),
                          Pr.corrected = c(correl.depth$pval.Bonferroni,
                                           correl.2d$pval.Bonferroni,
                                           correl.3d$pval.Bonferroni))
  n.class = nrow(df.correl)/3
  
  df.correl$distance <- rep(c("Depth|2D", "2D|Depth", "3D"), each = n.class)
  df.correl$distance <- factor(df.correl$distance, levels = c("Depth|2D", "2D|Depth", "3D"))
  df.correl$Signif <- ifelse(df.correl$Pr.corrected >= 0.05, "no", "yes")
  df.correl$group <- list.group[i]
  
  res.mantel.all <- rbind(res.mantel.all, df.correl)

}

write.csv(res.mantel.all,
          file.path(here::here(), "02_Results", new_dir, 
                    "Mantel_df_all.csv"))  

res.mantel.all$group <- factor(res.mantel.all$group,
                               levels = c("Spring2021", "Summer2021", "Fall2020"))
levels(res.mantel.all$group) <- c("Spring", "Summer", "Fall")

correl.plot <- res.mantel.all %>%
  filter(distance == "3D") %>%
  ggplot(., aes(x = class.index, y = Mantel.cor)) +
  geom_hline(yintercept = 0, color = "red") +
  geom_line() + 
  scale_x_continuous(labels = scales::label_comma()) +
  geom_point(aes(shape = Signif, size = 1)) +
  scale_shape_manual(values = c(1,19)) + 
  facet_grid(~ group, scales = "free_x") +
  xlab("Distance (m)") +
  ylab("Mantel correlation") +
  theme_bw() +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0)) +
  ggtitle("3D Mantel correlogram",
          subtitle = "(Fine scale samples, based on both temperature and salinity)")

ggsave(file.path(here::here(), "02_Results", new_dir, 
                 "Mantel_correl_all_temp_salinity.png"),
       correl.plot, 
       height = 2, width = 5, scale = 1.5)

## Par locus, car pas les mêmes sites -------------------------------------------
  
rm("metaData")

res.mantel.loci <- mat.or.vec(0,0)

list.group <- unique(metaData.all$group)
list.locus <- c("COI", "MiFishU", "16Schord")

for(j in 1:length(list.locus)){
  
  locus.name = list.locus[j]
  
  for(i in 1:length(list.group)){
    
    group.name <- list.group[i]
    
    metaData <- metaData.all %>%
      filter(group == group.name,
             Locus == locus.name) %>%
      distinct
    
    env.df <- metaData[, c("Longitude", "Latitude", "Depth_m")]
    rownames(env.df) <- metaData$ID_ADNe
    env.df <- env.df[rowSums(is.na(env.df)) == 0,]
    dist.depth <- vegdist(env.df$Depth_m, method = "euclidean")
    dist.2d <- distm(env.df[, c("Longitude", "Latitude")])
    dist.3d <- distance3D(env.df)

    
    dist.data <- vegdist(metaData %>% select(Salinity, Temperature), method = "euclidean")
    
    correl.depth <- mpmcorrelogram(xdis = dist.data, 
                                   geodis = dist.depth,
                                   zdis = dist.2d, 
                                   print = F)
    
    correl.2d <- mpmcorrelogram(xdis = dist.data, 
                                geodis = dist.2d,
                                zdis = dist.depth,
                                print = F)
    
    correl.3d <- mpmcorrelogram(xdis = dist.data, 
                                geodis = dist.3d,
                                print = F)
    
    df.correl <- data.frame(class.index = c(correl.depth$breaks[-1],
                                            correl.2d$breaks[-1],
                                            correl.3d$breaks[-1]),
                            Mantel.cor = c(correl.depth$rM,
                                           correl.2d$rM,
                                           correl.3d$rM),
                            Pr.Mantel = c(correl.depth$pvalues,
                                          correl.2d$pvalues,
                                          correl.3d$pvalues),
                            Pr.corrected = c(correl.depth$pval.Bonferroni,
                                             correl.2d$pval.Bonferroni,
                                             correl.3d$pval.Bonferroni))
    n.class = nrow(df.correl)/3
    
    df.correl$distance <- rep(c("Depth|2D", "2D|Depth", "3D"), each = n.class)
    df.correl$distance <- factor(df.correl$distance, levels = c("Depth|2D", "2D|Depth", "3D"))
    df.correl$Signif <- ifelse(df.correl$Pr.corrected >= 0.05, "no", "yes")
    df.correl$group <- group.name
    df.correl$Locus <- locus.name
    
    res.mantel.loci <- rbind(res.mantel.loci, df.correl)
    
  }
}


write.csv(res.mantel.loci,
          file.path(here::here(), "02_Results", new_dir, 
                    "Mantel_df_by_loci.csv"))  

res.mantel.loci$group <- factor(res.mantel.loci$group,
                                levels = c("Spring2021", "Summer2021", "Fall2020"))

list.plot <- list()

for(i in 1:length(list.locus)){
  
  correl.plot <- res.mantel.loci %>%
    filter(Locus == list.locus[i],
           distance == "3D") %>% 
    ggplot(., aes(x = class.index, y = Mantel.cor)) +
    geom_hline(yintercept = 0, color = "red") +
    geom_line() + 
    scale_x_continuous(labels = scales::label_comma()) +
    geom_point(aes(shape = Signif, size = 1)) +
    scale_shape_manual(values = c(1,19)) + 
    facet_grid(~ group, scales = "free_x") +
    xlab("Distance (m)") +
    ylab("Mantel correlation") +
    theme_bw() +
    theme(legend.position = "none",
          strip.background = element_blank(),
          strip.text = element_text(hjust = 0)) +
    ggtitle(paste("Godbout: stations with only", list.locus[i]),
            subtitle = "(based on both temmperature and salinity)")
  
  list.plot[[i]] <- correl.plot
  
}

correl.plot <- ggarrange(list.plot[[1]],
                         list.plot[[2]],
                         list.plot[[3]],
                         nrow = 3)

ggsave(file.path(here::here(), "02_Results", new_dir, 
                 "Mantel_correl_by_locus.png"),
       correl.plot, 
       height = 6, width = 5, scale = 1.5)
