# Info -------------------------------------------------------------------------

#
# Using G. Guenard package gslADNe to analyse eDNA metabarcoding data
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
# Using OTUs instead of ESVs (to avoid Mahalanobis distance typos)
# Test for biological effect compared to other factors
#
# CL
# 2023-09
#

# Library ----------------------------------------------------------------------
rm(list = ls())
gc()

library(gslADNe)
library(vegan)
library(mpmcorrelogram)

library(tidyverse)
library(seqinr)

library(ggpubr)
library(geosphere)

# library(mvpart)
# library(MVPARTwrap)

library(gridExtra)

# Creat new_dir for results  ---------------------------------------------------

new_dir <- "09_Biological_Replicate"
dir.create(file.path(here::here(), "02_Results", new_dir),
           showWarnings = FALSE)


# 1) Data transformation and RDA -----------------------------------------------

rm(list=setdiff(ls(), c("new_dir")))
gc()



## 1.1) Spatial analyses using Hellinger-Mahalanobis ---------------------------

list.locus <- c("COI", "MiFishU", "16Schord")


## Fine scale ------------------------------------------------------------------

### Godbout --------------------------------------------------------------------

project.name <- c("Godbout")

res.rda <- mat.or.vec(0,0)

for(i in 1:length(list.locus)){
  
  locus.name <- list.locus[i]
  
  OTUs.count <- read.csv(file.path(here::here(), "00_Data", "07_OTUs_3pc_tables",
                                   paste(project.name, locus.name, 
                                         "OTUs_count.csv", sep = "_")),
                         header = T, row.names = 1)
  
  OTUs.trans <- read.csv(file.path(here::here(), "02_Results", 
                                   "05c_Spatial_diversity_with_OTUs_3pc_NewFigures",
                                   "Transformed_data",
                                   paste("Hell_Maha", 
                                         "_", project.name, "_", locus.name, 
                                         ".csv", sep = "")),
                         row.names = 1)
  
  metaData <- read.csv(file.path(here::here(), "00_Data", "00_ESVs_tables_clean",
                                 paste(project.name, locus.name, 
                                       "metadata_ok.csv", sep = "_")),
                       header = T, row.names = 1)
  
  metaData$Transect <- sub("M.*", "M", metaData$ID_reception)
  metaData$Year <- sub("-.*", "", metaData$Date) %>% as.numeric
  metaData$Season <- ifelse(metaData$Date %in% c("2020-10-27", "2020-10-28"),
                            "Fall", 
                            ifelse(metaData$Date %in% c("2021-08-24", "2021-08-25"),
                                   "Summer",
                                   "Spring"))
  metaData$Station <- sub("-.*", "", metaData$ID_reception)
  metaData$biorep <- sub(".*-", "", metaData$ID_reception)
  
  nreads <- data.frame(ID_ADNe = rownames(OTUs.count),
                       nreads = rowSums(OTUs.count))
  
  metaData <- left_join(metaData, nreads)
  
  cv <- readRDS(file.path(here::here(), "02_Results", 
                          "05c_Spatial_diversity_with_OTUs_3pc_NewFigures",
                          "Transformed_data",
                          paste("Seq_CoVar", 
                                "_", project.name, "_", locus.name, 
                                ".rds", sep = "")))
  
  ## analysis by season 
  
  list.season <- c("Fall", "Summer", "Spring")
  
  for(j in 1:length(list.season)){
    
    meta.data <- subset(metaData, Season == list.season[j])
    count.trans <- OTUs.trans[meta.data$ID_ADNe, ]
    
    # sequencing depth effect 
    
    rda.seq.depth <- rda(count.trans ~ nreads, meta.data)
    res.seq.depth <- anova.cca(rda.seq.depth) %>% data.frame
    
    seq.depth.effect.p <- p.adjust(res.seq.depth$Pr..F.[1], 
                                   method = "bonferroni",
                                   n = 15)
    
    if(seq.depth.effect.p <= 0.05){
      
      rda.trans <- rda(count.trans ~ biorep + Station + Condition(nreads), 
                       meta.data, 
                       na.action = na.omit, scale = TRUE)
      
      res.temp <- anova.cca(rda.trans, by = "margin") %>% data.frame
      res.temp$R2 <- c(RsquareAdj(rda(count.trans ~ biorep + Condition(Station + nreads), meta.data))$r.squared,
                       RsquareAdj(rda(count.trans ~ Station + Condition(biorep + nreads), meta.data))$r.squared,
                       NA)
      res.temp$adj.R2 <- c(RsquareAdj(rda(count.trans ~ biorep + Condition(Station + nreads), meta.data))$adj.r.squared,
                           RsquareAdj(rda(count.trans ~ Station + Condition(biorep + nreads), meta.data))$adj.r.squared,
                           NA)
      res.temp$model <- "~ biorep + Station + Condition(nreads)"
      res.temp$factor <- rownames(res.temp)
      
      rda.trans1 <- rda(count.trans ~ biorep + Station + Condition(Transect + nreads), 
                        meta.data,
                        na.action = na.omit, scale = TRUE)
      
      res.temp1 <- anova.cca(rda.trans1, by = "margin") %>% data.frame
      res.temp1$R2 <- c(RsquareAdj(rda(count.trans ~ biorep + Condition(Transect + Station + nreads), meta.data))$r.squared,
                        RsquareAdj(rda(count.trans ~ Station + Condition(Transect + biorep + nreads), meta.data))$r.squared,
                        NA)
      res.temp1$adj.R2 <- c(RsquareAdj(rda(count.trans ~ biorep + Condition(Transect + Station + nreads), meta.data))$adj.r.squared,
                            RsquareAdj(rda(count.trans ~ Station + Condition(Transect + biorep + nreads), meta.data))$adj.r.squared,
                            NA)
      res.temp1$model <- "~ biorep + Station + Condition(Transect + nreads)"
      res.temp1$factor <- rownames(res.temp1)
      
      rda.trans2 <- rda(count.trans ~ biorep + Transect + Condition(nreads),
                        meta.data,
                        na.action = na.omit, scale = TRUE)
      
      res.temp2 <- anova.cca(rda.trans2, by = "margin") %>% data.frame
      res.temp2$R2 <- c(RsquareAdj(rda(count.trans ~ biorep + Condition(Transect + nreads), meta.data))$r.squared,
                        RsquareAdj(rda(count.trans ~ Transect + Condition(biorep + nreads), meta.data))$r.squared,
                        NA)
      res.temp2$adj.R2 <- c(RsquareAdj(rda(count.trans ~ biorep + Condition(Transect + nreads), meta.data))$adj.r.squared,
                            RsquareAdj(rda(count.trans ~ Transect + Condition(biorep + nreads), meta.data))$adj.r.squared,
                            NA)
      res.temp2$model <- "~ biorep + Transect + Condition(nreads)"
      res.temp2$factor <- rownames(res.temp2)
      
      res.temp$season <- list.season[j]
      res.temp1$season <- list.season[j]
      res.temp2$season <- list.season[j]
      
    }else{
      
      rda.trans <- rda(count.trans ~ biorep + Station, 
                       meta.data, 
                       na.action = na.omit, scale = TRUE)
      
      res.temp <- anova.cca(rda.trans, by = "margin") %>% data.frame
      res.temp$R2 <- c(RsquareAdj(rda(count.trans ~ biorep + Condition(Station), meta.data))$r.squared,
                       RsquareAdj(rda(count.trans ~ Station + Condition(biorep), meta.data))$r.squared,
                       NA)
      res.temp$adj.R2 <- c(RsquareAdj(rda(count.trans ~ biorep + Condition(Station), meta.data))$adj.r.squared,
                           RsquareAdj(rda(count.trans ~ Station + Condition(biorep), meta.data))$adj.r.squared,
                           NA)
      res.temp$model <- "~ biorep + Station"
      res.temp$factor <- rownames(res.temp)
      
      rda.trans1 <- rda(count.trans ~ biorep + Station + Condition(Transect), 
                        meta.data,
                        na.action = na.omit, scale = TRUE)
      
      res.temp1 <- anova.cca(rda.trans1, by = "margin") %>% data.frame
      res.temp1$R2 <- c(RsquareAdj(rda(count.trans ~ biorep + Condition(Transect + Station), meta.data))$r.squared,
                        RsquareAdj(rda(count.trans ~ Station + Condition(Transect + biorep), meta.data))$r.squared,
                        NA)
      res.temp1$adj.R2 <- c(RsquareAdj(rda(count.trans ~ biorep + Condition(Transect + Station), meta.data))$adj.r.squared,
                            RsquareAdj(rda(count.trans ~ Station + Condition(Transect + biorep), meta.data))$adj.r.squared,
                            NA)
      res.temp1$model <- "~ biorep + Station + Condition(Transect)"
      res.temp1$factor <- rownames(res.temp1)
      
      rda.trans2 <- rda(count.trans ~ biorep + Transect,
                        meta.data,
                        na.action = na.omit, scale = TRUE)
      
      res.temp2 <- anova.cca(rda.trans2, by = "margin") %>% data.frame
      res.temp2$R2 <- c(RsquareAdj(rda(count.trans ~ biorep + Condition(Transect), meta.data))$r.squared,
                        RsquareAdj(rda(count.trans ~ Transect + Condition(biorep), meta.data))$r.squared,
                        NA)
      res.temp2$adj.R2 <- c(RsquareAdj(rda(count.trans ~ biorep + Condition(Transect), meta.data))$adj.r.squared,
                            RsquareAdj(rda(count.trans ~ Transect + Condition(biorep), meta.data))$adj.r.squared,
                            NA)
      res.temp2$model <- "~ biorep + Transect"
      res.temp2$factor <- rownames(res.temp2)
      
      res.temp$season <- list.season[j]
      res.temp1$season <- list.season[j]
      res.temp2$season <- list.season[j]
      
      
    }

    res.all <- rbind(res.temp, res.temp1, res.temp2)
    res.all$locus <- locus.name
    
    res.rda <- rbind(res.rda, res.all)
    
  }
  

}
    
write.csv(res.rda,
          file.path(here::here(), "02_Results", new_dir, "Godbout_RDA_replicate_effect.csv"))
 

df.res <- subset(res.rda, factor != "Residual" & !model %in% c("~ biorep + Station + Condition(nreads)", "~ biorep + Station")) %>%
  mutate(factor = recode(factor, 
                         "biorep" = "Replicate"),
         model = recode(model, 
                        "~ biorep + Station + Condition(Transect)" = "Within a given transect",
                        "~ biorep + Station + Condition(Transect + nreads)" = "Within a given transect",
                        "~ biorep + Transect + Condition(nreads)" = "Among transects",
                        "~ biorep + Transect" = "Among transects"),
         padj = p.adjust(.$Pr..F., "bonferroni")) %>%
  mutate(Signif = ifelse(padj <= 0.05 & adj.R2 >= 0, "yes", "no"),
         locus = factor(locus, levels = c("COI", "MiFishU", "16Schord")),
         season = factor(season, levels = c("Spring", "Summer", "Fall")))

plot.res.among <- subset(df.res, model == "Among transects") %>%
  ggplot(aes(x = factor, y = adj.R2, col = Signif)) +
  geom_point(size = 7.5, stroke = 10, position = position_dodge2(w = 0.75), shape = "-") +
  facet_grid(locus ~ season, scale = "free_x") +
  scale_color_manual(name = expression(italic("P-value")),
                     values = c("black", "red"), label = c("> 0.05", "< 0.05")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0, size = 10)) +
  guides(fill=guide_legend(override.aes=list(colour=c(no="grey",yes="red")))) +
  ylab(expression(adjusted~R^2)) + xlab("") +
  ggtitle("Among transects")

plot.res.among

plot.res.within <- subset(df.res, model == "Within a given transect") %>%
  ggplot(aes(x = factor, y = adj.R2, col = Signif)) +
  geom_point(size = 7.5, stroke = 10, position = position_dodge2(w = 0.75), shape = "-") +
  facet_grid(locus ~ season, scale = "free_x") +
  scale_color_manual(name = expression(italic("P-value")),
                     values = c("black", "red"), label = c("> 0.05", "< 0.05")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0, size = 10)) +
  guides(fill=guide_legend(override.aes=list(colour=c(no="grey",yes="red")))) +
  ylab(expression(adjusted~R^2)) + xlab("") +
  ggtitle("Within transect")

plot.res.within

plot.res.all <- ggarrange(plot.res.within, plot.res.among,
                          ncol = 2,
                          labels = LETTERS, 
                          common.legend = TRUE, legend = "right") +
  theme(plot.background = element_rect(fill = "white"))
plot.res.all


ggsave(file.path(here::here(), "02_Results", new_dir, 
                   "Godbout_RDA_replicate_effect.png"),
         plot.res.all,
         width = 5, height = 2, scale = 1.75)
  
  