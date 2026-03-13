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
#
# CL
# 2023-06
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

new_dir <- "05c_Spatial_diversity_with_OTUs_3pc_NewFigures"
dir.create(file.path(here::here(), "02_Results", new_dir),
           showWarnings = FALSE)

dir_align <- file.path(here::here(), "02_Results", 
                       "05b_Spatial_diversity_with_OTUs_3pc",
                       "Aligned_OTUs_fasta")


# Using glsADNe 0.1-4 ----------------------------------------------------------
# step 1 to 3 in previous code 03a_gslADNe_Spacial_Diversity_with_OTUs_3pc


# 4) Data transformation and RDA -----------------------------------------------

rm(list=setdiff(ls(), c("new_dir", "dir_align")))
gc()


# Function - 3D dist
distance3D <- function (df) {
  planiDist <- distm(df[,1:2])
  depthDist <- dist(df[3], method = "euclidean") %>% as.matrix()
  dist3D <- sqrt(planiDist^2 + depthDist^2) %>% as.dist()
  return(dist3D)
}

## 4.1) Spatial analyses using Hellinger-Mahalanobis ---------------------------

## create new directory

new_dir2 <- "Result_Hellinger-Mahalanobis"
dir.create(file.path(here::here(), "02_Results", new_dir, new_dir2),
           showWarnings = FALSE)

list.locus <- c("COI", "MiFishU", "16Schord")

## Intermediate scale -----------------------------------------------------------------

### 1) BDA ---------------------------------------------------------------------

project.name <- "BDA"
# n.class = 15

for(i in 1:length(list.locus)){
  
  locus.name <- list.locus[i]
  
  file.align <- file.path(dir_align,
                          paste("Aligned_", project.name, "_", locus.name,
                                ".fasta", sep = ""))
  
  ## Data 
  
  OTUs.info <- read.csv(file.path(here::here(), "00_Data", "07_OTUs_3pc_tables",
                                  paste(project.name, locus.name, 
                                        "OTUs_infos.csv", sep = "_")),
                        header = T, row.names = 1)
  
  OTUs.count <- read.csv(file.path(here::here(), "00_Data", "07_OTUs_3pc_tables",
                                   paste(project.name, locus.name, 
                                         "OTUs_count.csv", sep = "_")),
                         header = T, row.names = 1)
  
  metaData <- read.csv(file.path(here::here(), "00_Data", "00_ESVs_tables_clean",
                                 paste(project.name, locus.name, 
                                       "metadata_ok.csv", sep = "_")),
                       header = T, row.names = 1) %>% 
    mutate(Site = ifelse(Site == "CrB" & Latitude >= 48.60, "CrA", Site))

  temp.nreads <- data.frame(ID_ADNe = row.names(OTUs.count),
                            nreads = rowSums(OTUs.count))
  
  metaData <- left_join(metaData, temp.nreads)
  
  nreads <- data.frame(OTUs = colnames(OTUs.count),
                       nreads = colSums(OTUs.count))
  
  OTUs.info <- left_join(OTUs.info, nreads)
  
  ## Pairwise genetic distance 
  
  dna_dist <- read.alignment(file.align,
                             format = "fasta") %>%
    dist.alignment(., gap = T)
  
  object <- dna_dist %>% as.matrix %>% .[rowSums(is.na(.)) == 0, colSums(is.na(.)) == 0] %>%
    as.dist %>% dna_PC
  
  cv <- dna_PCcov(object)
  
  ## Data transformation 
  
  otus.trans <- NULL
  
  otus.trans[[1]] <- OTUs.count[,colnames(OTUs.count) %in% names(cv$cm)] %>%
    decostand(method = "hellinger") %>% as.matrix
  
  otus.trans[[2]] <- OTUs.count[,colnames(OTUs.count) %in% names(cv$cm)] %>%
    decostand(method = "hellinger") %>% as.matrix %*% cv$tr_mat
  
  
  ## RDA 
  
  for(k in 1:2) {
    
    rda.seq.depth <- rda(otus.trans[[k]] ~ nreads, metaData)
    res.seq.depth <- data.frame(anova.cca(rda.seq.depth) %>% data.frame,
                                R2 = RsquareAdj(rda.seq.depth)$r.squared,
                                adj.R2 = RsquareAdj(rda.seq.depth)$adj.r.squared,
                                factor = "Sequencing depth",
                                locus = locus.name,
                                project = project.name)
    
    write.csv(res.seq.depth, 
              file.path(here::here(), "02_Results", new_dir, new_dir2,
                        paste("Sequencing_depth_effect_", 
                              ifelse(k==1, "1trans", "2trans"), 
                              "_", project.name, "_", locus.name, 
                              ".csv", sep = "")))
    
    seq.depth.effect.p <- p.adjust(res.seq.depth$Pr..F.[1], 
                                   method = "bonferroni",
                                   n = 15)
    
    if(seq.depth.effect.p <= 0.05){
      rda.trans <- rda(otus.trans[[k]] ~ Longitude + Latitude + Depth_m + Condition(nreads), metaData,
                       na.action = na.omit, scale = TRUE)
      ## Result pRDA
      
      res.trans.temp <- anova.cca(rda.trans, by = "margin") %>% data.frame
      
      res.trans.temp$R2 <- c(RsquareAdj(rda(otus.trans[[k]] ~ Longitude + Condition(Latitude + Depth_m + nreads), metaData, na.action = na.omit))$r.squared,
                             RsquareAdj(rda(otus.trans[[k]] ~ Latitude + Condition(Longitude + Depth_m + nreads), metaData, na.action = na.omit))$r.squared,
                             RsquareAdj(rda(otus.trans[[k]] ~ Depth_m + Condition(Longitude + Latitude + nreads), metaData, na.action = na.omit))$r.squared,
                             "NA")
      res.trans.temp$adj.R2 <- c(RsquareAdj(rda(otus.trans[[k]] ~ Longitude + Condition(Latitude + Depth_m + nreads), metaData, na.action = na.omit))$adj.r.squared,
                                 RsquareAdj(rda(otus.trans[[k]] ~ Latitude + Condition(Longitude + Depth_m + nreads), metaData, na.action = na.omit))$adj.r.squared,
                                 RsquareAdj(rda(otus.trans[[k]] ~ Depth_m + Condition(Longitude + Latitude + nreads), metaData, na.action = na.omit))$adj.r.squared,
                                 "NA")
      res.trans.temp$model <- "pRDA"
      res.trans.temp$transfo <- ifelse(k == 1, "Hellinger", "Hellinger-Mahalanobis")
      
      
      rda.longitude <- rda(otus.trans[[k]] ~ Longitude + Condition(nreads), metaData,
                           na.action = na.omit, scale = TRUE)
      rda.latitude <- rda(otus.trans[[k]] ~ Latitude + Condition(nreads), metaData,
                          na.action = na.omit, scale = TRUE)
      rda.depth <-rda(otus.trans[[k]] ~ Depth_m + Condition(nreads), metaData,
                      na.action = na.omit, scale = TRUE)
      
    }else{
      rda.trans <- rda(otus.trans[[k]] ~ Longitude + Latitude + Depth_m, metaData,
                       na.action = na.omit, scale = TRUE)
      
      ## Result pRDA
      
      res.trans.temp <- anova.cca(rda.trans, by = "margin") %>% data.frame
      
      res.trans.temp$R2 <- c(RsquareAdj(rda(otus.trans[[k]] ~ Longitude + Condition(Latitude + Depth_m), metaData, na.action = na.omit))$r.squared,
                             RsquareAdj(rda(otus.trans[[k]] ~ Latitude + Condition(Longitude + Depth_m), metaData, na.action = na.omit))$r.squared,
                             RsquareAdj(rda(otus.trans[[k]] ~ Depth_m + Condition(Longitude + Latitude), metaData, na.action = na.omit))$r.squared,
                             "NA")
      res.trans.temp$adj.R2 <- c(RsquareAdj(rda(otus.trans[[k]] ~ Longitude + Condition(Latitude + Depth_m), metaData, na.action = na.omit))$adj.r.squared,
                                 RsquareAdj(rda(otus.trans[[k]] ~ Latitude + Condition(Longitude + Depth_m), metaData, na.action = na.omit))$adj.r.squared,
                                 RsquareAdj(rda(otus.trans[[k]] ~ Depth_m + Condition(Longitude + Latitude), metaData, na.action = na.omit))$adj.r.squared,
                                 "NA")
      res.trans.temp$model <- "pRDA"
      res.trans.temp$transfo <- ifelse(k == 1, "Hellinger", "Hellinger-Mahalanobis")
      
      
      rda.longitude <- rda(otus.trans[[k]] ~ Longitude, metaData,
                           na.action = na.omit, scale = TRUE)
      rda.latitude <- rda(otus.trans[[k]] ~ Latitude, metaData,
                          na.action = na.omit, scale = TRUE)
      rda.depth <-rda(otus.trans[[k]] ~ Depth_m, metaData,
                      na.action = na.omit, scale = TRUE)
    }
    
    
    ## Complete model
    
    res.trans.all <- anova.cca(rda.trans) %>% data.frame
    
    res.trans.all$R2 <- c(RsquareAdj(rda.trans)$r.squared, 
                          "NA")
    res.trans.all$adj.R2 <- c(RsquareAdj(rda.trans)$adj.r.squared, 
                              "NA")
    res.trans.all$model <- "Complete"
    res.trans.all$transfo <- ifelse(k == 1, "Hellinger", "Hellinger-Mahalanobis")
    
    # each factor independently

    res.glob <- rbind(anova.cca(rda.longitude) %>% data.frame,
                      anova.cca(rda.latitude)%>% data.frame,
                      anova.cca(rda.depth)%>% data.frame)
    res.glob$R2 <- rbind(RsquareAdj(rda.longitude)$r.squared,
                         NA,
                         RsquareAdj(rda.latitude)$r.squared,
                         NA,
                         RsquareAdj(rda.depth)$r.squared,
                         NA)
    res.glob$adj.R2 <- rbind(RsquareAdj(rda.longitude)$adj.r.squared,
                         NA,
                         RsquareAdj(rda.latitude)$adj.r.squared,
                         NA,
                         RsquareAdj(rda.depth)$adj.r.squared,
                         NA)
    res.glob$model <- rep("global",6)
    res.glob$transfo <- ifelse(k == 1, "Hellinger", "Hellinger-Mahalanobis")
    rownames(res.glob) <- c("Longitude", "Residual", "Latitude", "Residual1", "Depth", "Residual3")
    
    
    
    ## Output RDA results
    
    res.rda <- rbind(res.trans.all, res.trans.temp, res.glob)
    
    write.csv(res.rda,
              file.path(here::here(), "02_Results", new_dir, new_dir2,
                        paste("Result_RDA_", 
                              ifelse(k==1, "1trans", "2trans"), 
                              "_", project.name, "_", locus.name, 
                              ".csv", sep = "")))
    
    
    ## Plot RDA 
    
    eigen.val <- rda.trans$CCA$eig
    rda1 <- round((eigen.val[1]/sum(eigen.val))*100,2)
    rda2 <- round((eigen.val[2]/sum(eigen.val))*100,2)
    
    percentVar.rda <- c(rda1, rda2)
    
    env.arrows <- data.frame(scores(rda.trans)$biplot)
    env.arrows$name <- rownames(env.arrows)
    
    RDA.scores <- data.frame(rda.trans$CCA$u,
                             ID_ADNe = rownames(scores(rda.trans$CCA$u, display="site"))) %>%
      left_join(., metaData) 
      
    write.csv(RDA.scores,
              file.path(here::here(), "02_Results", new_dir, new_dir2,
                        paste("RDA_scores_sites_", ifelse(k == 1, "1trans", "2trans"), 
                              "_", project.name, "_", 
                              locus.name, ".csv", sep = "")))
    
    rda.trans.plot <- RDA.scores %>%
      ggplot(aes(x = RDA1, y = RDA2, col = Site)) +
      geom_hline(yintercept = 0, linetype = "dashed", col = "grey") +
      geom_vline(xintercept = 0, linetype = "dashed", col = "grey") +
      geom_point(size = 4) +
      xlab(paste0("RDA1: ", percentVar.rda[1], "%")) +
      ylab(paste0("RDA2: ", percentVar.rda[2], "%")) +
      theme_bw() +
      scale_color_brewer(name = "Region", palette = "Set1") +
      theme(
        axis.text = element_text(color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_blank()
      ) +
      geom_segment(data = env.arrows,
                   aes(x = 0, y = 0, xend = RDA1, yend = RDA2), col = "black",
                   arrow = arrow(length = unit(0.5, "cm"))) +
      geom_text(data = env.arrows, aes(x = RDA1, y = RDA2, label = name, col = NULL), 
                hjust = 0.5, vjust = 1) +
      ggtitle(paste(project.name, locus.name), 
              subtitle = paste("Transformation:", ifelse(k == 1, "Hellinger", "Hellinger-Mahalanobis")))
    
    if(k == 1){
      
      sp.score.trans <- scores(rda.trans, display = "species") %>% data.frame
      sp.score.trans$OTUs <- rownames(sp.score.trans)
      
    } else {
      
      replot.rda <- replot.ordiplot(plot(rda.trans), cv$inv_tr_mat, type = "text")
      sp.score.trans <- replot.rda$species %>% data.frame()
      sp.score.trans$OTUs <- names(cv$cm)
      
    }
    
    sp.score.trans <- left_join(sp.score.trans, OTUs.info)
    sp.score.trans$Name <- ifelse(is.na(sp.score.trans$Taxon), 
                                  sp.score.trans$OTUs, 
                                  sp.score.trans$Taxon)
    
    write.csv(sp.score.trans,
              file.path(here::here(), "02_Results", new_dir, new_dir2,
                        paste("RDA_scores_species_", ifelse(k == 1, "1trans", "2trans"), 
                              "_", project.name, "_", 
                              locus.name, ".csv", sep = "")))
    
    sp.trans.plot <- ggplot(sp.score.trans, aes(x = RDA1, y = RDA2, label = Name)) +
      geom_point(col = "red", shape = 4) +
      ggrepel::geom_text_repel(size = 3) + 
      theme_bw() +
      ggtitle(paste(project.name, locus.name), 
              subtitle = paste("Transformation:", ifelse(k == 1, "Hellinger", "Hellinger-Mahalanobis")))
    
    plot.trans <- ggarrange(rda.trans.plot, sp.trans.plot,
                            ncol = 2)
    
    ggsave(file.path(here::here(), "02_Results", new_dir, new_dir2,
                     paste("RDA_plot_", ifelse(k == 1, "1trans", "2trans"), 
                           "_", project.name, "_", 
                           locus.name, ".png", sep = "")),
           plot.trans,
           height = 3, width = 7, scale = 1.75)
    
  }
  
  ## Spatial autocorrelation 
  
  env.df <- metaData[, c("Longitude", "Latitude", "Depth_m")]
  rownames(env.df) <- metaData$ID_ADNe
  env.df <- env.df[rowSums(is.na(env.df)) == 0,]
  dist.depth <- vegdist(env.df$Depth_m, method = "euclidean")
  dist.2d <- distm(env.df[, c("Longitude", "Latitude")])
  dist.3d <- distance3D(env.df)
  
  
  dist.data <- NULL
  dist.data[[1]] <- vegdist(otus.trans[[1]][rownames(env.df),], method = "euclidean")
  dist.data[[2]] <- vegdist(otus.trans[[2]][rownames(env.df),], method = "euclidean")
  
  
  for(k in 1:2){
    
    correl.depth <- mpmcorrelogram(xdis = dist.data[[k]], 
                                   geodis = dist.depth,
                                   zdis = dist.2d, 
                                   print = F)
    
    correl.2d <- mpmcorrelogram(xdis = dist.data[[k]], 
                                geodis = dist.2d,
                                zdis = dist.depth,
                                print = F)
    
    correl.3d <- mpmcorrelogram(xdis = dist.data[[k]], 
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
    
    write.csv(df.correl,
              file.path(here::here(), "02_Results", new_dir, new_dir2,
                        paste("Result_MantelCorrel_", 
                              ifelse(k==1, "1trans", "2trans"), 
                              "_", project.name, "_", locus.name, 
                              ".csv", sep = "")))
    
  }
  
}

### 2) Leim.VRoy.2019 ----------------------------------------------------------

project.name <- "PPO.Leim.VRoy.2019"


for(i in 1:length(list.locus)){
  
  locus.name <- list.locus[i]
  
  file.align <- file.path(dir_align,
                          paste("Aligned_", project.name, "_", locus.name,
                                ".fasta", sep = ""))
  
  ## Data 
  
  OTUs.info <- read.csv(file.path(here::here(), "00_Data", "07_OTUs_3pc_tables",
                                  paste(project.name, locus.name, 
                                        "OTUs_infos.csv", sep = "_")),
                        header = T, row.names = 1)
  
  OTUs.count <- read.csv(file.path(here::here(), "00_Data", "07_OTUs_3pc_tables",
                                   paste(project.name, locus.name, 
                                         "OTUs_count.csv", sep = "_")),
                         header = T, row.names = 1)
  
  metaData <- read.csv(file.path(here::here(), "00_Data", "00_ESVs_tables_clean",
                                 paste(project.name, locus.name, 
                                       "metadata_ok.csv", sep = "_")),
                       header = T, row.names = 1)
  
  temp.nreads <- data.frame(ID_ADNe = row.names(OTUs.count),
                            nreads = rowSums(OTUs.count))
  
  metaData <- left_join(metaData, temp.nreads)
  
  nreads <- data.frame(OTUs = colnames(OTUs.count),
                       nreads = colSums(OTUs.count))
  
  OTUs.info <- left_join(OTUs.info, nreads)
  
  ## Pairwise genetic distance 
  
  dna_dist <- read.alignment(file.align,
                             format = "fasta") %>%
    dist.alignment(., gap = T)
  
  object <- dna_dist %>% as.matrix %>% .[rowSums(is.na(.)) == 0, colSums(is.na(.)) == 0] %>%
    as.dist %>% dna_PC
  
  cv <- dna_PCcov(object)
  
  ## Data transformation 
  
  otus.trans <- NULL
  
  otus.trans[[1]] <- OTUs.count[,colnames(OTUs.count) %in% names(cv$cm)] %>%
    decostand(method = "hellinger") %>% as.matrix
  
  otus.trans[[2]] <- OTUs.count[,colnames(OTUs.count) %in% names(cv$cm)] %>%
    decostand(method = "hellinger") %>% as.matrix %*% cv$tr_mat
  
  
  ## RDA 
  
  for(k in 1:2) {
    
    rda.seq.depth <- rda(otus.trans[[k]] ~ nreads, metaData)
    res.seq.depth <- data.frame(anova.cca(rda.seq.depth) %>% data.frame,
                                R2 = RsquareAdj(rda.seq.depth)$r.squared,
                                adj.R2 = RsquareAdj(rda.seq.depth)$adj.r.squared,
                                factor = "Sequencing depth",
                                locus = locus.name,
                                project = project.name)
    
    write.csv(res.seq.depth, 
              file.path(here::here(), "02_Results", new_dir, new_dir2,
                        paste("Sequencing_depth_effect_", 
                              ifelse(k==1, "1trans", "2trans"), 
                              "_", project.name, "_", locus.name, 
                              ".csv", sep = "")))
    
    seq.depth.effect.p <- p.adjust(res.seq.depth$Pr..F.[1], 
                                   method = "bonferroni",
                                   n = 15)
    
    if(seq.depth.effect.p <= 0.05){
      rda.trans <- rda(otus.trans[[k]] ~ Longitude + Latitude + Depth_m +
                         Temperature + Salinity + Condition(nreads), metaData,
                       na.action = na.omit, scale = TRUE)
      
      ## Result pRDA
      
      res.trans.temp <- anova.cca(rda.trans, by = "margin") %>% data.frame
      
      res.trans.temp$R2 <- c(RsquareAdj(rda(otus.trans[[k]] ~ Longitude + Condition(Latitude + Depth_m + Temperature + Salinity + nreads), metaData, na.action = na.omit, scale = TRUE))$r.squared,
                             RsquareAdj(rda(otus.trans[[k]] ~ Latitude + Condition(Longitude + Depth_m + Temperature + Salinity + nreads), metaData, na.action = na.omit, scale = TRUE))$r.squared,
                             RsquareAdj(rda(otus.trans[[k]] ~ Depth_m + Condition(Longitude + Latitude + Temperature + Salinity + nreads), metaData, na.action = na.omit, scale = TRUE))$r.squared,
                             RsquareAdj(rda(otus.trans[[k]] ~ Temperature + Condition(Longitude + Latitude + Depth_m + Salinity + nreads), metaData, na.action = na.omit, scale = TRUE))$r.squared,
                             RsquareAdj(rda(otus.trans[[k]] ~ Salinity + Condition(Longitude + Latitude + Depth_m + Temperature + nreads), metaData, na.action = na.omit, scale = TRUE))$r.squared,
                             "NA")
      res.trans.temp$adj.R2 <- c(RsquareAdj(rda(otus.trans[[k]] ~ Longitude + Condition(Latitude + Depth_m + Temperature + Salinity + nreads), metaData, na.action = na.omit, scale = TRUE))$adj.r.squared,
                                 RsquareAdj(rda(otus.trans[[k]] ~ Latitude + Condition(Longitude + Depth_m + Temperature + Salinity + nreads), metaData, na.action = na.omit, scale = TRUE))$adj.r.squared,
                                 RsquareAdj(rda(otus.trans[[k]] ~ Depth_m + Condition(Longitude + Latitude + Temperature + Salinity + nreads), metaData, na.action = na.omit, scale = TRUE))$adj.r.squared,
                                 RsquareAdj(rda(otus.trans[[k]] ~ Temperature + Condition(Longitude + Latitude + Depth_m + Salinity + nreads), metaData, na.action = na.omit, scale = TRUE))$adj.r.squared,
                                 RsquareAdj(rda(otus.trans[[k]] ~ Salinity + Condition(Longitude + Latitude + Depth_m + Temperature + nreads), metaData, na.action = na.omit, scale = TRUE))$adj.r.squared,
                                 "NA")
      
      res.trans.temp$model <- "pRDA"
      res.trans.temp$transfo <- ifelse(k == 1, "Hellinger", "Hellinger-Mahalanobis")
      
      # each factor independently
      
      rda.longitude <- rda(otus.trans[[k]] ~ Longitude + Condition(nreads), metaData,
                           na.action = na.omit, scale = TRUE)
      
      rda.latitude <- rda(otus.trans[[k]] ~ Latitude + Condition(nreads), metaData,
                          na.action = na.omit, scale = TRUE)
      
      rda.depth <-rda(otus.trans[[k]] ~ Depth_m + Condition(nreads), metaData,
                      na.action = na.omit, scale = TRUE)
      
      rda.temperature <-rda(otus.trans[[k]] ~ Temperature + Condition(nreads), metaData,
                            na.action = na.omit, scale = TRUE)
      
      rda.salinity <-rda(otus.trans[[k]] ~ Salinity + Condition(nreads), metaData,
                         na.action = na.omit, scale = TRUE)
      
    }else{
      rda.trans <- rda(otus.trans[[k]] ~ Longitude + Latitude + Depth_m +
                         Temperature + Salinity, metaData,
                       na.action = na.omit, scale = TRUE)
      
      ## Result pRDA
      
      res.trans.temp <- anova.cca(rda.trans, by = "margin") %>% data.frame
      
      res.trans.temp$R2 <- c(RsquareAdj(rda(otus.trans[[k]] ~ Longitude + Condition(Latitude + Depth_m + Temperature + Salinity), metaData, na.action = na.omit, scale = TRUE))$r.squared,
                             RsquareAdj(rda(otus.trans[[k]] ~ Latitude + Condition(Longitude + Depth_m + Temperature + Salinity), metaData, na.action = na.omit, scale = TRUE))$r.squared,
                             RsquareAdj(rda(otus.trans[[k]] ~ Depth_m + Condition(Longitude + Latitude + Temperature + Salinity), metaData, na.action = na.omit, scale = TRUE))$r.squared,
                             RsquareAdj(rda(otus.trans[[k]] ~ Temperature + Condition(Longitude + Latitude + Depth_m + Salinity), metaData, na.action = na.omit, scale = TRUE))$r.squared,
                             RsquareAdj(rda(otus.trans[[k]] ~ Salinity + Condition(Longitude + Latitude + Depth_m + Temperature), metaData, na.action = na.omit, scale = TRUE))$r.squared,
                             "NA")
      res.trans.temp$adj.R2 <- c(RsquareAdj(rda(otus.trans[[k]] ~ Longitude + Condition(Latitude + Depth_m + Temperature + Salinity), metaData, na.action = na.omit, scale = TRUE))$adj.r.squared,
                                 RsquareAdj(rda(otus.trans[[k]] ~ Latitude + Condition(Longitude + Depth_m + Temperature + Salinity), metaData, na.action = na.omit, scale = TRUE))$adj.r.squared,
                                 RsquareAdj(rda(otus.trans[[k]] ~ Depth_m + Condition(Longitude + Latitude + Temperature + Salinity), metaData, na.action = na.omit, scale = TRUE))$adj.r.squared,
                                 RsquareAdj(rda(otus.trans[[k]] ~ Temperature + Condition(Longitude + Latitude + Depth_m + Salinity), metaData, na.action = na.omit, scale = TRUE))$adj.r.squared,
                                 RsquareAdj(rda(otus.trans[[k]] ~ Salinity + Condition(Longitude + Latitude + Depth_m + Temperature), metaData, na.action = na.omit, scale = TRUE))$adj.r.squared,
                                 "NA")
      
      res.trans.temp$model <- "pRDA"
      res.trans.temp$transfo <- ifelse(k == 1, "Hellinger", "Hellinger-Mahalanobis")
      
      # each factor independently
      
      rda.longitude <- rda(otus.trans[[k]] ~ Longitude, metaData,
                           na.action = na.omit, scale = TRUE)
      
      rda.latitude <- rda(otus.trans[[k]] ~ Latitude, metaData,
                          na.action = na.omit, scale = TRUE)
      
      rda.depth <-rda(otus.trans[[k]] ~ Depth_m, metaData,
                      na.action = na.omit, scale = TRUE)
      
      rda.temperature <-rda(otus.trans[[k]] ~ Temperature, metaData,
                            na.action = na.omit, scale = TRUE)
      
      rda.salinity <-rda(otus.trans[[k]] ~ Salinity, metaData,
                         na.action = na.omit, scale = TRUE)
      
    }
    
    
    
    
    ## Complete model
    
    res.trans.all <- anova.cca(rda.trans) %>% data.frame
    
    res.trans.all$R2 <- c(RsquareAdj(rda.trans)$r.squared, 
                          "NA")
    res.trans.all$adj.R2 <- c(RsquareAdj(rda.trans)$adj.r.squared, 
                              "NA")
    res.trans.all$model <- "Complete"
    res.trans.all$transfo <- ifelse(k == 1, "Hellinger", "Hellinger-Mahalanobis")
    
    # each factor independently
    
    
    res.glob <- rbind(anova.cca(rda.longitude) %>% data.frame,
                      anova.cca(rda.latitude)%>% data.frame,
                      anova.cca(rda.depth)%>% data.frame,
                      anova.cca(rda.temperature)%>% data.frame,
                      anova.cca(rda.salinity)%>% data.frame)
    res.glob$R2 <- rbind(RsquareAdj(rda.longitude)$r.squared,
                         NA,
                         RsquareAdj(rda.latitude)$r.squared,
                         NA,
                         RsquareAdj(rda.depth)$r.squared,
                         NA,
                         RsquareAdj(rda.temperature)$r.squared,
                         NA,
                         RsquareAdj(rda.salinity)$r.squared,
                         NA)
    res.glob$adj.R2 <- rbind(RsquareAdj(rda.longitude)$adj.r.squared,
                             NA,
                             RsquareAdj(rda.latitude)$adj.r.squared,
                             NA,
                             RsquareAdj(rda.depth)$adj.r.squared,
                             NA,
                             RsquareAdj(rda.temperature)$adj.r.squared,
                             NA,
                             RsquareAdj(rda.salinity)$adj.r.squared,
                             NA)
    res.glob$model <- rep("global", 10)
    res.glob$transfo <- ifelse(k == 1, "Hellinger", "Hellinger-Mahalanobis")
    rownames(res.glob) <- c("Longitude", "Residual", 
                            "Latitude", "Residual1", 
                            "Depth", "Residual2",
                            "Temperature", "Residual3",
                            "Salinity", "Residual4")
    
    ## Output RDA results
    
    res.rda <- rbind(res.trans.all, res.trans.temp, res.glob)
    
    write.csv(res.rda,
              file.path(here::here(), "02_Results", new_dir, new_dir2,
                        paste("Result_RDA_", 
                              ifelse(k==1, "1trans", "2trans"), 
                              "_", project.name, "_", locus.name, 
                              ".csv", sep = "")))
    
    
    ## Plot RDA 
    
    eigen.val <- rda.trans$CCA$eig
    rda1 <- round((eigen.val[1]/sum(eigen.val))*100,2)
    rda2 <- round((eigen.val[2]/sum(eigen.val))*100,2)
    
    percentVar.rda <- c(rda1, rda2)
    
    env.arrows <- data.frame(scores(rda.trans)$biplot)
    env.arrows$name <- rownames(env.arrows)
    
    RDA.scores <- data.frame(rda.trans$CCA$u,
                             ID_ADNe = rownames(scores(rda.trans$CCA$u, display="site"))) %>%
      left_join(., metaData) 
    
    write.csv(RDA.scores,
              file.path(here::here(), "02_Results", new_dir, new_dir2,
                        paste("RDA_scores_sites_", ifelse(k == 1, "1trans", "2trans"), 
                              "_", project.name, "_", 
                              locus.name, ".csv", sep = "")))
    
    rda.trans.plot <- RDA.scores %>%
      ggplot(aes(x = RDA1, y = RDA2, col = Site)) +
      geom_hline(yintercept = 0, linetype = "dashed", col = "grey") +
      geom_vline(xintercept = 0, linetype = "dashed", col = "grey") +
      geom_point(size = 4) +
      xlab(paste0("RDA1: ", percentVar.rda[1], "%")) +
      ylab(paste0("RDA2: ", percentVar.rda[2], "%")) +
      theme_bw() +
      scale_color_brewer(name = "Region", palette = "Set1") +
      theme(
        axis.text = element_text(color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_blank()
      ) +
      geom_segment(data = env.arrows,
                   aes(x = 0, y = 0, xend = RDA1, yend = RDA2), col = "black",
                   arrow = arrow(length = unit(0.5, "cm"))) +
      geom_text(data = env.arrows, aes(x = RDA1, y = RDA2, label = name, col = NULL), 
                hjust = 0.5, vjust = 1) +
      ggtitle(paste(project.name, locus.name), 
              subtitle = paste("Transformation:", ifelse(k == 1, "Hellinger", "Hellinger-Mahalanobis")))
    
    if(k == 1){
      
      sp.score.trans <- scores(rda.trans, display = "species") %>% data.frame
      sp.score.trans$OTUs <- rownames(sp.score.trans)
      
    } else {
      
      replot.rda <- replot.ordiplot(plot(rda.trans), cv$inv_tr_mat, type = "text")
      sp.score.trans <- replot.rda$species %>% data.frame()
      sp.score.trans$OTUs <- names(cv$cm)
      
    }
    
    sp.score.trans <- left_join(sp.score.trans, OTUs.info)
    sp.score.trans$Name <- ifelse(is.na(sp.score.trans$Taxon), 
                                  sp.score.trans$OTUs, 
                                  sp.score.trans$Taxon)
    
    write.csv(sp.score.trans,
              file.path(here::here(), "02_Results", new_dir, new_dir2,
                        paste("RDA_scores_species_", ifelse(k == 1, "1trans", "2trans"), 
                              "_", project.name, "_", 
                              locus.name, ".csv", sep = "")))
    
    sp.trans.plot <- ggplot(sp.score.trans, aes(x = RDA1, y = RDA2, label = Name)) +
      geom_point(col = "red", shape = 4) +
      ggrepel::geom_text_repel(size = 3) + 
      theme_bw() +
      ggtitle(paste(project.name, locus.name), 
              subtitle = paste("Transformation:", ifelse(k == 1, "Hellinger", "Hellinger-Mahalanobis")))
    
    plot.trans <- ggarrange(rda.trans.plot, sp.trans.plot,
                            ncol = 2)
    
    ggsave(file.path(here::here(), "02_Results", new_dir, new_dir2,
                     paste("RDA_plot_", ifelse(k == 1, "1trans", "2trans"), 
                           "_", project.name, "_", 
                           locus.name, ".png", sep = "")),
           plot.trans,
           height = 3, width = 7, scale = 1.75)
    
  }
  
  ## Spatial autocorrelation 
  
  env.df <- metaData[, c("Longitude", "Latitude", "Depth_m")]
  rownames(env.df) <- metaData$ID_ADNe
  env.df <- env.df[rowSums(is.na(env.df)) == 0,]
  dist.depth <- vegdist(env.df$Depth_m, method = "euclidean")
  dist.2d <- distm(env.df[, c("Longitude", "Latitude")])
  dist.3d <- distance3D(env.df)
  
  
  dist.data <- NULL
  dist.data[[1]] <- vegdist(otus.trans[[1]][rownames(env.df),], method = "euclidean")
  dist.data[[2]] <- vegdist(otus.trans[[2]][rownames(env.df),], method = "euclidean")
  
  
  for(k in 1:2){
    
    correl.depth <- mpmcorrelogram(xdis = dist.data[[k]], 
                                   geodis = dist.depth,
                                   zdis = dist.2d, 
                                   print = F)
    
    correl.2d <- mpmcorrelogram(xdis = dist.data[[k]], 
                                geodis = dist.2d,
                                zdis = dist.depth,
                                print = F)
    
    correl.3d <- mpmcorrelogram(xdis = dist.data[[k]], 
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
    
    write.csv(df.correl,
              file.path(here::here(), "02_Results", new_dir, new_dir2,
                        paste("Result_MantelCorrel_", 
                              ifelse(k==1, "1trans", "2trans"), 
                              "_", project.name, "_", locus.name, 
                              ".csv", sep = "")))
    
  }
  
}



## Fine scale ------------------------------------------------------------------

### 3) Godbout -----------------------------------------------------------------

project.name <- c("Godbout")


for(i in 1:length(list.locus)){
  
  locus.name <- list.locus[i]
  
  file.align <- file.path(dir_align,
                          paste("Aligned_", project.name, "_", locus.name,
                                ".fasta", sep = ""))
  
  ## Data 
  
  OTUs.info <- read.csv(file.path(here::here(), "00_Data", "07_OTUs_3pc_tables",
                                  paste(project.name, locus.name, 
                                        "OTUs_infos.csv", sep = "_")),
                        header = T, row.names = 1)
  
  OTUs.count <- read.csv(file.path(here::here(), "00_Data", "07_OTUs_3pc_tables",
                                   paste(project.name, locus.name, 
                                         "OTUs_count.csv", sep = "_")),
                         header = T, row.names = 1)
  
  metaData <- read.csv(file.path(here::here(), "00_Data", "00_ESVs_tables_clean",
                                 paste(project.name, locus.name, 
                                       "metadata_ok.csv", sep = "_")),
                       header = T, row.names = 1)
  
  temp.nreads <- data.frame(ID_ADNe = row.names(OTUs.count),
                            nreads = rowSums(OTUs.count))
  
  metaData <- left_join(metaData, temp.nreads)
  
  metaData$Station <- sub("M.*", "M", metaData$ID_reception)
  metaData$Year <- sub("-.*", "", metaData$Date) %>% as.numeric
  metaData$Season <- ifelse(metaData$Date %in% c("2020-10-27", "2020-10-28"),
                            "Fall", 
                            ifelse(metaData$Date %in% c("2021-08-24", "2021-08-25"),
                                   "Summer",
                                   "Spring"))
  metaData <- metaData %>% 
    mutate(Station = recode(Station, '2M' = '02M', '6M' = '06M'))
  
  
  nreads <- data.frame(OTUs = colnames(OTUs.count),
                       nreads = colSums(OTUs.count))
  
  OTUs.info <- left_join(OTUs.info, nreads)
  
  ## Pairwise genetic distance 
  
  dna_dist <- read.alignment(file.align,
                             format = "fasta") %>%
    dist.alignment(., gap = T)
  
  object <- dna_dist %>% as.matrix %>% .[rowSums(is.na(.)) == 0, colSums(is.na(.)) == 0] %>%
    as.dist %>% dna_PC
  
  cv <- dna_PCcov(object)
  
  ## Data transformation 
  
  otus.trans <- NULL
  
  otus.trans[[1]] <- OTUs.count[,colnames(OTUs.count) %in% names(cv$cm)] %>%
    decostand(method = "hellinger") %>% as.matrix
  
  otus.trans[[2]] <- OTUs.count[,colnames(OTUs.count) %in% names(cv$cm)] %>%
    decostand(method = "hellinger") %>% as.matrix %*% cv$tr_mat
  
  
  ## RDA 
  
  list.season <- paste(unique(metaData$Season))
  
  
  for(k in 1:2) {
    
    df.correl.all <- mat.or.vec(0,0)
    
    for(j in 1:length(list.season)){
      
      season.select <- list.season[j]
      meta.data <- subset(metaData, Season == season.select)
      
      env.df <- meta.data[, c("Longitude", "Latitude", "Depth_m")]
      rownames(env.df) <- meta.data$ID_ADNe
      env.df <- env.df[rowSums(is.na(env.df)) == 0,]
      dist.depth <- vegdist(env.df$Depth_m, method = "euclidean")
      dist.2d <- distm(env.df[, c("Longitude", "Latitude")])
      dist.3d <- distance3D(env.df)
      
      otus.trans.season <- otus.trans[[k]][rownames(env.df),]
      
      dist.data <- vegdist(otus.trans.season, method = "euclidean")
      
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
      df.correl$Season <- list.season[j]
      
      df.correl.all <- rbind(df.correl.all, df.correl)
      
      
      ## RDA by season ---------------------------------------------------------
      
      rda.seq.depth <- rda(otus.trans.season ~ nreads, meta.data)
      res.seq.depth <- data.frame(anova.cca(rda.seq.depth) %>% data.frame,
                                  R2 = RsquareAdj(rda.seq.depth)$r.squared,
                                  adj.R2 = RsquareAdj(rda.seq.depth)$adj.r.squared,
                                  factor = "Sequencing depth",
                                  locus = locus.name,
                                  project = project.name)
      
      write.csv(res.seq.depth,
                file.path(here::here(), "02_Results", new_dir, new_dir2,
                          paste("Sequencing_depth_effect_", 
                                ifelse(k==1, "1trans", "2trans"), 
                                "_", project.name, 
                                "_", season.select,
                                "_", locus.name, 
                                ".csv", sep = "")))
      
      seq.depth.effect.p <- p.adjust(res.seq.depth$Pr..F.[1], 
                                     method = "bonferroni",
                                     n = 15)
      
      if(seq.depth.effect.p <= 0.05){
        rda.trans <- rda(otus.trans.season ~ Longitude + Latitude + Depth_m +
                           Temperature + Salinity + Condition(nreads), meta.data,
                         na.action = na.omit, scale = TRUE)
        
        ## Result pRDA
        
        res.trans.temp <- anova.cca(rda.trans, by = "margin") %>% data.frame
        
        res.trans.temp$R2 <- c(RsquareAdj(rda(otus.trans.season ~ Longitude + Condition(Latitude + Depth_m + Temperature + Salinity + nreads), meta.data, na.action = na.omit, scale = TRUE))$r.squared,
                               RsquareAdj(rda(otus.trans.season ~ Latitude + Condition(Longitude + Depth_m + Temperature + Salinity + nreads), meta.data, na.action = na.omit, scale = TRUE))$r.squared,
                               RsquareAdj(rda(otus.trans.season ~ Depth_m + Condition(Longitude + Latitude + Temperature + Salinity + nreads), meta.data, na.action = na.omit, scale = TRUE))$r.squared,
                               RsquareAdj(rda(otus.trans.season ~ Temperature + Condition(Longitude + Latitude + Depth_m + Salinity + nreads), meta.data, na.action = na.omit, scale = TRUE))$r.squared,
                               RsquareAdj(rda(otus.trans.season ~ Salinity + Condition(Longitude + Latitude + Depth_m + Temperature + nreads), meta.data, na.action = na.omit, scale = TRUE))$r.squared,
                               "NA")
        res.trans.temp$adj.R2 <- c(RsquareAdj(rda(otus.trans.season ~ Longitude + Condition(Latitude + Depth_m + Temperature + Salinity + nreads), meta.data, na.action = na.omit, scale = TRUE))$adj.r.squared,
                                   RsquareAdj(rda(otus.trans.season ~ Latitude + Condition(Longitude + Depth_m + Temperature + Salinity + nreads), meta.data, na.action = na.omit, scale = TRUE))$adj.r.squared,
                                   RsquareAdj(rda(otus.trans.season ~ Depth_m + Condition(Longitude + Latitude + Temperature + Salinity + nreads), meta.data, na.action = na.omit, scale = TRUE))$adj.r.squared,
                                   RsquareAdj(rda(otus.trans.season ~ Temperature + Condition(Longitude + Latitude + Depth_m + Salinity + nreads), meta.data, na.action = na.omit, scale = TRUE))$adj.r.squared,
                                   RsquareAdj(rda(otus.trans.season ~ Salinity + Condition(Longitude + Latitude + Depth_m + Temperature + nreads), meta.data, na.action = na.omit, scale = TRUE))$adj.r.squared,
                                   "NA")
        
        res.trans.temp$model <- "pRDA"
        res.trans.temp$transfo <- ifelse(k == 1, "Hellinger", "Hellinger-Mahalanobis")
        
        # each factor independently
        
        rda.longitude <- rda(otus.trans.season ~ Longitude + Condition(nreads), meta.data,
                             na.action = na.omit, scale = TRUE)
        rda.latitude <- rda(otus.trans.season ~ Latitude + Condition(nreads), meta.data,
                            na.action = na.omit, scale = TRUE)
        rda.depth <-rda(otus.trans.season ~ Depth_m + Condition(nreads), meta.data,
                        na.action = na.omit, scale = TRUE)
        rda.temperature <-rda(otus.trans.season ~ Temperature + Condition(nreads), meta.data,
                              na.action = na.omit, scale = TRUE)
        rda.salinity <-rda(otus.trans.season ~ Salinity + Condition(nreads), meta.data,
                           na.action = na.omit, scale = TRUE)
      }else{
        
        rda.trans <- rda(otus.trans.season ~ Longitude + Latitude + Depth_m +
                           Temperature + Salinity, meta.data,
                         na.action = na.omit, scale = TRUE)
        
        ## Result pRDA
        
        res.trans.temp <- anova.cca(rda.trans, by = "margin") %>% data.frame
        
        res.trans.temp$R2 <- c(RsquareAdj(rda(otus.trans.season ~ Longitude + Condition(Latitude + Depth_m + Temperature + Salinity), meta.data, na.action = na.omit, scale = TRUE))$r.squared,
                               RsquareAdj(rda(otus.trans.season ~ Latitude + Condition(Longitude + Depth_m + Temperature + Salinity), meta.data, na.action = na.omit, scale = TRUE))$r.squared,
                               RsquareAdj(rda(otus.trans.season ~ Depth_m + Condition(Longitude + Latitude + Temperature + Salinity), meta.data, na.action = na.omit, scale = TRUE))$r.squared,
                               RsquareAdj(rda(otus.trans.season ~ Temperature + Condition(Longitude + Latitude + Depth_m + Salinity), meta.data, na.action = na.omit, scale = TRUE))$r.squared,
                               RsquareAdj(rda(otus.trans.season ~ Salinity + Condition(Longitude + Latitude + Depth_m + Temperature), meta.data, na.action = na.omit, scale = TRUE))$r.squared,
                               "NA")
        res.trans.temp$adj.R2 <- c(RsquareAdj(rda(otus.trans.season ~ Longitude + Condition(Latitude + Depth_m + Temperature + Salinity), meta.data, na.action = na.omit, scale = TRUE))$adj.r.squared,
                                   RsquareAdj(rda(otus.trans.season ~ Latitude + Condition(Longitude + Depth_m + Temperature + Salinity), meta.data, na.action = na.omit, scale = TRUE))$adj.r.squared,
                                   RsquareAdj(rda(otus.trans.season ~ Depth_m + Condition(Longitude + Latitude + Temperature + Salinity), meta.data, na.action = na.omit, scale = TRUE))$adj.r.squared,
                                   RsquareAdj(rda(otus.trans.season ~ Temperature + Condition(Longitude + Latitude + Depth_m + Salinity), meta.data, na.action = na.omit, scale = TRUE))$adj.r.squared,
                                   RsquareAdj(rda(otus.trans.season ~ Salinity + Condition(Longitude + Latitude + Depth_m + Temperature), meta.data, na.action = na.omit, scale = TRUE))$adj.r.squared,
                                   "NA")
        
        res.trans.temp$model <- "pRDA"
        res.trans.temp$transfo <- ifelse(k == 1, "Hellinger", "Hellinger-Mahalanobis")
        
        # each factor independently
        
        rda.longitude <- rda(otus.trans.season ~ Longitude, meta.data,
                             na.action = na.omit, scale = TRUE)
        rda.latitude <- rda(otus.trans.season ~ Latitude, meta.data,
                            na.action = na.omit, scale = TRUE)
        rda.depth <-rda(otus.trans.season ~ Depth_m, meta.data,
                        na.action = na.omit, scale = TRUE)
        rda.temperature <-rda(otus.trans.season ~ Temperature, meta.data,
                              na.action = na.omit, scale = TRUE)
        rda.salinity <-rda(otus.trans.season ~ Salinity, meta.data,
                           na.action = na.omit, scale = TRUE)
      }
      
      
      
      
      ## Complete model
      
      res.trans.all <- anova.cca(rda.trans) %>% data.frame
      
      res.trans.all$R2 <- c(RsquareAdj(rda.trans)$r.squared, 
                            "NA")
      res.trans.all$adj.R2 <- c(RsquareAdj(rda.trans)$adj.r.squared, 
                                "NA")
      res.trans.all$model <- "Complete"
      res.trans.all$transfo <- ifelse(k == 1, "Hellinger", "Hellinger-Mahalanobis")
      
      
      # each factor independently
      
      res.glob <- rbind(anova.cca(rda.longitude) %>% data.frame,
                        anova.cca(rda.latitude)%>% data.frame,
                        anova.cca(rda.depth)%>% data.frame,
                        anova.cca(rda.temperature)%>% data.frame,
                        anova.cca(rda.salinity)%>% data.frame)
      res.glob$R2 <- rbind(RsquareAdj(rda.longitude)$r.squared,
                           NA,
                           RsquareAdj(rda.latitude)$r.squared,
                           NA,
                           RsquareAdj(rda.depth)$r.squared,
                           NA,
                           RsquareAdj(rda.temperature)$r.squared,
                           NA,
                           RsquareAdj(rda.salinity)$r.squared,
                           NA)
      res.glob$adj.R2 <- rbind(RsquareAdj(rda.longitude)$adj.r.squared,
                               NA,
                               RsquareAdj(rda.latitude)$adj.r.squared,
                               NA,
                               RsquareAdj(rda.depth)$adj.r.squared,
                               NA,
                               RsquareAdj(rda.temperature)$adj.r.squared,
                               NA,
                               RsquareAdj(rda.salinity)$adj.r.squared,
                               NA)
      res.glob$model <- rep("global", 10)
      res.glob$transfo <- ifelse(k == 1, "Hellinger", "Hellinger-Mahalanobis")
      rownames(res.glob) <- c("Longitude", "Residual", 
                              "Latitude", "Residual1", 
                              "Depth_m", "Residual2",
                              "Temperature", "Residual3",
                              "Salinity", "Residual4")
      
      ## Output RDA results
      
      res.rda <- rbind(res.trans.all, res.trans.temp, res.glob)
      
      write.csv(res.rda,
                file.path(here::here(), "02_Results", new_dir, new_dir2,
                          paste("Result_RDA_", 
                                ifelse(k==1, "1trans", "2trans"), 
                                "_", project.name, 
                                "_", season.select,
                                "_", locus.name, 
                                ".csv", sep = "")))
      
      
      ## Plot RDA 
      
      eigen.val <- rda.trans$CCA$eig
      rda1 <- round((eigen.val[1]/sum(eigen.val))*100,2)
      rda2 <- round((eigen.val[2]/sum(eigen.val))*100,2)
      
      percentVar.rda <- c(rda1, rda2)
      
      env.arrows <- data.frame(scores(rda.trans)$biplot)
      env.arrows$name <- rownames(env.arrows)
      
      RDA.scores <- data.frame(rda.trans$CCA$u,
                               ID_ADNe = rownames(scores(rda.trans$CCA$u, display="site"))) %>%
        left_join(., metaData) 
      
      write.csv(RDA.scores,
                file.path(here::here(), "02_Results", new_dir, new_dir2,
                          paste("RDA_scores_sites_", ifelse(k == 1, "1trans", "2trans"), 
                                "_", project.name, 
                                "_", season.select,
                                "_", locus.name, ".csv", sep = "")))
      
      rda.trans.plot <- RDA.scores %>%
        ggplot(aes(x = RDA1, y = RDA2, col = Station)) +
        geom_hline(yintercept = 0, linetype = "dashed", col = "grey") +
        geom_vline(xintercept = 0, linetype = "dashed", col = "grey") +
        geom_point(size = 4) +
        xlab(paste0("RDA1: ", percentVar.rda[1], "%")) +
        ylab(paste0("RDA2: ", percentVar.rda[2], "%")) +
        theme_bw() +
        scale_color_manual(name = "Transect", 
                           values = c("black", RColorBrewer::brewer.pal(7, "Set1")[c(1:5,7)])) +
        theme(
          axis.text = element_text(color = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_blank()
        ) +
        geom_segment(data = env.arrows,
                     aes(x = 0, y = 0, xend = RDA1, yend = RDA2), col = "black",
                     arrow = arrow(length = unit(0.5, "cm"))) +
        geom_text(data = env.arrows, aes(x = RDA1, y = RDA2, label = name, col = NULL), 
                  hjust = 0.5, vjust = 1) +
        ggtitle(paste(project.name, season.select, locus.name), 
                subtitle = paste("Transformation:", ifelse(k == 1, "Hellinger", "Hellinger-Mahalanobis")))
      
      if(k == 1){
        
        sp.score.trans <- scores(rda.trans, display = "species") %>% data.frame
        sp.score.trans$OTUs <- rownames(sp.score.trans)
        
      } else {
        
        replot.rda <- replot.ordiplot(plot(rda.trans), cv$inv_tr_mat, 
                                      type = "text")
        sp.score.trans <- replot.rda$species %>% data.frame()
        sp.score.trans$OTUs <- names(cv$cm)
        
      }
      
      sp.score.trans <- left_join(sp.score.trans, OTUs.info)
      sp.score.trans$Name <- ifelse(is.na(sp.score.trans$Taxon), 
                                    sp.score.trans$OTUs, 
                                    sp.score.trans$Taxon)
      
      write.csv(sp.score.trans,
                file.path(here::here(), "02_Results", new_dir, new_dir2,
                          paste("RDA_scores_species_", ifelse(k == 1, "1trans", "2trans"), 
                                "_", project.name, 
                                "_", season.select,
                                "_", locus.name, ".csv", sep = "")))
      
      sp.trans.plot <- ggplot(sp.score.trans, aes(x = RDA1, y = RDA2, label = Name)) +
        geom_point(col = "red", shape = 4) +
        ggrepel::geom_text_repel(size = 3) + 
        theme_bw() +
        ggtitle(paste(project.name, season.select, locus.name), 
                subtitle = paste("Transformation:", ifelse(k == 1, "Hellinger", "Hellinger-Mahalanobis")))
      
      plot.trans <- ggarrange(rda.trans.plot, sp.trans.plot,
                              ncol = 2)
      
      ggsave(file.path(here::here(), "02_Results", new_dir, new_dir2,
                       paste("RDA_plot_", ifelse(k == 1, "1trans", "2trans"), 
                             "_", project.name, 
                             "_", season.select,
                             "_", locus.name, ".png", sep = "")),
             plot.trans,
             height = 3, width = 7, scale = 1.5)
      
      
    }
    
    write.csv(df.correl.all,
              file.path(here::here(), "02_Results", new_dir, new_dir2,
                        paste("Result_MantelCorrel_", 
                              ifelse(k==1, "1trans", "2trans"), 
                              "_", project.name, "_", locus.name, 
                              ".csv", sep = "")))
    
  }
  
  
}

# 
# # 5) Comparison of transformed vs non-transformed dataset ----------------------
# 
# rm(list=setdiff(ls(), c("new_dir", "new_dir_align")))
# gc()
# 
# new_dir3 <- "Comparison_Transf_vs_NonTransf"
# dir.create(file.path(here::here(), "02_Results", new_dir, new_dir3),
#            showWarnings = FALSE)
# 
# list.projects <- c("BDA", 
#                    "PPO.Leim.VRoy.2019", 
#                    "Godbout")
# 
# list.locus <- c("COI", "MiFishU", "16Schord")
# 
# 
# 
# for(i in 1:length(list.locus)){
#   
#   locus.name <- list.locus[i]
#   
#   for(j in 1:length(list.projects)){
#     
#     project.name <- list.projects[j]
#     
#     file.align <- file.path(here::here(), "02_Results", new_dir, new_dir_align,
#                             paste("Aligned_", project.name, "_", locus.name,
#                                   ".fasta", sep = ""))
#     
#     ## Data 
#     
#     OTUs.count <- read.csv(file.path(here::here(), "00_Data", "07_OTUs_3pc_tables",
#                                      paste(project.name, locus.name, 
#                                            "OTUs_count.csv", sep = "_")),
#                            header = T, row.names = 1)
#     
#     
#     ## Pairwise genetic distance 
#     
#     dna_dist <- read.alignment(file.align,
#                                format = "fasta") %>%
#       dist.alignment(., gap = T)
#     
#     object <- dna_dist %>% as.matrix %>% .[rowSums(is.na(.)) == 0, colSums(is.na(.)) == 0] %>%
#       as.dist %>% dna_PC
#     
#     cv <- dna_PCcov(object)
#     
#     ## Data transformation 
#     
#     otus.trans.hell<- OTUs.count[,colnames(OTUs.count) %in% names(cv$cm)] %>%
#       decostand(method = "hellinger") %>% as.matrix
#     
#     otus.trans.hellMahala <- OTUs.count[,colnames(OTUs.count) %in% names(cv$cm)] %>%
#       decostand(method = "hellinger") %>% as.matrix %*% cv$tr_mat
#     
#     ## Dist Euclidean
#     
#     eucli.hell <- vegdist(otus.trans.hell, method = "euclidean")
#     eucli.hellMahala <- vegdist(otus.trans.hellMahala, method = "euclidean")
#     
#     df <- data.frame(Hellinger = eucli.hell[lower.tri(eucli.hell)],
#                      Hellinger_Mahalanobis = eucli.hellMahala[lower.tri(eucli.hellMahala)])
#     
#     ## Mantel correlation
#     
#     res.mantel <- mantel(eucli.hell, eucli.hellMahala)
#     
#     plot.dist <- ggplot(df, aes(x = Hellinger, y = Hellinger_Mahalanobis)) +
#       geom_point() +
#       # stat_smooth(method = "lm") +
#       theme_bw() +
#       ggtitle(paste(project.name, locus.name),
#               subtitle = paste("r_mantel = ", round(res.mantel$statistic, 2),
#                                " (P = ", round(res.mantel$signif, 3), ")",
#                                sep = ""))
#     
#     ggsave(file.path(here::here(), "02_Results", new_dir, new_dir3,
#                      paste("Mantel_transfo_", 
#                            project.name, "_", locus.name,
#                            ".png", sep = "")),
#            plot.dist, width = 3, height = 3, scale = 1.5)
#     
#   }
#   
# }
# 

# 7) Table formatting of RDA results -------------------------------------------


rm(list=setdiff(ls(), c("new_dir", "dir_align")))
gc()

list.projects <- c("BDA", 
                   "PPO.Leim.VRoy.2019", 
                   "Godbout_Fall",
                   "Godbout_Spring", 
                   "Godbout_Summer")

list.locus <- c("COI", "MiFishU", "16Schord")

dir.res <- "Result_Hellinger-Mahalanobis"

df.rda <- mat.or.vec(0,0)

for(i in 1:length(list.locus)){
  
  locus.name <- list.locus[i]
  
  for(j in 1:length(list.projects)){
    
    project.name <- list.projects[j]
    
    res.H <- read.csv(file.path(here::here(), "02_Results", new_dir, dir.res, 
                                paste("Result_RDA_", 
                                      "1trans", 
                                      "_", project.name, "_", locus.name, 
                                      ".csv", sep = "")))
    
    res.HM <- read.csv(file.path(here::here(), "02_Results", new_dir, dir.res, 
                                 paste("Result_RDA_", 
                                       "2trans", 
                                       "_", project.name, "_", locus.name, 
                                       ".csv", sep = "")))
    
    res.temp <- rbind(res.H, res.HM)
    res.temp$locus <- locus.name
    res.temp$project <- project.name
    
    res.temp$X <- str_replace(res.temp$X, "Model", "Complete_model")
    
    df.rda <- rbind(df.rda, res.temp)
    
  }
  
}

df.rda$project <- factor(df.rda$project, levels = c("PPO.Leim.VRoy.2019", 
                                                    "BDA", 
                                                    "Godbout_Spring", 
                                                    "Godbout_Summer",
                                                    "Godbout_Fall"))

df.rda$X <- recode(df.rda$X, "Depth_m" = "Depth")
df.rda$project <- recode(df.rda$project, "Large" = "Broad")

df.rda$X <- gsub('[[:digit:]]+', '', df.rda$X)
df.rda <- df.rda %>% group_by(locus, project) %>% 
  mutate(padj = p.adjust(Pr..F., "bonferroni") %>% as.numeric)

write.csv(df.rda,
          file.path(here::here(), "02_Results", new_dir,
                    "Result_RDA_all.csv"))


levels(df.rda$project) <- c("Broad\n(Summer)", 
                            "Intermediate\n(Summer)", 
                            "Fine\n(Spring)", "Fine\n(Summer)", "Fine\n(Fall)")

df.rda$X <- recode(df.rda$X, "Depth_m" = "Depth")

plot.res.rda <- df.rda %>% filter(!str_detect(X, "Residual")) %>%
  mutate(Signif = ifelse(padj <= 0.05 & adj.R2 >= 0, "yes", "no"),
         locus = factor(locus, levels = c("COI", "MiFishU", "16Schord"))) %>% 
  ggplot(aes(x = X, y = adj.R2, shape = transfo, fill = Signif)) +
  geom_point(size = 3.5, position = position_dodge2(w = 0.75)) +
  facet_grid(model * locus ~ project , scales = "free", space = "free_x") +
  scale_shape_manual(name = "Data tranformation", values = c(21, 24)) +
  scale_fill_manual(name = "P value", values = c("grey", "red"), label = c("> 0.05", "< 0.05")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0, size = 12)) +
  guides(fill=guide_legend(override.aes=list(colour=c(no="grey",yes="red")))) +
  ylab("adjusted R2") + xlab("Complete model and marginal effect") +
  ggtitle("RDA results")
plot.res.rda 

ggsave(file.path(here::here(), "02_Results", new_dir, 
                 "Result_RDA_All.png"),
       plot.res.rda, width = 7, height = 5, scale = 1.5)


plot.complete.mod <- subset(df.rda, transfo == "Hellinger-Mahalanobis" & model == "Complete") %>% 
  filter(!str_detect(X, "Residual")) %>%
  mutate(Signif = ifelse(padj <= 0.05 & adj.R2 >= 0, "yes", "no"),
         locus = factor(locus, levels = c("COI", "MiFishU", "16Schord"))) %>% 
  ggplot(aes(x = locus, y = adj.R2, shape = transfo, col = Signif)) +
  geom_point(size = 7.5, position = position_dodge2(w = 0.75), stroke = 10) +
  facet_grid("" ~ project, scales = "free_x", space = "free_x") +
  scale_shape_manual(name = "Data tranformation", values = c("-", "-")) +
  scale_colour_manual(name = expression(italic("P-value")), values = c("black", "red"), label = c("> 0.05", "< 0.05")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0)) +
  guides(colour=guide_legend(override.aes=list(size = c(no = 5,yes = 5),
                                               shape = c(no = "-", yes = "-"), 
                                               colour=c(no="black",yes="red"))),
         shape = "none") +
  ylab(expression("adjusted R"^2)) + 
  xlab("Locus") +
  ggtitle("Complete model")
plot.complete.mod


plot.global.mod <- subset(df.rda, transfo == "Hellinger-Mahalanobis" & model == "global") %>% 
  filter(!str_detect(X, "Residual")) %>%
  mutate(Signif = ifelse(padj <= 0.05 & adj.R2 >= 0, "yes", "no"),
         locus = factor(locus, levels = c("COI", "MiFishU", "16Schord"))) %>% 
  ggplot(aes(x = X, y = adj.R2, shape = transfo, col = Signif)) +
  geom_point(size = 7.5, position = position_dodge2(w = 0.75), stroke = 10) +
  facet_grid(locus ~ project, scales = "free_x") +
  scale_shape_manual(name = "Data tranformation", values = c("-", "-")) +
  scale_colour_manual(name = expression(italic("P-value")), values = c("black", "red"), label = c("> 0.05", "< 0.05")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0)) +
  guides(colour=guide_legend(override.aes=list(size = c(no = 5,yes = 5),
                                               shape = c(no = "-", yes = "-"), 
                                               colour=c(no="black",yes="red"))),
         shape = "none") +
  ylab(expression("adjusted R"^2)) + 
  xlab("Explanatory variable") +
  ggtitle("Global effect")
plot.global.mod

plot.marginal.mod <- subset(df.rda, transfo == "Hellinger-Mahalanobis" & model == "pRDA") %>% 
  filter(!str_detect(X, "Residual")) %>%
  mutate(Signif = ifelse(padj <= 0.05 & adj.R2 >= 0, "yes", "no"),
         locus = factor(locus, levels = c("COI", "MiFishU", "16Schord"))) %>% 
  ggplot(aes(x = X, y = adj.R2, shape = transfo, col = Signif)) +
  geom_point(size = 7.5, position = position_dodge2(w = 0.75), stroke = 10) +
  facet_grid(locus ~ project, scales = "free_x") +
  scale_shape_manual(name = "Data tranformation", values = c("-", "-")) +
  scale_colour_manual(name = expression(italic("P-value")), values = c("black", "red"), label = c("> 0.05", "< 0.05")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0)) +
  guides(colour=guide_legend(override.aes=list(size = c(no = 5,yes = 5),
                                               shape = c(no = "-", yes = "-"), 
                                               colour=c(no="black",yes="red"))),
         shape = "none") +
  ylab(expression("adjusted R"^2)) + 
  xlab("Explanatory variable") +
  ggtitle("Marginal effect")
plot.marginal.mod

fig.rda <- ggarrange(plot.complete.mod,
                     plot.global.mod,
                     plot.marginal.mod,
                     nrow = 3, 
                     labels = c("A.", 
                                "B.", 
                                "C."),
                     hjust = -1, vjust = 1.5,
                     heights = c(2, 3.5, 3.5), 
                     common.legend = F, legend = "right") +
  bgcolor("white") 
fig.rda

ggsave(file.path(here::here(), "02_Results", new_dir, 
                 "Result_RDA_OTUs_HM_only.png"),
       fig.rda, width = 5.5, height = 7, scale = 1.5)

### RDA sequencing depth --------------------------------------------------------
rm(list=setdiff(ls(), c("new_dir", "dir_align")))
gc()

list.projects <- c("BDA", 
                   "PPO.Leim.VRoy.2019", 
                   "Godbout_Fall",
                   "Godbout_Spring", 
                   "Godbout_Summer")

list.locus <- c("COI", "MiFishU", "16Schord")

dir.res <- "Result_Hellinger-Mahalanobis"

df.seq.depth <- mat.or.vec(0,0)

for(i in 1:length(list.locus)){
  
  locus.name <- list.locus[i]
  
  for(j in 1:length(list.projects)){
    
    project.name <- list.projects[j]
    
    res.temp <- read.csv(file.path(here::here(), "02_Results", new_dir, dir.res, 
                                   paste("Sequencing_depth_effect_", 
                                         "2trans", 
                                         "_", project.name, "_", locus.name, 
                                         ".csv", sep = "")))
    
    res.temp <- subset(res.temp, select = c("X", "Df", "Variance","F","Pr..F.","R2","adj.R2"))
    res.temp$factor <- "sequencing_depth"
    res.temp$locus <- locus.name
    res.temp$project <- project.name
    
    res.temp$X <- str_replace(res.temp$X, "Model", "Complete_model")
    
    df.seq.depth <- rbind(df.seq.depth, res.temp)
    
  }
  
}


df.seq.depth$project <- recode(df.seq.depth$project, 
                               "PPO.Leim.VRoy.2019" = "Broad\n(Summer)", 
                               "BDA" = "Intermediate\n(Summer)", 
                               "Godbout_Spring" = "Fine\n(Spring)", 
                               "Godbout_Summer" = "Fine\n(Summer)",
                               "Godbout_Fall" = "Fine\n(Fall)")

df.seq.depth$project <- factor(df.seq.depth$project, 
                               levels = c("Broad\n(Summer)", 
                                          "Intermediate\n(Summer)", 
                                          "Fine\n(Spring)", 
                                          "Fine\n(Summer)",
                                          "Fine\n(Fall)"))


df.seq.depth$padj <- p.adjust(df.seq.depth$Pr..F., "holm")

plot.res.seq.def <- subset(df.seq.depth, X == "Complete_model") %>%
  mutate(Signif = ifelse(.$padj <= 0.05 & adj.R2 >= 0, "yes", "no"),
         locus = factor(locus, levels = c("COI", "MiFishU", "16Schord"))) %>% 
  ggplot(aes(x = project, y = adj.R2, col = Signif)) +
  geom_point(size = 7.5, position = position_dodge2(w = 0.75), stroke = 10, shape = "-") +
  facet_wrap(~locus, scales = "free") +
  scale_colour_manual(name = expression(italic("P-value")), values = c("black", "red"), label = c("> 0.05", "< 0.05")) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0, size = 11)) +
  guides(colour=guide_legend(override.aes=list(size = c(no = 5,yes = 5),
                                               shape = c(no = "-", yes = "-"), 
                                               colour=c(no="black",yes="red"))),
         shape = "none") +
  ylab(expression("adjusted R"^2)) + 
  xlab("") +
  scale_y_continuous(limits = c(-0.03, 0.5))
plot.res.seq.def

### Reads taxonomic assignment -------------------------------------------------

list.projects <- c("BDA", 
                   "PPO.Leim.VRoy.2019", 
                   "Godbout")

list.locus <- c("COI", "MiFishU", "16Schord")

df.otus <- mat.or.vec(0,0)
df.meta <- mat.or.vec(0,0)

for(i in 1:length(list.projects)){
  
  project.name <- list.projects[i]
  
  for(j in 1:length(list.locus)){
    
    locus.name <- list.locus[j]
    
    data.MOTUs <- read.csv(file.path(here::here(), "00_Data", "07_OTUs_3pc_tables",
                                     paste(project.name, locus.name, 
                                           "OTUs_infos.csv", sep = "_")),
                           header = T, row.names = 1) %>%
      select(OTUs, sequence, Taxon, Levels, species, genus,
             family, order, class, phylum, kingdom)
    
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
                                          "Spring2021"))))
    OTUs.count <- read.csv(file.path(here::here(), "00_Data", "07_OTUs_3pc_tables",
                                     paste(project.name, locus.name, 
                                           "OTUs_count.csv", sep = "_")),
                           header = T, row.names = 1)
    
    temp.nreads <- data.frame(ID_ADNe = row.names(OTUs.count),
                              nreads = rowSums(OTUs.count))
    
    metaData <- left_join(metaData, temp.nreads)
    
    
    metaData <- metaData %>% mutate(Site = ifelse(Site == "CrB" & Latitude >= 48.60, "CrA", Site))
    
    data.MOTUs$locus <- locus.name
    data.MOTUs$project <- project.name
    
    data.OTUs <- read.csv(file.path(here::here(), "00_Data", "07_OTUs_3pc_tables",
                                    paste(project.name, locus.name, 
                                          "OTUs_count.csv", sep = "_")),
                          header = T, row.names = 1)
    
    # change OTU names by sequences
    
    names(data.OTUs)[match(data.MOTUs$OTUs, names(data.OTUs))] = data.MOTUs$sequence
    
    data.OTUs$ID_ADNe <- rownames(data.OTUs)
    data.OTUs <- left_join(data.OTUs, metaData[,c("ID_ADNe", "group")])
    
    data.OTUs$group <- paste(project.name, data.OTUs$group, sep ="_") 
    
    data.OTU.temp <- select(data.OTUs, -ID_ADNe)
    
    nreads <- data.OTU.temp %>% reshape::melt() %>% group_by(group, variable) %>%
      summarise(nreads = sum(value))
    names(nreads) <- c("group", "sequence", "nreads")
    
    motus <- left_join(data.MOTUs, nreads)
    
    df.otus <- rbind(df.otus, motus)
    
    df.meta <- rbind(df.meta, metaData)
    
  }
}


## seq depth by samples --------------------------------------------------------

plot.seq.depth <- df.meta %>%
  mutate(Project = recode(Project, 
                          "BDA" = "Intermediate\n(Summer)",
                          "PPO.Leim.VRoy.2019" = "Broad\n(Summer)",
                          "PPO.Godbout.2021" = "Fine",
                          "PPO.KMcGregor.2020" = "Fine\n(Fall)")) %>%
  mutate(Project = ifelse(Project == "Fine", 
                          ifelse(Date %in% c("2021-08-24", "2021-08-25"),
                                 "Fine\n(Summer)",
                                 "Fine\n(Spring)"), 
                          Project)) %>%
  mutate(Project = factor(Project, levels = c("Broad\n(Summer)", 
                                              "Intermediate\n(Summer)", 
                                              "Fine\n(Spring)",
                                              "Fine\n(Summer)",
                                              "Fine\n(Fall)")),
         Locus = factor(Locus, levels = c("COI", "MiFishU", "16Schord"))) %>%
  ggplot(aes(x = Project, y = nreads, fill = Project)) +
  geom_jitter(
    aes(color = Project), alpha = 0.5,
    position = position_jitter(0.2),
    size = 1.2
  ) +
  stat_summary(
    fun.data="mean_sdl",  fun.args = list(mult=1), 
    geom = "pointrange",  size = 0.4, 
  )+
  scale_color_manual(values = c("gray18", "gray28", "gray38", "gray38", "gray38"),
                     guide = guide_legend(override.aes = list(alpha = 0))) +
  scale_y_continuous(labels = scales::comma) +
  facet_wrap(~ Locus, scales = "free") +
  theme_bw() +
  xlab("") + ylab("Number of reads")  +
  theme(axis.text.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0, size = 11),
        legend.title = element_text(color = "transparent"),
        legend.text = element_text(color = "transparent"))


plot.seq.depth 



## otus ------------------------------------------------------------------------
df.otus$group <- factor(df.otus$group, 
                        levels = c("PPO.Leim.VRoy.2019_Forestville", "PPO.Leim.VRoy.2019_Colombier",
                                   "PPO.Leim.VRoy.2019_Betsiamites", "PPO.Leim.VRoy.2019_Manicouagan",
                                   "PPO.Leim.VRoy.2019_Baie-Comeau", "PPO.Leim.VRoy.2019_Godbout",
                                   "BDA_CrA", "BDA_CrB", "BDA_CrC",
                                   "Godbout_Fall2020", "Godbout_Spring2021", "Godbout_Summer2021"))


### By scale of sampling area ---------------------------------------------------

assigned.read.ok <- df.otus %>% 
  mutate(locus = factor(locus, levels = c("COI", "MiFishU", "16Schord")),
         assigned = factor(ifelse(is.na(phylum), "No", "Yes"), levels = c("Yes", "No")),
         group = paste0(map_chr(str_split(group, "_"), 2))) %>%
  mutate(scale = recode(group, 
                        "Forestville" = "Broad\n(Summer)",
                        "Colombier" = "Broad\n(Summer)",
                        "Betsiamites" = "Broad\n(Summer)",
                        "Manicouagan" = "Broad\n(Summer)",
                        "Baie-Comeau" = "Broad\n(Summer)",
                        "Godbout"= "Broad\n(Summer)",
                        "CrA" = "Intermediate\n(Summer)",
                        "CrB" = "Intermediate\n(Summer)", 
                        "CrC" = "Intermediate\n(Summer)",
                        "Spring2021" = "Fine\n(Spring)", 
                        "Summer2021" = "Fine\n(Summer)", 
                        "Fall2020" = "Fine\n(Fall)")) %>%
  mutate(scale = factor(scale, levels = c("Broad\n(Summer)", 
                                          "Intermediate\n(Summer)", 
                                          "Fine\n(Spring)",
                                          "Fine\n(Summer)",
                                          "Fine\n(Fall)"))) %>%
  group_by(assigned, locus, scale) %>%
  summarise(N = sum(nreads), All = "All") %>% ungroup() %>%
  group_by(scale, locus) %>%
  mutate(Total = sum(N),
         Percent = N/Total, 
         Lab = paste0(round(100*Percent,0),'%')) %>%
  ggplot(aes(x = scale, y = N, fill = assigned, width=.6)) +
  geom_bar(stat = "identity", position = position_stack()) +
  geom_text(aes(label = Lab), position = position_stack(vjust = .5), size = 2) +
  # geom_text(aes(y = Total, label = Total), vjust = -0.25, size = 3) +
  scale_fill_manual(values = c("red", "grey"),
                    name = "Taxonomic\nassignment") +
  theme_bw() +
  xlab("Sampling scale") + ylab("Number of reads") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0, size = 11)) +
  ggh4x::facet_grid2(~ locus, scales = "free", space = "free_x",
                     independent = "y") +
  scale_y_continuous(labels = scales::comma)

fig.seq.depth <- ggarrange(plot.seq.depth, 
                           plot.res.seq.def,
                           assigned.read.ok,
                           nrow = 3, align = "v",
                           heights = c(1.5, 1.5, 2),
                           labels = c("A. Sequencing depth per sample", 
                                      "B. Sequencing depth effect",
                                      "C. Taxonomic assignement"), 
                           hjust = 0, vjust = 0.5,
                           font.label = list(color = "black", size = 12)) +
  theme(plot.margin = margin(t = 0.5, l = 0.2, unit = "cm"),
        plot.background = element_rect(fill = "white"))

fig.seq.depth 


ggsave(file.path(here::here(), "02_Results", new_dir, "Fig_Sequencing_depth.png"),
       fig.seq.depth,
       width = 5.5, height = 4, scale = 1.75)




# Plot Mantel correlogram -----------------------------------------------------

rm(list=setdiff(ls(), c("new_dir", "dir_align")))
gc()


### For Large and intermediate scale
list.projects <- c("PPO.Leim.VRoy.2019", "BDA", "Godbout")
list.scale <- c("Broad\n(Summer)", "Intermediate\n(Summer)", "Fine")

df.correl.res <- mat.or.vec(0,0)


for(i in 1:length(list.projects)){
  
  project.name <- list.projects[i]
  scale.name <- list.scale[i]
  
  df.coi <- read.csv(file.path(here::here(), "02_Results", new_dir,
                               "Result_Hellinger-Mahalanobis",
                               paste("Result_MantelCorrel_2trans_",
                                     project.name, 
                                     "_COI.csv",sep = "")), 
                     row.names = 1)
  df.coi$locus <- "COI"
  
  df.mifish <- read.csv(file.path(here::here(), "02_Results", new_dir,
                                  "Result_Hellinger-Mahalanobis",
                                  paste("Result_MantelCorrel_2trans_",
                                        project.name, 
                                        "_MiFishU.csv",sep = "")), 
                        row.names = 1)
  df.mifish$locus <- "MiFishU"
  
  df.16s <- read.csv(file.path(here::here(), "02_Results", new_dir,
                               "Result_Hellinger-Mahalanobis",
                               paste("Result_MantelCorrel_2trans_",
                                     project.name, 
                                     "_16Schord.csv",sep = "")), 
                     row.names = 1)
  df.16s$locus <- "16Schord"
  
  df.correl <- rbind(df.coi, df.mifish, df.16s)
  
  df.correl$locus <- factor(df.correl$locus,
                            levels = c("COI", "MiFishU", "16Schord"))
  df.correl$distance <- factor(df.correl$distance, 
                               levels = c("Depth|2D", "2D|Depth", "3D"))
  df.correl$scale <- if(i == 3){
    paste(list.scale[i], "\n(", df.correl$Season, ")", sep = "")
  }else{
    list.scale[i]
  }
    

  df.correl.res <- if(i == 3) {
    rbind(df.correl.res, subset(df.correl, select = -c(Season)))
  }else{
    rbind(df.correl.res, df.correl)
  }
    
}


df.correl.res$scale <- factor(df.correl.res$scale,
                              levels = c("Broad\n(Summer)",
                                         "Intermediate\n(Summer)",
                                         "Fine\n(Spring)",
                                         "Fine\n(Summer)",
                                         "Fine\n(Fall)"))


correl.plot <- subset(df.correl.res, distance == "3D") %>%
  ggplot(., aes(x = class.index, y = Mantel.cor)) +
  geom_hline(yintercept = 0, color = "red") +
  geom_line() + 
  scale_x_continuous(labels = scales::label_comma()) +
  geom_point(aes(shape = Signif), size = 2) +
  scale_shape_manual(values = c(1,19)) + 
  facet_grid(locus ~ scale , scales = "free_x") +
  labs(y = expression("Mantel statistic (r"["M"]*")"),
       x = "Distance classes (m)") +
  theme_bw() +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0),
        axis.text.x = element_text(angle = 45,
                                   hjust = 1, vjust = 1)) +
  ggtitle("3D Mantel correlogram")

correl.plot

ggsave(file.path(here::here(), "02_Results", new_dir, 
                 "Mantel_correl_3D.png"),
       correl.plot, 
       height = 3, width = 5.5, scale = 1.5)

correl.plot.depth <- subset(df.correl.res, distance == "Depth|2D") %>%
  ggplot(., aes(x = class.index, y = Mantel.cor)) +
  geom_hline(yintercept = 0, color = "red") +
  geom_line() + 
  scale_x_continuous(labels = scales::label_comma()) +
  geom_point(aes(shape = Signif), size = 2) +
  scale_shape_manual(values = c(1,19)) + 
  facet_grid(locus ~ scale, scales = "free_x") +
  labs(y = expression("Mantel statistic (r"["M"]*")"),
       x = "Distance classes (m)") +
  theme_bw() +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0),
        axis.text.x = element_text(angle = 45,
                                   hjust = 1, vjust = 1)) +
  ggtitle("Depth|2D")

correl.plot.2d <- subset(df.correl.res, distance == "2D|Depth") %>%
  ggplot(., aes(x = class.index, y = Mantel.cor)) +
  geom_hline(yintercept = 0, color = "red") +
  geom_line() + 
  scale_x_continuous(labels = scales::label_comma()) +
  geom_point(aes(shape = Signif), size = 2) +
  scale_shape_manual(values = c(1,19)) + 
  facet_grid(locus ~ scale, scales = "free_x") +
  labs(y = expression("Mantel statistic (r"["M"]*")"),
       x = "Distance classes (m)") +
  theme_bw() +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0),
        axis.text.x = element_text(angle = 45,
                                   hjust = 1, vjust = 1)) +
  ggtitle("2D|Depth")

correl.plotbis <- ggarrange(correl.plot.depth, correl.plot.2d,
                            nrow = 2,
                            labels = LETTERS,
                            hjust = -1, vjust = 1.5)

ggsave(file.path(here::here(), "02_Results", new_dir, 
                 "Mantel_correl_2D-Depth.png"),
       correl.plotbis, 
       height = 5.5, width = 5.5, scale = 1.5)

