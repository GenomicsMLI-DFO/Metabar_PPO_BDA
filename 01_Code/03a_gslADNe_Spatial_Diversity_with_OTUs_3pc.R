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

new_dir <- "05b_Spatial_diversity_with_OTUs_3pc"
dir.create(file.path(here::here(), "02_Results", new_dir),
           showWarnings = FALSE)

new_dir_align <- "Aligned_OTUs_fasta"
dir.create(file.path(here::here(), "02_Results", new_dir, new_dir_align),
           showWarnings = FALSE)


# Using glsADNe 0.1-4 ----------------------------------------------------------


# 0) Number of samples per subproject ------------------------------------------

## initial 
temp <- read.csv(file.path(here::here(), 
                           "00_Data", "00_FileInfos",
                           "metaData_2019_to_2021.csv"))
temp <- subset(temp, Project %in% c("BDA", "PPO.Leim.VRoy.2019",
                                        "PPO.KMcGregor.2020",
                                        "PPO.Godbout.2021")) %>%
  mutate(group = ifelse(Project %in% c("BDA", "PPO.Leim.VRoy.2019"), 
                Project,
                ifelse(Project %in% c("PPO.KMcGregor.2020"),
                       "Fall2020", 
                       ifelse(Date %in% c("2021-08-24", "2021-08-25"),
                              "Summer2021",
                              "Spring2021"))))

table(temp$group, temp$Locus)

## number of reads 

cutadapt.res1 <- read.csv("/media/genobiwan/Storage/Projets/Metabar_Biodiversite_2021/00_Data/02a_Cutadapt/log/Cutadapt_Stats.csv")
cutadapt.res2 <- read.csv("/media/genobiwan/Extra_Storage/Projets/Metabar_Marin2021/00_Data/02a_Cutadapt/log/Cutadapt_Stats.csv")
cutadapt.res <- rbind(cutadapt.res1, cutadapt.res2)


temp2 <- left_join(temp, cutadapt.res, by = c("ID_ADNe" = "ID_labo", "Locus" = "Loci"))

read.init <- temp2 %>% group_by(Locus)%>%
  summarise(N.init = sum(Raw))
data.frame(read.init)


## after filatration

list.projects <- c("BDA", 
                   "PPO.Leim.VRoy.2019", 
                   "Godbout")

list.locus <- c("COI", "MiFishU", "16Schord")

metaData.all <- mat.or.vec(0,0)

for(i in 1:length(list.projects)){
 
   project.name = list.projects[i]
  
  for(j in 1:length(list.locus)){
    
    locus.name = list.locus[j]
    
    metaData <- read.csv(file.path(here::here(), "00_Data", "00_ESVs_tables_clean",
                                   paste(project.name, locus.name, "metadata_ok.csv",
                                         sep = "_")),
                         header = T)
    metaData$Project <- project.name
    metaData <- metaData %>%
      mutate(group = ifelse(Project %in% c("BDA", "PPO.Leim.VRoy.2019"), 
                            Project,
                            ifelse(Date %in% c("2020-10-27", "2020-10-28"),
                                   "Fall2020", 
                                   ifelse(Date %in% c("2021-08-24", "2021-08-25"),
                                          "Summer2021",
                                          "Spring2021"))))
    
    esv <- read.csv(file.path(here::here(), "00_Data", "00_ESVs_tables_clean",
                              paste(project.name, locus.name, "ESVs_ok.csv",
                                    sep = "_")),
                    header = T, row.names = 1) %>% rowSums() %>% data.frame
    colnames(esv) <- "Nreads"
    esv$ID_ADNe <- rownames(esv)
    
    metaData <- left_join(metaData, esv)
    
    metaData.all <- rbind(metaData.all, metaData)
    
  }
}

read.OTU <- metaData.all %>% group_by(Locus)%>%
  summarise(Nreads = sum(Nreads))
data.frame(read.OTU)

table(metaData.all$group, metaData.all$Locus)

# 1) Align sequences by locus and by project -----------------------------------

for(i in 1:length(list.locus)){
  
  locus.name <- list.locus[i]
  
  for(j in 1:length(list.projects)){
    
    project.name <- list.projects[j]
    
    ## Data --------------------------------------------------------------------

    data.MOTUs <- read.csv(file.path(here::here(), "00_Data", "07_OTUs_3pc_tables",
                                     paste(project.name, locus.name,
                                           "OTUs_infos.csv", sep = "_")),
                           header = T, row.names = 1)
    ## Genetic distance --------------------------------------------------------
    
    # fasta file
    x <- data.MOTUs$sequence %>% strsplit("")
    names(x) <- data.MOTUs$OTUs
    write.dna(x, file = file.path(here::here(), "02_Results", new_dir, new_dir_align,
                                  paste(project.name, "_", locus.name,
                                        ".fasta", sep = "")),
              format = "fasta", nbcol =  -1)
    rm(x)
    
    # alignment
    file.fasta <- file.path(here::here(), "02_Results", new_dir, new_dir_align,
                            paste(project.name, "_", locus.name,
                                  ".fasta", sep = ""))
    
    file.align <- file.path(here::here(), "02_Results", new_dir, new_dir_align,
                            paste("Aligned_", project.name, "_", locus.name,
                                  ".fasta", sep = ""))
    
    cmd <- paste("--auto", file.fasta, ">", file.align, sep = " ")
    system2("mafft", cmd)
    
    
  }
}



# 2) ESVs pairwise distance and PC scores (sequences) --------------------------


for(i in 1:length(list.locus)){
  
  locus.name <- list.locus[i]
  
  for(j in 1:length(list.projects)){
    
    project.name <- list.projects[j]  
    
    file.align <- file.path(here::here(), "02_Results", new_dir, new_dir_align,
                            paste("Aligned_", project.name, "_", locus.name,
                                  ".fasta", sep = ""))
    
    ## Data --------------------------------------------------------------------

    data.MOTUs <- read.csv(file.path(here::here(), "00_Data", "07_OTUs_3pc_tables",
                                     paste(project.name, locus.name, 
                                           "OTUs_infos.csv", sep = "_")),
                           header = T, row.names = 1)
    
    
    data.OTUs <- read.csv(file.path(here::here(), "00_Data", "07_OTUs_3pc_tables",
                                    paste(project.name, locus.name, 
                                          "OTUs_count.csv", sep = "_")),
                          header = T, row.names = 1)
    
    data.ESVs <- read.csv(file.path(here::here(), "00_Data", "00_ESVs_tables_clean",
                                    paste(project.name, locus.name, 
                                          "ESVs_ok.csv", sep = "_")),
                          header = T, row.names = 1)
    
    
    nreads <- data.frame(OTUs = colnames(data.OTUs),
                         nreads = colSums(data.OTUs))
    
    data.MOTUs <- left_join(data.MOTUs, nreads)
    
    ## Pairwise genetic distance -----------------------------------------------
    
    dna_dist <- read.alignment(file.align,
                               format = "fasta") %>%
      dist.alignment(., gap = T)
    
    
    
    ## Sequences diversity -----------------------------------------------------
    
    
    temp.align <- dna_dist %>% as.matrix %>% .[rowSums(is.na(.)) == 0, colSums(is.na(.)) == 0]
    if(nrow(temp.align) == 0){
      object.align <- matrix(nrow = 0, ncol = 0)
    } else {
      object.align <- temp.align %>% as.dist %>% dna_PC##(a=0.5, r = 0.7)
    }
    
    
    if(nrow(temp.align) == 0){
      seq.pcoa.align <- NULL
    } else {
      seq.pcoa.align <- data.frame(OTUs = names(object.align$cm), 
                                   object.align$u) %>%
        left_join(., data.MOTUs) %>%
        ggplot() +
        geom_point(aes(x = PC_1, y = PC_2, col = order, size = log1p(nreads)),
                   alpha = 0.7) +
        #   scale_color_manual(values = rainbow(n = 388)) +
        theme_classic() +
        ggtitle("OTUs - taxon assignment",
                subtitle = paste(object.align$u %>% nrow, " OTUs from ",
                                 ncol(data.ESVs), " ESVs, ",
                                 locus.name, ", ", project.name, sep = ""))
    }
    
    
    ggsave(file.path(here::here(), "02_Results", new_dir, new_dir_align,
                     paste("Fig_PC_sequences", 
                           "_", project.name, "_", locus.name, ".png", sep = "")),
           seq.pcoa.align,
           width = 15, height = 8, units = "cm",
           scale = 1.5)
  }
}


# 3) Most abundant taxon -------------------------------------------------------


## 3.1) Histogram of OTUs assigned to a taxon ----------------------------------

rm(list=setdiff(ls(), c("new_dir", "new_dir_align")))
gc()


list.projects <- c("BDA", 
                   "PPO.Leim.VRoy.2019", 
                   "Godbout")

list.locus <- c("COI", "MiFishU", "16Schord")

df.otus <- mat.or.vec(0,0)

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
    
  }
}

df.otus$group <- factor(df.otus$group, 
                        levels = c("PPO.Leim.VRoy.2019_Forestville", "PPO.Leim.VRoy.2019_Colombier",
                                   "PPO.Leim.VRoy.2019_Betsiamites", "PPO.Leim.VRoy.2019_Manicouagan",
                                   "PPO.Leim.VRoy.2019_Baie-Comeau", "PPO.Leim.VRoy.2019_Godbout",
                                   "BDA_CrA", "BDA_CrB", "BDA_CrC",
                                   "Godbout_Fall2020", "Godbout_Spring2021", "Godbout_Summer2021"))

taxon.reads.ok <- subset(df.otus, !locus == "COI") %>% 
  group_by(phylum, group, locus, project) %>%
  mutate(project = factor(project, levels = c("BDA", "PPO.Leim.VRoy.2019", "Godbout"))) %>%
  summarise(N = sum(nreads), All = "All") %>% ungroup() %>%
  group_by(group, locus) %>%
  mutate(Total = sum(N),
         Percent = N/Total, 
         Lab = paste0(round(100*Percent,0),'%')) %>%
  ggplot(aes(x = group, y = N, fill = phylum)) +
  geom_bar(stat = "identity", position = position_stack()) +
  geom_text(aes(label = Lab), position = position_stack(vjust = .5), size = 2) +
  # geom_text(aes(y = Total, label = Total), vjust = -0.25, size = 3) +
  scale_fill_brewer(palette = "Set2", na.value = "grey",
                    name = "") +
  theme_bw() +
  ggtitle("OTUs taxonomic assignment") +
  xlab("") + ylab("Number of reads") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~ locus*project, scales = "free") +
  scale_y_continuous(labels=function(x) format(x, big.mark = ",", decimal.mark = ".", scientific = FALSE))
taxon.reads.ok 

ggsave(file.path(here::here(), "02_Results", new_dir, "Nb_assign_OTUs_MiFish_16Schord.png"),
       taxon.reads.ok,
       width = 5, height = 6, scale = 1.75)


coi.plot <- subset(df.otus, locus == "COI") %>% 
  group_by(phylum, group, locus, project) %>%
  mutate(project = factor(project, levels = c("BDA", "PPO.Leim.VRoy.2019", "Godbout"))) %>%
  summarise(N = sum(nreads), All = "All") %>% ungroup() %>%
  group_by(group, locus) %>%
  mutate(Total = sum(N),
         Percent = N/Total, 
         Lab = paste0(round(100*Percent,0),'%')) %>%
  ggplot(aes(x = group, y = N, fill = phylum)) +
  geom_bar(stat = "identity", position = position_stack()) +
  geom_text(aes(label = Lab), position = position_stack(vjust = .5), size = 2) +
  # geom_text(aes(y = Total, label = Total), vjust = -0.25, size = 3) +
  # scale_fill_brewer(palette = "Set2", na.value = "grey",
  #                   name = "") +
  theme_bw() +
  ggtitle("OTUs taxonomic assignment") +
  xlab("") + ylab("Number of reads") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~ locus*project, scale = "free") +
  scale_y_continuous(##trans = "log10", 
    labels=function(x) format(x, big.mark = ",", decimal.mark = ".", scientific = FALSE))
coi.plot


ggsave(file.path(here::here(), "02_Results", new_dir, "Nb_assign_OTUs_COI.png"),
       coi.plot,
       width = 5, height = 3, scale = 1.75)

### By scale of sampling area ---------------------------------------------------

assigned.read.ok <- df.otus %>% 
  mutate(locus = factor(locus, levels = c("COI", "MiFishU", "16Schord")),
         assigned = factor(ifelse(is.na(phylum), "No", "Yes"), levels = c("Yes", "No")),
         group = paste0(map_chr(str_split(group, "_"), 2)),
         scale = ifelse(project == "PPO.Leim.VRoy.2019", "Large",
                        ifelse(project == "BDA", "Intermediate", "Fine"))) %>%
  mutate(group = factor(group, levels = c("Forestville", "Colombier", "Betsiamites", "Manicouagan", "Baie-Comeau", "Godbout",
                                          "CrA", "CrB", "CrC",
                                          "Fall2020", "Spring2021", "Summer2021")),
         scale = factor(scale, levels = c("Large", "Intermediate", "Fine"))) %>%
  group_by(assigned, group, locus, scale) %>%
  summarise(N = sum(nreads), All = "All") %>% ungroup() %>%
  group_by(scale, group, locus) %>%
  mutate(Total = sum(N),
         Percent = N/Total, 
         Lab = paste0(round(100*Percent,0),'%')) %>%
  ggplot(aes(x = group, y = N, fill = assigned)) +
  geom_bar(stat = "identity", position = position_stack()) +
  geom_text(aes(label = Lab), position = position_stack(vjust = .5), size = 2) +
  # geom_text(aes(y = Total, label = Total), vjust = -0.25, size = 3) +
  scale_fill_manual(values = c("red", "grey"),
                    name = "Taxonomic\nassignment") +
  theme_bw() +
  xlab("Sampling region or season") + ylab("Number of reads") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0, size = 12)) +
  ggh4x::facet_grid2(locus ~ scale, scales = "free", space = "free_x",
                     independent = "y") +
  scale_y_continuous(labels=function(x) format(x, big.mark = ",", decimal.mark = ".", scientific = FALSE))
assigned.read.ok 

ggsave(file.path(here::here(), "02_Results", new_dir, "Nb_assign_OTUs_all.png"),
       assigned.read.ok,
       width = 5, height = 3, scale = 1.75)


## 3.2) 10 most abundant taxon ----------------------------------------------------

graph.10.taxon <- df.otus %>% 
  mutate(Taxon = ifelse(is.na(Taxon), "Unknown", Taxon),
         locus = factor(locus, levels = c("COI", "MiFishU", "16Schord")),
         project = factor(project, levels = c("BDA", "PPO.Leim.VRoy.2019", "Godbout"))) %>%
  filter(Taxon != "Unknown") %>% # Change and remove unknown taxon
  group_by(Taxon, group, locus, project, phylum, class) %>%
  summarise(Nreads = sum(nreads)) %>% # Group by taxon per group and sum
  filter(Nreads > 0) %>%
  group_by(locus, project, phylum, class, group) %>%
  arrange(desc(Nreads)) %>% 
  group_by(group, locus) %>%
  slice_head(n = 10) %>% 
  ggplot(aes(x = group, y = Taxon, fill = Nreads)) + 
  geom_bin2d() +
  scale_fill_distiller(palette = "Spectral", trans = "log10", na.value = "white") +
  facet_grid(phylum ~ locus * project, scale = "free", space = "free") +
  theme_bw() +
  labs(title = "10 most abundant taxon per group") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
        axis.text.y = element_text(face = "italic", size = 8, vjust = 0.5),
        axis.ticks.y = element_blank(),
        strip.text.y = element_text(angle = 0, size = 8), 
        panel.spacing.y=unit(0.1, "lines"),
        legend.position = "right",
        legend.text = element_text(size = 8),
        panel.grid = element_blank())

graph.10.taxon 

ggsave(file.path(here::here(), "02_Results", new_dir, "10_most_abudant_taxon.png"),
       graph.10.taxon,
       width = 5, height = 6, scale = 2)

graph.5.taxon <- df.otus %>% 
  mutate(Taxon = ifelse(is.na(Taxon), "Unknown", Taxon),
         locus = factor(locus, levels = c("COI", "MiFishU", "16Schord")),
         project = factor(project, levels = c("BDA", "PPO.Leim.VRoy.2019", "Godbout"))) %>%
  filter(Taxon != "Unknown") %>% # Change and remove unknown taxon
  group_by(Taxon, group, locus, project, phylum, class) %>%
  summarise(Nreads = sum(nreads)) %>% # Group by taxon per group and sum
  filter(Nreads > 0) %>%
  group_by(locus, project, phylum, class, group) %>%
  arrange(desc(Nreads)) %>% 
  group_by(group, locus) %>%
  slice_head(n = 5) %>% 
  ggplot(aes(x = group, y = Taxon, fill = Nreads)) + 
  geom_bin2d() +
  scale_fill_distiller(palette = "Spectral", trans = "log10", na.value = "white") +
  facet_grid(phylum ~ locus * project, scale = "free", space = "free") +
  theme_bw() +
  labs(title = "5 most abundant taxon per group") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
        axis.text.y = element_text(face = "italic", size = 8, vjust = 0.5),
        axis.ticks.y = element_blank(),
        strip.text.y = element_text(angle = 0, size = 8), 
        panel.spacing.y=unit(0.1, "lines"),
        legend.position = "right",
        legend.text = element_text(size = 8),
        panel.grid = element_blank())

graph.5.taxon 

ggsave(file.path(here::here(), "02_Results", new_dir, "05_most_abudant_taxon.png"),
       graph.5.taxon,
       width = 5, height = 5, scale = 2)


df.otus$project <- factor(df.otus$project, levels = c("PPO.Leim.VRoy.2019", "BDA", "Godbout"))
levels(df.otus$project) <- c("Large", "Intermediate", "Fine")

levels(df.otus$group) <- c("Forestville", "Colombier",
                           "Betsiamites", "Manicouagan",
                           "Baie-Comeau", "Godbout",
                           "CrA", "CrB", "CrC",
                           "Fall2020", "Spring2021", "Summer2021")


graph.all.taxon <- df.otus %>% 
  mutate(Taxon = ifelse(is.na(Taxon), "Unknown", Taxon),
         locus = factor(locus, levels = c("COI", "MiFishU", "16Schord"))) %>%
  group_by(Taxon, group, locus, project, phylum, class) %>%
  summarise(Nreads = sum(nreads)) %>% # Group by taxon per group and sum
  filter(Nreads > 0) %>%
  ggplot(aes(x = group, y = Taxon, fill = Nreads)) + 
  geom_bin2d() +
  scale_fill_distiller(palette = "Spectral", trans = "log10", na.value = "white") +
  facet_grid(phylum ~ locus * project, scale = "free", space = "free") +
  theme_bw() +
  labs(title = "Taxon abundance per group") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
        axis.text.y = element_text(face = "italic", size = 8, vjust = 0.5),
        axis.ticks.y = element_blank(),
        strip.text.y = element_text(angle = 0, size = 8), 
        panel.spacing.y=unit(0.1, "lines"),
        legend.position = "right",
        legend.text = element_text(size = 8),
        panel.grid = element_blank())

graph.all.taxon 

ggsave(file.path(here::here(), "02_Results", new_dir, "All_taxon.png"),
       graph.all.taxon,
       width = 6, height = 8, scale = 2)




# 4) Data transformation and RDA -----------------------------------------------

rm(list=setdiff(ls(), c("new_dir", "new_dir_align")))
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
  
  file.align <- file.path(here::here(), "02_Results", new_dir, new_dir_align,
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
  
  
  ## ORDISTEP
  

  ## RDA 
  
  for(k in 1:2) {
    
    rda.trans <- rda(otus.trans[[k]] ~ Longitude + Latitude + Depth_m, metaData,
                     na.action = na.omit)
 
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
    
    
    ## Complete model
    
    res.trans.all <- anova.cca(rda.trans) %>% data.frame
    
    res.trans.all$R2 <- c(RsquareAdj(rda.trans)$r.squared, 
                          "NA")
    res.trans.all$adj.R2 <- c(RsquareAdj(rda.trans)$adj.r.squared, 
                              "NA")
    res.trans.all$model <- "Complete"
    res.trans.all$transfo <- ifelse(k == 1, "Hellinger", "Hellinger-Mahalanobis")
    
    ## Output RDA results
    
    res.rda <- rbind(res.trans.all, res.trans.temp)
    
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
    
    rda.trans.plot <- data.frame(rda.trans$CCA$u,
                                 ID_ADNe = rownames(scores(rda.trans$CCA$u, display="site"))) %>%
      left_join(., metaData) %>%
      ggplot(aes(x = RDA1, y = RDA2, col = Site)) +
      geom_hline(yintercept = 0, linetype = "dashed", col = "grey") +
      geom_vline(xintercept = 0, linetype = "dashed", col = "grey") +
      geom_point(size = 4) +
      xlab(paste0("RDA1: ", percentVar.rda[1], "%")) +
      ylab(paste0("RDA2: ", percentVar.rda[2], "%")) +
      theme_bw() +
      scale_color_brewer(palette = "Set1") +
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
           height = 3, width = 6, scale = 1.75)
    
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
    
    
    correl.plot <- df.correl %>%
      ggplot(., aes(x = class.index, y = Mantel.cor)) +
      geom_hline(yintercept = 0, color = "red") +
      geom_line() + 
      scale_x_continuous(labels = scales::label_comma()) +
      geom_point(aes(shape = Signif, size = 1)) +
      scale_shape_manual(values = c(1,19)) + 
      facet_grid(~ distance, scales = "free_x") +
      xlab("Distance (m)") +
      ylab("Mantel correlation") +
      theme_bw() +
      theme(legend.position = "none",
            strip.background = element_blank(),
            strip.text = element_text(hjust = 0)) +
      ggtitle(paste(project.name, locus.name),
              subtitle = paste("Transformation:", ifelse(k==1, "Hellinger", "Hellinger-Mahalanobis"), sep = " "))
    
    ggsave(file.path(here::here(), "02_Results", new_dir, new_dir2,
                     paste("Mantel_correl_", ifelse(k==1, "1trans", "2trans"), "_", project.name, "_", locus.name, ".png", sep = "")),
           correl.plot, 
           height = 2, width = 6, scale = 1.5)
    
  }
  
}

### 2) Leim.VRoy.2019 ----------------------------------------------------------

project.name <- "PPO.Leim.VRoy.2019"


for(i in 1:length(list.locus)){
  
  locus.name <- list.locus[i]
  
  file.align <- file.path(here::here(), "02_Results", new_dir, new_dir_align,
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
    
    rda.trans <- rda(otus.trans[[k]] ~ Longitude + Latitude + Depth_m +
                       Temperature + Salinity, metaData,
                     na.action = na.omit)
    
    ## Result pRDA
    
    res.trans.temp <- anova.cca(rda.trans, by = "margin") %>% data.frame
    
    res.trans.temp$R2 <- c(RsquareAdj(rda(otus.trans[[k]] ~ Longitude + Condition(Latitude + Depth_m + Temperature + Salinity), metaData, na.action = na.omit))$r.squared,
                           RsquareAdj(rda(otus.trans[[k]] ~ Latitude + Condition(Longitude + Depth_m + Temperature + Salinity), metaData, na.action = na.omit))$r.squared,
                           RsquareAdj(rda(otus.trans[[k]] ~ Depth_m + Condition(Longitude + Latitude + Temperature + Salinity), metaData, na.action = na.omit))$r.squared,
                           RsquareAdj(rda(otus.trans[[k]] ~ Temperature + Condition(Longitude + Latitude + Depth_m + Salinity), metaData, na.action = na.omit))$r.squared,
                           RsquareAdj(rda(otus.trans[[k]] ~ Salinity + Condition(Longitude + Latitude + Depth_m + Temperature), metaData, na.action = na.omit))$r.squared,
                           "NA")
    res.trans.temp$adj.R2 <- c(RsquareAdj(rda(otus.trans[[k]] ~ Longitude + Condition(Latitude + Depth_m + Temperature + Salinity), metaData, na.action = na.omit))$adj.r.squared,
                               RsquareAdj(rda(otus.trans[[k]] ~ Latitude + Condition(Longitude + Depth_m + Temperature + Salinity), metaData, na.action = na.omit))$adj.r.squared,
                               RsquareAdj(rda(otus.trans[[k]] ~ Depth_m + Condition(Longitude + Latitude + Temperature + Salinity), metaData, na.action = na.omit))$adj.r.squared,
                               RsquareAdj(rda(otus.trans[[k]] ~ Temperature + Condition(Longitude + Latitude + Depth_m + Salinity), metaData, na.action = na.omit))$adj.r.squared,
                               RsquareAdj(rda(otus.trans[[k]] ~ Salinity + Condition(Longitude + Latitude + Depth_m + Temperature), metaData, na.action = na.omit))$adj.r.squared,
                               "NA")
    
    res.trans.temp$model <- "pRDA"
    res.trans.temp$transfo <- ifelse(k == 1, "Hellinger", "Hellinger-Mahalanobis")
    
    
    ## Complete model
    
    res.trans.all <- anova.cca(rda.trans) %>% data.frame
    
    res.trans.all$R2 <- c(RsquareAdj(rda.trans)$r.squared, 
                          "NA")
    res.trans.all$adj.R2 <- c(RsquareAdj(rda.trans)$adj.r.squared, 
                              "NA")
    res.trans.all$model <- "Complete"
    res.trans.all$transfo <- ifelse(k == 1, "Hellinger", "Hellinger-Mahalanobis")
    
    ## Output RDA results
    
    res.rda <- rbind(res.trans.all, res.trans.temp)
    
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
    
    rda.trans.plot <- data.frame(rda.trans$CCA$u,
                                 ID_ADNe = rownames(scores(rda.trans$CCA$u, display="site"))) %>%
      left_join(., metaData) %>%
      ggplot(aes(x = RDA1, y = RDA2, col = Site)) +
      geom_hline(yintercept = 0, linetype = "dashed", col = "grey") +
      geom_vline(xintercept = 0, linetype = "dashed", col = "grey") +
      geom_point(size = 4) +
      xlab(paste0("RDA1: ", percentVar.rda[1], "%")) +
      ylab(paste0("RDA2: ", percentVar.rda[2], "%")) +
      theme_bw() +
      scale_color_brewer(palette = "Set1") +
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
           height = 3, width = 6, scale = 1.75)
    
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
    
    correl.plot <- df.correl %>%
      ggplot(., aes(x = class.index, y = Mantel.cor)) +
      geom_hline(yintercept = 0, color = "red") +
      geom_line() + 
      scale_x_continuous(labels = scales::label_comma()) +
      geom_point(aes(shape = Signif, size = 1)) +
      scale_shape_manual(values = c(1,19)) + 
      facet_grid(~ distance, scales = "free_x") +
      xlab("Distance (m)") +
      ylab("Mantel correlation") +
      theme_bw() +
      theme(legend.position = "none",
            strip.background = element_blank(),
            strip.text = element_text(hjust = 0)) +
      ggtitle(paste(project.name, locus.name),
              subtitle = paste("Transformation:", ifelse(k==1, "Hellinger", "Hellinger-Mahalanobis"), sep = " "))
    
    ggsave(file.path(here::here(), "02_Results", new_dir, new_dir2,
                     paste("Mantel_correl_", ifelse(k==1, "1trans", "2trans"), "_", project.name, "_", locus.name, ".png", sep = "")),
           correl.plot, 
           height = 2, width = 6, scale = 1.5)
    
  }
  
}



## Fine scale ------------------------------------------------------------------

### 3) Godbout -----------------------------------------------------------------

project.name <- c("Godbout")


for(i in 1:length(list.locus)){
  
  locus.name <- list.locus[i]
  
  file.align <- file.path(here::here(), "02_Results", new_dir, new_dir_align,
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
  
  metaData$Station <- sub("M.*", "M", metaData$ID_reception)
  metaData$Year <- sub("-.*", "", metaData$Date) %>% as.numeric
  metaData$Season <- ifelse(metaData$Date %in% c("2020-10-27", "2020-10-28"),
                            "Fall", 
                            ifelse(metaData$Date %in% c("2021-08-24", "2021-08-25"),
                                   "Summer",
                                   "Spring"))
  
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
    
    rda.trans <- rda(otus.trans[[k]] ~ Longitude + Latitude + Depth_m +
                       Temperature + Season + Salinity, metaData,
                     na.action = na.omit)
    
    ## Result pRDA
    
    res.trans.temp <- anova.cca(rda.trans, by = "margin") %>% data.frame
    
    res.trans.temp$R2 <- c(RsquareAdj(rda(otus.trans[[k]] ~ Longitude + Condition(Latitude + Depth_m + Season + Temperature + Salinity), metaData, na.action = na.omit))$r.squared,
                           RsquareAdj(rda(otus.trans[[k]] ~ Latitude + Condition(Longitude + Depth_m + Season + Temperature + Salinity), metaData, na.action = na.omit))$r.squared,
                           RsquareAdj(rda(otus.trans[[k]] ~ Depth_m + Condition(Longitude + Latitude + Season + Temperature + Salinity), metaData, na.action = na.omit))$r.squared,
                           RsquareAdj(rda(otus.trans[[k]] ~ Season + Condition(Longitude + Latitude + Depth_m + Temperature + Salinity), metaData, na.action = na.omit))$r.squared,
                           RsquareAdj(rda(otus.trans[[k]] ~ Temperature + Condition(Longitude + Latitude + Depth_m + Season + Salinity), metaData, na.action = na.omit))$r.squared,
                           RsquareAdj(rda(otus.trans[[k]] ~ Salinity + Condition(Longitude + Latitude + Depth_m + Season + Temperature), metaData, na.action = na.omit))$r.squared,
                           "NA")
    res.trans.temp$adj.R2 <- c(RsquareAdj(rda(otus.trans[[k]] ~ Longitude + Condition(Latitude + Depth_m + Season + Temperature + Salinity), metaData, na.action = na.omit))$adj.r.squared,
                               RsquareAdj(rda(otus.trans[[k]] ~ Latitude + Condition(Longitude + Depth_m + Season + Temperature + Salinity), metaData, na.action = na.omit))$adj.r.squared,
                               RsquareAdj(rda(otus.trans[[k]] ~ Depth_m + Condition(Longitude + Latitude + Season + Temperature + Salinity), metaData, na.action = na.omit))$adj.r.squared,
                               RsquareAdj(rda(otus.trans[[k]] ~ Season + Condition(Longitude + Latitude + Depth_m + Temperature + Salinity), metaData, na.action = na.omit))$adj.r.squared,
                               RsquareAdj(rda(otus.trans[[k]] ~ Temperature + Condition(Longitude + Latitude + Depth_m + Season + Salinity), metaData, na.action = na.omit))$adj.r.squared,
                               RsquareAdj(rda(otus.trans[[k]] ~ Salinity + Condition(Longitude + Latitude + Depth_m + Season + Temperature), metaData, na.action = na.omit))$adj.r.squared,
                               "NA")
    
    res.trans.temp$model <- "pRDA"
    res.trans.temp$transfo <- ifelse(k == 1, "Hellinger", "Hellinger-Mahalanobis")
    
    
    ## Complete model
    
    res.trans.all <- anova.cca(rda.trans) %>% data.frame
    
    res.trans.all$R2 <- c(RsquareAdj(rda.trans)$r.squared, 
                          "NA")
    res.trans.all$adj.R2 <- c(RsquareAdj(rda.trans)$adj.r.squared, 
                              "NA")
    res.trans.all$model <- "Complete"
    res.trans.all$transfo <- ifelse(k == 1, "Hellinger", "Hellinger-Mahalanobis")
    
    ## Output RDA results
    
    res.rda <- rbind(res.trans.all, res.trans.temp)
    
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
    
    rda.trans.plot <- data.frame(rda.trans$CCA$u,
                                 ID_ADNe = rownames(scores(rda.trans$CCA$u, display="site"))) %>%
      left_join(., metaData) %>%
      ggplot(aes(x = RDA1, y = RDA2, col = Season)) +
      geom_hline(yintercept = 0, linetype = "dashed", col = "grey") +
      geom_vline(xintercept = 0, linetype = "dashed", col = "grey") +
      geom_point(size = 4) +
      xlab(paste0("RDA1: ", percentVar.rda[1], "%")) +
      ylab(paste0("RDA2: ", percentVar.rda[2], "%")) +
      theme_bw() +
      scale_color_brewer(palette = "Set1") +
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
           height = 3, width = 6, scale = 1.75)
    
  }
  
  ## Spatial autocorrelation 
  
  
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

        rda.trans <- rda(otus.trans.season ~ Longitude + Latitude + Depth_m +
                           Temperature + Salinity, meta.data,
                         na.action = na.omit)
        
        ## Result pRDA
        
        res.trans.temp <- anova.cca(rda.trans, by = "margin") %>% data.frame
        
        res.trans.temp$R2 <- c(RsquareAdj(rda(otus.trans[[k]] ~ Longitude + Condition(Latitude + Depth_m + Temperature + Salinity), metaData, na.action = na.omit))$r.squared,
                               RsquareAdj(rda(otus.trans[[k]] ~ Latitude + Condition(Longitude + Depth_m + Temperature + Salinity), metaData, na.action = na.omit))$r.squared,
                               RsquareAdj(rda(otus.trans[[k]] ~ Depth_m + Condition(Longitude + Latitude + Temperature + Salinity), metaData, na.action = na.omit))$r.squared,
                               RsquareAdj(rda(otus.trans[[k]] ~ Temperature + Condition(Longitude + Latitude + Depth_m + Salinity), metaData, na.action = na.omit))$r.squared,
                               RsquareAdj(rda(otus.trans[[k]] ~ Salinity + Condition(Longitude + Latitude + Depth_m + Temperature), metaData, na.action = na.omit))$r.squared,
                               "NA")
        res.trans.temp$adj.R2 <- c(RsquareAdj(rda(otus.trans[[k]] ~ Longitude + Condition(Latitude + Depth_m + Temperature + Salinity), metaData, na.action = na.omit))$adj.r.squared,
                                   RsquareAdj(rda(otus.trans[[k]] ~ Latitude + Condition(Longitude + Depth_m + Temperature + Salinity), metaData, na.action = na.omit))$adj.r.squared,
                                   RsquareAdj(rda(otus.trans[[k]] ~ Depth_m + Condition(Longitude + Latitude + Temperature + Salinity), metaData, na.action = na.omit))$adj.r.squared,
                                   RsquareAdj(rda(otus.trans[[k]] ~ Temperature + Condition(Longitude + Latitude + Depth_m + Salinity), metaData, na.action = na.omit))$adj.r.squared,
                                   RsquareAdj(rda(otus.trans[[k]] ~ Salinity + Condition(Longitude + Latitude + Depth_m + Temperature), metaData, na.action = na.omit))$adj.r.squared,
                                   "NA")
        
        res.trans.temp$model <- "pRDA"
        res.trans.temp$transfo <- ifelse(k == 1, "Hellinger", "Hellinger-Mahalanobis")
        
        
        ## Complete model
        
        res.trans.all <- anova.cca(rda.trans) %>% data.frame
        
        res.trans.all$R2 <- c(RsquareAdj(rda.trans)$r.squared, 
                              "NA")
        res.trans.all$adj.R2 <- c(RsquareAdj(rda.trans)$adj.r.squared, 
                                  "NA")
        res.trans.all$model <- "Complete"
        res.trans.all$transfo <- ifelse(k == 1, "Hellinger", "Hellinger-Mahalanobis")
        
        ## Output RDA results
        
        res.rda <- rbind(res.trans.all, res.trans.temp)
        
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
        
        rda.trans.plot <- data.frame(rda.trans$CCA$u,
                                     ID_ADNe = rownames(scores(rda.trans$CCA$u, display="site"))) %>%
          left_join(., meta.data) %>%
          ggplot(aes(x = RDA1, y = RDA2, col = Station)) +
          geom_hline(yintercept = 0, linetype = "dashed", col = "grey") +
          geom_vline(xintercept = 0, linetype = "dashed", col = "grey") +
          geom_point(size = 4) +
          xlab(paste0("RDA1: ", percentVar.rda[1], "%")) +
          ylab(paste0("RDA2: ", percentVar.rda[2], "%")) +
          theme_bw() +
          scale_color_manual(values = c("black", RColorBrewer::brewer.pal(7, "Set1")[c(1:5,7)])) +
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
               height = 3, width = 6, scale = 1.5)

      
    }
    
    write.csv(df.correl.all,
              file.path(here::here(), "02_Results", new_dir, new_dir2,
                        paste("Result_MantelCorrel_", 
                              ifelse(k==1, "1trans", "2trans"), 
                              "_", project.name, "_", locus.name, 
                              ".csv", sep = "")))
    
    correl.plot <- df.correl.all %>%
      ggplot(., aes(x = class.index, y = Mantel.cor, col = Season)) +
      geom_hline(yintercept = 0, color = "red") +
      geom_line() + 
      scale_x_continuous(labels = scales::label_comma()) +
      geom_point(aes(shape = Signif, size = 1)) +
      scale_shape_manual(values = c(1,19)) + 
      scale_color_brewer(palette = "Set1") +
      facet_grid(Season ~ distance, scales = "free_x") +
      xlab("Distance (m)") +
      ylab("Mantel correlation") +
      theme_bw() +
      theme(legend.position = "none",
            strip.background = element_blank(),
            strip.text = element_text(hjust = 0)) +
      ggtitle(paste(project.name, locus.name),
              subtitle = paste("Transformation:", ifelse(k==1, "Hellinger", "Hellinger-Mahalanobis"), sep = " "))
    
    ggsave(file.path(here::here(), "02_Results", new_dir, new_dir2,
                     paste("Mantel_correl_", ifelse(k==1, "1trans", "2trans"), "_", project.name, "_", locus.name, ".png", sep = "")),
           correl.plot, 
           height = 3, width = 6, scale = 1.5)
    
  }
  
  
}


# 5) Comparison of transformed vs non-transformed dataset ----------------------

rm(list=setdiff(ls(), c("new_dir", "new_dir_align")))
gc()

new_dir3 <- "Comparison_Transf_vs_NonTransf"
dir.create(file.path(here::here(), "02_Results", new_dir, new_dir3),
           showWarnings = FALSE)

list.projects <- c("BDA", 
                   "PPO.Leim.VRoy.2019", 
                   "Godbout")

list.locus <- c("COI", "MiFishU", "16Schord")



for(i in 1:length(list.locus)){
  
  locus.name <- list.locus[i]
  
  for(j in 1:length(list.projects)){
    
    project.name <- list.projects[j]
    
    file.align <- file.path(here::here(), "02_Results", new_dir, new_dir_align,
                            paste("Aligned_", project.name, "_", locus.name,
                                  ".fasta", sep = ""))
    
    ## Data 
    
    OTUs.count <- read.csv(file.path(here::here(), "00_Data", "07_OTUs_3pc_tables",
                                     paste(project.name, locus.name, 
                                           "OTUs_count.csv", sep = "_")),
                           header = T, row.names = 1)
    
    
    ## Pairwise genetic distance 
    
    dna_dist <- read.alignment(file.align,
                               format = "fasta") %>%
      dist.alignment(., gap = T)
    
    object <- dna_dist %>% as.matrix %>% .[rowSums(is.na(.)) == 0, colSums(is.na(.)) == 0] %>%
      as.dist %>% dna_PC
    
    cv <- dna_PCcov(object)
    
    ## Data transformation 
    
    otus.trans.hell<- OTUs.count[,colnames(OTUs.count) %in% names(cv$cm)] %>%
      decostand(method = "hellinger") %>% as.matrix
    
    otus.trans.hellMahala <- OTUs.count[,colnames(OTUs.count) %in% names(cv$cm)] %>%
      decostand(method = "hellinger") %>% as.matrix %*% cv$tr_mat
    
    ## Dist Euclidean
    
    eucli.hell <- vegdist(otus.trans.hell, method = "euclidean")
    eucli.hellMahala <- vegdist(otus.trans.hellMahala, method = "euclidean")
    
    df <- data.frame(Hellinger = eucli.hell[lower.tri(eucli.hell)],
                     Hellinger_Mahalanobis = eucli.hellMahala[lower.tri(eucli.hellMahala)])
    
    ## Mantel correlation
    
    res.mantel <- mantel(eucli.hell, eucli.hellMahala)
    
    plot.dist <- ggplot(df, aes(x = Hellinger, y = Hellinger_Mahalanobis)) +
      geom_point() +
      # stat_smooth(method = "lm") +
      theme_bw() +
      ggtitle(paste(project.name, locus.name),
              subtitle = paste("r_mantel = ", round(res.mantel$statistic, 2),
                               " (P = ", round(res.mantel$signif, 3), ")",
                               sep = ""))
    
    ggsave(file.path(here::here(), "02_Results", new_dir, new_dir3,
                     paste("Mantel_transfo_", 
                           project.name, "_", locus.name,
                           ".png", sep = "")),
           plot.dist, width = 3, height = 3, scale = 1.5)
    
  }
  
}


# 7) Table formatting of RDA results -------------------------------------------


rm(list=setdiff(ls(), c("new_dir", "new_dir_align")))
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
                                                    "Godbout_Fall",
                                                    "Godbout_Spring", 
                                                    "Godbout_Summer"))

write.csv(df.rda, 
          file.path(here::here(), "02_Results", new_dir, 
                    "Result_RDA_all.csv"))


levels(df.rda$project) <- c("Large", 
                            "Intermediate", 
                            "Fine\n(Fall)", "Fine\n(Spring)", "Fine\n(Summer)")


plot.res.rda <- subset(df.rda, !X %in% c("Residual", "Residual1")) %>%
  mutate(Signif = ifelse(.$Pr..F. <= 0.05 & adj.R2 >= 0, "yes", "no"),
         locus = factor(locus, levels = c("COI", "MiFishU", "16Schord"))) %>% 
  ggplot(aes(x = X, y = adj.R2, shape = transfo, fill = Signif)) +
  geom_point(size = 3.5, position = position_dodge2(w = 0.75)) +
  facet_grid(locus ~ project, scales = "free_x", space = "free_x") +
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

plot.res.rda <- subset(df.rda, transfo == "Hellinger-Mahalanobis" & !X %in% c("Residual", "Residual1")) %>%
  mutate(Signif = ifelse(.$Pr..F. <= 0.05 & adj.R2 >= 0, "yes", "no"),
         locus = factor(locus, levels = c("COI", "MiFishU", "16Schord"))) %>% 
  ggplot(aes(x = X, y = adj.R2, shape = transfo, col = Signif)) +
  geom_point(size = 7.5, position = position_dodge2(w = 0.75), stroke = 10) +
  facet_grid(locus ~ project, scales = "free_x", space = "free_x") +
  scale_shape_manual(name = "Data tranformation", values = "-") +
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
  xlab("Complete model and marginal effect") 
plot.res.rda 

ggsave(file.path(here::here(), "02_Results", new_dir, 
                 "Result_RDA_OTUs_HM_only.png"),
       plot.res.rda, width = 5, height = 3, scale = 1.75)



## compare OTUs vs ESVs --------------------------------------------------------


rm(list=setdiff(ls(), c("new_dir", "new_dir_align")))
gc()

list.projects <- c("BDA", 
                   "PPO.Leim.VRoy.2019", 
                   "Godbout")

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
    
    res.esv <- read.csv(file.path(here::here(), "02_Results", "05_Spatial_diversity", "Result_Hellinger", 
                                  paste("RDA_results_", 
                                        project.name, "_", locus.name, 
                                        ".csv", sep = ""))) %>% select(!c("Transfo", "Project_Locus"))
    
    res.esv$transfo <- "Hellinger-Mahalanobis"
    res.esv$locus <- locus.name
    res.esv$project <- project.name
    res.esv$method <- "ESVs"
    
    res.temp <- rbind(res.H, res.HM) %>% select(!c("model"))
    res.temp$locus <- locus.name
    res.temp$project <- project.name
    res.temp$method <- "OTUs"
    
    res.temp2 <- rbind(res.esv, res.temp)
    res.temp2$X <- str_replace(res.temp2$X, "Model", "Complete_model")
    res.temp2$method <- paste(res.temp2$transfo, "on", res.temp2$method)
    
    df.rda <- rbind(df.rda, res.temp2)
    
  }
  
}



plot.res.rda <- subset(df.rda, !X %in% c("Residual", "Residual1")) %>%
  mutate(Signif = ifelse(.$Pr..F. >= 0.05, "no", "yes"),
         project = factor(project, levels = c("BDA", "PPO.Leim.VRoy.2019", "Godbout")),
         locus = factor(locus, levels = c("COI", "MiFishU", "16Schord"))) %>% 
  ggplot(aes(x = X, y = R2, shape = method, fill = Signif)) +
  geom_point(size = 3.5, position = position_dodge2(w = 0.75)) +
  facet_grid(locus ~ project, scales = "free_x", space = "free_x") +
  scale_shape_manual(name = "Data tranformation", values = c(21, 24, 23)) +
  scale_fill_manual(name = "P value", values = c("grey", "red"), label = c("> 0.05", "< 0.05")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  guides(fill=guide_legend(override.aes=list(colour=c(no="grey",yes="red")))) +
  ylab("R2") + xlab("Complete model and marginal effect") +
  ggtitle("RDA results")
plot.res.rda 

ggsave(file.path(here::here(), "02_Results", new_dir, 
                 "Result_RDA_OTUs_vs_ESVs.png"),
       plot.res.rda, width = 7, height = 4, scale = 1.5)


# Plot Mantel correlogram -----------------------------------------------------

rm(list=setdiff(ls(), c("new_dir", "new_dir_align")))
gc()


### For Large and intermediate scale
list.projects <- c("PPO.Leim.VRoy.2019", "BDA")
list.scale <- c("Large scale", "Intermediate scale")

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
   
   
   correl.plot <- df.correl %>%
     ggplot(., aes(x = class.index, y = Mantel.cor)) +
     geom_hline(yintercept = 0, color = "red") +
     geom_line() + 
     scale_x_continuous(labels = scales::label_comma()) +
     geom_point(aes(shape = Signif, size = 1)) +
     scale_shape_manual(values = c(1,19)) + 
     facet_grid(locus ~ distance, scales = "free_x") +
     xlab("Distance (m)") +
     ylab("Mantel correlation") +
     theme_bw() +
     theme(legend.position = "none",
           strip.background = element_blank(),
           strip.text = element_text(hjust = 0)) +
     ggtitle(scale.name)
   
   correl.plot
   
   ggsave(file.path(here::here(), "02_Results", new_dir, 
                    paste("Mantel_correl_", scale.name, ".png", sep = "")),
          correl.plot, 
          height = 4, width = 6.5, scale = 1.2)
   
}

### for fine scale

df.coi <- read.csv(file.path(here::here(), "02_Results", new_dir,
                             "Result_Hellinger-Mahalanobis",
                             paste("Result_MantelCorrel_2trans_",
                                   "Godbout", 
                                   "_COI.csv",sep = "")), 
                   row.names = 1)
df.coi$locus <- "COI"

df.mifish <- read.csv(file.path(here::here(), "02_Results", new_dir,
                                "Result_Hellinger-Mahalanobis",
                                paste("Result_MantelCorrel_2trans_",
                                      "Godbout", 
                                      "_MiFishU.csv",sep = "")), 
                      row.names = 1)
df.mifish$locus <- "MiFishU"

df.16s <- read.csv(file.path(here::here(), "02_Results", new_dir,
                             "Result_Hellinger-Mahalanobis",
                             paste("Result_MantelCorrel_2trans_",
                                   "Godbout", 
                                   "_16Schord.csv",sep = "")), 
                   row.names = 1)
df.16s$locus <- "16Schord"

df.correl <- rbind(df.coi, df.mifish, df.16s)

df.correl$locus <- factor(df.correl$locus,
                          levels = c("COI", "MiFishU", "16Schord"))
df.correl$distance <- factor(df.correl$distance, 
                             levels = c("Depth|2D", "2D|Depth", "3D"))


correl.plot <- df.correl %>%
  ggplot(., aes(x = class.index, y = Mantel.cor, col = Season)) +
  geom_hline(yintercept = 0, color = "red") +
  geom_line() + 
  scale_x_continuous(labels = scales::label_comma()) +
  geom_point(aes(shape = Signif, size = 1)) +
  scale_shape_manual(values = c(1,19)) + 
  scale_color_brewer(palette = "Set1") +
  facet_grid(locus*Season ~ distance, scales = "free_x") +
  xlab("Distance (m)") +
  ylab("Mantel correlation") +
  theme_bw() +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0)) +
  ggtitle("Fine scale")

correl.plot

ggsave(file.path(here::here(), "02_Results", new_dir, 
                 paste("Mantel_correl_Fine.png", sep = "")),
       correl.plot, 
       height = 7, width = 6.5, scale = 1.2)






