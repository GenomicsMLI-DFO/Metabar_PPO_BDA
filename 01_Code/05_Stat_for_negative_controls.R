# Info -------------------------------------------------------------------------
#
# Statistics for the different negative control
#
# Subset ESVs tables for 
# Locus:
#     - MiFishU
#     - 16Schord
#     - COI
# 
# for - BDA 2020
#     - VRoy 2019
#     - Godbout 2020 and 2021
#
# CL
# 2023-11
#

# Library ----------------------------------------------------------------------
rm(list = ls())
gc()

library(readxl)
library(tidyverse)
library(ggridges)
library(ggpubr)
library(ggVennDiagram)

# Create new_dir for results  ---------------------------------------------------

new_dir <- "10_Negative_control"
dir.create(file.path(here::here(), "02_Results", new_dir),
           showWarnings = FALSE)

new_dir2 <- "Tables"
dir.create(file.path(here::here(), "02_Results", new_dir, new_dir2),
           showWarnings = FALSE)

# List Data --------------------------------------------------------------------

list.projects <- c("BDA", "PPO.Leim.VRoy.2019", "PPO.KMcGregor.2020", "PPO.Godbout.2021")

list.locus <- c("COI", "MiFishU", "16Schord")

ESV.table.control.long <- data.frame(stringsAsFactors = F)

for(i in 1:length(list.projects)) {
  
  project.name <- list.projects[i]
  
  for(j in 1:length(list.locus)) {
    
    locus.name <- list.locus[j]
    
    esv.tab <- read.csv(file.path(here::here(), "00_Data", "00_ESV_postTagJump",
                                  paste0("ESVtab.postTagjump_",
                                         locus.name, "_", project.name,
                                         ".csv")),
                        row.names = 1)
    
    seq.info <- read.csv(file.path(here::here(), "00_Data", "00_ESV_postTagJump",
                                   paste0("Samples.Metabarinfo.postTagjump_",
                                          locus.name, "_", project.name,
                                          ".csv"))) %>%
      subset(., Type_echantillon != "ECH")
    
    motus.infos.all <- read.csv(file.path(here::here(), "00_Data", "00_ESV_postTagJump",
                                          paste0("MOTUs.Metabarinfo.postTagjump_",
                                                 locus.name, "_", project.name,
                                                 ".csv")),
                                row.names = 1)
    
    
    esv.control <- esv.tab[seq.info$sample_id,]
    meta.data <- subset(seq.info, sample_id %in% rownames(esv.control))
    motus.info <- subset(motus.infos.all, QueryAccVer %in% colnames(esv.control))
    
    
    write.csv(esv.control,
              file.path(here::here(), "02_Results", new_dir, new_dir2,
                        paste0("ESV_postTagjump_",
                               locus.name, "_", project.name,
                               ".csv")))
    write.csv(meta.data,
              file.path(here::here(), "02_Results", new_dir, new_dir2,
                        paste0("MetaData_postTagjump_",
                               locus.name, "_", project.name,
                               ".csv")))
    write.csv(motus.info,
              file.path(here::here(), "02_Results", new_dir, new_dir2,
                        paste0("MOTUs.Metabarinfo_postTagjump_",
                               locus.name, "_", project.name,
                               ".csv")))
    
    esv.control$sample_id <- rownames(esv.control)
    
    ESV.table.control.long.int <- esv.control %>%
      pivot_longer(-sample_id, names_to = "ID", values_to = "Nreads") %>%
      mutate(Loci = locus.name) %>%
      left_join(seq.info %>% select(sample_id, ID_projet, Loci, Type_echantillon, Station, Site_echantillonnage) %>%  distinct()) %>%
      left_join(motus.info %>% select(-Loci), by = c("ID" = "QueryAccVer")) %>%
      mutate(Project = project.name)
    
    ESV.table.control.long <- bind_rows(ESV.table.control.long, ESV.table.control.long.int)
    
    
  }
  
}


# # Plots negative control -------------------------------------------------------
write.csv(ESV.table.control.long,
          file.path(here::here(), "02_Results", new_dir, new_dir2,
                    paste0("ESV_table_control_long",
                           ".csv")))


# rm(list = ls())
# gc()


# ESV.table.control.long <- read.csv("ESV_table_control_long.csv")
ESV.table.control.long <- ESV.table.control.long %>% 
  mutate(Sampling.scale = ifelse(Project == "BDA", "Intermediate", ifelse(Project == "PPO.Leim.VRoy.2019", "Broad", "Fine")))

ESV.table.control.long$Type_echantillon <- recode(ESV.table.control.long$Type_echantillon,
                                                  "SNC_FNC" = "FNC",
                                                  "MNC" = "PNC",
                                                  "MPC" = "PPC")
ESV.table.control.long$Type_echantillon <- factor(ESV.table.control.long$Type_echantillon,
                                                  levels = c("FNC", "ENC", "PNC", "NTC", "PPC"))

# Légende
# ECH = Échantillon
# SNC = Contrôle négatif de terrain (field negative control)
# FNC = Contrôle négatif de filtration (filtration negative control)
# ENC = Contrôle négatif d’extraction (extraction negative control)
# PNC/NTC = Contrôle negatif de PCR (PCR negatif control)
# PPC = Contrôle positif de PCR (PCR positive control)

list.locus <- c("COI", "MiFishU", "16Schord")
list.plots <- list()

for(i in 1:length(list.locus)){
  
  locus.name <- list.locus[i]
  
  resi.int <- subset(ESV.table.control.long, Loci == locus.name & Type_echantillon != "PPC") %>% 
    mutate(Taxon= ifelse(is.na(Taxon), "Unassigned", Taxon),
           Category = ifelse(not_a_max_conta == "TRUE", "Good MOTUs", "Contaminants")) %>% 
    group_by(Category, Loci, sample_id, Type_echantillon, Taxon, Sampling.scale ) %>% 
    summarise(Nreads = sum(Nreads)) %>% 
    mutate(Taxon_loci_Cat = paste(Taxon, Loci, Category),
           Sampling.scale = factor(Sampling.scale, levels = c("Broad", "Intermediate", "Fine"))) 
  
  focus.MOTUs <- resi.int %>% filter(Nreads > 0) %>% pull(Taxon_loci_Cat) %>% unique()
  
  
  graph.resi <- subset(resi.int, Taxon_loci_Cat %in% focus.MOTUs) %>%
    ggplot(aes(fill = Nreads , x = sample_id, y = Taxon)) +
    labs(x= "", y = "") + 
    geom_bin2d(color = "darkgray", aes(group = Nreads))+
    scale_fill_distiller(trans = "log10",
                         palette = "Spectral",
                         na.value = "white"#,
                         #breaks = c(1, 10, 100, 1000, 10000, 100000, 1000000,10000000), labels = c("1", "10", "100", "1,000", "10,000", "100,000", "1,000,000", "10,000,000")
    ) +
    theme_minimal()+
    ggh4x::facet_nested(~ Sampling.scale + Type_echantillon, scale = "free", space = "free") + 
    ggtitle(locus.name) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          strip.background = element_rect(fill = "white"),
          strip.text.y = element_text(angle = 0),
          strip.text.x = element_text(angle = 0),
          panel.spacing = unit(0, "in"),
          legend.text = element_text(angle = 0, hjust = 1),
          legend.title = element_text(angle = 0, vjust = 0.9),
          legend.position = "bottom")
  
  list.plots[[i]] <- graph.resi
  
}

fig.scale <- ggarrange(list.plots[[1]],
                       list.plots[[2]], 
                       list.plots[[3]], ncol = 1, nrow = 3,
                       common.legend = T, legend = "right",
                       labels = LETTERS, align = "hv",
                       heights = c(.9, 1.8, 1.4)) +
  theme(plot.background = element_rect(fill = "white"))
fig.scale

ggsave(file.path(here::here(), "02_Results", new_dir,
                 "negative_control.png"),
       width = 6, height = 6, scale = 2)


test <- resi.int %>% group_by(Sampling.scale, Type_echantillon) %>% summarise(n = length(unique(sample_id)))
test


# Number of OTUs and ESVs per locus --------------------------------------------

rm(list = ls())

list.locus <- c("COI", "MiFishU", "16Schord")
list.project <- c("BDA", "Godbout", "PPO.Leim.VRoy.2019")

df.ESVs <- mat.or.vec(0,0)
df.OTUs <- mat.or.vec(0,0)

for(i in 1:length(list.project)){
  
  project.name <- list.project[i]
  
  for(j in 1:length(list.locus)){
    
    locus.name <- list.locus[j]
    
    df.temp <- read.csv(file.path(here::here(), "00_Data", "07_OTUs_3pc_tables",
                                  paste(project.name, locus.name, "OTUs_all.csv", 
                                        sep = "_"))) %>%
      select(Name, OTUs, sequence)
    
    df.temp$locus <- locus.name
    df.temp$project <- project.name
    
    df.ESVs <- rbind(df.ESVs, df.temp)
    
    
    df.temp2 <- read.csv(file.path(here::here(), "00_Data", "07_OTUs_3pc_tables",
                                  paste(project.name, locus.name, "OTUs_infos.csv", 
                                        sep = "_"))) %>%
      select(OTUs, sequence)
    
    df.temp2$locus <- locus.name
    df.temp2$project <- project.name
    
    df.OTUs <- rbind(df.OTUs, df.temp2)
    

  }
}


esv <- df.ESVs %>% group_by(locus) %>% summarise(N_ESVs = unique(sequence) %>% length)
otu <- df.OTUs %>% group_by(locus) %>% summarise(N_ESVs = unique(sequence) %>% length)

