# Info -------------------------------------------------------------------------

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
# 2023-03
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

new_dir <- "00_ESVs_tables_clean"
dir.create(file.path(here::here(), "00_Data", new_dir),
           showWarnings = FALSE)

# List Data --------------------------------------------------------------------

list.projects <- c("BDA", "PPO.Leim.VRoy.2019", "PPO.KMcGregor.2020", "PPO.Godbout.2021")

list.locus <- c("COI", "MiFishU", "16Schord")

# ## additional subsamples for G. Guénard
# 
# list.projects <- c( "PPO.PO.DLevesque.2019",
#                     "PPO.Leim.MJRoux.2019",
#                     "PPO.Kildir.DLevesque.2020",
#                     "PPO.Leim.MJRoux.2020")
# list.locus <- c("MiFishU")

## Sample and MOTUs infos ------------------------------------------------------

sample.info <- mat.or.vec(0,0)
MOTUs.info <- mat.or.vec(0,0)

for(i in 1:length(list.locus)){
  
  locus <- list.locus[i]
  
  for(j in 1:length(list.projects)){
    
    project <- list.projects[j]
    
    files.samples <- unlist(paste("Samples.Metabarinfo.corrected_", 
                                  locus, "_",
                                  project,
                                  ".csv",
                                  sep = ""))
  
    temp <- read.csv(file.path(here::here(), "00_Data", "00_Samples_infos",
                               files.samples))
    temp$locus <- locus
    
    sample.info <- rbind(sample.info, temp[, c("sample_id", "project", "locus", "seqdepth_ok", "artefact_type",
                                               "nb_reads_postmetabaR", "nb_motus_postmetabaR")])
    
    rm(temp)
    
    
    file.MOTUs <- unlist(paste("MOTUs.Metabarinfo.corrected_", 
                               locus, "_",
                               project,
                               ".csv",
                               sep = ""))
    
    temp <- read.csv(file.path(here::here(), "00_Data", "00_MOTUS_infos",
                               file.MOTUs)) %>% data.frame
    temp$locus <- locus
    temp$project <- project
    temp$new_ESV <- paste(project, temp$QueryAccVer, sep = "_")
    temp$Seq_length <- str_length(temp$sequence)
    
    MOTUs.info <- rbind(MOTUs.info, temp)
    
    rm(temp)
  }}
  
    
## Find remaining adaptors (Nextera transposase):
df.temp <- MOTUs.info[,c("sequence", "new_ESV")]
df.temp$adapt_f <- grepl("CTGTCTCTTATACACATCT", df.temp$seq, ignore.case = TRUE)
df.temp$adapt_r <- grepl("AGATGTGTATAAGAGACAG", df.temp$seq, ignore.case = TRUE)
df.temp$trans1_f <- grepl("GTCTCGTGGGCTCGG", df.temp$seq, ignore.case = TRUE)
df.temp$trans1_r <- grepl("CCGAGCCCACGAGAC", df.temp$seq, ignore.case = TRUE)
df.temp$trans2_f <- grepl("TCGTCGGCAGCGTC", df.temp$seq, ignore.case = TRUE)
df.temp$trans2_r <- grepl("GACGCTGCCGACGA", df.temp$seq, ignore.case = TRUE)
df.temp$fld_f <- grepl("ACACTGACGACATGGTTCTACA", df.temp$seq, ignore.case = TRUE)
df.temp$fld_r <- grepl("TACGGTAGCAGAGACTTGGTCT", df.temp$seq, ignore.case = TRUE)

df.temp$Adapt <- apply(df.temp[,c("adapt_f", "adapt_r",
                                  "trans1_f", "trans1_r",
                                  "trans2_f", "trans2_r",
                                  "fld_f", "fld_r")],
                       1, sum)
MOTUs.info$Adapt <- ifelse(df.temp$Adapt == 0, "FALSE", "TRUE") %>% 
  factor(levels = c("TRUE", "FALSE"))

MOTUs.info <- MOTUs.info %>%
  mutate(Selected = ifelse(Adapt == "TRUE", 
                           "No_Adapters", 
                           ifelse(locus != "MiFishU",
                                  "Yes_Good",
                                  ifelse(Seq_length > 180,
                                         "No_long_seq", 
                                         "Yes_Good")))
                )

# remove possible contamination ESVs
MOTUs.info <- subset(MOTUs.info, !Taxon %in% c("Alburnus alburnus", "Bos", 
                                              "Candidatus Pelagibacter", 
                                              "Canis lupus", 
                                              "Catostomus commersonii", 
                                              "Chrosomus neogaeus", 
                                              "Columba livia", 
                                              "Formosa", "Gallus", 
                                              "Hominidae", "Homo sapiens", 
                                              "Lepus americanus", "Mus", 
                                              "Odocoileus virginianus", 
                                              "Phasianidae", 
                                              "Planktomarina temperata", 
                                              "Pseudoalteromonas", 
                                              "Pseudoalteromonas prydzensis", 
                                              "Rhynchobrunnera orthospora", 
                                              "Sus", "Sus scrofa", 
                                              "Synechococcus", 
                                              "Xiphinema brevicolle"))
           
## Metadata --------------------------------------------------------------------

metaData.all <- read.csv(file.path(here::here(), 
                                  "00_Data", "00_FileInfos",
                                  "metaData_2019_to_2021.csv"))

## ESVs table ------------------------------------------------------------------

## Subset only ESVs with no adaptors and good length

for(i in 1:length(list.locus)){
  
  locus.name <- list.locus[i]
  
  for(j in 1:length(list.projects)){
    
    project.name = list.projects[j]
    
    sample.ok <- subset(sample.info,
                        project == project.name &
                          locus == locus.name &
                          artefact_type %in% c("Not artefactual", "Low sequencing depth"))
    
    motus.ok <- subset(MOTUs.info,
                       locus == locus.name & 
                         project == project.name & 
                         Selected == "Yes_Good")

    esv.temp <- read.csv(file.path(here::here(), 
                                   "00_Data",
                                   "00_ESVs_Corrected",
                                   paste("ESVtab.corrected_", 
                                         locus.name, 
                                         "_",
                                         project.name,
                                         ".csv",
                                         sep = "")),
                                   header = T, row.names = 1)
      
      list.esv.ok <- paste(motus.ok$QueryAccVer)
      list.sample.ok <- paste(sample.ok$sample_id)
      
      count.esv.ok <- esv.temp[list.sample.ok, list.esv.ok] %>% .[rowSums(.) > 0, colSums(.) > 0]
      
      meta.ok <-  metaData.all %>% subset(., ID_ADNe %in% rownames(count.esv.ok) & Locus == locus.name)
      
      count.esv.ok <- count.esv.ok[match(meta.ok$ID_ADNe, rownames(count.esv.ok)),]
      
      sample.ok <- subset(sample.ok, sample_id %in% rownames(count.esv.ok))
      
      motus.ok <- subset(motus.ok, QueryAccVer %in% colnames(count.esv.ok))
      
      # export tables
      
      write.csv(count.esv.ok,
                file.path(here::here(), "00_Data", new_dir,
                          paste(project.name, locus.name, "ESVs_ok.csv", sep = "_")))
      
      write.csv(motus.ok,
                file.path(here::here(), "00_Data", new_dir,
                          paste(project.name, locus.name, "MOTUs_ok.csv", sep = "_")))
      
      write.csv(meta.ok,
                file.path(here::here(), "00_Data", new_dir,
                          paste(project.name, locus.name, "metadata_ok.csv", sep = "_")))
      
      write.csv(sample.ok,
                file.path(here::here(), "00_Data", new_dir,
                          paste(project.name, locus.name, "samples_ok.csv", sep = "_")))
      
  }
}


# Merging Godbout 2020 and 2021


for(i in 1:length(list.locus)){
  
  locus.name <- list.locus[i]
  
  esv.2020 <- read.csv(file.path(here::here(), "00_Data", "00_ESVs_tables_clean",
                                 paste("PPO.KMcGregor.2020", locus.name, 
                                       "ESVs_ok.csv", sep = "_")),
                       header = T, row.names = 1)
  
  esv.2021 <- read.csv(file.path(here::here(), "00_Data", "00_ESVs_tables_clean",
                                 paste("PPO.Godbout.2021", locus.name, 
                                       "ESVs_ok.csv", sep = "_")),
                       header = T, row.names = 1)
  
  
  MOTUs.2020 <- read.csv(file.path(here::here(), "00_Data", "00_ESVs_tables_clean",
                                   paste("PPO.KMcGregor.2020", locus.name, 
                                         "MOTUs_ok.csv", sep = "_")),
                         header = T, row.names = 1)

  MOTUs.2021 <- read.csv(file.path(here::here(), "00_Data", "00_ESVs_tables_clean",
                                   paste("PPO.Godbout.2021", locus.name, 
                                         "MOTUs_ok.csv", sep = "_")),
                         header = T, row.names = 1)
 
  meta.2020 <- read.csv(file.path(here::here(), "00_Data", "00_ESVs_tables_clean",
                                  paste("PPO.KMcGregor.2020", locus.name, 
                                        "metadata_ok.csv", sep = "_")),
                        header = T, row.names = 1)
  
  meta.2021 <- read.csv(file.path(here::here(), "00_Data", "00_ESVs_tables_clean",
                                  paste("PPO.Godbout.2021", locus.name, 
                                        "metadata_ok.csv", sep = "_")),
                        header = T, row.names = 1)
  
  
  sample.2020 <- read.csv(file.path(here::here(), "00_Data", "00_ESVs_tables_clean",
                                  paste("PPO.KMcGregor.2020", locus.name, 
                                        "samples_ok.csv", sep = "_")),
                        header = T, row.names = 1)
  
  sample.2021 <- read.csv(file.path(here::here(), "00_Data", "00_ESVs_tables_clean",
                                  paste("PPO.Godbout.2021", locus.name, 
                                        "samples_ok.csv", sep = "_")),
                        header = T, row.names = 1)
  
  
  # change ESV names by sequences
  
  names(esv.2020)[match(MOTUs.2020$QueryAccVer, names(esv.2020))] = MOTUs.2020$sequence
  names(esv.2021)[match(MOTUs.2021$QueryAccVer, names(esv.2021))] = MOTUs.2021$sequence
  
  esv.2020$Project <- "Godbout.2020"
  esv.2021$Project <- "Godbout.2021"
  
  esv.2020$Samples <- rownames(esv.2020)
  esv.2021$Samples <- rownames(esv.2021)
  
  esv.all <- full_join(esv.2020, esv.2021)
  rownames(esv.all) <- esv.all$Samples
  
  MOTUs.all <- full_join(MOTUs.2020, MOTUs.2021, by = c("sequence", "Taxon", "Levels", "species", "genus", "family",
                                                        "order", "class", "phylum", "kingdom", "locus", "Seq_length", "Adapt", "Selected"))
  
  MOTUs.all$QueryAccVer <- paste("Godbout_newESV", 
                                 locus.name, 
                                 seq(1, nrow(MOTUs.all)), sep = "_")
  
  esv.all <- esv.all[, !colnames(esv.all) %in% c("Project", "Samples")]
  esv.all[is.na(esv.all)] <- 0
  names(esv.all)[match(names(esv.all), MOTUs.all$sequence)] = MOTUs.all$QueryAccVer
  
  meta.all <- full_join(meta.2020, meta.2021)
  
  sample.all <- full_join(sample.2020, sample.2021)
  
  # export tables
  
  write.csv(sample.all,
            file.path(here::here(), "00_Data", new_dir,
                      paste("Godbout", locus.name, "samples_ok.csv", sep = "_")))
  
  write.csv(meta.all,
            file.path(here::here(), "00_Data", new_dir,
                      paste("Godbout", locus.name, "metadata_ok.csv", sep = "_")))
  
  write.csv(MOTUs.all,
            file.path(here::here(), "00_Data", new_dir,
                      paste("Godbout", locus.name, "MOTUs_ok.csv", sep = "_")))
  
  write.csv(esv.all,
            file.path(here::here(), "00_Data", new_dir,
                      paste("Godbout", locus.name, "ESVs_ok.csv", sep = "_")))
}





