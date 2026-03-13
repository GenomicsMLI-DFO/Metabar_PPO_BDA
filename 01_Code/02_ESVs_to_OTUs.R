# Info --------------------------------------------------------------------

#
# Cluster ESVs to OTUs
# Using kmer R package
#
# CL
# 2023-06-14
#

# Library -----------------------------------------------------------------
rm(list = ls())
gc()

library(ape)
library(kmer)
library(tidyverse)
library(Biostrings)

# Creat new_dir for results  ---------------------------------------------------

# new_dir <- "07_OTUs_tables" ## ===> with threshold of 5%
new_dir <- "07_OTUs_3pc_tables" ## ===> with threshold of 3%
dir.create(file.path(here::here(), "00_Data", new_dir),
           showWarnings = FALSE)


# Clean tables are already available, outputs of 00_Fromat_ESVs_tales script

list.projects <- c("BDA", 
                   "PPO.Leim.VRoy.2019", 
                   "Godbout")

list.locus <- c("COI", "MiFishU", "16Schord")


# Cluster ESVs by locus and by project -----------------------------------------

for(i in 1:length(list.locus)){
  
  locus.name <- list.locus[i]
  
  for(j in 1:length(list.projects)){
    
    ## 1) Read data ------------------------------------------------------------
    
    project.name <- list.projects[j]
    
    MOTUs.info <- read.csv(file.path(here::here(), "00_Data", "00_ESVs_tables_clean",
                                     paste(project.name, locus.name, 
                                           "MOTUs_ok.csv", sep = "_")),
                           header = T, row.names = 1)
    
    data.ESVs <- read.csv(file.path(here::here(), "00_Data", "00_ESVs_tables_clean",
                                    paste(project.name, locus.name, 
                                          "ESVs_ok.csv", sep = "_")),
                          header = T, row.names = 1)
    # 
    # metaData <- read.csv(file.path(here::here(), "00_Data", "00_ESVs_tables_clean",
    #                                paste(project.name, locus.name, 
    #                                      "metadata_ok.csv", sep = "_")),
    #                      header = T, row.names = 1)
    
    ## 2) Create OTUs ----------------------------------------------------------
    
    ### DNA sequences for all ESVs
    
    seq.test <- DNAStringSet(MOTUs.info$sequence)
    names(seq.test) <- MOTUs.info$QueryAccVer
    
    seq.obj1 <- as.DNAbin(seq.test)
    
    test.OTUs <- otu(seq.obj1, threshold = 0.97, nstart = 40, method = "farthest")
    nlevels(as.factor(test.OTUs))
    
    
    OTUs.all <- data.frame(Name = gsub("[*]", "", names(test.OTUs)),
                           OTUs = test.OTUs,
                           Central = grepl("[*]", names(test.OTUs)))
    OTUs.all.seq <- left_join(OTUs.all, MOTUs.info, by = c("Name" = "QueryAccVer"))
    

    write.csv(OTUs.all.seq, 
              file.path(here::here(), "00_Data", new_dir,
                        paste(project.name, locus.name, 
                              "OTUs_all.csv", sep = "_")))
    
    ### Creat table OTUs infos 
    
    list.otus <- unique(OTUs.all.seq$OTUs)
    
    OTUs.info <- mat.or.vec(0,0)
    
    for(k in 1:length(list.otus)){
      
      temp1 <- subset(OTUs.all.seq, OTUs == list.otus[k])
      
      temp2 <- subset(temp1, Central == "TRUE")
      
      temp2 <- if(ncol(temp2 == 0)){
        temp1[1,]
      }
      
      temp3 <- temp1[!is.na(temp1$Taxon),]
      taxon.temp <- if(nrow(temp3) == 0){
        temp2[,c("Taxon", "Levels", "species", "genus", "family", "order", "class", "phylum", "kingdom")]
      } else {
        temp3[1,c("Taxon", "Levels", "species", "genus", "family", "order", "class", "phylum", "kingdom")]
      }
                          
      sequence <- temp2[,c(2:4)]
      
      taxon <- if(is.na(temp2$Taxon)){
        taxon.temp
      } else {
        temp2[1,c("Taxon", "Levels", "species", "genus", "family", "order", "class", "phylum", "kingdom")]
      }
                    
      otus <- cbind(sequence, taxon)
      
      OTUs.info <- rbind(OTUs.info, otus)
      
    }
    
    OTUs.info$OTUs <- paste("OTU", OTUs.info$OTUs, sep = "_")
    
    write.csv(OTUs.info, 
              file.path(here::here(), "00_Data", new_dir,
                        paste(project.name, locus.name, 
                              "OTUs_infos.csv", sep = "_")))
    

    ## 3) New count table based on OTUs ----------------------------------------

    
    new.otus <- mat.or.vec(nrow(data.ESVs),0)
    
    list.otus <- unique(OTUs.all.seq$OTUs)
    
    for(k in 1:length(list.otus)) {
      OTU <- list.otus[k]
      list.esv.temp <- subset(OTUs.all.seq, OTUs == OTU)$Name
      esv.temp <- subset(data.ESVs, select = list.esv.temp) %>% apply(., 1, sum) 
      esv.temp <- data.frame(OTU = esv.temp)
      colnames(esv.temp) <- paste("OTU", OTU, sep = "_")
      new.otus <- cbind(new.otus, esv.temp)
      
      rm(esv.temp)
    }
    
    write.csv(new.otus, 
              file.path(here::here(), "00_Data", new_dir,
                        paste(project.name, locus.name, "OTUs_count.csv", sep = "_")))
    
  }
  
}

