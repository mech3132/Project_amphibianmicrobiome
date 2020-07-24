#!/bin/bash Rscript

#### Loading data #####
library(tidyverse)
load("3_5sp_mapping_otu_downstream/mf_alt_filt_final.RData")
load("3_5sp_mapping_otu_downstream/otu_filt.RData")
load("3_5sp_mapping_otu_downstream/taxonomy.RData")
load("4_Bayesian_models/all_p.RData")

load("5_random_forest_LOO/importance_infect_count_wsp_LOO.RData")
load("5_random_forest_LOO/importance_infect_PA_wsp_LOO.RData")
load("5_random_forest_LOO/importance_infect_onlyp_wsp_LOO.RData")

load("5_random_forest_LOO/importance_PABD_count_wsp_LOO.RData")
load("5_random_forest_LOO/importance_PABD_PA_wsp_LOO.RData")
load("5_random_forest_LOO/importance_PABD_onlyp_wsp_LOO.RData")

source("./R_code/ANCOM_updated_code_MYC.R")

#### Which OTUs are important for which species? ####

# Get list of OTU taxonomies; make unique and set as rownames
allOTUs <- make.unique(names = as.character(taxonomy[match(otu_filt$`#OTU ID`, taxonomy$Sequence),"taxonomy"]))
allSeqs <- otu_filt$`#OTU ID`
# rownames(otu_filt) <- allOTUs
rownames(otu_filt) <- allSeqs

# Mapping file adjust
mf_alt_filt_final$prepost <- factor(mf_alt_filt_final$prepost)
mf_alt_filt_final$species <- factor(mf_alt_filt_final$species)

# Collapse taxonomy by level
taxonomy_expanded <- taxonomy %>% separate(taxonomy, into=c("Kingdom","Phylum","Class","Order","Family","Genus","Species"), sep="; ", remove = FALSE) %>%
  mutate(Kingdom = ifelse(is.na(Kingdom), "Unassigned",Kingdom)
         , Kingdom = ifelse(Kingdom == "k__", "Unassigned",Kingdom)
         , Phylum = ifelse(is.na(Phylum), "Unassigned",ifelse(Phylum == "p__", Kingdom, Phylum))
         , Class = ifelse(is.na(Class), "Unassigned",ifelse(Class == "c__", Phylum, Class))
         , Order = ifelse(is.na(Order), "Unassigned",ifelse(Order == "o__", Class, Order))
         , Family = ifelse(is.na(Family), "Unassigned",ifelse(Family == "f__", Order, Family))
         , Genus = ifelse(is.na(Genus), "Unassigned",ifelse(Genus == "g__", Family, Genus))
         , Species = ifelse(is.na(Species), "Unassigned",ifelse(Species == "s__", Genus, Species))
  ) 

# Combine OTU table and MF table
OTU_t <- t(otu_filt[,-which(colnames(otu_filt)=="#OTU ID")]) %>% 
  as.data.frame() %>%
  rownames_to_column(var="SampleID")
mfotu <- OTU_t %>% 
  left_join(mf_alt_filt_final) %>%
  select(one_of(allSeqs), species, Bd_exposure, time, PABD, Bd_load,prepost, SampleID) %>%
  gather(-c(species, SampleID, Bd_exposure, time,PABD, Bd_load,prepost), key=Sequence, value=count) %>%
  left_join(taxonomy_expanded) 

mfotu %>% 
    group_by(SampleID, species, Class, prepost) %>% summarize(reads=sum(count)) %>% ungroup() %>%
    group_by(SampleID, species,prepost) %>% mutate(RelAbund = reads/sum(reads)) %>% 
    ggplot(aes(x=SampleID, y=RelAbund, fill=Class)) + 
    geom_bar(stat = "identity") +
    facet_wrap(prepost~species, nrow=2, scales="free", drop = TRUE) + theme(axis.text.x = element_text(angle=90))
  

# Find out ANCOM results for this; pre-post
OTU_ancom <- OTU_t %>% rename(Sample.ID = SampleID)
MF_ancom <- mf_alt_filt_final %>% rename(Sample.ID= SampleID)
allW <- data.frame()
for ( sp in unique(MF_ancom$species) ) {
  MF_temp <- MF_ancom %>% filter(species == sp)
  OTU_temp <- OTU_ancom %>% filter(Sample.ID %in% MF_temp$Sample.ID)
  ancom_temp <- ANCOM.main.myc(OTUdat = OTU_temp, Vardat = MF_temp,
                               adjusted = F,
                               repeated = F,
                               main.var = "prepost",
                               adj.formula = NULL,
                               repeat.var = NULL,
                               longitudinal = F,
                               random.formula = NULL,
                               multcorr = 2,
                               sig = 0.05,
                               prev.cut = 0.9)
  tempW <- ancom_temp$W.taxa %>% mutate(change = sign(Post_relAbund-Pre_relAbund)) %>%
    select(otu.names,W_stat, W_f, change, detected_wcrit, detected_0.9, effSize_statcrit, effSize_0.9) %>%
    mutate(species = sp)
  
  allW <- rbind(allW, tempW)
  
}

allW_withtax <- allW %>% rename(Sequence = otu.names) %>% left_join(taxonomy_expanded)

allW_withtax %>%ggplot(aes(x=W_f, y=W_stat, col = factor(inhibitory))) + geom_point() +
  facet_wrap(.~ species)

## All of them
ancom_all <- ANCOM.main.myc(OTUdat = OTU_ancom, Vardat = MF_ancom,
                            adjusted = F,
                            repeated = F,
                            main.var = "prepost",
                            adj.formula = NULL,
                            repeat.var = NULL,
                            longitudinal = F,
                            random.formula = NULL,
                            multcorr = 2,
                            sig = 0.05,
                            prev.cut = 0.8)

ancom_all$Plot.volcano


#### Which OTUs are important for predicting infection risk? ####

#### Which OTUs are important for predicting infection load? ####

#### Which OTUs change after infection with infection presence/absence? ####

#### Which OTUs change after infection with infection load? ####


