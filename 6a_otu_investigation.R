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

# make folder
dir.create("6a_otu_investigation")

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
View(OTU_t)
# make relabund
OTU_t_relabund <- OTU_t %>% gather(-SampleID, key=OTU, value=reads) %>% 
  group_by(SampleID) %>%
  mutate(relAbund = reads/sum(reads)) %>%
  select(SampleID, OTU, relAbund)%>%
  spread(key=OTU, value=relAbund)
View(OTU_t_relabund)

mfotu <- OTU_t %>% 
  left_join(mf_alt_filt_final) %>%
  select(one_of(allSeqs), species, Bd_exposure, time, PABD, Bd_load,prepost, indivID, SampleID) %>%
  gather(-c(species, SampleID, Bd_exposure, time,PABD, Bd_load,indivID, prepost), key=Sequence, value=count) %>%
  left_join(taxonomy_expanded) 

mfotu_relabund <- OTU_t_relabund %>% 
  left_join(mf_alt_filt_final) %>%
  select(one_of(allSeqs), species, Bd_exposure, time, PABD, Bd_load,prepost, indivID, SampleID) %>%
  gather(-c(species, SampleID, Bd_exposure, time,PABD, Bd_load,indivID, prepost), key=Sequence, value=relAbund) %>%
  left_join(taxonomy_expanded) 

# mfotu %>% 
#     group_by(SampleID, species, Class, prepost) %>% summarize(reads=sum(count)) %>% ungroup() %>%
#     group_by(SampleID, species,prepost) %>% mutate(RelAbund = reads/sum(reads)) %>% 
#     ggplot(aes(x=SampleID, y=RelAbund, fill=Class)) + 
#     geom_bar(stat = "identity") +
#     facet_wrap(prepost~species, nrow=2, scales="free", drop = TRUE) + theme(axis.text.x = element_text(angle=90))
#   

# Find out ANCOM results for this; pre-post
MF_ancom <- mf_alt_filt_final %>% rename(Sample.ID= SampleID) %>% filter(Bd_exposure == "Bd-exposed") %>% mutate(PABD_fact = ifelse(PABD==1, "Infected","NotInfected"))
OTU_ancom <- OTU_t %>% rename(Sample.ID = SampleID) %>% filter(Sample.ID %in% MF_ancom$Sample.ID)
allW <- data.frame()
if ( !file.exists("6a_otu_investigation/allW_withtax.RData") ) {
  for ( sp in unique(MF_ancom$species) ) {
    sp = unique(MF_ancom$species)[3]
    MF_temp <- MF_ancom %>% filter(species == sp)
    if (length(unique(MF_temp$PABD))<2) {
      next 
    } else {
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
                                   prev.cut = 0.6)
      tempW <- ancom_temp$W.taxa %>% mutate(change = sign(Post_relAbund-Pre_relAbund)) %>%
        select(otu.names,W_stat, W_f, change, detected_wcrit, detected_0.9, effSize_statcrit, effSize_0.9) %>%
        mutate(species = sp)
      # ancom_temp <- ANCOM.main.myc(OTUdat = OTU_temp, Vardat = MF_temp,
      #                              adjusted = F,
      #                              repeated = F,
      #                              main.var = "PABD_fact",
      #                              adj.formula = NULL,
      #                              repeat.var = NULL,
      #                              longitudinal = F,
      #                              random.formula = NULL,
      #                              multcorr = 2,
      #                              sig = 0.05,
      #                              prev.cut = 0.6)
      # 
      # tempW <- ancom_temp$W.taxa %>% mutate(change = sign(Infected_relAbund-NotInfected_relAbund)) %>%
      #   select(otu.names,W_stat, W_f, change, detected_wcrit, detected_0.9, effSize_statcrit, effSize_0.9) %>%
      #   mutate(species = sp)
      
      allW <- rbind(allW, tempW)
    }
  }
  allW_withtax <- allW %>% rename(Sequence = otu.names) %>% left_join(taxonomy_expanded)
  
  save(allW_withtax, file="./6a_otu_investigation/allW_withtax.RData")
} else {
  load("6a_otu_investigation/allW_withtax.RData")
}


allW_withtax %>%ggplot(aes(x=W_f, y=W_stat, col = factor(inhibitory))) + geom_point() +
  facet_wrap(.~ species)



## Also do one with PABD

MF_ancom_postbd <- MF_ancom %>% filter(prepost == "Post")
OTU_ancom_postbd <- OTU_ancom %>% filter(Sample.ID %in% MF_ancom_postbd$Sample.ID)
ancom_postPABD <- ANCOM.main.myc(OTUdat = OTU_ancom_postbd, Vardat = MF_ancom_postbd,
                             adjusted = F,
                             repeated = F,
                             main.var = "PABD_fact",
                             adj.formula = NULL,
                             repeat.var = NULL,
                             longitudinal = F,
                             random.formula = NULL,
                             multcorr = 2,
                             sig = 0.05,
                             prev.cut = 0.6)
ancom_postPABD$W.taxa %>% filter(W_stat>20) %>% rename(Sequence = otu.names) %>%left_join(taxonomy_expanded) %>%
  mutate(change = sign(Infected_relAbund-NotInfected_relAbund))



#### Which OTUs are important for predicting infection risk? ####

# Which individuals become infected?
infectedIndiv <- mf_alt_filt_final %>% filter(PABD == 1) %>% select(indivID) %>% pull() %>% unique()
# What is max bd load of each infected individual?
infectedIndiv_maxBdLoad <- mf_alt_filt_final %>% filter(indivID %in% infectedIndiv) %>% group_by(indivID) %>% summarize(maxBdLoad = max(Bd_load))

# Summarize by individual
mfotu_byIndiv_pre <- mfotu_relabund %>% filter(prepost == "Pre", Bd_exposure == "Bd-exposed") %>% group_by(species, indivID, inhibitory, Sequence, taxonomy, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>% 
  summarize(aveRelAbund = mean(relAbund)) %>%
  ungroup()

# Separate mf and otu table
OTU_ancom_byindiv_pre <- mfotu_byIndiv_pre %>% 
  select(indivID, Sequence, aveRelAbund) %>% spread(key=Sequence, value=aveRelAbund) %>%
  rename(Sample.ID = indivID)
MF_ancom_byindiv_pre <- mfotu %>% filter(indivID %in% OTU_ancom_byindiv_pre$Sample.ID) %>% group_by(indivID) %>% summarize(PABD_future = ifelse(sum(PABD)>0, 1,0)) %>%
  rename(Sample.ID = indivID)

# Run ANCOM
ancom_predictPABD <- ANCOM.main.myc(OTUdat = OTU_ancom_byindiv_pre, Vardat = MF_ancom_byindiv_pre,
                                 adjusted = F,
                                 repeated = F,
                                 main.var = "PABD_future",
                                 adj.formula = NULL,
                                 repeat.var = NULL,
                                 longitudinal = F,
                                 random.formula = NULL,
                                 multcorr = 2,
                                 sig = 0.05,
                                 prev.cut = 0.6)





ggplot(aes(x=prepost, y=aveRelAbund, fill=Sequence, lty=factor(inhibitory))) + geom_bar(stat = "identity", show.legend = FALSE) +
  facet_grid(indivID~ Bd_exposure)

# Plot correlation between each OTU and the infection risk, colored by species.



  

#### Which OTUs are important for predicting infection load? ####

# Which individuals become infected?
exposedIndiv <- mf_alt_filt_final %>% filter(Bd_exposure=="Bd-exposed") %>% select(indivID) %>% pull() %>% unique()
# What is max bd load of each infected individual?
exposedIndiv_maxBdLoad <- mf_alt_filt_final %>% filter(indivID %in% exposedIndiv) %>% group_by(indivID) %>% summarize(maxBdLoad = max(Bd_load))

# Summarize by individual
mfotu_byIndiv_pre <- mfotu_relabund %>% filter(prepost == "Pre", Bd_exposure == "Bd-exposed") %>% group_by(species, indivID, inhibitory, Sequence, taxonomy, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>% 
  summarize(aveRelAbund = mean(relAbund)) %>%
  ungroup()

# Separate mf and otu table

MFOTU_byindiv_pre <- mfotu_byIndiv_pre %>% 
  select(indivID, Sequence, aveRelAbund) %>% spread(key=Sequence, value=aveRelAbund) %>%
  full_join(exposedIndiv_maxBdLoad)

# Which OTUs have zero counts
allOTUs <- colnames(MFOTU_byindiv_pre)[!(colnames(MFOTU_byindiv_pre) %in% c("indivID","maxBdLoad"))]
colSums_allOTUs <- MFOTU_byindiv_pre[,allOTUs] %>% colSums() 
allOTUs_tokeep <- allOTUs[which(colSums_allOTUs>0)]

sigcorr_otus <- c()
sigcorr_otus_neg <- c()
for ( tempOTU in allOTUs_tokeep) {
  tempData <- MFOTU_byindiv_pre[,c(tempOTU,"maxBdLoad")] %>% rename(otu = paste0(tempOTU))
  tempResult <- cor.test(tempData$otu, tempData$maxBdLoad)
  if ( tempResult$p.value < 0.05) {
    sigcorr_otus <- c(sigcorr_otus, tempOTU)
    if ( tempResult$statistic < 0 ) {
      sigcorr_otus_neg <- c(sigcorr_otus_neg, tempOTU)
    }
  }
  # tempResult <- anova(lm(maxBdLoad ~ otu, data = tempData))
  # if (tempResult$`Pr(>F)`[1] < 0.05) {
  #   sigcorr_otus <- c(sigcorr_otus, tempOTU)
  # }
}

dir.create("6a_otu_investigation/Bdload_corr")
for ( o in sigcorr_otus) {
  tempTab <- MFOTU_byindiv_pre[,c(o, "maxBdLoad")] %>% 
    rename(otu = paste0(o))
  ggsave(file = paste0("6a_otu_investigation/Bdload_corr/",o,".pdf")
       ,ggplot(tempTab) + geom_point(aes(x=otu, y=maxBdLoad))
       )
  
}

# The one neg corr
tempTab <- MFOTU_byindiv_pre[,c("TACGTAGGGTGCGAGCGTTGTCCGGAATTACTGGGCGTAAAGAGCTCGTAGGTGGTTTGTCGCGTCGTCTGTGAAATTCCGGGGCTTAAC", "maxBdLoad")] %>% 
  rename(otu = paste0("TACGTAGGGTGCGAGCGTTGTCCGGAATTACTGGGCGTAAAGAGCTCGTAGGTGGTTTGTCGCGTCGTCTGTGAAATTCCGGGGCTTAAC"))
ggplot(tempTab) + geom_point(aes(x=otu, y=maxBdLoad))
taxonomy_expanded %>% filter(Sequence %in% sigcorr_otus_neg)

# Some other interesting ones
pos_corr <- c("TACAGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGTTTGTTAAGTTGGATGTGAAATCCCCGGGCTCAAC"
  , "TACAGAGGGTGCGAGCGTTAATCGGATTTACTGGGCGTAAAGCGTACGTAGGCGGCTTTTTAAGTCGGATGTGAAATCCCTGAGCTTAAC"
  , "TACGAAGGGGGCTAGCGTTGCTCGGAATCACTGGGCGTAAAGGGCGCGTAGGCGGACTCTTAAGTCGGGGGTGAAAGCCCAGGGCTCAAC"
  , "TACGAAGGGGGCTAGCGTTGCTCGGAATTACTGGGCGTAAAGGGCGCGTAGGCGGATCGTTAAGTCAGAGGTGAAATCCCGGAGCTCAAC"
  , "TACGAAGGGGGCTAGCGTTGTTCGGAATTACTGGGCGTAAAGCGCACGTAGGCGGATATTTAAGTCAGGGGTGAAATCCCAGAGCTCAAC"
  , "TACGGAGGATCCAAGCGTTATCCGGAATCATTGGGTTTAAAGGGTCCGTAGGCGGTTTAATAAGTCAGTGGTGAAAGCCCATCGCTCAAC")
taxonomy_expanded %>% filter(Sequence %in% pos_corr)


#### Which OTUs change after infection with infection presence/absence? ####

#### Which OTUs change after infection with infection load? ####


