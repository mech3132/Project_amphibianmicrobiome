# Load packages
library(tidyverse)
library(rstanarm)
library(car) #Anova
library(vegan) # for permanova and bray curtis
library(gridExtra)
library(betareg) # for beta distr
library(lmtest) # for beta analysis
library(MASS) # for isoMDS and stepAIC
library(ade4)
library(reshape2)
library(lme4)
dir.create("./6_Statistics")

# Mapping files
load("./3_5sp_mapping_otu_downstream/mf_alt_filt_final.RData")
# OTU table of inhibitory bacteria
load("./3_5sp_mapping_otu_downstream/otu_filt.RData")

# Distance matrices
load("./3_5sp_mapping_otu_downstream/braycurtis_filt.RData")
load("./3_5sp_mapping_otu_downstream/unweighted_unifrac_filt.RData")
load("./3_5sp_mapping_otu_downstream/weighted_unifrac_filt.RData")

# Previous analyses summaries
load("./4_Bayesian_models/all_p_con.RData")
load("./4_Bayesian_models/all_p.RData")
load("./4_Bayesian_models/all_p_pred.RData")
# load("./5_logratio_tests/lr_con_inhib.RData")
# load("./5_logratio_tests/lr_treat_inhib.RData")

# Random forest LOO
load("./5_random_forest_LOO/compare_infect_onlyp_wsp_LOO.RData")
load("./5_random_forest_LOO/compare_infect_onlyp_nosp_LOO.RData")
load("./5_random_forest_LOO/compare_infect_count_wsp_LOO.RData")
load("./5_random_forest_LOO/compare_infect_count_nosp_LOO.RData")
load("./5_random_forest_LOO/compare_infect_PA_wsp_LOO.RData")
load("./5_random_forest_LOO/compare_infect_PA_nosp_LOO.RData")

load("./5_random_forest_LOO/compare_PABD_onlyp_wsp_LOO.RData")
load("./5_random_forest_LOO/compare_PABD_onlyp_nosp_LOO.RData")
load("./5_random_forest_LOO/compare_PABD_count_wsp_LOO.RData")
load("./5_random_forest_LOO/compare_PABD_count_nosp_LOO.RData")
load("./5_random_forest_LOO/compare_PABD_PA_wsp_LOO.RData")
load("./5_random_forest_LOO/compare_PABD_PA_nosp_LOO.RData")

load("./5_random_forest_LOO/importance_infect_onlyp_wsp_LOO.RData")
load("./5_random_forest_LOO/importance_infect_onlyp_nosp_LOO.RData")
load("./5_random_forest_LOO/importance_infect_count_wsp_LOO.RData")
load("./5_random_forest_LOO/importance_infect_count_nosp_LOO.RData")
load("./5_random_forest_LOO/importance_infect_PA_wsp_LOO.RData")
load("./5_random_forest_LOO/importance_infect_PA_nosp_LOO.RData")

load("./5_random_forest_LOO/importance_PABD_onlyp_wsp_LOO.RData")
load("./5_random_forest_LOO/importance_PABD_onlyp_nosp_LOO.RData")
load("./5_random_forest_LOO/importance_PABD_count_wsp_LOO.RData")
load("./5_random_forest_LOO/importance_PABD_count_nosp_LOO.RData")
load("./5_random_forest_LOO/importance_PABD_PA_wsp_LOO.RData")
load("./5_random_forest_LOO/importance_PABD_PA_nosp_LOO.RData")

load("./5_random_forest_LOO/MSE_infect_onlyp_wsp_LOO.RData")
load("./5_random_forest_LOO/MSE_infect_onlyp_nosp_LOO.RData")
load("./5_random_forest_LOO/MSE_infect_count_wsp_LOO.RData")
load("./5_random_forest_LOO/MSE_infect_count_nosp_LOO.RData")
load("./5_random_forest_LOO/MSE_infect_PA_wsp_LOO.RData")
load("./5_random_forest_LOO/MSE_infect_PA_nosp_LOO.RData")





# add a species column and PABD column
all_p <- all_p %>%
    mutate(PABD=ifelse(infect>0,1,0)) %>%
    rename(Bd_load=infect) %>%
    separate(indivID, into=c("species","indiv"), remove=FALSE) %>%
    mutate(species=factor(species, levels=c("Anbo","Rhma","Osse","Raca","Rapi"))) 
all_p_pred <- all_p_pred %>%
    mutate(PABD=ifelse(infect>0,1,0)) %>%
    rename(Bd_load=infect) %>%
    separate(indivID, into=c("species","indiv"), remove=FALSE)
all_p_pred_infectonly <- all_p_pred %>%
    filter(PABD>0) %>%
    mutate(log_Bd_load=log(Bd_load +1))

### Filter mapping file in mf_con and mf_treat
mf_con <- mf_alt_filt_final %>%
    filter(Bd_exposure=="Control")
mf_treat <- mf_alt_filt_final %>%
    filter(Bd_exposure=="Bd-exposed")
### Filter dm into con and treat
braycurtis_filt_con <- braycurtis_filt[mf_con$SampleID,mf_con$SampleID]
braycurtis_filt_treat <- braycurtis_filt[mf_treat$SampleID,mf_treat$SampleID]

unweighted_unifrac_filt_con <- unweighted_unifrac_filt[mf_con$SampleID,mf_con$SampleID]
unweighted_unifrac_filt_treat <- unweighted_unifrac_filt[mf_treat$SampleID,mf_treat$SampleID]

weighted_unifrac_filt_con <- unweighted_unifrac_filt[mf_con$SampleID,mf_con$SampleID]
weighted_unifrac_filt_treat <- unweighted_unifrac_filt[mf_treat$SampleID,mf_treat$SampleID]

# Test to see if different alpha and beta metrics show the same result
options(repr.plot.height=5, repr.plot.width=8)
gg_mf <- ggplot(data=mf_alt_filt_final)            
grid.arrange(gg_mf + geom_point(aes(x=observed_otus, y=chao1))
             , gg_mf + geom_point(aes(x=observed_otus, y=faith_pd))
             , gg_mf + geom_point(aes(x=observed_otus, y=shannon))
             , gg_mf + geom_point(aes(x=chao1, y=faith_pd))
             , gg_mf + geom_point(aes(x=chao1, y=shannon))
             , gg_mf + geom_point(aes(x=faith_pd, y=shannon))  
             , nrow=2)

## Correlate beta diversity metrics with each other
mantel.rtest(as.dist(braycurtis_filt)
, as.dist(unweighted_unifrac_filt)
            , nrepet = 99)
mantel.rtest(as.dist(braycurtis_filt)
, as.dist(weighted_unifrac_filt)
            , nrepet = 99)
mantel.rtest(as.dist(unweighted_unifrac_filt)
, as.dist(weighted_unifrac_filt)
            , nrepet = 99)

# Test to see if different beta metrics show the same result
options(repr.plot.height=10, repr.plot.width=8)
grid.arrange(gg_mf + geom_point(aes(x=disper_braycurtis, y=disper_unweighted_unifrac))
            , gg_mf + geom_point(aes(x=dist_braycurtis, y=dist_unweighted_unifrac))

             , gg_mf + geom_point(aes(x=disper_braycurtis, y=disper_weighted_unifrac))
            , gg_mf + geom_point(aes(x=dist_braycurtis, y=dist_weighted_unifrac))

             , gg_mf + geom_point(aes(x=disper_unweighted_unifrac, y=disper_weighted_unifrac))
             , gg_mf + geom_point(aes(x=dist_unweighted_unifrac, y=dist_weighted_unifrac))
             , nrow=3
            )

#### CURSORY GLANCE AT DATA ####
gg_NMDS <- mf_con %>%
    ggplot(aes(x=NMDS1_unweighted_unifrac, y=NMDS2_unweighted_unifrac)) +
    geom_point(aes(col=species), cex=3, show.legend = FALSE)
gg_infect <- mf_treat  %>%
    ggplot(aes(x=species, y=Bd_load)) +
    geom_point(aes(col=species), cex=3, position = position_jitter(width=0.1, height=0.05), show.legend = FALSE)

temp1a <-  mf_con %>%
    dplyr::select(species, shannon) %>%
    mutate(metric="Shannon_diversity") %>%
    rename(value=shannon)
temp1b <-  mf_con %>%
    dplyr::select(species, chao1) %>%
    mutate(metric="Chao1_richness") %>%
    rename(value=chao1)
temp2 <- mf_con %>%
    dplyr::select(species, inhibRich) %>%
    mutate(metric="Inhibitory_OTU_Richness")%>%
    rename(value=inhibRich)
temp3 <- mf_con %>%
    dplyr::select(species, percInhib) %>%
    mutate(metric="Percent_Inhibitory")%>%
    rename(value=percInhib)
temp4 <- mf_con %>%
    dplyr::select(species, disper_unweighted_unifrac) %>%
    mutate(metric="Dispersion_from_centroid_log")%>%
    rename(value=disper_unweighted_unifrac)
temp5 <- mf_con %>%
    dplyr::select(species, dist_unweighted_unifrac) %>%
    mutate(metric="Distance_from_previous_timepoint")%>%
    rename(value=dist_unweighted_unifrac)


gg_all <- rbind(temp1a,temp1b,temp2,temp3,temp4, temp5) %>%
    rename(Species=species) %>%
    mutate(Metric = gsub("_"," ",metric, fixed=TRUE)) %>%
    mutate(Metric = factor(Metric, levels=c("Shannon diversity","Chao1 richness","Dispersion from centroid log", "Distance from previous timepoint","Inhibitory OTU Richness","Percent Inhibitory"))) %>%
    ggplot(aes(x=Species, y=value)) +
    geom_boxplot() +
    geom_point(aes(col=Species), position = position_jitter(width=0.1, height=0), alpha=1/3)+
    facet_wrap(Metric~., scales = "free",nrow=3) +
    ylab("")+
    xlab("Species") 
lay <- rbind(c(1,2,2),
             c(3,2,2))

#+ fig.height=12, fig.width=10
options(repr.plot.height=8, repr.plot.width=12)
grid.arrange(gg_NMDS, gg_all, gg_infect, layout_matrix = lay)

print("Braycurtis")
adonis_composition_bc_timexspecies_con <- adonis2(dist(braycurtis_filt_con) ~ species*time + indivID, data=mf_con)
adonis_composition_bc_timexspecies_con
write.table(adonis_composition_bc_timexspecies_con
            , file="./6_Statistics/adonis_composition_bc_timexspecies_con.txt"
            , quote=FALSE, sep="\t", col.names=NA, row.names=TRUE)
print("Unweighted Unifrac")
adonis_composition_uwu_timexspecies_con <- adonis2(dist(unweighted_unifrac_filt_con) ~ species*time + indivID, data=mf_con)
adonis_composition_uwu_timexspecies_con
write.table(adonis_composition_uwu_timexspecies_con
            , file="./6_Statistics/adonis_composition_uwu_timexspecies_con.txt"
            , quote=FALSE, sep="\t", col.names=NA, row.names=TRUE)
print("Weighted Unifrac")
adonis_composition_wu_timexspecies_con <- adonis2(dist(weighted_unifrac_filt_con) ~ species*time + indivID, data=mf_con)
adonis_composition_wu_timexspecies_con
write.table(adonis_composition_wu_timexspecies_con
            , file="./6_Statistics/adonis_composition_wu_timexspecies_con.txt"
            , quote=FALSE, sep="\t", col.names=NA, row.names=TRUE)


print("Braycurtis")
adonis_composition_bc_timexspeciesxpabd_treat <- adonis2(dist(braycurtis_filt_treat) ~ species*time*PABD + indivID, data=mf_treat)
adonis_composition_bc_timexspeciesxpabd_treat
write.table(adonis_composition_bc_timexspeciesxpabd_treat
            , file="./6_Statistics/adonis_composition_bc_timexspeciesxpabd_treat.txt"
            , quote=FALSE, sep="\t", col.names=NA, row.names=TRUE)

print("Unweighted Unifrac")
adonis_composition_uwu_timexspeciesxpabd_treat <- adonis2(dist(unweighted_unifrac_filt_treat) ~ species*time*PABD + indivID, data=mf_treat)
adonis_composition_uwu_timexspeciesxpabd_treat
write.table(adonis_composition_uwu_timexspeciesxpabd_treat
            , file="./6_Statistics/adonis_composition_uwu_timexspeciesxpabd_treat.txt"
            , quote=FALSE, sep="\t", col.names=NA, row.names=TRUE)

print("Weighted Unifrac")
adonis_composition_wu_timexspeciesxpabd_treat <- adonis2(dist(weighted_unifrac_filt_treat) ~ species*time*PABD + indivID, data=mf_treat)
adonis_composition_wu_timexspeciesxpabd_treat
write.table(adonis_composition_wu_timexspeciesxpabd_treat
            , file="./6_Statistics/adonis_composition_wu_timexspeciesxpabd_treat.txt"
            , quote=FALSE, sep="\t", col.names=NA, row.names=TRUE)

### All data, see if there is an effect of exposure/bd on overall composition

print("Braycurtis")
adonis_composition_bc_timexspeciesxpabdxtreatment <- adonis2(dist(braycurtis_filt) ~ species + time + indivID + prepost:Bd_exposure + prepost:Bd_exposure:PABD, data=mf_alt_filt_final)
adonis_composition_bc_timexspeciesxpabdxtreatment
write.table(adonis_composition_bc_timexspeciesxpabdxtreatment
            , file="./6_Statistics/adonis_composition_bc_timexspeciesxpabdxtreatment.txt"
            , quote=FALSE, sep="\t", col.names=NA, row.names=TRUE)

print("Unweighted Unifrac")
adonis_composition_uwu_timexspeciesxpabdxtreatment <- adonis2(dist(unweighted_unifrac_filt) ~ species + time + indivID + prepost:Bd_exposure + prepost:Bd_exposure:PABD, data=mf_alt_filt_final)
adonis_composition_uwu_timexspeciesxpabdxtreatment
write.table(adonis_composition_uwu_timexspeciesxpabdxtreatment
            , file="./6_Statistics/adonis_composition_uwu_timexspeciesxpabdxtreatment.txt"
            , quote=FALSE, sep="\t", col.names=NA, row.names=TRUE)

print("Weighted Unifrac")
adonis_composition_wu_timexspeciesxpabdxtreatment <- adonis2(dist(weighted_unifrac_filt) ~ species + time + indivID + prepost:Bd_exposure + prepost:Bd_exposure:PABD, data=mf_alt_filt_final)
adonis_composition_wu_timexspeciesxpabdxtreatment
write.table(adonis_composition_wu_timexspeciesxpabdxtreatment
            , file="./6_Statistics/adonis_composition_wu_timexspeciesxpabdxtreatment.txt"
            , quote=FALSE, sep="\t", col.names=NA, row.names=TRUE)

### All data, see if there is an effect of exposure/bd on overall composition
print("Braycurtis")
adonis_composition_bc_timexspeciesxbdloadxtreatment <- adonis2(dist(braycurtis_filt) ~ species + time + indivID + prepost:Bd_exposure + prepost:Bd_exposure:Bd_load, data=mf_alt_filt_final)
adonis_composition_bc_timexspeciesxbdloadxtreatment
write.table(adonis_composition_bc_timexspeciesxbdloadxtreatment
            , file="./6_Statistics/adonis_composition_bc_timexspeciesxbdloadxtreatment.txt"
            , quote=FALSE, sep="\t", col.names=NA, row.names=TRUE)

print("Unweighted Unifrac")
adonis_composition_uwu_timexspeciesxbdloadxtreatment <- adonis2(dist(unweighted_unifrac_filt) ~ species + time + indivID + prepost:Bd_exposure + prepost:Bd_exposure:Bd_load, data=mf_alt_filt_final)
adonis_composition_uwu_timexspeciesxbdloadxtreatment
write.table(adonis_composition_uwu_timexspeciesxbdloadxtreatment
            , file="./6_Statistics/adonis_composition_uwu_timexspeciesxbdloadxtreatment.txt"
            , quote=FALSE, sep="\t", col.names=NA, row.names=TRUE)

print("Weighted Unifrac")
adonis_composition_wu_timexspeciesxbdloadxtreatment <- adonis2(dist(weighted_unifrac_filt) ~ species + time + indivID + prepost:Bd_exposure + prepost:Bd_exposure:Bd_load, data=mf_alt_filt_final)
adonis_composition_wu_timexspeciesxbdloadxtreatment
write.table(adonis_composition_wu_timexspeciesxbdloadxtreatment
            , file="./6_Statistics/adonis_composition_wu_timexspeciesxbdloadxtreatment.txt"
            , quote=FALSE, sep="\t", col.names=NA, row.names=TRUE)

### Observed_otus
# Type I ANOVA to see if interactions present
Anova(lmer(log(observed_otus) ~ species*time + (1|indivID), data=mf_con), type = 2)
# Type III ANOVA (regardless of interactions)
anova_richness_logobsotu_speciesxtime_con <- Anova(lmer(log(observed_otus) ~ species*time + (1|indivID), data=mf_con), type=3)
anova_richness_logobsotu_speciesxtime_con

write.table(anova_richness_logobsotu_speciesxtime_con
            , file="./6_Statistics/anova_richness_logobsotu_speciesxtime_con.txt"
            , quote=FALSE, sep="\t", col.names=NA, row.names=TRUE)

### Chao1
# Type I ANOVA to see if interactions present
Anova(lmer(log(chao1) ~ species*time + (1|indivID), data=mf_con), type=2)
# Type III ANOVA (regardless of interactions)
anova_richness_logchao1_speciesxtime_con <- Anova(lmer(log(chao1) ~ species*time + (1|indivID), data=mf_con), type=3)
anova_richness_logchao1_speciesxtime_con

write.table(anova_richness_logchao1_speciesxtime_con
            , file="./6_Statistics/anova_richness_logchao1_speciesxtime_con.txt"
            , quote=FALSE, sep="\t", col.names=NA, row.names=TRUE)

### Faith's PD
# Type I ANOVA to see if interactions present
Anova(lmer(log(faith_pd) ~ species*time + (1|indivID), data=mf_con), type=2)
# Type III ANOVA (regardless of interactions)
anova_richness_logfaithpd_speciesxtime_con <- Anova(lmer(log(faith_pd) ~ species*time+ (1|indivID), data=mf_con), type=3)
anova_richness_logfaithpd_speciesxtime_con

write.table(anova_richness_logfaithpd_speciesxtime_con
            , file="./6_Statistics/anova_richness_logfaithpd_speciesxtime_con.txt"
            , quote=FALSE, sep="\t", col.names=NA, row.names=TRUE)

### Shannon
# Type I ANOVA to see if interactions present
Anova(lmer(shannon ~ species*time+ (1|indivID), data=mf_con), type=2)
# Type III ANOVA (regardless of interactions)
anova_richness_shannon_speciesxtime_con <- Anova(lmer(shannon ~ species*time+ (1|indivID), data=mf_con), type=3)
anova_richness_shannon_speciesxtime_con

write.table(anova_richness_shannon_speciesxtime_con
            , file="./6_Statistics/anova_richness_shannon_speciesxtime_con.txt"
            , quote=FALSE, sep="\t", col.names=NA, row.names=TRUE)

### Observed otus
# Type I ANOVA to see if interactions present
Anova(lmer(log(observed_otus) ~ species*time+ (1|indivID), data=mf_treat), type=2)
# Type III ANOVA (regardless of interactions)
anova_richness_logobsotu_speciesxtime_treat <- Anova(lmer(log(observed_otus) ~ species*time+ (1|indivID), data=mf_treat), type=3)
anova_richness_logobsotu_speciesxtime_treat

write.table(anova_richness_logobsotu_speciesxtime_treat
            , file="./6_Statistics/anova_richness_logobsotu_speciesxtime_treat.txt"
            , quote=FALSE, sep="\t", col.names=NA, row.names=TRUE)

### Chao1
# Type I ANOVA to see if interactions present
Anova(lmer(log(chao1) ~ species*time+ (1|indivID), data=mf_treat), type=2)
# Type III ANOVA when interactions ARE present
anova_richness_logchao1_speciesxtime_treat <- Anova(lmer(log(chao1) ~ species*time + (1|indivID), data=mf_treat), type=3)
anova_richness_logchao1_speciesxtime_treat

write.table(anova_richness_logchao1_speciesxtime_treat
            , file="./6_Statistics/anova_richness_logchao1_speciesxtime_treat.txt"
            , quote=FALSE, sep="\t", col.names=NA, row.names=TRUE)

### Faith's PD
# Type I ANOVA to see if interactions present
Anova(lmer(log(faith_pd) ~ species*time+ (1|indivID), data=mf_treat), type=2)
# Type III ANOVA when interactions ARE present
anova_richness_logfaithpd_speciesxtime_treat <- Anova(lmer(log(faith_pd) ~ species*time+ (1|indivID), data=mf_treat), type=3)
anova_richness_logfaithpd_speciesxtime_treat

write.table(anova_richness_logfaithpd_speciesxtime_treat
            , file="./6_Statistics/anova_richness_logfaithpd_speciesxtime_treat.txt"
            , quote=FALSE, sep="\t", col.names=NA, row.names=TRUE)

### Shannon
# Type I ANOVA to see if interactions present
Anova(lmer(shannon ~ species*time+ (1|indivID), data=mf_treat), type=2)
# Type III ANOVA when interactions ARE present
anova_richness_shannon_speciesxtime_treat <- Anova(lmer(shannon ~ species*time+ (1|indivID), data=mf_treat), type=3)
anova_richness_shannon_speciesxtime_treat

write.table(anova_richness_logobsotu_speciesxtime_treat
            , file="./6_Statistics/anova_richness_shannon_speciesxtime_treat.txt"
            , quote=FALSE, sep="\t", col.names=NA, row.names=TRUE)

#### STOP HERE ####

### Braycurtis
# Type I ANOVA to see if interactions present
anova(lm(log(disper_braycurtis) ~ species*time, data=mf_con))
# Type III ANOVA (regardless of interactions)
anova_disper_braycurtis_speciesxtime_con <- Anova(lm(log(disper_braycurtis) ~ species*time, data=mf_con), type=3)

anova_disper_braycurtis_speciesxtime_con

write.table(anova_disper_braycurtis_speciesxtime_con
            , file="./6_Statistics/anova_disper_braycurtis_speciesxtime_con.txt"
            , quote=FALSE, sep="\t", col.names=NA, row.names=TRUE)

### Unweighted Unifrac
# Type I ANOVA to see if interactions present
anova(lm(log(disper_unweighted_unifrac) ~ species*time, data=mf_con))
# Type III ANOVA (regardless of interactions)
anova_disper_unweighted_unifrac_speciesxtime_con <- Anova(lm(log(disper_unweighted_unifrac) ~ species*time, data=mf_con), type=3)

anova_disper_unweighted_unifrac_speciesxtime_con

write.table(anova_disper_unweighted_unifrac_speciesxtime_con
            , file="./6_Statistics/anova_disper_unweighted_unifrac_speciesxtime_con.txt"
            , quote=FALSE, sep="\t", col.names=NA, row.names=TRUE)

### Unweighted Unifrac
# Type I ANOVA to see if interactions present
anova(lm(log(disper_weighted_unifrac) ~ species*time, data=mf_con))
# Type III ANOVA (regardless of interactions)
anova_disper_weighted_unifrac_speciesxtime_con <- Anova(lm(log(disper_weighted_unifrac) ~ species*time, data=mf_con), type=3)

anova_disper_weighted_unifrac_speciesxtime_con

write.table(anova_disper_weighted_unifrac_speciesxtime_con
            , file="./6_Statistics/anova_disper_weighted_unifrac_speciesxtime_con.txt"
            , quote=FALSE, sep="\t", col.names=NA, row.names=TRUE)




### Braycurtis
# Type I ANOVA to see if interactions present
anova(lm(log(disper_braycurtis) ~ species*time, data=mf_treat))
# Type III ANOVA (regardless of interactions)
anova_disper_braycurtis_speciesxtime_treat <- Anova(lm(log(disper_braycurtis) ~ species*time, data=mf_treat), type=3)

anova_disper_braycurtis_speciesxtime_treat

write.table(anova_disper_braycurtis_speciesxtime_treat
            , file="./6_Statistics/anova_disper_braycurtis_speciesxtime_treat.txt"
            , quote=FALSE, sep="\t", col.names=NA, row.names=TRUE)

### Unweighted Unifrac
# Type I ANOVA to see if interactions present
anova(lm(log(disper_unweighted_unifrac) ~ species*time, data=mf_treat))
# Type III ANOVA (regardless of interactions)
anova_disper_unweighted_unifrac_speciesxtime_treat <- Anova(lm(log(disper_unweighted_unifrac) ~ species*time, data=mf_treat), type=3)

anova_disper_unweighted_unifrac_speciesxtime_treat

write.table(anova_disper_unweighted_unifrac_speciesxtime_treat
            , file="./6_Statistics/anova_disper_unweighted_unifrac_speciesxtime_treat.txt"
            , quote=FALSE, sep="\t", col.names=NA, row.names=TRUE)

### Unweighted Unifrac
# Type I ANOVA to see if interactions present
anova(lm(log(disper_weighted_unifrac) ~ species*time, data=mf_treat))
# Type III ANOVA (regardless of interactions)
anova_disper_weighted_unifrac_speciesxtime_treat <- Anova(lm(log(disper_weighted_unifrac) ~ species*time, data=mf_treat), type=3)

anova_disper_weighted_unifrac_speciesxtime_treat

write.table(anova_disper_weighted_unifrac_speciesxtime_treat
            , file="./6_Statistics/anova_disper_weighted_unifrac_speciesxtime_treat.txt"
            , quote=FALSE, sep="\t", col.names=NA, row.names=TRUE)




### Braycurtis
# Type I ANOVA to see if interactions present
Anova(glm(dist_braycurtis ~ species*time, data=mf_con, family=Gamma(link="identity")), type=3)
# Type III ANOVA (regardless of interactions)
anova_dist_braycurtis_speciesxtime_con <- Anova(glm(dist_braycurtis ~ species*time, data=mf_con, family=Gamma(link="identity")), type=3)

anova_dist_braycurtis_speciesxtime_con

write.table(anova_dist_braycurtis_speciesxtime_con
            , file="./6_Statistics/anova_dist_braycurtis_speciesxtime_con.txt"
            , quote=FALSE, sep="\t", col.names=NA, row.names=TRUE)

### Unweighted Unifrac
# Type I ANOVA to see if interactions present
Anova(glm(dist_unweighted_unifrac ~ species*time, data=mf_con, family=Gamma(link="identity")), type=3)
# Type III ANOVA (regardless of interactions)
anova_dist_unweighted_unifrac_speciesxtime_con <- Anova(glm(dist_unweighted_unifrac ~ species*time, data=mf_con, family=Gamma(link="identity")), type=3)

anova_dist_unweighted_unifrac_speciesxtime_con

write.table(anova_dist_unweighted_unifrac_speciesxtime_con
            , file="./6_Statistics/anova_dist_unweighted_unifrac_speciesxtime_con.txt"
            , quote=FALSE, sep="\t", col.names=NA, row.names=TRUE)


### Weighted Unifrac
# Type I ANOVA to see if interactions present
Anova(glm(dist_weighted_unifrac ~ species*time, data=mf_con, family=Gamma(link="identity")), type=3)
# Type III ANOVA (regardless of interactions)
anova_dist_weighted_unifrac_speciesxtime_con <- Anova(glm(dist_weighted_unifrac ~ species*time, data=mf_con, family=Gamma(link="identity")), type=3)

anova_dist_weighted_unifrac_speciesxtime_con

write.table(anova_dist_weighted_unifrac_speciesxtime_con
            , file="./6_Statistics/anova_dist_weighted_unifrac_speciesxtime_con.txt"
            , quote=FALSE, sep="\t", col.names=NA, row.names=TRUE)





### Braycurtis
# Type I ANOVA to see if interactions present
Anova(glm(dist_braycurtis ~ species*time, data=mf_treat, family=Gamma(link="identity")), type=3)
# Type III ANOVA (regardless of interactions)
anova_dist_braycurtis_speciesxtime_treat <- Anova(glm(dist_braycurtis ~ species*time, data=mf_treat, family=Gamma(link="identity")), type=3)

anova_dist_braycurtis_speciesxtime_treat

write.table(anova_dist_braycurtis_speciesxtime_treat
            , file="./6_Statistics/anova_dist_braycurtis_speciesxtime_treat.txt"
            , quote=FALSE, sep="\t", col.names=NA, row.names=TRUE)

### Unweighted Unifrac
# Type I ANOVA to see if interactions present
Anova(glm(dist_unweighted_unifrac ~ species*time, data=mf_treat, family=Gamma(link="identity")), type=3)
# Type II ANOVA when interactions ARE NOT present
anova_dist_unweighted_unifrac_speciesxtime_treat <- Anova(glm(dist_unweighted_unifrac ~ species*time, data=mf_treat, family=Gamma(link="identity")), type=3)

anova_dist_unweighted_unifrac_speciesxtime_treat

write.table(anova_dist_unweighted_unifrac_speciesxtime_treat
            , file="./6_Statistics/anova_dist_unweighted_unifrac_speciesxtime_treat.txt"
            , quote=FALSE, sep="\t", col.names=NA, row.names=TRUE)


### Weighted Unifrac
# Type I ANOVA to see if interactions present
Anova(glm(dist_weighted_unifrac ~ species*time, data=mf_treat, family=Gamma(link="identity")), type=3)
# Type III ANOVA (regardless of interactions)
anova_dist_weighted_unifrac_speciesxtime_treat <- Anova(glm(dist_weighted_unifrac ~ species*time, data=mf_treat, family=Gamma(link="identity")), type=3)

anova_dist_weighted_unifrac_speciesxtime_treat

write.table(anova_dist_weighted_unifrac_speciesxtime_treat
            , file="./6_Statistics/anova_dist_weighted_unifrac_speciesxtime_treat.txt"
            , quote=FALSE, sep="\t", col.names=NA, row.names=TRUE)





# Type I ANOVA to see if interactions present
Anova(glm(inhibRich ~ species*time, data=mf_con, family=poisson(link="identity")), type=3)
# Type III ANOVA (regardless of interactions)
anova_inhibRich_speciesxtime_con <- Anova(glm(inhibRich ~ species*time, data=mf_con, family=poisson(link="identity")), type=3)

anova_inhibRich_speciesxtime_con

write.table(anova_inhibRich_speciesxtime_con
            , file="./6_Statistics/anova_inhibRich_speciesxtime_con.txt"
            , quote=FALSE, sep="\t", col.names=NA, row.names=TRUE)


# Type I ANOVA to see if interactions present
Anova(glm(inhibRich ~ species*time, data=mf_treat, family=poisson(link="identity")), type=3)
# Type III ANOVA (regardless of interactions)
anova_inhibRich_speciesxtime_treat <- Anova(glm(inhibRich ~ species*time, data=mf_treat, family=poisson(link="identity")), type=3)

anova_inhibRich_speciesxtime_treat

write.table(anova_inhibRich_speciesxtime_treat
            , file="./6_Statistics/anova_inhibRich_speciesxtime_treat.txt"
            , quote=FALSE, sep="\t", col.names=NA, row.names=TRUE)


# Type I ANOVA to see if interactions present
Anova(glm(percInhib ~ species*time, data=mf_con, family=Gamma), type=3)
# Type III ANOVA (regardless of interactions)
anova_percInhib_speciesxtime_con <- Anova(glm(percInhib ~ species*time, data=mf_con, family=Gamma), type=3)

anova_percInhib_speciesxtime_con

write.table(anova_percInhib_speciesxtime_con
            , file="./6_Statistics/anova_percInhib_speciesxtime_con.txt"
            , quote=FALSE, sep="\t", col.names=NA, row.names=TRUE)


# Type I ANOVA to see if interactions present
Anova(glm(percInhib ~ species*time, data=mf_treat, family=Gamma), type=3)
# Type III ANOVA (regardless of interactions)
anova_percInhib_speciesxtime_treat <- Anova(glm(percInhib ~ species*time, data=mf_treat, family=), type=3)

anova_percInhib_speciesxtime_treat

write.table(anova_percInhib_speciesxtime_treat
            , file="./6_Statistics/anova_percInhib_speciesxtime_treat.txt"
            , quote=FALSE, sep="\t", col.names=NA, row.names=TRUE)


get_stats <- function(stat_test, variable) {
    indep_var <- rownames(stat_test)
stat_col <- colnames(stat_test)
if ("LR Chisq" %in% stat_col) {
    stat_val <- stat_test[variable,"LR Chisq"]
    df_val <- stat_test[variable,"Df"]
    stat_phrase <- paste0("LR Chisq(",df_val,")=",round(stat_val,3))
    
    p_val <- stat_test[variable, "Pr(>Chisq)"]
    test <- "LR Chisq"
} else if ("F value" %in% stat_col) {
    stat_val <- stat_test[variable,"F value"]
    df_val <- stat_test[variable,"Df"]
    df_res <- stat_test["Residuals","Df"]
    stat_phrase <- paste0("F(",df_val,",",df_res,")=",round(stat_val,3))
    
    p_val <- stat_test[variable, "Pr(>F)"]
    test <- "ANOVA"

} else {
    print("NO TEST MATCHING")
}
if ( p_val < 0.001 ) {
    p_phrase <- "p<0.001"
} else {
    p_phrase <- paste0("p=",round(p_val,3))

}

output<- list()
output$p.phrase <- p_phrase
output$stat_phrase <- stat_phrase
output$all <- paste0(p_phrase,", ", stat_phrase)
output$test <- test
return(output)
}



stats_table <- data.frame("Microbial Community Trait"=c(
                                                        "Alpha diversity"
                                                       , "Alpha diversity"
                                                       , "Alpha diversity"
                                                       , "Alpha diversity"
                                                       , "Dispersion"
                                                       , "Dispersion"
                                                       , "Dispersion"
                                                       , "Distance"
                                                       , "Distance"
                                                       , "Distance"
                                                       , "Inhibitory Richness"
                                                       , "Percent Inhibitory"
                                                       )
                          , "Metric"=c(
                                       "Observed otus"
                                      , "Chao1"
                                      , "Faith's PD"
                                      , "Shannon"
                                      , "Bray-curtis"
                                      , "Unweighted Unifrac"
                                      , "Weighted Unifrac"
                                      , "Bray-curtis"
                                      , "Unweighted Unifrac"
                                      , "Weighted Unifrac"
                                      , "Count"
                                      , "Proportion"
                                      )
                         , "Test"=NA
                         , "Control_species"=NA
                         , "Control_time"=NA
                         , "Control_speciesxtime"=NA
                         , "Treatment_species"=NA
                         , "Treatment_time"=NA
                         , "Treatment_speciesxtime"=NA)
# List of all tests
all_tests <- c("anova_richness_logobsotu_speciesxtime"
,"anova_richness_logchao1_speciesxtime"
,"anova_richness_logfaithpd_speciesxtime"
,"anova_richness_shannon_speciesxtime"
,"anova_disper_braycurtis_speciesxtime"
,"anova_disper_unweighted_unifrac_speciesxtime"
,"anova_disper_weighted_unifrac_speciesxtime"
,"anova_dist_braycurtis_speciesxtime"
,"anova_dist_unweighted_unifrac_speciesxtime"
,"anova_dist_weighted_unifrac_speciesxtime"
,"anova_inhibRich_speciesxtime"
,"anova_percInhib_speciesxtime")
# Make these the rownames
rownames(stats_table) <- all_tests




for ( t in all_tests ) {
    for ( group in c("con","treat") ) {
        full_group <- ifelse (group=="con","Control_","Treatment_")
        for ( variable in c("species","time","speciesxtime")) {
            if ( variable == "speciesxtime") {
                varget = "species:time"
            } else {
                varget = variable
            }
            temp_stats <- get_stats(stat_test = get(paste0(t,"_",group)), varget)
            stats_table[t,"Test"] <- temp_stats$test
            stats_table[t,paste0(full_group,variable)] <- temp_stats$all
        }
    }
}

stats_table
write.table(stats_table, file="./6_Statistics/general_patterns_stats.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep="\t")


# anova_disper_braycurtis_speciesxtime_treat

# First, find out if any metrics are correlated to each other
all_metrics <- c("p_observed_otus"
 , "p_chao1"
 , "p_faith_pd"
 , "p_shannon"
 , "p_disper_braycurtis"
 , "p_disper_unweighted_unifrac"
 , "p_disper_weighted_unifrac"
 , "p_dist_braycurtis"
 , "p_dist_unweighted_unifrac"
 , "p_dist_weighted_unifrac"
 , "p_percInhib"
 , "p_inhibRich")

all_p_ponly <- all_p %>%
dplyr::select(one_of(all_metrics))

all_corr <- cor(all_p_ponly)
all_corr


options(repr.plot.height=4, repr.plot.width=5)
melt(all_corr) %>%
ggplot() +
geom_tile(aes(x=Var1, y=Var2, fill=value)) +
theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_gradient2(low="red", mid="white", high="green")


# Looping through all metrics

df_AIC_PABD <- data.frame(metric=all_metrics, AIC=NA, stringsAsFactors = FALSE)
for ( m in all_metrics ) {
    temp_rename <- all_p %>%
    rename(dep=paste0(m)) 
    temp_glm <- glm(PABD ~ species*dep, data=temp_rename, family=binomial)
    df_AIC_PABD[match(m,all_metrics),"AIC"] <- AIC(temp_glm)
}
df_AIC_PABD <- rbind(df_AIC_PABD, c("species_only",AIC(glm(PABD ~ species, data=temp_rename, family=binomial))))
df_AIC_PABD[order(df_AIC_PABD$AIC),]


stepAIC(glm(PABD ~ species*p_inhibRich*p_disper_unweighted_unifrac, data=all_p, family=binomial))


df_AIC_Bdload <- data.frame(metric=all_metrics, AIC=NA,stringsAsFactors = FALSE)
for ( m in all_metrics ) {
    temp_rename <- all_p %>%
    rename(dep=paste0(m)) %>%
    filter(PABD>0)
    temp_glm <- glm(Bd_load+1 ~ species*dep, data=temp_rename, family=gaussian(link="log"))
    df_AIC_Bdload[match(m,all_metrics),"AIC"] <- AIC(temp_glm)
}
df_AIC_Bdload <- rbind(df_AIC_Bdload, c("species_only",AIC(glm(Bd_load+1 ~ species, data=temp_rename, family=gaussian(link="log")))))

df_AIC_Bdload[order(df_AIC_Bdload$AIC),]

stepAIC(lm(log(Bd_load+1) ~ species*p_dist_weighted_unifrac*p_percInhib, data=temp_rename))


# Quick and dirty

anova_PABD_observed_otus <- Anova(glm(PABD ~ species*p_observed_otus, data=all_p, family=binomial))
anova_PABD_chao1 <- Anova(glm(PABD ~ species*p_chao1, data=all_p, family=binomial))
anova_PABD_faith_pd <- Anova(glm(PABD ~ species*p_faith_pd, data=all_p, family=binomial))
anova_PABD_shannon <- Anova(glm(PABD ~ species*p_shannon, data=all_p, family=binomial))

anova_PABD_dist_braycurtis <- Anova(glm(PABD ~ species*p_dist_braycurtis, data=all_p, family=binomial))
anova_PABD_dist_unweighted_unifrac <- Anova(glm(PABD ~ species*p_dist_unweighted_unifrac, data=all_p, family=binomial))
anova_PABD_dist_weighted_unifrac <- Anova(glm(PABD ~ species*p_dist_weighted_unifrac, data=all_p, family=binomial))


anova_PABD_disper_braycurtis <- Anova(glm(PABD ~ species*p_disper_braycurtis, data=all_p, family=binomial))
anova_PABD_disper_unweighted_unifrac <- Anova(glm(PABD ~ species*p_disper_unweighted_unifrac, data=all_p, family=binomial))
anova_PABD_disper_weighted_unifrac <- Anova(glm(PABD ~ species*p_disper_weighted_unifrac, data=all_p, family=binomial)) ### DIFFERENT


anova_PABD_inhibRich <- Anova(glm(PABD ~ species*p_inhibRich, data=all_p, family=binomial))
anova_PABD_percInhib <- Anova(glm(PABD ~ species*p_percInhib, data=all_p, family=binomial))


# Print to see results
anova_PABD_observed_otus
anova_PABD_chao1
anova_PABD_faith_pd
anova_PABD_shannon

anova_PABD_dist_braycurtis
anova_PABD_dist_unweighted_unifrac
anova_PABD_dist_weighted_unifrac

anova_PABD_disper_braycurtis
anova_PABD_disper_unweighted_unifrac
anova_PABD_disper_weighted_unifrac

anova_PABD_inhibRich
anova_PABD_percInhib

## Get individual correlations
species_list <- as.vector(all_p %>% dplyr::select(species) %>% pull() %>% unique())
metrics <- colnames(all_p)[grep("^p_",colnames(all_p))]
all_anova_p_PABD <- list()
for ( m in metrics) {
  all_anova_p_PABD[[m]] <- list()
  for ( sp in species_list) {
    all_p_temp <- all_p %>% filter(species ==sp)
    form <- formula(paste0("PABD ~ ",(m)))
    # test <- summary(glm(form, data=all_p_temp, family=binomial))
    # test$coefficients
    all_anova_p_PABD[[m]][[sp]] <- (summary(glm(form, data=all_p_temp, family=binomial))$coefficients)
  }
}

# Get sample size
nrow(all_p)

# Quick and dirty
all_p_noinfect <- all_p %>%
filter(PABD>0)

anova_Bdload_observed_otus <- Anova(lm(log(Bd_load+1) ~ species*p_observed_otus, data=all_p_noinfect))
anova_Bdload_chao1 <- Anova(lm(log(Bd_load+1) ~ species*p_chao1, data=all_p_noinfect))
anova_Bdload_faith_pd <- Anova(lm(log(Bd_load+1) ~ species*p_faith_pd, data=all_p_noinfect))
anova_Bdload_shannon <- Anova(lm(log(Bd_load+1) ~ species*p_shannon, data=all_p_noinfect))

anova_Bdload_dist_braycurtis <- Anova(lm(log(Bd_load+1) ~ species*p_dist_braycurtis, data=all_p_noinfect))
anova_Bdload_dist_unweighted_unifrac <- Anova(lm(log(Bd_load+1) ~ species*p_dist_unweighted_unifrac, data=all_p_noinfect))
anova_Bdload_dist_weighted_unifrac <- Anova(lm(log(Bd_load+1) ~ species*p_dist_weighted_unifrac, data=all_p_noinfect))


anova_Bdload_disper_braycurtis <- Anova(lm(log(Bd_load+1) ~ species*p_disper_braycurtis, data=all_p_noinfect))
anova_Bdload_disper_unweighted_unifrac <- Anova(lm(log(Bd_load+1) ~ species*p_disper_unweighted_unifrac, data=all_p_noinfect))
anova_Bdload_disper_weighted_unifrac <- Anova(lm(log(Bd_load+1) ~ species*p_disper_weighted_unifrac, data=all_p_noinfect))


anova_Bdload_inhibRich <- Anova(lm(log(Bd_load+1) ~ species*p_inhibRich, data=all_p_noinfect))
anova_Bdload_percInhib <- Anova(lm(log(Bd_load+1) ~ species*p_percInhib, data=all_p_noinfect))


# Print to see it
anova_Bdload_observed_otus
anova_Bdload_chao1
anova_Bdload_faith_pd
anova_Bdload_shannon

anova_Bdload_dist_braycurtis
anova_Bdload_dist_unweighted_unifrac
anova_Bdload_dist_weighted_unifrac

anova_Bdload_disper_braycurtis
anova_Bdload_disper_unweighted_unifrac
anova_Bdload_disper_weighted_unifrac

anova_Bdload_inhibRich
anova_Bdload_percInhib

## Get individual correlations
species_list <- as.vector(all_p %>% dplyr::select(species) %>% pull() %>% unique())
metrics <- colnames(all_p)[grep("^p_",colnames(all_p))]
all_anova_p_Bd <- list()
for ( m in metrics) {
  all_anova_p_Bd[[m]] <- list()
  for ( sp in species_list) {
    all_p_temp <- all_p_noinfect %>% filter(species == sp)
    if (nrow(all_p_temp)>0) {
      form <- formula(paste0("log(Bd_load+1) ~ ",(m)))
      # test <- summary(glm(form, data=all_p_temp, family=binomial))
      # test$coefficients
      all_anova_p_Bd[[m]][[sp]] <- (summary(lm(form, data=all_p_temp))$coefficients)
    }
    
  }
}

# Get sample size
nrow(all_p_noinfect)

cor.test(all_p$p_inhibRich, all_p$p_disper_unweighted_unifrac)
cor.test(all_p$p_percInhib, all_p$p_dist_weighted_unifrac)

# What is the variance explained by inhibitory richness?

summary(lm(log(Bd_load+1) ~ species*p_percInhib, data=all_p_noinfect))
anova(lm(log(Bd_load+1) ~ species*p_percInhib, data=all_p_noinfect))


# Filter out zeros
all_p_pred_infectonly <- all_p_pred %>%
filter(PABD>0)

# Presence/absence
# WITH zeros
print("Alpha diversity")
anova_observed_otus_PABD <- Anova(lm(p_observed_otus ~ species*PABD, data=all_p_pred), type=2)
anova_chao1_PABD <- Anova(lm(p_chao1 ~ species*PABD, data=all_p_pred), type=2)
anova_faith_pd_PABD <- Anova(lm(p_faith_pd ~ species*PABD, data=all_p_pred), type=2)
anova_shannon_PABD <- Anova(lm(p_shannon ~ species*PABD, data=all_p_pred), type=2)

print("Dist")
anova_dist_braycurtis_PABD <- Anova(lm(p_dist_braycurtis ~ species*PABD, data=all_p_pred), type=2)
anova_dist_unweighted_unifrac_PABD <- Anova(lm(p_dist_unweighted_unifrac ~ species*PABD, data=all_p_pred), type=2)
anova_dist_weighted_unifrac_PABD <- Anova(lm(p_dist_weighted_unifrac ~ species*PABD, data=all_p_pred), type=2)


print("Disper")
anova_disper_braycurtis_PABD <- Anova(lm(p_disper_braycurtis ~ species*PABD, data=all_p_pred), type=2)
anova_disper_unweighted_unifrac_PABD <- Anova(lm(p_disper_unweighted_unifrac ~ species*PABD, data=all_p_pred), type=2)
anova_disper_weighted_unifrac_PABD <- Anova(lm(p_disper_weighted_unifrac ~ species*PABD, data=all_p_pred), type=2)


print("inhibitory richness")
anova_inhibRich_PABD <- Anova(lm(p_inhibRich ~ species*PABD, data=all_p_pred), type=2)
print("inhibitory percent")
anova_percInhib_PABD <- Anova(lm(p_percInhib ~ species*PABD, data=all_p_pred), type=2)


# Bd Load
# NO zeros

print("Alpha diversity")
anova_observed_otus_Bdload <- Anova(lm(p_observed_otus ~ species*Bd_load, data=all_p_pred_infectonly), type=2)
anova_chao1_Bdload <- Anova(lm(p_chao1 ~ species*Bd_load, data=all_p_pred_infectonly), type=2)
anova_faith_pd_Bdload <- Anova(lm(p_faith_pd ~ species*Bd_load, data=all_p_pred_infectonly), type=2)
anova_shannon_Bdload <- Anova(lm(p_shannon ~ species*Bd_load, data=all_p_pred_infectonly), type=2)

print("Dist")
anova_dist_braycurtis_Bdload <- Anova(lm(p_dist_braycurtis ~ species*Bd_load, data=all_p_pred_infectonly), type=2)
anova_dist_unweighted_unifrac_Bdload <- Anova(lm(p_dist_unweighted_unifrac ~ species*Bd_load, data=all_p_pred_infectonly), type=2)
anova_dist_weighted_unifrac_Bdload <- Anova(lm(p_dist_weighted_unifrac ~ species*Bd_load, data=all_p_pred_infectonly), type=2)


print("Disper")
anova_disper_braycurtis_Bdload <- Anova(lm(p_disper_braycurtis ~ species*Bd_load, data=all_p_pred_infectonly), type=2)
anova_disper_unweighted_unifrac_Bdload <- Anova(lm(p_disper_unweighted_unifrac ~ species*Bd_load, data=all_p_pred_infectonly), type=2)
anova_disper_weighted_unifrac_Bdload <- Anova(lm(p_disper_weighted_unifrac ~ species*Bd_load, data=all_p_pred_infectonly), type=2)


print("inhibitory richness")
anova_inhibRich_Bdload <- Anova(lm(p_inhibRich ~ species*Bd_load, data=all_p_pred_infectonly), type=2)
print("inhibitory percent")
anova_percInhib_Bdload <- Anova(lm(p_percInhib ~ species*Bd_load, data=all_p_pred_infectonly), type=2)


stats_table_causeffect <- data.frame("Microbial Community Trait"=c(
                                                        "Alpha diversity"
                                                       , "Alpha diversity"
                                                       , "Alpha diversity"
                                                       , "Alpha diversity"
                                                       , "Dispersion"
                                                       , "Dispersion"
                                                       , "Dispersion"
                                                       , "Distance"
                                                       , "Distance"
                                                       , "Distance"
                                                       , "Inhibitory Richness"
                                                       , "Percent Inhibitory"
                                                       )
                          , "Metric"=c(
                                       "Observed otus"
                                      , "Chao1"
                                      , "Faith's PD"
                                      , "Shannon"
                                      , "Bray-curtis"
                                      , "Unweighted Unifrac"
                                      , "Weighted Unifrac"
                                      , "Bray-curtis"
                                      , "Unweighted Unifrac"
                                      , "Weighted Unifrac"
                                      , "Count"
                                      , "Proportion"
                                      )
                         , "Effect_on_PABD"=NA
                         , "Effect_on_Bd_load"=NA
                         , "Effect_of_PABD"=NA
                         , "Effect_of_Bd_load"=NA
                         )
# List of all tests
all_test_p <- c("observed_otus"
,"chao1"
,"faith_pd"
,"shannon"
,"disper_braycurtis"
,"disper_unweighted_unifrac"
,"disper_weighted_unifrac"
,"dist_braycurtis"
,"dist_unweighted_unifrac"
,"dist_weighted_unifrac"
,"inhibRich"
,"percInhib")
# Make these the rownames
rownames(stats_table_causeffect) <- all_test_p


# Copy this table to look at interactive effects separately
stats_table_causeffect_interactiononly <- stats_table_causeffect

# Extracting main effects
for ( t in all_test_p ) {
    for ( type in c("cause","effect")) {
        for ( test in c("PABD","Bdload") ) {
            if ( test=="Bdload") {
                varget<-"Bd_load"
            } else {
                varget<-test
            }
            # get fullgroup name
            if ( type=="cause") {
                full_group <- paste0("Effect_on_",varget)
            } else {
                full_group <- paste0("Effect_of_",varget)
            }
             
            if ( type=="effect") {
                temp_test <- get_stats(stat_test=get(paste0("anova_",t,"_",test)), variable=varget)

            } else if ( type=="cause" ) {
                temp_test <- get_stats(stat_test=get(paste0("anova_",test,"_",t)), variable=paste0("p_",t))
            }
        
        stats_table_causeffect[t,paste0(full_group)] <- temp_test$all
        
        }
    }
    
}

stats_table_causeffect
write.table(stats_table_causeffect, file="./6_Statistics/cause_effect_stats.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep="\t")


# Extracting interaction effects
for ( t in all_test_p ) {
    for ( type in c("cause","effect")) {
        for ( test in c("PABD","Bdload") ) {
            if ( test=="Bdload") {
                varget<-"Bd_load"
            } else {
                varget<-paste0(test)
            }
            # get fullgroup name
            if ( type=="cause") {
                full_group <- paste0("Effect_on_",varget)
            } else {
                full_group <- paste0("Effect_of_",varget)
            }
             
            if ( type=="effect") {
                temp_test <- get_stats(stat_test=get(paste0("anova_",t,"_",test)), variable=paste0("species:",varget))

            } else if ( type=="cause" ) {
                temp_test <- get_stats(stat_test=get(paste0("anova_",test,"_",t)), variable=paste0("species:p_",t))
            }
        
        stats_table_causeffect_interactiononly[t,paste0(full_group)] <- temp_test$all
        
        }
    }
    
}

stats_table_causeffect_interactiononly
write.table(stats_table_causeffect_interactiononly, file="./6_Statistics/cause_effect_interactions_stats.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep="\t")


# Bd Load
# NO zeros

print("Alpha diversity")
Anova(lm(p_observed_otus ~ species*Bd_load, data=all_p_pred), type=2)
Anova(lm(p_chao1 ~ species*Bd_load, data=all_p_pred), type=2)
Anova(lm(p_faith_pd ~ species*Bd_load, data=all_p_pred), type=2)
Anova(lm(p_shannon ~ species*Bd_load, data=all_p_pred), type=2)

print("Dist")
Anova(lm(p_dist_braycurtis ~ species*Bd_load, data=all_p_pred), type=2)
Anova(lm(p_dist_unweighted_unifrac ~ species*Bd_load, data=all_p_pred), type=2)
 Anova(lm(p_dist_weighted_unifrac ~ species*Bd_load, data=all_p_pred), type=2)


print("Disper")
Anova(lm(p_disper_braycurtis ~ species*Bd_load, data=all_p_pred), type=2)
Anova(lm(p_disper_unweighted_unifrac ~ species*Bd_load, data=all_p_pred), type=2)
Anova(lm(p_disper_weighted_unifrac ~ species*Bd_load, data=all_p_pred), type=2)


print("inhibitory richness")
Anova(lm(p_inhibRich ~ species*Bd_load, data=all_p_pred), type=2)
print("inhibitory percent")
Anova(lm(p_percInhib ~ species*Bd_load, data=all_p_pred), type=2)


Anova(glm(PABD ~ species+p_inhibRich, data=all_p, family=binomial), type=3)


options(repr.plot.height=3, repr.plot.width=10)
gg_PABD_inhibRich <-all_p %>%
ggplot(aes(y=PABD)) +
       geom_point(aes(x=p_inhibRich, col=species), show.legend=FALSE)
gg_PABD_disper_unweighted_unifrac <- all_p %>%
ggplot(aes(y=PABD)) +
       geom_point(aes(x=p_disper_unweighted_unifrac,col=species))
grid.arrange(gg_PABD_inhibRich,gg_PABD_disper_unweighted_unifrac, nrow=1)


pearson_inhibRich_observed_otus <- cor.test(all_p$p_inhibRich, all_p$p_observed_otus, method="pearson")
pearson_inhibRich_chao1 <- cor.test(all_p$p_inhibRich, all_p$p_chao1, method="pearson")
pearson_inhibRich_faithPD <- cor.test(all_p$p_inhibRich, all_p$p_faith_pd, method="pearson")
pearson_inhibRich_shannon <- cor.test(all_p$p_inhibRich, all_p$p_shannon, method="pearson")
pearson_inhibRich_percInhib <- cor.test(all_p$p_inhibRich, all_p$p_percInhib, method="pearson")

pearson_inhibRich_observed_otus
pearson_inhibRich_chao1
pearson_inhibRich_faithPD
pearson_inhibRich_shannon
pearson_inhibRich_percInhib

capture.output(pearson_inhibRich_observed_otus
            , file="./6_Statistics/pearson_inhibRich_observed_otus.txt"
            )
capture.output(pearson_inhibRich_chao1
            , file="./6_Statistics/pearson_inhibRich_chao1.txt"
            )
capture.output(pearson_inhibRich_faithPD
            , file="./6_Statistics/pearson_inhibRich_faithPD.txt"
            )
capture.output(pearson_inhibRich_shannon
            , file="./6_Statistics/pearson_inhibRich_shannon.txt"
            )
capture.output(pearson_inhibRich_percInhib
            , file="./6_Statistics/pearson_inhibRich_percInhib.txt"
            )


pearson_inhibRich_disper_UWU <- cor.test(all_p$p_inhibRich, all_p$p_disper_unweighted_unifrac, method="pearson")

pearson_inhibRich_disper_UWU
capture.output(pearson_inhibRich_disper_UWU
            , file="./6_Statistics/pearson_inhibRich_disper_UWU.txt")

options(repr.plot.height=3, repr.plot.width=10)
gg_dist_Bd_load <- all_p %>%
filter(PABD>0) %>%
ggplot(aes(y=log(Bd_load+1))) +
       geom_point(aes(x=p_dist_weighted_unifrac, col=species))

gg_percInhib_Bd_load <- all_p %>%
filter(PABD>0) %>%
ggplot(aes(y=log(Bd_load+1))) +
       geom_point(aes(x=p_percInhib, col=species))
grid.arrange(gg_dist_Bd_load,gg_percInhib_Bd_load, nrow=1)

# Effect of traits on risk of infection (PABD)
options(repr.plot.height=8, repr.plot.width=10)
all_p %>%
dplyr::select(indivID, species, p_observed_otus, p_chao1, p_shannon, p_faith_pd
       , p_disper_braycurtis, p_disper_unweighted_unifrac, p_disper_weighted_unifrac
      , p_dist_braycurtis, p_dist_unweighted_unifrac, p_dist_weighted_unifrac
      , p_percInhib, p_inhibRich, Bd_load, PABD) %>%
gather(p_observed_otus, p_chao1, p_shannon, p_faith_pd
       , p_disper_braycurtis, p_disper_unweighted_unifrac, p_disper_weighted_unifrac
      , p_dist_braycurtis, p_dist_unweighted_unifrac, p_dist_weighted_unifrac
      , p_percInhib, p_inhibRich, key=Metric, value=Value) %>%
ggplot() +
geom_point(aes(x=Value, y=PABD, col=species)) +
facet_wrap(~Metric, ncol=3)

# Effect of traits on future Bd load
options(repr.plot.height=10, repr.plot.width=10)
all_p_noinfect %>%
dplyr::select(indivID, species, p_observed_otus, p_chao1, p_shannon, p_faith_pd
       , p_disper_braycurtis, p_disper_unweighted_unifrac, p_disper_weighted_unifrac
      , p_dist_braycurtis, p_dist_unweighted_unifrac, p_dist_weighted_unifrac
      , p_percInhib, p_inhibRich, Bd_load, PABD) %>%
gather(p_observed_otus, p_chao1, p_shannon, p_faith_pd
       , p_disper_braycurtis, p_disper_unweighted_unifrac, p_disper_weighted_unifrac
      , p_dist_braycurtis, p_dist_unweighted_unifrac, p_dist_weighted_unifrac
      , p_percInhib, p_inhibRich, key=Metric, value=Value) %>%
ggplot() +
geom_point(aes(x=Value, y=Bd_load, col=species)) +
facet_wrap(~Metric, ncol=3)

# Effect of PABD on traits
options(repr.plot.height=12, repr.plot.width=10)
all_p_pred %>%
dplyr::select(indivID, species, p_observed_otus, p_chao1, p_shannon, p_faith_pd
       , p_disper_braycurtis, p_disper_unweighted_unifrac, p_disper_weighted_unifrac
      , p_dist_braycurtis, p_dist_unweighted_unifrac, p_dist_weighted_unifrac
      , p_percInhib, p_inhibRich, Bd_load, PABD) %>%
gather(p_observed_otus, p_chao1, p_shannon, p_faith_pd
       , p_disper_braycurtis, p_disper_unweighted_unifrac, p_disper_weighted_unifrac
      , p_dist_braycurtis, p_dist_unweighted_unifrac, p_dist_weighted_unifrac
      , p_percInhib, p_inhibRich, key=Metric, value=Value) %>%
ggplot() +
geom_point(aes(x=PABD, y=Value, col=species)) +
facet_wrap(~Metric, ncol=3)

# Effect of Bd load on traits
options(repr.plot.height=12, repr.plot.width=10)
all_p_pred_infectonly %>%
dplyr::select(indivID, species, p_observed_otus, p_chao1, p_shannon, p_faith_pd
       , p_disper_braycurtis, p_disper_unweighted_unifrac, p_disper_weighted_unifrac
      , p_dist_braycurtis, p_dist_unweighted_unifrac, p_dist_weighted_unifrac
      , p_percInhib, p_inhibRich, Bd_load, PABD) %>%
gather(p_observed_otus, p_chao1, p_shannon, p_faith_pd
       , p_disper_braycurtis, p_disper_unweighted_unifrac, p_disper_weighted_unifrac
      , p_dist_braycurtis, p_dist_unweighted_unifrac, p_dist_weighted_unifrac
      , p_percInhib, p_inhibRich, key=Metric, value=Value) %>%
ggplot() +
geom_point(aes(x=Bd_load, y=Value, col=species)) +
facet_wrap(~Metric, ncol=3)

print("observed_otus")
anova(lm(p_observed_otus ~ Bd_load, data=all_p_pred_infectonly))
print("chao1")
anova(lm(p_chao1 ~ Bd_load, data=all_p_pred_infectonly))
print("shannon")
anova(lm(p_shannon ~ Bd_load, data=all_p_pred_infectonly))
print("faith pd")
anova(lm(p_faith_pd ~ Bd_load, data=all_p_pred_infectonly))
print("dist_bc")
anova(lm(p_dist_braycurtis ~ Bd_load, data=all_p_pred_infectonly))
print("dist_uwu")
anova(lm(p_dist_unweighted_unifrac ~ Bd_load, data=all_p_pred_infectonly))
print("dist_wu")
anova(lm(p_dist_weighted_unifrac ~ Bd_load, data=all_p_pred_infectonly))
print("disper_bc")
anova(lm(p_disper_braycurtis ~ Bd_load, data=all_p_pred_infectonly))
print("disper_uwu")
anova(lm(p_disper_unweighted_unifrac ~ Bd_load, data=all_p_pred_infectonly))
print("disper_wu")
anova(lm(p_disper_weighted_unifrac ~ Bd_load, data=all_p_pred_infectonly))
print("percinhib")
anova(lm(p_percInhib ~ Bd_load, data=all_p_pred_infectonly))
print("inhibRich")
anova(lm(p_inhibRich ~ Bd_load, data=all_p_pred_infectonly))




adonis2(weighted_unifrac_filt ~ species + time + species:time + Bd_exposure:prepost + species:Bd_exposure:prepost + Bd_exposure:prepost:Bd_load, data=mf_alt_filt_final, by="term")
adonis2(unweighted_unifrac_filt ~ species + time + species:time + Bd_exposure:prepost + species:Bd_exposure:prepost + Bd_exposure:prepost:Bd_load, data=mf_alt_filt_final, by="term")
adonis2(braycurtis_filt ~ species + time + species:time + Bd_exposure:prepost + species:Bd_exposure:prepost + Bd_exposure:prepost:Bd_load, data=mf_alt_filt_final, by="term")


adonis2(weighted_unifrac_filt ~ species + time + species:time + Bd_exposure:prepost + species:Bd_exposure:prepost + Bd_exposure:prepost:PABD, data=mf_alt_filt_final, by="term")
adonis2(unweighted_unifrac_filt ~ species + time + species:time + Bd_exposure:prepost + species:Bd_exposure:prepost + Bd_exposure:prepost:PABD, data=mf_alt_filt_final, by="term")
adonis2(braycurtis_filt ~ species + time + species:time + Bd_exposure:prepost + species:Bd_exposure:prepost + Bd_exposure:prepost:PABD, data=mf_alt_filt_final, by="term")


errorRate_test_all <- rbind(cbind(compare_PABD_onlyp_nosp_LOO, species=FALSE), cbind(compare_PABD_onlyp_wsp_LOO, species=TRUE)
                            , cbind(compare_PABD_count_nosp_LOO, species=FALSE), cbind(compare_PABD_count_wsp_LOO, species=TRUE)
                            ,cbind(compare_PABD_PA_nosp_LOO, species=FALSE), cbind(compare_PABD_PA_wsp_LOO, species=TRUE)) %>%
  mutate(Training_data = ifelse(otu_type=="NA", "Community traits",ifelse(otu_type=="count","ASV counts","ASV prevalence")))

errorRate_test_all_summary <- errorRate_test_all %>%
  group_by(species, Training_data) %>%
  summarize(correct=mean(Test_pred==Test_obs))
errorRate_test_all_summary

MSE_test_all <- rbind(cbind(compare_infect_onlyp_nosp_LOO, species="No species predictor"), cbind(compare_infect_onlyp_wsp_LOO, species="With species as predictor")
                      , cbind(compare_infect_count_nosp_LOO, species="No species predictor"), cbind(compare_infect_count_wsp_LOO, species="With species as predictor")
                      ,cbind(compare_infect_PA_nosp_LOO, species="No species predictor"), cbind(compare_infect_PA_wsp_LOO, species="With species as predictor")) %>%
  group_by(species, otu_type) %>%
  summarize(MSE=mean((Test_obs-Test_pred)^2)) %>%
  unite(otu_type, species, col=group, remove=FALSE) %>%
  mutate(Training_data = ifelse(otu_type=="NA", "Community traits",ifelse(otu_type=="count","ASV counts","ASV prevalence"))) %>%
  dplyr::select(Training_data, species, MSE)

MSE_test_all

## Testing to see if MSE differs between training and test set


community_withsp <- MSE_infect_onlyp_wsp_LOO$error-MSE_test_all%>%filter(Training_data=="Community traits", species=="With species as predictor")%>%pull(MSE)
community_nosp <- MSE_infect_onlyp_nosp_LOO$error-MSE_test_all%>%filter(Training_data=="Community traits", species=="No species predictor")%>%pull(MSE)
count_withsp <- MSE_infect_count_wsp_LOO$error-MSE_test_all%>%filter(Training_data=="ASV counts", species=="With species as predictor")%>%pull(MSE)
count_nosp <- MSE_infect_count_nosp_LOO$error-MSE_test_all%>%filter(Training_data=="ASV counts", species=="No species predictor")%>%pull(MSE)
prevalence_wsp <- MSE_infect_PA_wsp_LOO$error-MSE_test_all%>%filter(Training_data=="ASV prevalence", species=="With species as predictor")%>%pull(MSE)
prevalence_nosp <- MSE_infect_PA_nosp_LOO$error-MSE_test_all%>%filter(Training_data=="ASV prevalence", species=="No species predictor")%>%pull(MSE)
t.test(community_withsp)
t.test(community_nosp)
t.test(count_withsp)
t.test(count_nosp)
t.test(prevalence_wsp)
t.test(prevalence_nosp)

### Looking at importance; PABD

importance_PABD_onlyp_wsp_LOO %>% group_by(taxonomy) %>%
  summarize(meanMDA=mean(MeanDecreaseAccuracy), sdMDA=sd(MeanDecreaseAccuracy)) %>%
  arrange(-meanMDA)

### Looking at importance; infect
importance_infect_onlyp_wsp_LOO %>% group_by(taxonomy) %>%
  summarize(meanIncMSE=mean(X.IncMSE), sdIncMSE=sd(X.IncMSE)) %>%
  arrange(-meanIncMSE)




