# Load packages
library(tidyverse) # for data manipulation
library(rstanarm) # For bayesian estimates of alpha and beta diversity
library(gridExtra) # For arranging ggplots
library(mgcv) # For beta distribution (beta diversity)
library(RColorBrewer) # colors for barplots
library(grid) # for text grobs in gridExtra
library(vegan) # for adonis
library(MASS) # for fitdistr
library(car) # for Anova
library(betareg) # for beta regression
library(lme4) # for lmer
source("msc_codebits_Baysian.R")
# Read in mapping files and OTU tables

load("./3_5sp_mapping_otu_downstream/mf_alt_filt_final.RData")
load("./3_5sp_mapping_otu_downstream/otu_filt.RData")
load("./3_5sp_mapping_otu_downstream/otu_filt_inhibOnly.RData")
load("./3_5sp_mapping_otu_downstream/braycurtis_filt.Rdata")
load("./3_5sp_mapping_otu_downstream/unweighted_unifrac_filt.Rdata")
load("./3_5sp_mapping_otu_downstream/weighted_unifrac_filt.Rdata")
load("./3_5sp_mapping_otu_downstream/alpha_metrics.Rdata")
load("./3_5sp_mapping_otu_downstream/beta_metrics.Rdata")

mf_con <- mf_alt_filt_final %>%
filter(Bd_exposure == "Control")

mf_treat <- mf_alt_filt_final %>%
filter(Bd_exposure == "Bd-exposed")

dir.create("./4_Bayesian_models")

#### PLOTTING EXP DESIGN ####
options(repr.plot.width=5, repr.plot.height=10)
mf_alt_filt_final %>%
    mutate(Contaminated = factor(ifelse(orig_contam ==1, "!Contaminated",NA), levels=c("!Contaminated"))
          , Bd_logload = (Bd_load)) %>%
    ggplot(aes(x=time, y=indiv)) +
    geom_line(aes(group=indivID, col=Bd_exposure)) +
    geom_point(aes(group=indivID,bg=Bd_logload), cex=4, pch=21)+
    scale_color_manual(values=c("black","blue","orange")) +
    scale_fill_gradient(low = "white", high = "red") +
    geom_vline(aes(xintercept=5.5), col="orange")+
    geom_point(aes(group=indivID, col=Contaminated), cex=1, pch=19)+ ## NEW LINE
    facet_wrap(~species, nrow=5) +
    xlab("Time") +
    ylab("Individual Toad")


#### PLOTTING BETA COMPOSITION PLOTS ####
# What do different species look like? (Control only)
options(repr.plot.width=12, repr.plot.height=12)
gg_control_bc_sp <- mf_con %>%
    ggplot(aes(x=NMDS1_braycurtis, y=NMDS2_braycurtis)) +
    geom_point(aes(col=species), cex=3) +
    ggtitle(label = paste0("Bray-Curtis (Stress:",as.character(round(mean(mf_alt_filt_final$NMDS_stress_braycurtis),0)/100),")"))

gg_control_uwu_sp <- mf_con %>%
    ggplot(aes(x=NMDS1_unweighted_unifrac, y=NMDS2_unweighted_unifrac)) +
    geom_point(aes(col=species), cex=3) +
    ggtitle(label = paste0("Unweighted Unifrac (Stress:",as.character(round(mean(mf_alt_filt_final$NMDS_stress_unweighted_unifrac),0)/100),")"))

gg_control_wu_sp <- mf_con %>%
    ggplot(aes(x=NMDS1_weighted_unifrac, y=NMDS2_weighted_unifrac)) +
    geom_point(aes(col=species), cex=3) +
    ggtitle(label = paste0("Weighted Unifrac (Stress:",as.character(round(mean(mf_alt_filt_final$NMDS_stress_weighted_unifrac),0)/100),")"))

gg_control_bc_time <- mf_con %>%
    ggplot(aes(x=NMDS1_braycurtis, y=NMDS2_braycurtis)) +
    geom_point(aes(col=time), cex=3)
gg_control_uwu_time <- mf_con %>%
    ggplot(aes(x=NMDS1_unweighted_unifrac, y=NMDS2_unweighted_unifrac)) +
    geom_point(aes(col=time), cex=3)
gg_control_wu_time <- mf_con %>%
    ggplot(aes(x=NMDS1_weighted_unifrac, y=NMDS2_weighted_unifrac)) +
    geom_point(aes(col=time), cex=3)

grid.arrange( gg_control_bc_sp,gg_control_bc_time
             , gg_control_uwu_sp,gg_control_uwu_time
             , gg_control_wu_sp,gg_control_wu_time 
            , nrow=3)


options(repr.plot.height=12, repr.plot.width=6) 
gg_all_bc <- mf_alt_filt_final %>%
    ggplot(aes(x=NMDS1_braycurtis, y=NMDS2_braycurtis)) +
    geom_point(aes(col=species,pch=Bd_exposure ), cex=2) +
    ggtitle(label = paste0("Bray-Curtis (Stress:",as.character(round(mean(mf_alt_filt_final$NMDS_stress_braycurtis),0)/100),")"))
gg_all_uwu <- mf_alt_filt_final %>%
    ggplot(aes(x=NMDS1_unweighted_unifrac, y=NMDS2_unweighted_unifrac)) +
    geom_point(aes(col=species,pch=Bd_exposure ), cex=2) +
    ggtitle(label = paste0("Unweighted Unifrac (Stress:",as.character(round(mean(mf_alt_filt_final$NMDS_stress_unweighted_unifrac),0)/100),")"))
gg_all_wu <- mf_alt_filt_final %>%
    ggplot(aes(x=NMDS1_weighted_unifrac, y=NMDS2_weighted_unifrac)) +
    geom_point(aes(col=species,pch=Bd_exposure ), cex=2) +
    ggtitle(label = paste0("Weighted Unifrac (Stress:",as.character(round(mean(mf_alt_filt_final$NMDS_stress_weighted_unifrac),0)/100),")"))
grid.arrange(gg_all_bc,gg_all_uwu, gg_all_wu
            , nrow=3)



# Check if marginal has an effect
adonis2(braycurtis_filt ~ species:time, data=mf_alt_filt_final, by="margin")
# Check main effects
adonis2(braycurtis_filt ~ species + time + species:time, data=mf_alt_filt_final)
# Include Bd_exposure
adonis2(braycurtis_filt ~ species*time*Bd_exposure, data=mf_alt_filt_final)

# Check if marginal has an effect
adonis2(unweighted_unifrac_filt ~ species:time, data=mf_alt_filt_final, by="margin")
# Check main effects
adonis2(unweighted_unifrac_filt ~ species + time + species:time, data=mf_alt_filt_final)
# Include Bd_exposure
adonis2(unweighted_unifrac_filt ~ species*time*Bd_exposure, data=mf_alt_filt_final)

# Check if marginal has an effect
adonis2(weighted_unifrac_filt ~ species:time, data=mf_alt_filt_final, by="margin")
# Check main effects
adonis2(weighted_unifrac_filt ~ species + time + species:time, data=mf_alt_filt_final)
# Include Bd_exposure
adonis2(weighted_unifrac_filt ~ species*time*Bd_exposure, data=mf_alt_filt_final)

options(repr.plot.height=5, repr.plot.width=12)
# Bray curtis
mf_alt_filt_final %>%
    mutate(Exposure=factor(prepost, levels=c("Pre","Post"))) %>%
    ggplot(aes(x=NMDS1_braycurtis,y=NMDS2_braycurtis)) +
    geom_point(aes(bg=Bd_exposure, col=Exposure), cex=2, alpha=0.8, pch=21) +
    facet_grid(Bd_exposure ~ species) +
    scale_color_manual(values=c("white","black")) +
    scale_fill_manual(values=c("blue","red"))
# Unweighted Unifrac
mf_alt_filt_final %>%
    mutate(Exposure=factor(prepost, levels=c("Pre","Post"))) %>%
    ggplot(aes(x=NMDS1_unweighted_unifrac,y=NMDS2_unweighted_unifrac)) +
    geom_point(aes(bg=Bd_exposure, col=Exposure), cex=2, alpha=0.8, pch=21) +
    facet_grid(Bd_exposure ~ species) +
    scale_color_manual(values=c("white","black")) +
    scale_fill_manual(values=c("blue","red"))
# Weighted Unifrac
mf_alt_filt_final %>%
    mutate(Exposure=factor(prepost, levels=c("Pre","Post"))) %>%
    ggplot(aes(x=NMDS1_weighted_unifrac,y=NMDS2_weighted_unifrac)) +
    geom_point(aes(bg=Bd_exposure, col=Exposure), cex=2, alpha=0.8, pch=21) +
    facet_grid(Bd_exposure ~ species) +
    scale_color_manual(values=c("white","black")) +
    scale_fill_manual(values=c("blue","red"))

options(repr.plot.height=5, repr.plot.width=12)
# Bray curtis
mf_alt_filt_final %>%
    mutate(Infection_status=factor(ifelse(PABD==0,"Not infected","Infected"), levels=c("Not infected","Infected"))) %>%
    ggplot(aes(x=NMDS1_braycurtis,y=NMDS2_braycurtis)) +
    geom_point(aes(bg=Bd_exposure, col=Infection_status), cex=2, alpha=0.8, pch=21) +
    facet_grid(Bd_exposure ~ species) +
    scale_color_manual(values=c("white","black")) +
    scale_fill_manual(values=c("blue","red"))
# Unweighted Unifrac
mf_alt_filt_final %>%
    mutate(Infection_status=factor(ifelse(PABD==0,"Not infected","Infected"), levels=c("Not infected","Infected"))) %>%
    ggplot(aes(x=NMDS1_unweighted_unifrac,y=NMDS2_unweighted_unifrac)) +
    geom_point(aes(bg=Bd_exposure, col=Infection_status), cex=2, alpha=0.8, pch=21) +
    facet_grid(Bd_exposure ~ species) +
    scale_color_manual(values=c("white","black")) +
    scale_fill_manual(values=c("blue","red"))
# Weighted Unifrac
mf_alt_filt_final %>%
    mutate(Infection_status=factor(ifelse(PABD==0,"Not infected","Infected"), levels=c("Not infected","Infected"))) %>%
    ggplot(aes(x=NMDS1_weighted_unifrac,y=NMDS2_weighted_unifrac)) +
    geom_point(aes(bg=Bd_exposure, col=Infection_status), cex=2, alpha=0.8, pch=21) +
    facet_grid(Bd_exposure ~ species) +
    scale_color_manual(values=c("white","black")) +
    scale_fill_manual(values=c("blue","red"))

options(repr.plot.height=3, repr.plot.width=12)
# Bray curtis
mf_alt_filt_final %>%
    mutate(Infection_status=factor(ifelse(PABD==0,"Not infected","Infected"), levels=c("Not infected","Infected"))) %>%
    filter(prepost=="Post") %>%
    ggplot(aes(x=NMDS1_braycurtis, y=NMDS2_braycurtis)) +
    geom_point(aes(bg=time, col=Infection_status), cex=3, pch=21) +
    scale_color_manual(values=c("white","red"))+
    facet_wrap(~species, nrow=1)
# Unweighted Unifrac
mf_alt_filt_final %>%
    mutate(Infection_status=factor(ifelse(PABD==0,"Not infected","Infected"), levels=c("Not infected","Infected"))) %>%
    filter(prepost=="Post") %>%
    ggplot(aes(x=NMDS1_unweighted_unifrac, y=NMDS2_unweighted_unifrac)) +
    geom_point(aes(bg=time, col=Infection_status), cex=3, pch=21) +
    scale_color_manual(values=c("white","red"))+
    facet_wrap(~species, nrow=1)
# Weighted Unifrac
mf_alt_filt_final %>%
    mutate(Infection_status=factor(ifelse(PABD==0,"Not infected","Infected"), levels=c("Not infected","Infected"))) %>%
    filter(prepost=="Post") %>%
    ggplot(aes(x=NMDS1_weighted_unifrac, y=NMDS2_weighted_unifrac)) +
    geom_point(aes(bg=time, col=Infection_status), cex=3, pch=21) +
    scale_color_manual(values=c("white","red"))+
    facet_wrap(~species, nrow=1)

##### Plotting alpha diversity ######
# Fit normal and lognormal distributions to each of the alpha diversity values

# Observed OTUs
    x.fit.norm <- seq(min(mf_con$observed_otus)-sd(mf_con$observed_otus)
                     , max(mf_con$observed_otus)+sd(mf_con$observed_otus)
                     , length.out = 100)
    # Fit normal distribution
    param.norm <- fitdistr(mf_con$observed_otus, densfun="normal")
    y.pred.norm <- dnorm(x.fit.norm, mean = param.norm$estimate[1], sd = param.norm$estimate[2])
    # Fit a lognormal distribution
    param.lnorm <- fitdistr(mf_con$observed_otus, densfun="lognormal")
    y.pred.lnorm <- dlnorm(x.fit.norm, meanlog=param.lnorm$estimate[1], sdlog = param.lnorm$estimate[2])
    
    ggplot(data=mf_con, aes(x=observed_otus)) +
    geom_histogram(aes(y=..density..), bins=20) +
    geom_line(data=data.frame(x=x.fit.norm, y=y.pred.norm), aes(x=x, y=y), col="red") +
    geom_line(data=data.frame(x=x.fit.norm, y=y.pred.lnorm), aes(x=x, y=y), col="blue")

# Chao1
    x.fit.norm <- seq(min(mf_con$chao1)-sd(mf_con$chao1)
                     , max(mf_con$chao1)+sd(mf_con$chao1)
                     , length.out = 100)
    # Fit normal distribution
    param.norm <- fitdistr(mf_con$chao1, densfun="normal")
    y.pred.norm <- dnorm(x.fit.norm, mean = param.norm$estimate[1], sd = param.norm$estimate[2])
    # Fit a lognormal distribution
    param.lnorm <- fitdistr(mf_con$chao1, densfun="lognormal")
    y.pred.lnorm <- dlnorm(x.fit.norm, meanlog=param.lnorm$estimate[1], sdlog = param.lnorm$estimate[2])
    
    ggplot(data=mf_con, aes(x=chao1)) +
    geom_histogram(aes(y=..density..), bins=20) +
    geom_line(data=data.frame(x=x.fit.norm, y=y.pred.norm), aes(x=x, y=y), col="red") +
    geom_line(data=data.frame(x=x.fit.norm, y=y.pred.lnorm), aes(x=x, y=y), col="blue")

# Shannon
    x.fit.norm <- seq(min(mf_con$shannon)-sd(mf_con$shannon)
                     , max(mf_con$shannon)+sd(mf_con$shannon)
                     , length.out = 100)
    # Fit normal distribution
    param.norm <- fitdistr(mf_con$shannon, densfun="normal")
    y.pred.norm <- dnorm(x.fit.norm, mean = param.norm$estimate[1], sd = param.norm$estimate[2])
    # Fit a lognormal distribution
    param.lnorm <- fitdistr(mf_con$shannon, densfun="lognormal")
    y.pred.lnorm <- dlnorm(x.fit.norm, meanlog=param.lnorm$estimate[1], sdlog = param.lnorm$estimate[2])
    
    ggplot(data=mf_con, aes(x=shannon)) +
    geom_histogram(aes(y=..density..), bins=20) +
    geom_line(data=data.frame(x=x.fit.norm, y=y.pred.norm), aes(x=x, y=y), col="red") +
    geom_line(data=data.frame(x=x.fit.norm, y=y.pred.lnorm), aes(x=x, y=y), col="blue")

# Faith's PD
    x.fit.norm <- seq(min(mf_con$faith_pd)-sd(mf_con$faith_pd)
                     , max(mf_con$faith_pd)+sd(mf_con$faith_pd)
                     , length.out = 100)
    # Fit normal distribution
    param.norm <- fitdistr(mf_con$faith_pd, densfun="normal")
    y.pred.norm <- dnorm(x.fit.norm, mean = param.norm$estimate[1], sd = param.norm$estimate[2])
    # Fit a lognormal distribution
    param.lnorm <- fitdistr(mf_con$faith_pd, densfun="lognormal")
    y.pred.lnorm <- dlnorm(x.fit.norm, meanlog=param.lnorm$estimate[1], sdlog = param.lnorm$estimate[2])
    # Fit a gamma distribution
    param.gamma <- fitdistr(mf_con$faith_pd, densfun="gamma")
    y.pred.gamma <- dgamma(x.fit.norm, shape=param.gamma$estimate[1], rate = param.gamma$estimate[2])
    
    ggplot(data=mf_con, aes(x=faith_pd)) +
    geom_histogram(aes(y=..density..), bins=20) +
    geom_line(data=data.frame(x=x.fit.norm, y=y.pred.norm), aes(x=x, y=y), col="red") +
    geom_line(data=data.frame(x=x.fit.norm, y=y.pred.lnorm), aes(x=x, y=y), col="blue") +
    geom_line(data=data.frame(x=x.fit.norm, y=y.pred.gamma), aes(x=x, y=y), col="green")

# Check to see if turnover is changing with time significantly
options(repr.plot.height=6, repr.plot.width=10)
gg_divtime_con <- mf_con %>%
    filter(!is.na(observed_otus)) %>%
    mutate(PABD=factor(PABD)) %>%
    ggplot(aes(x=time, y=observed_otus)) + 
    geom_line(aes(group=indivID)) +
    geom_point(aes(col=PABD)) +
    scale_color_manual(values=c("blue","red"))+
    facet_grid(~species)
gg_divtime_treat <- mf_treat %>%
    filter(!is.na(observed_otus)) %>%
    mutate(PABD=factor(PABD)) %>%
    ggplot(aes(x=time, y=observed_otus)) + 
    geom_line(aes(group=indivID)) +
    geom_point(aes(group=indivID, col=PABD)) +
    scale_color_manual(values=c("blue","red"))+
    geom_vline(aes(xintercept=5.5)) +
    facet_grid(~species)
grid.arrange(gg_divtime_con, gg_divtime_treat, nrow=2)


# Check to see if turnover is changing with time significantly
options(repr.plot.height=6, repr.plot.width=10)
gg_divtime_con <- mf_con %>%
    filter(!is.na(chao1)) %>%
    mutate(PABD=factor(PABD)) %>%
    ggplot(aes(x=time, y=chao1)) + 
    geom_line(aes(group=indivID)) +
    geom_point(aes(col=PABD)) +
    scale_color_manual(values=c("blue","red"))+
    facet_grid(~species)
gg_divtime_treat <- mf_treat %>%
    filter(!is.na(chao1)) %>%
    mutate(PABD=factor(PABD)) %>%
    ggplot(aes(x=time, y=chao1)) + 
    geom_line(aes(group=indivID)) +
    geom_point(aes(group=indivID, col=PABD)) +
    scale_color_manual(values=c("blue","red"))+
    geom_vline(aes(xintercept=5.5)) +
    facet_grid(~species)
grid.arrange(gg_divtime_con, gg_divtime_treat, nrow=2)


# Check to see if turnover is changing with time significantly
options(repr.plot.height=6, repr.plot.width=10)
gg_divtime_con <- mf_con %>%
    filter(!is.na(shannon)) %>%
    mutate(PABD=factor(PABD)) %>%
    ggplot(aes(x=time, y=shannon)) + 
    geom_line(aes(group=indivID)) +
    geom_point(aes(col=PABD)) +
    scale_color_manual(values=c("blue","red"))+
    facet_grid(~species)
gg_divtime_treat <- mf_treat %>%
    filter(!is.na(shannon)) %>%
    mutate(PABD=factor(PABD)) %>%
    ggplot(aes(x=time, y=shannon)) + 
    geom_line(aes(group=indivID)) +
    geom_point(aes(group=indivID, col=PABD)) +
    scale_color_manual(values=c("blue","red"))+
    geom_vline(aes(xintercept=5.5)) +
    facet_grid(~species)
grid.arrange(gg_divtime_con, gg_divtime_treat, nrow=2)


# Check to see if turnover is changing with time significantly
options(repr.plot.height=6, repr.plot.width=10)
gg_divtime_con <- mf_con %>%
    filter(!is.na(faith_pd)) %>%
    mutate(PABD=factor(PABD)) %>%
    ggplot(aes(x=time, y=faith_pd)) + 
    geom_line(aes(group=indivID)) +
    geom_point(aes(col=PABD)) +
    scale_color_manual(values=c("blue","red"))+
    facet_grid(~species)
gg_divtime_treat <- mf_treat %>%
    filter(!is.na(faith_pd)) %>%
    mutate(PABD=factor(PABD)) %>%
    ggplot(aes(x=time, y=faith_pd)) + 
    geom_line(aes(group=indivID)) +
    geom_point(aes(group=indivID, col=PABD)) +
    scale_color_manual(values=c("blue","red"))+
    geom_vline(aes(xintercept=5.5)) +
    facet_grid(~species)
grid.arrange(gg_divtime_con, gg_divtime_treat, nrow=2)

####### Alpha diversity metrics test ######
# Loop through all alpha diversity metrics
for ( a in alpha_metrics ) {
    print("-------------")
    print(a)
    assign(paste0("anova_",a), anova_log_lm_3way(mf=mf_alt_filt_final, dep=paste0(a), indep1="species", indep2="time", indep3="Bd_exposure"))
}
# Now with indivID
for ( a in alpha_metrics ) {
    print("-------------")
    print(a)
    assign(paste0("anova_lmer_",a), anova_log_lmer_3way(mf=mf_alt_filt_final, dep=paste0(a), indep1="species", indep2="time", indep3="Bd_exposure"))
}

# separated controls and Bd_exposures
for ( a in alpha_metrics ) {
    print("-------")
    print(a)
    print("CONTROL")
    assign(paste0("anova_", a, "con"), anova_log_lm_2way(mf=mf_con, dep = paste0(a), indep1="species", indep2="time"))
    print("Bd_exposure")
    assign(paste0("anova_", a, "treat"), anova_log_lm_2way(mf=mf_treat, dep = paste0(a), indep1="species", indep2="time"))
}
# separated controls and Bd_exposures, with indivID
for ( a in alpha_metrics ) {
    print("-------")
    print(a)
    print("CONTROL")
    assign(paste0("anova_lmer_", a, "con"), anova_log_lmer_2way(mf=mf_con, dep = paste0(a), indep1="species", indep2="time"))
    print("Bd_exposure")
    assign(paste0("anova_lmer_", a, "treat"), anova_log_lmer_2way(mf=mf_treat, dep = paste0(a), indep1="species", indep2="time"))
}

#### Plotting dispersion ####
# Bray-curtis
    x.fit.norm <- seq(min(mf_con$disper_braycurtis)-sd(mf_con$disper_braycurtis)
                     , max(mf_con$disper_braycurtis)+sd(mf_con$disper_braycurtis)
                     , length.out = 100)
    # Fit normal distribution
    param.norm <- fitdistr(mf_con$disper_braycurtis, densfun="normal")
    y.pred.norm <- dnorm(x.fit.norm, mean = param.norm$estimate[1], sd = param.norm$estimate[2])
    # Fit a lognormal distribution
    param.lnorm <- fitdistr(mf_con$disper_braycurtis, densfun="lognormal")
    y.pred.lnorm <- dlnorm(x.fit.norm, meanlog=param.lnorm$estimate[1], sdlog = param.lnorm$estimate[2])
    # Fit a gamma
    param.gamma <- fitdistr(mf_con$disper_braycurtis, densfun="gamma")
    y.pred.gamma <- dgamma(x.fit.norm, shape=param.gamma$estimate[1], rate = param.gamma$estimate[2])
    
    ggplot(data=mf_con, aes(x=disper_braycurtis)) +
    geom_histogram(aes(y=..density..), bins=20) +
    geom_line(data=data.frame(x=x.fit.norm, y=y.pred.norm), aes(x=x, y=y), col="red") +
    geom_line(data=data.frame(x=x.fit.norm, y=y.pred.lnorm), aes(x=x, y=y), col="blue") +
    geom_line(data=data.frame(x=x.fit.norm, y=y.pred.gamma), aes(x=x, y=y), col="green")

# Unweighted Unifrac
    x.fit.norm <- seq(min(mf_con$disper_unweighted_unifrac)-sd(mf_con$disper_unweighted_unifrac)
                     , max(mf_con$disper_unweighted_unifrac)+sd(mf_con$disper_unweighted_unifrac)
                     , length.out = 100)
    # Fit normal distribution
    param.norm <- fitdistr(mf_con$disper_unweighted_unifrac, densfun="normal")
    y.pred.norm <- dnorm(x.fit.norm, mean = param.norm$estimate[1], sd = param.norm$estimate[2])
    # Fit a lognormal distribution
    param.lnorm <- fitdistr(mf_con$disper_unweighted_unifrac, densfun="lognormal")
    y.pred.lnorm <- dlnorm(x.fit.norm, meanlog=param.lnorm$estimate[1], sdlog = param.lnorm$estimate[2])
    # Fit a gamma
    param.gamma <- fitdistr(mf_con$disper_unweighted_unifrac, densfun="gamma")
    y.pred.gamma <- dgamma(x.fit.norm, shape=param.gamma$estimate[1], rate = param.gamma$estimate[2])
    
    ggplot(data=mf_con, aes(x=disper_unweighted_unifrac)) +
    geom_histogram(aes(y=..density..), bins=20) +
    geom_line(data=data.frame(x=x.fit.norm, y=y.pred.norm), aes(x=x, y=y), col="red") +
    geom_line(data=data.frame(x=x.fit.norm, y=y.pred.lnorm), aes(x=x, y=y), col="blue") +
    geom_line(data=data.frame(x=x.fit.norm, y=y.pred.gamma), aes(x=x, y=y), col="green")

# Weighted Unifrac
    x.fit.norm <- seq(min(mf_con$disper_weighted_unifrac)-sd(mf_con$disper_weighted_unifrac)
                     , max(mf_con$disper_weighted_unifrac)+sd(mf_con$disper_weighted_unifrac)
                     , length.out = 100)
    # Fit normal distribution
    param.norm <- fitdistr(mf_con$disper_weighted_unifrac, densfun="normal")
    y.pred.norm <- dnorm(x.fit.norm, mean = param.norm$estimate[1], sd = param.norm$estimate[2])
    # Fit a lognormal distribution
    param.lnorm <- fitdistr(mf_con$disper_weighted_unifrac, densfun="lognormal")
    y.pred.lnorm <- dlnorm(x.fit.norm, meanlog=param.lnorm$estimate[1], sdlog = param.lnorm$estimate[2])
    # Fit a gamma
    param.gamma <- fitdistr(mf_con$disper_weighted_unifrac, densfun="gamma")
    y.pred.gamma <- dgamma(x.fit.norm, shape=param.gamma$estimate[1], rate = param.gamma$estimate[2])
    
    ggplot(data=mf_con, aes(x=disper_weighted_unifrac)) +
    geom_histogram(aes(y=..density..), bins=20) +
    geom_line(data=data.frame(x=x.fit.norm, y=y.pred.norm), aes(x=x, y=y), col="red") +
    geom_line(data=data.frame(x=x.fit.norm, y=y.pred.lnorm), aes(x=x, y=y), col="blue") +
    geom_line(data=data.frame(x=x.fit.norm, y=y.pred.gamma), aes(x=x, y=y), col="green")

## Braycurtis
# Check to see if turnover is changing with time significantly
gg_dispertime_con <- mf_con %>%
    filter(!is.na(disper_braycurtis)) %>%
    mutate(PABD=factor(PABD)) %>%
    ggplot(aes(x=time, y=disper_braycurtis)) + 
    geom_line(aes(group=indivID)) +
    geom_point(aes(group=indivID, col=PABD)) +
    scale_color_manual(values=c("blue","red"))+
    facet_grid(~species)
gg_dispertime_treat <- mf_treat %>%
    filter(!is.na(disper_braycurtis)) %>%
    mutate(PABD=factor(PABD)) %>%
    ggplot(aes(x=time, y=disper_braycurtis)) + 
    geom_line(aes(group=indivID)) +
    geom_point(aes(group=indivID, col=PABD)) +
    scale_color_manual(values=c("blue","red"))+
    geom_vline(aes(xintercept=5.5))+
    facet_grid(~species)
grid.arrange(gg_dispertime_con, gg_dispertime_treat, nrow=2)



## Unweighted Unifrac
# Check to see if turnover is changing with time significantly
gg_dispertime_con <- mf_con %>%
    filter(!is.na(disper_unweighted_unifrac)) %>%
    mutate(PABD=factor(PABD)) %>%
    ggplot(aes(x=time, y=disper_unweighted_unifrac)) + 
    geom_line(aes(group=indivID)) +
    geom_point(aes(group=indivID, col=PABD)) +
    scale_color_manual(values=c("blue","red"))+
    facet_grid(~species)
gg_dispertime_treat <- mf_treat %>%
    filter(!is.na(disper_unweighted_unifrac)) %>%
    mutate(PABD=factor(PABD)) %>%
    ggplot(aes(x=time, y=disper_unweighted_unifrac)) + 
    geom_line(aes(group=indivID)) +
    geom_point(aes(group=indivID, col=PABD)) +
    scale_color_manual(values=c("blue","red"))+
    geom_vline(aes(xintercept=5.5))+
    facet_grid(~species)
grid.arrange(gg_dispertime_con, gg_dispertime_treat, nrow=2)


## Weighted Unifrac
# Check to see if turnover is changing with time significantly
gg_dispertime_con <- mf_con %>%
    filter(!is.na(disper_weighted_unifrac)) %>%
    mutate(PABD=factor(PABD)) %>%
    ggplot(aes(x=time, y=disper_weighted_unifrac)) + 
    geom_line(aes(group=indivID)) +
    geom_point(aes(group=indivID, col=PABD)) +
    scale_color_manual(values=c("blue","red"))+
    facet_grid(~species)
gg_dispertime_treat <- mf_treat %>%
    filter(!is.na(disper_weighted_unifrac)) %>%
    mutate(PABD=factor(PABD)) %>%
    ggplot(aes(x=time, y=disper_weighted_unifrac)) + 
    geom_line(aes(group=indivID)) +
    geom_point(aes(group=indivID, col=PABD)) +
    scale_color_manual(values=c("blue","red"))+
    geom_vline(aes(xintercept=5.5))+
    facet_grid(~species)
grid.arrange(gg_dispertime_con, gg_dispertime_treat, nrow=2)

####### Dispersion metrics test ######
# try regular 3-way iterative test
for ( b in beta_metrics ) {
    print("--------")
    print(b)
    assign(paste0("anova_disper_",b), anova_log_lm_3way(mf=mf_alt_filt_final, dep = paste0("disper_",b), indep1="species", indep2="time", indep3="Bd_exposure"))
}
# with indivID
for ( b in beta_metrics ) {
    print("--------")
    print(b)
    assign(paste0("anova_disper_lmer_",b), anova_log_lmer_3way(mf=mf_alt_filt_final, dep = paste0("disper_",b), indep1="species", indep2="time", indep3="Bd_exposure"))
}


# Manual for weighted unifrac
Anova(lm(log(disper_weighted_unifrac) ~ species*time*Bd_exposure, data = mf_alt_filt_final), type=3)
Anova(lmer(log(disper_weighted_unifrac) ~ species*time*Bd_exposure + (1|indivID), data = mf_alt_filt_final), type=3)

for ( b in beta_metrics ) {
    print("-------")
    print(b)
    print("CONTROL")
    assign(paste0("anova_disper_",b,"_con"),anova_log_lm_2way(mf=mf_con, dep = paste0("disper_",b), indep1="species", indep2="time"))
    print("Bd_exposure")
    assign(paste0("anova_disper_",b,"_treat"),anova_log_lm_2way(mf=mf_treat, dep = paste0("disper_",b), indep1="species", indep2="time"))

}
# with indivID
for ( b in beta_metrics ) {
    print("-------")
    print(b)
    print("CONTROL")
    assign(paste0("anova_disper_lmer_",b,"_con"),anova_log_lmer_2way(mf=mf_con, dep = paste0("disper_",b), indep1="species", indep2="time"))
    print("Bd_exposure")
    assign(paste0("anova_disper_lmer_",b,"_treat"),anova_log_lmer_2way(mf=mf_treat, dep = paste0("disper_",b), indep1="species", indep2="time"))
    
}

#### Plotting distance ####

# Bray-curtis
    # eliminate NAs
    beta_values <- mf_con$dist_braycurtis[!is.na(mf_con$dist_braycurtis)]
    x.fit.norm <- seq(min(beta_values)-sd(beta_values)
                     , max(beta_values)+sd(beta_values)
                     , length.out = 100)
    # Fit normal distribution
    param.norm <- fitdistr(beta_values, densfun="normal")
    y.pred.norm <- dnorm(x.fit.norm, mean = param.norm$estimate[1], sd = param.norm$estimate[2])
    # Fit a lognormal distribution
    param.lnorm <- fitdistr(beta_values, densfun="lognormal")
    y.pred.lnorm <- dlnorm(x.fit.norm, meanlog=param.lnorm$estimate[1], sdlog = param.lnorm$estimate[2])
    # Fit a gamma
    param.gamma <- fitdistr(beta_values, densfun="gamma")
    y.pred.gamma <- dgamma(x.fit.norm, shape=param.gamma$estimate[1], rate = param.gamma$estimate[2])
    # Fit a beta
    param.beta <- fitdistr(beta_values, densfun="beta", start=list(shape1=6,shape2=6))
    y.pred.beta <- dbeta(x.fit.norm, shape1=param.beta$estimate[1], shape2 = param.beta$estimate[2])
    
    
    ggplot(data=mf_con, aes(x=dist_braycurtis)) +
    geom_histogram(aes(y=..density..), bins=20) +
    geom_line(data=data.frame(x=x.fit.norm, y=y.pred.norm), aes(x=x, y=y), col="red") +
    geom_line(data=data.frame(x=x.fit.norm, y=y.pred.lnorm), aes(x=x, y=y), col="blue") +
    geom_line(data=data.frame(x=x.fit.norm, y=y.pred.gamma), aes(x=x, y=y), col="green") +
    geom_line(data=data.frame(x=x.fit.norm, y=y.pred.beta), aes(x=x, y=y), col="purple") 



# Unweighted Unifrac
    # eliminate NAs
    beta_values <- mf_con$dist_unweighted_unifrac[!is.na(mf_con$dist_unweighted_unifrac)]
    x.fit.norm <- seq(min(beta_values)-sd(beta_values)
                     , max(beta_values)+sd(beta_values)
                     , length.out = 100)
    # Fit normal distribution
    param.norm <- fitdistr(beta_values, densfun="normal")
    y.pred.norm <- dnorm(x.fit.norm, mean = param.norm$estimate[1], sd = param.norm$estimate[2])
    # Fit a lognormal distribution
    param.lnorm <- fitdistr(beta_values, densfun="lognormal")
    y.pred.lnorm <- dlnorm(x.fit.norm, meanlog=param.lnorm$estimate[1], sdlog = param.lnorm$estimate[2])
    # Fit a gamma
    param.gamma <- fitdistr(beta_values, densfun="gamma")
    y.pred.gamma <- dgamma(x.fit.norm, shape=param.gamma$estimate[1], rate = param.gamma$estimate[2])
    # Fit a beta
    param.beta <- fitdistr(beta_values, densfun="beta", start=list(shape1=6,shape2=6))
    y.pred.beta <- dbeta(x.fit.norm, shape1=param.beta$estimate[1], shape2 = param.beta$estimate[2])
    
    ggplot(data=mf_con, aes(x=dist_unweighted_unifrac)) +
    geom_histogram(aes(y=..density..), bins=20) +
    geom_line(data=data.frame(x=x.fit.norm, y=y.pred.norm), aes(x=x, y=y), col="red") +
    geom_line(data=data.frame(x=x.fit.norm, y=y.pred.lnorm), aes(x=x, y=y), col="blue") +
    geom_line(data=data.frame(x=x.fit.norm, y=y.pred.gamma), aes(x=x, y=y), col="green") +
    geom_line(data=data.frame(x=x.fit.norm, y=y.pred.beta), aes(x=x, y=y), col="purple")

# Weighted Unifrac
    # eliminate NAs
    beta_values <- mf_con$dist_weighted_unifrac[!is.na(mf_con$dist_weighted_unifrac)]
    x.fit.norm <- seq(min(beta_values)-sd(beta_values)
                     , max(beta_values)+sd(beta_values)
                     , length.out = 100)
    # Fit normal distribution
    param.norm <- fitdistr(beta_values, densfun="normal")
    y.pred.norm <- dnorm(x.fit.norm, mean = param.norm$estimate[1], sd = param.norm$estimate[2])
    # Fit a lognormal distribution
    param.lnorm <- fitdistr(beta_values, densfun="lognormal")
    y.pred.lnorm <- dlnorm(x.fit.norm, meanlog=param.lnorm$estimate[1], sdlog = param.lnorm$estimate[2])
    # Fit a gamma
    param.gamma <- fitdistr(beta_values, densfun="gamma")
    y.pred.gamma <- dgamma(x.fit.norm, shape=param.gamma$estimate[1], rate = param.gamma$estimate[2])
    # Fit a beta
    param.beta <- fitdistr(beta_values, densfun="beta", start=list(shape1=6,shape2=6))
    y.pred.beta <- dbeta(x.fit.norm, shape1=param.beta$estimate[1], shape2 = param.beta$estimate[2])

    
    ggplot(data=mf_con, aes(x=dist_weighted_unifrac)) +
    geom_histogram(aes(y=..density..), bins=20) +
    geom_line(data=data.frame(x=x.fit.norm, y=y.pred.norm), aes(x=x, y=y), col="red") +
    geom_line(data=data.frame(x=x.fit.norm, y=y.pred.lnorm), aes(x=x, y=y), col="blue") +
    geom_line(data=data.frame(x=x.fit.norm, y=y.pred.gamma), aes(x=x, y=y), col="green") + 
    geom_line(data=data.frame(x=x.fit.norm, y=y.pred.beta), aes(x=x, y=y), col="purple")

## Braycurtis
# Check to see if turnover is changing with time significantly
gg_disttime_con <- mf_con %>%
    filter(!is.na(dist_braycurtis)) %>%
    mutate(PABD=factor(PABD)) %>%
    ggplot(aes(x=time, y=dist_braycurtis)) + 
    geom_line(aes(group=indivID)) +
    geom_point(aes(group=indivID, col=PABD)) +
    scale_color_manual(values=c("blue","red"))+
    facet_grid(~species)
gg_disttime_treat <- mf_treat %>%
    filter(!is.na(dist_braycurtis)) %>%
    mutate(PABD=factor(PABD)) %>%
    ggplot(aes(x=time, y=dist_braycurtis)) + 
    geom_line(aes(group=indivID)) +
    geom_point(aes(group=indivID, col=PABD)) +
    scale_color_manual(values=c("blue","red"))+
    geom_vline(aes(xintercept=5.5))+
    facet_grid(~species)
grid.arrange(gg_disttime_con, gg_disttime_treat, nrow=2)

## Unweighted Unifrac
# Check to see if turnover is changing with time significantly
gg_disttime_con <- mf_con %>%
    filter(!is.na(dist_unweighted_unifrac)) %>%
    mutate(PABD=factor(PABD)) %>%
    ggplot(aes(x=time, y=dist_unweighted_unifrac)) + 
    geom_line(aes(group=indivID)) +
    geom_point(aes(group=indivID, col=PABD)) +
    scale_color_manual(values=c("blue","red"))+
    facet_grid(~species)
gg_disttime_treat <- mf_treat %>%
    filter(!is.na(dist_unweighted_unifrac)) %>%
    mutate(PABD=factor(PABD)) %>%
    ggplot(aes(x=time, y=dist_unweighted_unifrac)) + 
    geom_line(aes(group=indivID)) +
    geom_point(aes(group=indivID, col=PABD)) +
    scale_color_manual(values=c("blue","red"))+
    geom_vline(aes(xintercept=5.5))+
    facet_grid(~species)
grid.arrange(gg_disttime_con, gg_disttime_treat, nrow=2)

## Weighted Unifrac
# Check to see if turnover is changing with time significantly
gg_disttime_con <- mf_con %>%
    filter(!is.na(dist_weighted_unifrac)) %>%
    mutate(PABD=factor(PABD)) %>%
    ggplot(aes(x=time, y=dist_weighted_unifrac)) + 
    geom_line(aes(group=indivID)) +
    geom_point(aes(group=indivID, col=PABD)) +
    scale_color_manual(values=c("blue","red"))+
    facet_grid(~species)
gg_disttime_treat <- mf_treat %>%
    filter(!is.na(dist_weighted_unifrac)) %>%
    mutate(PABD=factor(PABD)) %>%
    ggplot(aes(x=time, y=dist_weighted_unifrac)) + 
    geom_line(aes(group=indivID)) +
    geom_point(aes(group=indivID, col=PABD)) +
    scale_color_manual(values=c("blue","red"))+
    geom_vline(aes(xintercept=5.5))+
    facet_grid(~species)
grid.arrange(gg_disttime_con, gg_disttime_treat, nrow=2)

####### Distance metrics test ######

# No indivID versions for these, because betareg doesn't have an easy mixed effects model
for ( b in beta_metrics) {
    print("--------")
    print(b)
    assign(paste0("anova_dist_",b), anova_betareg_3way(mf=mf_alt_filt_final, dep = paste0("dist_",b), indep1="species", indep2="time", indep3="Bd_exposure"))
}


for ( b in beta_metrics ) {
    print("-------")
    print(b)
    print("CONTROL")
    assign(paste0("anova_dist_",b,"_con"),anova_betareg_2way(mf=mf_con, dep = paste0("dist_",b), indep1="species", indep2="time"))
    print("Bd_exposure")
    assign(paste0("anova_dist_",b,"_treat"),anova_betareg_2way(mf=mf_treat, dep = paste0("dist_",b), indep1="species", indep2="time"))

}

#### Plotting inhib prop #####
# Inhibitory proportion
    # eliminate NAs
    x.fit.norm <- seq(min(mf_con$percInhib)-sd(mf_con$percInhib)
                     , max(mf_con$percInhib)+sd(mf_con$percInhib)
                     , length.out = 100)
    # Fit normal distribution
    param.norm <- fitdistr(mf_con$percInhib, densfun="normal")
    y.pred.norm <- dnorm(x.fit.norm, mean = param.norm$estimate[1], sd = param.norm$estimate[2])
    # Fit a lognormal distribution
    param.lnorm <- fitdistr(mf_con$percInhib, densfun="lognormal")
    y.pred.lnorm <- dlnorm(x.fit.norm, meanlog=param.lnorm$estimate[1], sdlog = param.lnorm$estimate[2])
    # Fit a gamma
    param.gamma <- fitdistr(mf_con$percInhib, densfun="gamma")
    y.pred.gamma <- dgamma(x.fit.norm, shape=param.gamma$estimate[1], rate = param.gamma$estimate[2])
    # Fit a beta
    param.beta <- fitdistr(mf_con$percInhib, densfun="beta", start=list(shape1=6,shape2=6))
    y.pred.beta <- dbeta(x.fit.norm, shape1=param.beta$estimate[1], shape2 = param.beta$estimate[2])

    
    ggplot(data=mf_con, aes(x=percInhib)) +
    geom_histogram(aes(y=..density..), bins=20) +
    geom_line(data=data.frame(x=x.fit.norm, y=y.pred.norm), aes(x=x, y=y), col="red") +
    geom_line(data=data.frame(x=x.fit.norm, y=y.pred.lnorm), aes(x=x, y=y), col="blue") +
    geom_line(data=data.frame(x=x.fit.norm, y=y.pred.gamma), aes(x=x, y=y), col="green") + 
    geom_line(data=data.frame(x=x.fit.norm, y=y.pred.beta), aes(x=x, y=y), col="purple")

# Check to see if turnover is changing with time significantly
gg_percInhibtime_con <- mf_con %>%
    filter(!is.na(percInhib)) %>%
    mutate(PABD=factor(PABD)) %>%
    ggplot(aes(x=time, y=percInhib)) + 
    geom_line(aes(group=indivID)) +
    geom_point(aes(group=indivID, col=PABD)) +
    scale_color_manual(values=c("blue","red"))+
    facet_grid(~species)
gg_percInhibtime_treat <- mf_treat %>%
    filter(!is.na(percInhib)) %>%
    mutate(PABD=factor(PABD)) %>%
    ggplot(aes(x=time, y=percInhib)) + 
    geom_line(aes(group=indivID)) +
    geom_point(aes(group=indivID, col=PABD)) +
    scale_color_manual(values=c("blue","red"))+
    geom_vline(aes(xintercept=5.5))+
    facet_grid(~species)
grid.arrange(gg_percInhibtime_con, gg_percInhibtime_treat, nrow=2)

#### Inhib prop metrics test #####
# Check ANOVAs to see if statistical change in percent inhibitory. I think the beta distribution above looks the best.
anova_percInhib <- anova_betareg_3way(mf=mf_alt_filt_final, dep="percInhib", indep1="species", indep2="time", indep3="Bd_exposure")


# Separate control and treat
print("CONTROL")
anova_percInhib_con <- anova_betareg_2way_percInhib(mf=mf_con, dep = "percInhib", indep1="species", indep2="time")
print("Bd_exposure")
anova_percInhib_treat <- anova_betareg_2way_percInhib(mf=mf_treat, dep = "percInhib", indep1="species", indep2="time")

##### Plotting inhibitory richness #####
# Inhibitory richness
    # eliminate NAs
    inhibRich_dat <- mf_con$inhibRich
    x.fit.norm <- seq(min(inhibRich_dat)-sd(inhibRich_dat)
                     , max(inhibRich_dat)+sd(inhibRich_dat)
                     , length.out = 100)
    # Fit normal distribution
    param.norm <- fitdistr(inhibRich_dat, densfun="normal")
    y.pred.norm <- dnorm(x.fit.norm, mean = param.norm$estimate[1], sd = param.norm$estimate[2])
    # Fit a lognormal distribution
    param.lnorm <- fitdistr(inhibRich_dat, densfun="lognormal")
    y.pred.lnorm <- dlnorm(x.fit.norm, meanlog=param.lnorm$estimate[1], sdlog = param.lnorm$estimate[2])
    # Fit a gamma
    param.gamma <- fitdistr(inhibRich_dat, densfun="gamma")
    y.pred.gamma <- dgamma(x.fit.norm, shape=param.gamma$estimate[1], rate = param.gamma$estimate[2])
    # Fit a poisson
    param.pois <- fitdistr(inhibRich_dat, densfun="poisson")
    y.pred.pois <- dpois(round(x.fit.norm), lambda=param.pois$estimate[1])

    
    ggplot(data=mf_con, aes(x=inhibRich)) +
    geom_histogram(aes(y=..density..), bins=20) +
    geom_line(data=data.frame(x=x.fit.norm, y=y.pred.norm), aes(x=x, y=y), col="red") +
    geom_line(data=data.frame(x=x.fit.norm, y=y.pred.lnorm), aes(x=x, y=y), col="blue") +
    geom_line(data=data.frame(x=x.fit.norm, y=y.pred.gamma), aes(x=x, y=y), col="green") + 
    geom_line(data=data.frame(x=x.fit.norm, y=y.pred.pois), aes(x=x, y=y), col="purple")

# Check to see if turnover is changing with time significantly
gg_inhibRichtime_con <- mf_con %>%
    filter(!is.na(inhibRich)) %>%
    mutate(PABD=factor(PABD)) %>%
    ggplot(aes(x=time, y=inhibRich)) + 
    geom_line(aes(group=indivID)) +
    geom_point(aes(group=indivID, col=PABD)) +
    scale_color_manual(values=c("blue","red"))+
    facet_grid(~species)
gg_inhibRichtime_treat <- mf_treat %>%
    filter(!is.na(inhibRich)) %>%
    mutate(PABD=factor(PABD)) %>%
    ggplot(aes(x=time, y=inhibRich)) + 
    geom_line(aes(group=indivID)) +
    geom_point(aes(group=indivID, col=PABD)) +
    scale_color_manual(values=c("blue","red"))+
    geom_vline(aes(xintercept=5.5))+
    facet_grid(~species)
grid.arrange(gg_inhibRichtime_con, gg_inhibRichtime_treat, nrow=2)

##### Inhib Rich metrics test #####

# Check ANOVAs to see if statistical change in percent inhibitory. I think the beta distribution above looks the best.
anova_inhibRich <- anova_betareg_3way(mf=mf_alt_filt_final, dep="percInhib", indep1="species", indep2="time", indep3="Bd_exposure")

print("CONTROL")
anova_inhibRich_con <- anova_pois_2way(mf=mf_con, dep="inhibRich", indep1="species", indep2="time")
print("Bd_exposure")
anova_inhibRich_treat <- anova_pois_2way(mf=mf_treat, dep="inhibRich", indep1="species", indep2="time")


# Re-run everything?
RERUN_RICH = FALSE
RERUN_DIST = FALSE
RERUN_DISP = FALSE
RERUN_PERCINHIB = FALSE
RERUN_INHIBRICH = FALSE

##### RICHNESS (observed otus) (I) #######
anova_observed_otus
anova_lmer_observed_otus

if ( RERUN_RICH ) {

  lmer_log_observed_otus <- stan_lmer(log(observed_otus) ~ -1 + species + (1|indivID), data=mf_con
                              , prior = normal(0, 10, autoscale = TRUE)
                              #, family = gaussian(link="log")
                              , seed = 988735
                              , adapt_delta = 0.999

  )
  # anova_observed_otuscon
  lmer_log_observed_otus_wtime <- stan_lmer(log(observed_otus) ~ -1 + species*time + (1|indivID), data=mf_con
                                      , prior = normal(0, 10, autoscale = TRUE)
                                      #, family = gaussian(link="log")
                                      , seed = 988735
                                      , adapt_delta = 0.999
                                      
  )
  save(lmer_log_observed_otus, file="./4_Bayesian_models/lmer_log_observed_otus.RData")
  save(lmer_log_observed_otus_wtime, file="./4_Bayesian_models/lmer_log_observed_otus_wtime.RData")
} else {
  load("./4_Bayesian_models/lmer_log_observed_otus.RData")
    load("./4_Bayesian_models/lmer_log_observed_otus_wtime.RData")
}
prior_summary(lmer_log_observed_otus)
prior_summary(lmer_log_observed_otus_wtime)

observed_otus_processed <- process_glmer(g_lmer = lmer_log_observed_otus
                      , dep= "observed_otus"
                      , name_dep="log_observed_otus"
                      , transform_raw_func="log"
                      , intercept_present=F
                      , fit_distr="normal"
                      , mf_con=mf_con
                      , mf_treat=mf_treat)
observed_otus_processed$gg_model_Distribution_of_all_values
observed_otus_processed$gg_p
observed_otus_processed$gg_ExpectedDistribution_and_Bd_exposed
observed_otus_processed$gg_ExpectedDistribution_controls

observed_otus_processed_wtime <- process_glmer_withtime(g_lmer = lmer_log_observed_otus_wtime
                                         , dep= "observed_otus"
                                         , name_dep="log_observed_otus"
                                         , transform_raw_func="log"
                                         , intercept_present=F
                                         , fit_distr="normal"
                                         , mf_con=mf_con
                                         , mf_treat=mf_treat
                                         , time_factor=T
                                         , time_inter = T)
observed_otus_processed_wtime$gg_model_Distribution_of_all_values
observed_otus_processed_wtime$gg_p
observed_otus_processed_wtime$gg_ExpectedDistribution_and_Bd_exposed
observed_otus_processed_wtime$gg_ExpectedDistribution_controls



##### RICHNESS (Chao1 otus) (I) #######

if ( RERUN_RICH ) {
  
  lmer_log_chao1 <- stan_lmer(log(chao1) ~ -1 + species + (1|indivID), data=mf_con
                                      , prior = normal(0, 10, autoscale = TRUE)
                                      #, family = gaussian(link="log")
                                      , seed = 5793482
                                      , adapt_delta = 0.999
                                      
  )
  # anova_chao1con
  lmer_log_chao1_wtime <- stan_lmer(log(chao1) ~ -1 + species*time + (1|indivID), data=mf_con
                              , prior = normal(0, 10, autoscale = TRUE)
                              #, family = gaussian(link="log")
                              , seed = 5793482
                              , adapt_delta = 0.999
                              
  )
  save(lmer_log_chao1, file="./4_Bayesian_models/lmer_log_chao1.RData")
  save(lmer_log_chao1_wtime, file="./4_Bayesian_models/lmer_log_chao1_wtime.RData")
} else {
  load("./4_Bayesian_models/lmer_log_chao1.RData")
    load("./4_Bayesian_models/lmer_log_chao1_wtime.RData")
    
}
prior_summary(lmer_log_chao1)
prior_summary(lmer_log_chao1_wtime)


chao1_processed <- process_glmer(g_lmer = lmer_log_chao1
              , dep= "chao1"
              , name_dep="log_chao1"
              , transform_raw_func="log"
              , intercept_present=F
              , fit_distr="normal"
              , mf_con=mf_con
              , mf_treat=mf_treat)
chao1_processed$gg_model_Distribution_of_all_values
chao1_processed$gg_p
chao1_processed$gg_ExpectedDistribution_and_Bd_exposed
chao1_processed$gg_ExpectedDistribution_controls

chao1_processed_wtime <- process_glmer_withtime(g_lmer = lmer_log_chao1_wtime
                                 , dep= "chao1"
                                 , name_dep="log_chao1"
                                 , transform_raw_func="log"
                                 , intercept_present=F
                                 , fit_distr="normal"
                                 , mf_con=mf_con
                                 , mf_treat=mf_treat
                                 , time_factor = T
                                 , time_inter = T)
chao1_processed_wtime$gg_model_Distribution_of_all_values
chao1_processed_wtime$gg_p
chao1_processed_wtime$gg_ExpectedDistribution_and_Bd_exposed
chao1_processed_wtime$gg_ExpectedDistribution_controls

##### RICHNESS (Shannon) (I) #######

if ( RERUN_RICH ) {
  
  lmer_shannon <- stan_lmer(shannon ~ -1 + species + (1|indivID), data=mf_con
                              , prior = normal(0, 10, autoscale = TRUE)
                              , seed = 5793482
                              , adapt_delta = 0.999
                              
  )
  # anova_shannoncon
  lmer_shannon_wtime <- stan_lmer(shannon ~ -1 + species*time + (1|indivID), data=mf_con
                            , prior = normal(0, 10, autoscale = TRUE)
                            , seed = 5793482
                            , adapt_delta = 0.999
                            
  )
  save(lmer_shannon, file="./4_Bayesian_models/lmer_shannon.RData")
  save(lmer_shannon_wtime, file="./4_Bayesian_models/lmer_shannon_wtime.RData")
} else {
    load("./4_Bayesian_models/lmer_shannon.RData")
    load("./4_Bayesian_models/lmer_shannon_wtime.RData")
}
prior_summary(lmer_shannon)
prior_summary(lmer_shannon_wtime)

shannon_processed <- process_glmer(g_lmer = lmer_shannon
                                 , dep= "shannon"
                                 , name_dep="shannon"
                                 , transform_raw_func="None"
                                 , intercept_present=F
                                 , fit_distr="normal"
                                 , mf_con=mf_con
                                 , mf_treat=mf_treat)
shannon_processed$gg_model_Distribution_of_all_values
shannon_processed$gg_p
shannon_processed$gg_ExpectedDistribution_and_Bd_exposed
shannon_processed$gg_ExpectedDistribution_controls

shannon_processed_wtime <- process_glmer_withtime(g_lmer = lmer_shannon_wtime
                                   , dep= "shannon"
                                   , name_dep="shannon"
                                   , transform_raw_func="None"
                                   , intercept_present=F
                                   , fit_distr="normal"
                                   , mf_con=mf_con
                                   , mf_treat=mf_treat
                                   , time_inter=T
                                   , time_factor = T)
shannon_processed_wtime$gg_model_Distribution_of_all_values
shannon_processed_wtime$gg_p
shannon_processed_wtime$gg_ExpectedDistribution_and_Bd_exposed
shannon_processed_wtime$gg_ExpectedDistribution_controls

##### RICHNESS (Faith's PD) (I) #######

if ( RERUN_RICH ) {
  
  glmer_faith_pd <- stan_glmer(faith_pd ~ species + (1|indivID), data=mf_con
                            , prior = normal(0, 10, autoscale = TRUE)
                            , seed = 5793482
                            , family=Gamma(link="identity")
                            , adapt_delta = 0.999

  )
  # anova_faith_pdcon
  glmer_faith_pd_wtime <- stan_glmer(faith_pd ~ species*time + (1|indivID), data=mf_con
                               , prior = normal(0, 10, autoscale = TRUE)
                               , seed = 5793482
                               , family=Gamma(link="identity")
                               , adapt_delta = 0.999
                               
  )
  save(glmer_faith_pd, file="./4_Bayesian_models/glmer_faith_pd.RData")
  save(glmer_faith_pd_wtime, file="./4_Bayesian_models/glmer_faith_pd_wtime.RData")
} else {
  load("./4_Bayesian_models/glmer_faith_pd.RData")
    load("./4_Bayesian_models/glmer_faith_pd_wtime.RData")
}
prior_summary(glmer_faith_pd)
prior_summary(glmer_faith_pd_wtime)

faith_pd_processed <- process_glmer(g_lmer = glmer_faith_pd
                                   , dep= "faith_pd"
                                   , name_dep="faith_pd"
                                   , transform_raw_func="None"
                                   , intercept_present=T
                                   , fit_distr="Gamma"
                                   , mf_con = mf_con
                                   , mf_treat=mf_treat
                                   )
faith_pd_processed$gg_model_Distribution_of_all_values
faith_pd_processed$gg_p
faith_pd_processed$gg_ExpectedDistribution_and_Bd_exposed
faith_pd_processed$gg_ExpectedDistribution_controls

faith_pd_processed_wtime <- process_glmer_withtime(g_lmer = glmer_faith_pd_wtime
                                    , dep= "faith_pd"
                                    , name_dep="faith_pd"
                                    , transform_raw_func="None"
                                    , intercept_present=T
                                    , fit_distr="Gamma"
                                    , mf_con = mf_con
                                    , mf_treat=mf_treat
                                    , time_inter = T
                                    , time_factor = T
)
faith_pd_processed_wtime$gg_model_Distribution_of_all_values
faith_pd_processed_wtime$gg_p
faith_pd_processed_wtime$gg_ExpectedDistribution_and_Bd_exposed
faith_pd_processed_wtime$gg_ExpectedDistribution_controls

##### Dispersion (BC) (I) #######

if ( RERUN_DISP ) {
  
  lmer_disper_braycurtis <- stan_lmer(log(disper_braycurtis) ~ -1 + species + (1|indivID) + time
                             , data=mf_con
                             , prior_intercept = normal(location = 0,scale = 5, autoscale = TRUE)
                             , prior = normal(location=0, scale=5, autoscale=TRUE)
                             , seed= 29473
  )
  
  lmer_disper_braycurtis_wtime <- stan_lmer(log(disper_braycurtis) ~ -1 + species*time + (1|indivID)
                                      , data=mf_con
                                      , prior_intercept = normal(location = 0,scale = 5, autoscale = TRUE)
                                      , prior = normal(location=0, scale=5, autoscale=TRUE)
                                      , seed= 29473
  )
  save(lmer_disper_braycurtis, file="./4_Bayesian_models/lmer_disper_braycurtis.RData")
  save(lmer_disper_braycurtis_wtime, file="./4_Bayesian_models/lmer_disper_braycurtis_wtime.RData")
} else {
  load("./4_Bayesian_models/lmer_disper_braycurtis.RData")
    load("./4_Bayesian_models/lmer_disper_braycurtis_wtime.RData")
    
}
prior_summary(lmer_disper_braycurtis)
prior_summary(lmer_disper_braycurtis_wtime)

disper_braycurtis_processed <- process_glmer(g_lmer = lmer_disper_braycurtis
                                    , dep= "disper_braycurtis"
                                    , name_dep="log_disper_braycurtis"
                                    , transform_raw_func="log"
                                    , intercept_present=F
                                    , time_factor = T
                                    , fit_distr="normal"
                                    , mf_con=mf_con
                                    , mf_treat=mf_treat
)
disper_braycurtis_processed$gg_model_Distribution_of_all_values
disper_braycurtis_processed$gg_p
disper_braycurtis_processed$gg_ExpectedDistribution_and_Bd_exposed
disper_braycurtis_processed$gg_ExpectedDistribution_controls


disper_braycurtis_processed_wtime <- process_glmer_withtime(g_lmer = lmer_disper_braycurtis_wtime
                                             , dep= "disper_braycurtis"
                                             , name_dep="log_disper_braycurtis"
                                             , transform_raw_func="log"
                                             , intercept_present=F
                                             , time_factor = T
                                             , fit_distr="normal"
                                             , mf_con=mf_con
                                             , mf_treat=mf_treat
                                             , time_inter = T
)
disper_braycurtis_processed_wtime$gg_model_Distribution_of_all_values
disper_braycurtis_processed_wtime$gg_p
disper_braycurtis_processed_wtime$gg_ExpectedDistribution_and_Bd_exposed
disper_braycurtis_processed_wtime$gg_ExpectedDistribution_controls

##### Dispersion (UWU) (I) #######

if ( RERUN_DISP ) {
  
  lmer_disper_unweighted_unifrac <- stan_lmer(log(disper_unweighted_unifrac) ~ -1 + species + (1|indivID) + time
                                      , data=mf_con
                                      , prior_intercept = normal(location = 0,scale = 5, autoscale = TRUE)
                                      , prior = normal(location=0, scale=5, autoscale=TRUE)
                                      , seed= 29473
  )
  
  lmer_disper_unweighted_unifrac_wtime <- stan_lmer(log(disper_unweighted_unifrac) ~ -1 + species*time + (1|indivID) 
                                              , data=mf_con
                                              , prior_intercept = normal(location = 0,scale = 5, autoscale = TRUE)
                                              , prior = normal(location=0, scale=5, autoscale=TRUE)
                                              , seed= 29473
  )
  save(lmer_disper_unweighted_unifrac, file="./4_Bayesian_models/lmer_disper_unweighted_unifrac.RData")
  save(lmer_disper_unweighted_unifrac_wtime, file="./4_Bayesian_models/lmer_disper_unweighted_unifrac_wtime.RData")
} else {
    load("./4_Bayesian_models/lmer_disper_unweighted_unifrac.RData")
    load("./4_Bayesian_models/lmer_disper_unweighted_unifrac_wtime.RData")
}
prior_summary(lmer_disper_unweighted_unifrac)
prior_summary(lmer_disper_unweighted_unifrac_wtime)

disper_unweighted_unifrac_processed <- process_glmer(g_lmer = lmer_disper_unweighted_unifrac
                                             , dep= "disper_unweighted_unifrac"
                                             , name_dep="log_disper_unweighted_unifrac"
                                             , transform_raw_func="log"
                                             , intercept_present=F
                                             , time_factor = T
                                             , fit_distr="normal"
                                             , mf_con=mf_con
                                             , mf_treat=mf_treat
)
disper_unweighted_unifrac_processed$gg_model_Distribution_of_all_values
disper_unweighted_unifrac_processed$gg_p
disper_unweighted_unifrac_processed$gg_ExpectedDistribution_and_Bd_exposed
disper_unweighted_unifrac_processed$gg_ExpectedDistribution_controls

disper_unweighted_unifrac_processed_wtime <- process_glmer_withtime(g_lmer = lmer_disper_unweighted_unifrac_wtime
                                                     , dep= "disper_unweighted_unifrac"
                                                     , name_dep="log_disper_unweighted_unifrac"
                                                     , transform_raw_func="log"
                                                     , intercept_present=F
                                                     , time_factor = T
                                                     , fit_distr="normal"
                                                     , mf_con=mf_con
                                                     , mf_treat=mf_treat
                                                     , time_inter = T
)
disper_unweighted_unifrac_processed_wtime$gg_model_Distribution_of_all_values
disper_unweighted_unifrac_processed_wtime$gg_p
disper_unweighted_unifrac_processed_wtime$gg_ExpectedDistribution_and_Bd_exposed
disper_unweighted_unifrac_processed_wtime$gg_ExpectedDistribution_controls

##### Dispersion (WU) (I) #######

if ( RERUN_DISP ) {
  
  lmer_disper_weighted_unifrac <- stan_lmer(log(disper_weighted_unifrac) ~ -1 + species + (1|indivID) + time
                                              , data=mf_con
                                              , prior_intercept = normal(location = 0,scale = 5, autoscale = TRUE)
                                              , prior = normal(location=0, scale=5, autoscale=TRUE)
                                              , seed= 29473
  )

  lmer_disper_weighted_unifrac_wtime <- stan_lmer(log(disper_weighted_unifrac) ~ -1 + species*time + (1|indivID) 
                                            , data=mf_con
                                            , prior_intercept = normal(location = 0,scale = 5, autoscale = TRUE)
                                            , prior = normal(location=0, scale=5, autoscale=TRUE)
                                            , seed= 29473
  )
  save(lmer_disper_weighted_unifrac, file="./4_Bayesian_models/lmer_disper_weighted_unifrac.RData")
  save(lmer_disper_weighted_unifrac_wtime, file="./4_Bayesian_models/lmer_disper_weighted_unifrac_wtime.RData")
} else {
    load("./4_Bayesian_models/lmer_disper_weighted_unifrac.RData")
    load("./4_Bayesian_models/lmer_disper_weighted_unifrac_wtime.RData")
}
prior_summary(lmer_disper_weighted_unifrac)
prior_summary(lmer_disper_weighted_unifrac_wtime)

disper_weighted_unifrac_processed <- process_glmer(g_lmer = lmer_disper_weighted_unifrac
                                                     , dep= "disper_weighted_unifrac"
                                                     , name_dep="log_disper_weighted_unifrac"
                                                     , transform_raw_func="log"
                                                     , intercept_present=F
                                                     , time_factor = T
                                                     , fit_distr="normal"
                                                     , mf_con=mf_con
                                                     , mf_treat=mf_treat
)
disper_weighted_unifrac_processed$gg_model_Distribution_of_all_values
disper_weighted_unifrac_processed$gg_p
disper_weighted_unifrac_processed$gg_ExpectedDistribution_and_Bd_exposed
disper_weighted_unifrac_processed$gg_ExpectedDistribution_controls

disper_weighted_unifrac_processed_wtime <- process_glmer_withtime(g_lmer = lmer_disper_weighted_unifrac_wtime
                                                   , dep= "disper_weighted_unifrac"
                                                   , name_dep="log_disper_weighted_unifrac"
                                                   , transform_raw_func="log"
                                                   , intercept_present=F
                                                   , time_factor = T
                                                   , fit_distr="normal"
                                                   , mf_con=mf_con
                                                   , mf_treat=mf_treat
                                                   , time_inter = T
)
disper_weighted_unifrac_processed_wtime$gg_model_Distribution_of_all_values
disper_weighted_unifrac_processed_wtime$gg_p
disper_weighted_unifrac_processed_wtime$gg_ExpectedDistribution_and_Bd_exposed
disper_weighted_unifrac_processed_wtime$gg_ExpectedDistribution_controls

##### Distance (BC) (I) #######

if ( RERUN_DIST ) {
  
  glmer_dist_braycurtis <- stan_glmer(dist_braycurtis ~ -1 + species + (1|indivID)
                           , data=mf_con
                           , family =mgcv::betar
                           , prior_intercept = normal(location = 0.5,scale = 2.5, autoscale = TRUE)
                           , prior = normal(location=0.5, scale=2.5, autoscale=TRUE)
                           , seed= 623445
  )
  
  glmer_dist_braycurtis_wtime <- stan_glmer(dist_braycurtis ~ -1 + species*time + (1|indivID)
                                      , data=mf_con
                                      , family =mgcv::betar
                                      , prior_intercept = normal(location = 0.5,scale = 2.5, autoscale = TRUE)
                                      , prior = normal(location=0.5, scale=2.5, autoscale=TRUE)
                                      , seed= 623445
  )
  save(glmer_dist_braycurtis, file="./4_Bayesian_models/glmer_dist_braycurtis.RData")
  save(glmer_dist_braycurtis_wtime, file="./4_Bayesian_models/glmer_dist_braycurtis_wtime.RData")
} else {
    load("./4_Bayesian_models/glmer_dist_braycurtis.RData")
    load("./4_Bayesian_models/glmer_dist_braycurtis_wtime.RData")
}
prior_summary(glmer_dist_braycurtis)
prior_summary(glmer_dist_braycurtis_wtime)

dist_braycurtis_processed <- process_glmer(g_lmer = glmer_dist_braycurtis
                                    , dep= "dist_braycurtis"
                                    , name_dep="dist_braycurtis"
                                    , transform_raw_func="beta"
                                    , intercept_present=F
                                    , fit_distr="Beta"
                                    , mf_con=mf_con
                                    , mf_treat=mf_treat
)

dist_braycurtis_processed$gg_model_Distribution_of_all_values
dist_braycurtis_processed$gg_p
dist_braycurtis_processed$gg_ExpectedDistribution_and_Bd_exposed
dist_braycurtis_processed$gg_ExpectedDistribution_controls

dist_braycurtis_processed_wtime <- process_glmer_withtime(g_lmer = glmer_dist_braycurtis_wtime
                                           , dep= "dist_braycurtis"
                                           , name_dep="dist_braycurtis"
                                           , transform_raw_func="beta"
                                           , intercept_present=F
                                           , fit_distr="Beta"
                                           , mf_con=mf_con
                                           , mf_treat=mf_treat
                                           , time_factor=T
                                           , time_inter = T
                                           
)

dist_braycurtis_processed_wtime$gg_model_Distribution_of_all_values
dist_braycurtis_processed_wtime$gg_p
dist_braycurtis_processed_wtime$gg_ExpectedDistribution_and_Bd_exposed
dist_braycurtis_processed_wtime$gg_ExpectedDistribution_controls

##### Distance (UWU) (I) #######

if ( RERUN_DIST ) {
  
  glmer_dist_unweighted_unifrac <- stan_glmer(dist_unweighted_unifrac ~ -1 + species + (1|indivID)
                                      , data=mf_con
                                      , family =mgcv::betar
                                      , prior_intercept = normal(location = 0.5,scale = 2.5, autoscale = TRUE)
                                      , prior = normal(location=0.5, scale=2.5, autoscale=TRUE)
                                      , seed= 623445
  )
  
  glmer_dist_unweighted_unifrac_wtime <- stan_glmer(dist_unweighted_unifrac ~ -1 + species*time + (1|indivID)
                                              , data=mf_con
                                              , family =mgcv::betar
                                              , prior_intercept = normal(location = 0.5,scale = 2.5, autoscale = TRUE)
                                              , prior = normal(location=0.5, scale=2.5, autoscale=TRUE)
                                              , seed= 623445
  )
  save(glmer_dist_unweighted_unifrac, file="./4_Bayesian_models/glmer_dist_unweighted_unifrac.RData")
  save(glmer_dist_unweighted_unifrac_wtime, file="./4_Bayesian_models/glmer_dist_unweighted_unifrac_wtime.RData")
} else {
    load("./4_Bayesian_models/glmer_dist_unweighted_unifrac.RData")
    load("./4_Bayesian_models/glmer_dist_unweighted_unifrac_wtime.RData")
}
prior_summary(glmer_dist_unweighted_unifrac)
prior_summary(glmer_dist_unweighted_unifrac_wtime)

dist_unweighted_unifrac_processed <- process_glmer(g_lmer = glmer_dist_unweighted_unifrac
                                           , dep= "dist_unweighted_unifrac"
                                           , name_dep="dist_unweighted_unifrac"
                                           , transform_raw_func="beta"
                                           , intercept_present=F
                                           , fit_distr="Beta"
                                           , mf_con=mf_con
                                           , mf_treat=mf_treat
)

dist_unweighted_unifrac_processed$gg_model_Distribution_of_all_values
dist_unweighted_unifrac_processed$gg_p
dist_unweighted_unifrac_processed$gg_ExpectedDistribution_and_Bd_exposed
dist_unweighted_unifrac_processed$gg_ExpectedDistribution_controls

dist_unweighted_unifrac_processed_wtime <- process_glmer_withtime(g_lmer = glmer_dist_unweighted_unifrac_wtime
                                                   , dep= "dist_unweighted_unifrac"
                                                   , name_dep="dist_unweighted_unifrac"
                                                   , transform_raw_func="beta"
                                                   , intercept_present=F
                                                   , fit_distr="Beta"
                                                   , mf_con=mf_con
                                                   , mf_treat=mf_treat
                                                   , time_factor = T
                                                   , time_inter = T
                                                   
)

dist_unweighted_unifrac_processed_wtime$gg_model_Distribution_of_all_values
dist_unweighted_unifrac_processed_wtime$gg_p
dist_unweighted_unifrac_processed_wtime$gg_ExpectedDistribution_and_Bd_exposed
dist_unweighted_unifrac_processed_wtime$gg_ExpectedDistribution_controls

##### Distance (WU) (I) #######

if ( RERUN_DIST ) {
  
  glmer_dist_weighted_unifrac <- stan_glmer(dist_weighted_unifrac ~ -1 + species + (1|indivID)
                                              , data=mf_con
                                              , family =mgcv::betar
                                              , prior_intercept = normal(location = 0.5,scale = 2.5, autoscale = TRUE)
                                              , prior = normal(location=0.5, scale=2.5, autoscale=TRUE)
                                              , seed= 623445
  )
  # anova_dist_weighted_unifrac_con
  glmer_dist_weighted_unifrac_wtime <- stan_glmer(dist_weighted_unifrac ~ -1 + species*time + (1|indivID)
                                            , data=mf_con
                                            , family =mgcv::betar
                                            , prior_intercept = normal(location = 0.5,scale = 2.5, autoscale = TRUE)
                                            , prior = normal(location=0.5, scale=2.5, autoscale=TRUE)
                                            , seed= 623445
  )
  save(glmer_dist_weighted_unifrac, file="./4_Bayesian_models/glmer_dist_weighted_unifrac.RData")
  save(glmer_dist_weighted_unifrac_wtime, file="./4_Bayesian_models/glmer_dist_weighted_unifrac_wtime.RData")
} else {
    load("./4_Bayesian_models/glmer_dist_weighted_unifrac.RData")
    load("./4_Bayesian_models/glmer_dist_weighted_unifrac_wtime.RData")
}
prior_summary(glmer_dist_weighted_unifrac)
prior_summary(glmer_dist_weighted_unifrac_wtime)

dist_weighted_unifrac_processed <- process_glmer(g_lmer = glmer_dist_weighted_unifrac
                                                   , dep= "dist_weighted_unifrac"
                                                   , name_dep="dist_weighted_unifrac"
                                                   , transform_raw_func="beta"
                                                   , intercept_present=F
                                                   , fit_distr="Beta"
                                                   , mf_con=mf_con
                                                   , mf_treat=mf_treat
)

dist_weighted_unifrac_processed$gg_model_Distribution_of_all_values
dist_weighted_unifrac_processed$gg_p
dist_weighted_unifrac_processed$gg_ExpectedDistribution_and_Bd_exposed
dist_weighted_unifrac_processed$gg_ExpectedDistribution_controls

dist_weighted_unifrac_processed_wtime <- process_glmer_withtime(g_lmer = glmer_dist_weighted_unifrac_wtime
                                                 , dep= "dist_weighted_unifrac"
                                                 , name_dep="dist_weighted_unifrac"
                                                 , transform_raw_func="beta"
                                                 , intercept_present=F
                                                 , fit_distr="Beta"
                                                 , mf_con=mf_con
                                                 , mf_treat=mf_treat
                                                 , time_factor = T
                                                 , time_inter = T
)

dist_weighted_unifrac_processed_wtime$gg_model_Distribution_of_all_values
dist_weighted_unifrac_processed_wtime$gg_p
dist_weighted_unifrac_processed_wtime$gg_ExpectedDistribution_and_Bd_exposed
dist_weighted_unifrac_processed_wtime$gg_ExpectedDistribution_controls

##### Percent Inhibitory (I) #######


if ( RERUN_PERCINHIB ) {
  
  glmer_percInhib <- stan_glmer(percInhib ~ -1 + species + (1|indivID)
                                      , data=mf_con
                                      , family =mgcv::betar
                                      , prior_intercept = normal(location = 0.5,scale = 2.5, autoscale = TRUE)
                                      , prior = normal(location=0.5, scale=2.5, autoscale=TRUE)
                                      , seed= 59283
  )
  # anova_percInhib_con
  glmer_percInhib_wtime <- stan_glmer(percInhib ~ -1 + species*time + (1|indivID)
                                , data=mf_con
                                , family =mgcv::betar
                                , prior_intercept = normal(location = 0.5,scale = 2.5, autoscale = TRUE)
                                , prior = normal(location=0.5, scale=2.5, autoscale=TRUE)
                                , seed= 59283
  )
  save(glmer_percInhib, file="./4_Bayesian_models/glmer_percInhib.RData")
  save(glmer_percInhib_wtime, file="./4_Bayesian_models/glmer_percInhib_wtime.RData")
} else {
    load("./4_Bayesian_models/glmer_percInhib.RData")
    load("./4_Bayesian_models/glmer_percInhib_wtime.RData")
}
prior_summary(glmer_percInhib)
prior_summary(glmer_percInhib_wtime)

percInhib_processed <- process_glmer(g_lmer = glmer_percInhib
                                           , dep= "percInhib"
                                           , name_dep="percInhib"
                                           , transform_raw_func="beta"
                                           , intercept_present=F
                                           , fit_distr="Beta"
                                           , mf_con=mf_con
                                           , mf_treat=mf_treat
)

percInhib_processed$gg_model_Distribution_of_all_values
percInhib_processed$gg_p
percInhib_processed$gg_ExpectedDistribution_and_Bd_exposed
percInhib_processed$gg_ExpectedDistribution_controls


percInhib_processed_wtime <- process_glmer_withtime(g_lmer = glmer_percInhib_wtime
                                     , dep= "percInhib"
                                     , name_dep="percInhib"
                                     , transform_raw_func="beta"
                                     , intercept_present=F
                                     , fit_distr="Beta"
                                     , mf_con=mf_con
                                     , mf_treat=mf_treat
                                     , time_inter = T
                                     , time_factor=T
)

percInhib_processed_wtime$gg_model_Distribution_of_all_values
percInhib_processed_wtime$gg_p
percInhib_processed_wtime$gg_ExpectedDistribution_and_Bd_exposed
percInhib_processed_wtime$gg_ExpectedDistribution_controls

##### Inhibitory Richness (I) #######

if ( RERUN_INHIBRICH ) {
  # Here, I include a random variable of sample because of over-dispersion from poisson
  glmer_inhibRich <- stan_glmer(inhibRich ~ species + (1|indivID) + (1|SampleID) + time, data=mf_con
                                    , prior = normal(0, 10, autoscale = TRUE)
                                    , family= poisson(link="identity")
                                    , seed = 5423409)
  # anova_inhibRich_con
  glmer_inhibRich_wtime <- stan_glmer(inhibRich ~ species*time + (1|indivID) + (1|SampleID), data=mf_con
                                , prior = normal(0, 10, autoscale = TRUE)
                                , family= poisson(link="identity")
                                , seed = 5423409)
  save(glmer_inhibRich, file="./4_Bayesian_models/glmer_inhibRich.RData")
  save(glmer_inhibRich_wtime, file="./4_Bayesian_models/glmer_inhibRich_wtime.RData")
} else {
    load("./4_Bayesian_models/glmer_inhibRich.RData")
    load("./4_Bayesian_models/glmer_inhibRich_wtime.RData")
}
prior_summary(glmer_inhibRich)
prior_summary(glmer_inhibRich_wtime)

inhibRich_processed <- process_glmer(g_lmer = glmer_inhibRich
                                     , dep= "inhibRich"
                                     , name_dep="inhibRich"
                                     , transform_raw_func="None"
                                     , intercept_present=T
                                     , fit_distr="Poisson"
                                     , mf_con=mf_con
                                     , mf_treat=mf_treat
                                     , time_factor=T
)

inhibRich_processed$gg_model_Distribution_of_all_values
inhibRich_processed$gg_p
inhibRich_processed$gg_ExpectedDistribution_and_Bd_exposed
inhibRich_processed$gg_ExpectedDistribution_controls

inhibRich_processed_wtime <- process_glmer_withtime(g_lmer = glmer_inhibRich_wtime
                                     , dep= "inhibRich"
                                     , name_dep="inhibRich"
                                     , transform_raw_func="None"
                                     , intercept_present=T
                                     , fit_distr="Poisson"
                                     , mf_con=mf_con
                                     , mf_treat=mf_treat
                                     , time_factor=T
                                     , time_inter=T
)

inhibRich_processed_wtime$gg_model_Distribution_of_all_values
inhibRich_processed_wtime$gg_p
inhibRich_processed_wtime$gg_ExpectedDistribution_and_Bd_exposed
inhibRich_processed_wtime$gg_ExpectedDistribution_controls

# Combine all p values and save work
all_tabs <- c("observed_otus"
, "chao1"
, "shannon"
, "faith_pd"
, "disper_braycurtis"
, "disper_unweighted_unifrac"
, "disper_weighted_unifrac"
, "dist_braycurtis"
, "dist_unweighted_unifrac"
, "dist_weighted_unifrac"
, "percInhib"
, "inhibRich")
# Change names of columns so they're different
for ( tab in all_tabs) {
    temp_tab <- get(paste0(tab,"_processed_wtime"))$all_p
    name1 <- paste0("exp_",tab)
    name2 <- paste0("p_",tab)
    colnames(temp_tab) <- c("indivID",name1, name2,"infect")
    assign(paste0("all_p_",tab), temp_tab)
}
# combine all_p's
all_p <- data.frame(indivID=observed_otus_processed_wtime$all_p$indivID, infect=observed_otus_processed_wtime$all_p$infect)
for ( tab in all_tabs ) {
    all_p <- get(paste0("all_p_",tab)) %>%
    dplyr::select(indivID, paste0("exp_",tab), paste0("p_",tab)) %>%
    full_join(all_p, by="indivID")
}
save(all_p, file="./4_Bayesian_models/all_p.RData")


# Combine all p values from con and save work
# Change names of columns so they're different

for ( tab in all_tabs) {
    temp_tab <- get(paste0(tab,"_processed_wtime"))$Control_individuals
    name1 <- paste0("exp_",tab)
    name2 <- paste0("p_",tab)
    colnames(temp_tab) <- c("indivID","species","indiv",name1, name2)
    assign(paste0("all_p_con",tab), temp_tab)
}
# combine all_p's
all_p_con <- data.frame(indivID=observed_otus_processed_wtime$Control_individuals$indivID)
for ( tab in all_tabs ) {
    all_p_con <- get(paste0("all_p_con",tab)) %>%
    dplyr::select(indivID, paste0("exp_",tab), paste0("p_",tab)) %>%
    full_join(all_p_con, by="indivID")
}
save(all_p_con, file="./4_Bayesian_models/all_p_con.RData")

## Finally, combine all_p and all_p_con
all_p_temp <- all_p %>%
  mutate(Bd_exposure="Bd-exposed")
all_p_combined <- all_p_con %>%
  mutate(infect=0, Bd_exposure="Control") %>%
  rbind(all_p_temp)

save(all_p_combined, file="./4_Bayesian_models/all_p_combined.RData")

# Make a mf with pre-exposure treatment individuals
mf_all_noinfect <- mf_alt_filt_final %>%
filter(!(Bd_exposure=="Bd-exposed"&prepost=="Post"))



##### RICHNESS (observed otus) (II) #######
# There was no effect of time
if ( RERUN_RICH ) {
  lmer_log_observed_otus_all <- stan_lmer(log(observed_otus) ~ -1 + species + (1|indivID), data=mf_all_noinfect
                                      , prior = normal(0, 10, autoscale = TRUE)
                                      #, family = gaussian(link="log")
                                      , seed = 298473
                                      , adapt_delta = 0.999
                                      
  )
  
  lmer_log_observed_otus_all_wtime <- stan_lmer(log(observed_otus) ~ -1 + species*time + (1|indivID), data=mf_all_noinfect
                                          , prior = normal(0, 10, autoscale = TRUE)
                                          #, family = gaussian(link="log")
                                          , seed = 298473
                                          , adapt_delta = 0.999
                                          
  )
  save(lmer_log_observed_otus_all, file="./4_Bayesian_models/lmer_log_observed_otus_all.RData")
  save(lmer_log_observed_otus_all_wtime, file="./4_Bayesian_models/lmer_log_observed_otus_all_wtime.RData")
} else {
    load("./4_Bayesian_models/lmer_log_observed_otus_all.RData")
    load("./4_Bayesian_models/lmer_log_observed_otus_all_wtime.RData")
}
prior_summary(lmer_log_observed_otus_all)
prior_summary(lmer_log_observed_otus_all_wtime)


observed_otus_processed_all <- process_glmer_all(g_lmer = lmer_log_observed_otus_all
                                         , dep= "observed_otus"
                                         , name_dep="log_observed_otus"
                                         , transform_raw_func="log"
                                         , intercept_present=F
                                         , fit_distr="normal"
                                         , mf_con=mf_con
                                         , mf_treat=mf_treat)
observed_otus_processed_all$gg_model_Distribution_of_all_values
observed_otus_processed_all$gg_p
observed_otus_processed_all$gg_ExpectedDistribution_and_Bd_exposed
observed_otus_processed_all$gg_ExpectedDistribution_controls

observed_otus_processed_all_wtime <- process_glmer_all_wtime(g_lmer = lmer_log_observed_otus_all_wtime
                                                 , dep= "observed_otus"
                                                 , name_dep="log_observed_otus"
                                                 , transform_raw_func="log"
                                                 , intercept_present=F
                                                 , fit_distr="normal"
                                                 , mf_con=mf_con
                                                 , mf_treat=mf_treat
                                                 , time_factor = T
                                                 , time_inter = T)
observed_otus_processed_all_wtime$gg_model_Distribution_of_all_values
observed_otus_processed_all_wtime$gg_p
observed_otus_processed_all_wtime$gg_ExpectedDistribution_and_Bd_exposed
observed_otus_processed_all_wtime$gg_ExpectedDistribution_controls


##### RICHNESS (Chao1) (II) #######

if ( RERUN_RICH ) {
  
  lmer_log_chao1_all <- stan_lmer(log(chao1) ~ -1 + species + (1|indivID), data=mf_all_noinfect
                              , prior = normal(0, 10, autoscale = TRUE)
                              #, family = gaussian(link="log")
                              , seed = 5793482
                              , adapt_delta = 0.999

  )
  
  lmer_log_chao1_all_wtime <- stan_lmer(log(chao1) ~ -1 + species*time + (1|indivID), data=mf_all_noinfect
                                  , prior = normal(0, 10, autoscale = TRUE)
                                  #, family = gaussian(link="log")
                                  , seed = 5793482
                                  , adapt_delta = 0.999
                                  
  )
  save(lmer_log_chao1_all, file="./4_Bayesian_models/lmer_log_chao1_all.RData")
  save(lmer_log_chao1_all_wtime, file="./4_Bayesian_models/lmer_log_chao1_all_wtime.RData")
} else {
    load("./4_Bayesian_models/lmer_log_chao1_all.RData")
    load("./4_Bayesian_models/lmer_log_chao1_all_wtime.RData")
}
prior_summary(lmer_log_chao1_all)
prior_summary(lmer_log_chao1_all_wtime)


chao1_processed_all <- process_glmer_all(g_lmer = lmer_log_chao1_all
                                 , dep= "chao1"
                                 , name_dep="log_chao1"
                                 , transform_raw_func="log"
                                 , intercept_present=F
                                 , fit_distr="normal"
                                 , mf_con=mf_con
                                 , mf_treat=mf_treat)
chao1_processed_all$gg_model_Distribution_of_all_values
chao1_processed_all$gg_p
chao1_processed_all$gg_ExpectedDistribution_and_Bd_exposed
chao1_processed_all$gg_ExpectedDistribution_controls

chao1_processed_all_wtime <- process_glmer_all_wtime(g_lmer = lmer_log_chao1_all_wtime
                                         , dep= "chao1"
                                         , name_dep="log_chao1"
                                         , transform_raw_func="log"
                                         , intercept_present=F
                                         , fit_distr="normal"
                                         , mf_con=mf_con
                                         , mf_treat=mf_treat
                                         , time_factor = T
                                         , time_inter = T)
chao1_processed_all_wtime$gg_model_Distribution_of_all_values
chao1_processed_all_wtime$gg_p
chao1_processed_all_wtime$gg_ExpectedDistribution_and_Bd_exposed
chao1_processed_all_wtime$gg_ExpectedDistribution_controls


##### DIVERSITY (Shannon) (II) #######

if ( RERUN_RICH ) {
  
  lmer_shannon_all <- stan_lmer(shannon ~ -1 + species + (1|indivID), data=mf_all_noinfect
                            , prior = normal(0, 10, autoscale = TRUE)
                            , seed = 5793482
                            , adapt_delta = 0.999
                            
  )
  lmer_shannon_all_wtime <- stan_lmer(shannon ~ -1 + species*time + (1|indivID), data=mf_all_noinfect
                                , prior = normal(0, 10, autoscale = TRUE)
                                , seed = 5793482
                                , adapt_delta = 0.999
                                
  )
  save(lmer_shannon_all, file="./4_Bayesian_models/lmer_shannon_all.RData")
  save(lmer_shannon_all_wtime, file="./4_Bayesian_models/lmer_shannon_all_wtime.RData")
} else {
    load("./4_Bayesian_models/lmer_shannon_all.RData")
    load("./4_Bayesian_models/lmer_shannon_all_wtime.RData")
}
prior_summary(lmer_shannon_all)
prior_summary(lmer_shannon_all_wtime)


shannon_processed_all <- process_glmer_all(g_lmer = lmer_shannon_all
                                   , dep= "shannon"
                                   , name_dep="shannon"
                                   , transform_raw_func="None"
                                   , intercept_present=F
                                   , fit_distr="normal"
                                   , mf_con=mf_con
                                   , mf_treat=mf_treat)
shannon_processed_all$gg_model_Distribution_of_all_values
shannon_processed_all$gg_p
shannon_processed_all$gg_ExpectedDistribution_and_Bd_exposed
shannon_processed_all$gg_ExpectedDistribution_controls


shannon_processed_all_wtime <- process_glmer_all_wtime(g_lmer = lmer_shannon_all_wtime
                                           , dep= "shannon"
                                           , name_dep="shannon"
                                           , transform_raw_func="None"
                                           , intercept_present=F
                                           , fit_distr="normal"
                                           , mf_con=mf_con
                                           , mf_treat=mf_treat
                                           , time_factor = T
                                           , time_inter = T)
shannon_processed_all_wtime$gg_model_Distribution_of_all_values
shannon_processed_all_wtime$gg_p
shannon_processed_all_wtime$gg_ExpectedDistribution_and_Bd_exposed
shannon_processed_all_wtime$gg_ExpectedDistribution_controls

##### DIVERSITY (Faith's PD) (II) #######

if ( RERUN_RICH ) {
  
  glmer_faith_pd_all <- stan_glmer(faith_pd ~ species + (1|indivID), data=mf_all_noinfect
                               , prior = normal(0, 10, autoscale = TRUE)
                               , seed = 5793482
                               , family=Gamma(link="identity")
                               , adapt_delta = 0.999
                               
  )
  
  glmer_faith_pd_all_wtime <- stan_glmer(faith_pd ~ species*time + (1|indivID), data=mf_all_noinfect
                                   , prior = normal(0, 10, autoscale = TRUE)
                                   , seed = 5793482
                                   , family=Gamma(link="identity")
                                   , adapt_delta = 0.999
                                   
  )
  save(glmer_faith_pd_all, file="./4_Bayesian_models/glmer_faith_pd_all.RData")
  save(glmer_faith_pd_all_wtime, file="./4_Bayesian_models/glmer_faith_pd_all_wtime.RData")
} else {
    load("./4_Bayesian_models/glmer_faith_pd_all.RData")
    load("./4_Bayesian_models/glmer_faith_pd_all_wtime.RData")
}
prior_summary(glmer_faith_pd_all)
prior_summary(glmer_faith_pd_all_wtime)

faith_pd_processed_all <- process_glmer_all(g_lmer = glmer_faith_pd_all
                                    , dep= "faith_pd"
                                    , name_dep="faith_pd"
                                    , transform_raw_func="None"
                                    , intercept_present=T
                                    , fit_distr="Gamma"
                                    , mf_con = mf_con
                                    , mf_treat=mf_treat
)

faith_pd_processed_all$gg_model_Distribution_of_all_values
faith_pd_processed_all$gg_p
faith_pd_processed_all$gg_ExpectedDistribution_and_Bd_exposed
faith_pd_processed_all$gg_ExpectedDistribution_controls

faith_pd_processed_all_wtime <- process_glmer_all_wtime(g_lmer = glmer_faith_pd_all_wtime
                                            , dep= "faith_pd"
                                            , name_dep="faith_pd"
                                            , transform_raw_func="None"
                                            , intercept_present=T
                                            , fit_distr="Gamma"
                                            , mf_con = mf_con
                                            , mf_treat=mf_treat
                                            , time_factor=T
                                            , time_inter = T
)

faith_pd_processed_all_wtime$gg_model_Distribution_of_all_values
faith_pd_processed_all_wtime$gg_p
faith_pd_processed_all_wtime$gg_ExpectedDistribution_and_Bd_exposed
faith_pd_processed_all_wtime$gg_ExpectedDistribution_controls

##### DISPERSION (Bray-curtis) (II) #######
if ( RERUN_DISP ) {
  
  lmer_disper_braycurtis_all <- stan_lmer(log(disper_braycurtis) ~ -1 + species + (1|indivID) + time
                                      , data=mf_all_noinfect
                                      , prior_intercept = normal(location = 0,scale = 5, autoscale = TRUE)
                                      , prior = normal(location=0, scale=5, autoscale=TRUE)
                                      , seed= 29473
  )
  lmer_disper_braycurtis_all_wtime <- stan_lmer(log(disper_braycurtis) ~ -1 + species*time + (1|indivID)
                                          , data=mf_all_noinfect
                                          , prior_intercept = normal(location = 0,scale = 5, autoscale = TRUE)
                                          , prior = normal(location=0, scale=5, autoscale=TRUE)
                                          , seed= 29473
  )
  save(lmer_disper_braycurtis_all, file="./4_Bayesian_models/lmer_disper_braycurtis_all.RData")
  save(lmer_disper_braycurtis_all_wtime, file="./4_Bayesian_models/lmer_disper_braycurtis_all_wtime.RData")
} else {
    load("./4_Bayesian_models/lmer_disper_braycurtis_all.RData")
    load("./4_Bayesian_models/lmer_disper_braycurtis_all_wtime.RData")
}
prior_summary(lmer_disper_braycurtis_all)
prior_summary(lmer_disper_braycurtis_all_wtime)

disper_braycurtis_processed_all <- process_glmer_all(g_lmer = lmer_disper_braycurtis_all
                                             , dep= "disper_braycurtis"
                                             , name_dep="log_disper_braycurtis"
                                             , transform_raw_func="log"
                                             , intercept_present=F
                                             , time_factor = T
                                             , fit_distr="normal"
                                             , mf_con=mf_con
                                             , mf_treat=mf_treat
)
disper_braycurtis_processed_all$gg_model_Distribution_of_all_values
disper_braycurtis_processed_all$gg_p
disper_braycurtis_processed_all$gg_ExpectedDistribution_and_Bd_exposed
disper_braycurtis_processed_all$gg_ExpectedDistribution_controls

disper_braycurtis_processed_all_wtime <- process_glmer_all_wtime(g_lmer = lmer_disper_braycurtis_all_wtime
                                                     , dep= "disper_braycurtis"
                                                     , name_dep="log_disper_braycurtis"
                                                     , transform_raw_func="log"
                                                     , intercept_present=F
                                                     , fit_distr="normal"
                                                     , mf_con=mf_con
                                                     , mf_treat=mf_treat
                                                     , time_inter = T
                                                     , time_factor = T
)
disper_braycurtis_processed_all_wtime$gg_model_Distribution_of_all_values
disper_braycurtis_processed_all_wtime$gg_p
disper_braycurtis_processed_all_wtime$gg_ExpectedDistribution_and_Bd_exposed
disper_braycurtis_processed_all_wtime$gg_ExpectedDistribution_controls


##### DISPERSION (unweighted Unifrac) (II) #######
if ( RERUN_DISP ) {
  
  lmer_disper_unweighted_unifrac_all <- stan_lmer(log(disper_unweighted_unifrac) ~ -1 + species + (1|indivID) + time
                                              , data=mf_all_noinfect
                                              , prior_intercept = normal(location = 0,scale = 5, autoscale = TRUE)
                                              , prior = normal(location=0, scale=5, autoscale=TRUE)
                                              , seed= 29473
  )
  
  lmer_disper_unweighted_unifrac_all_wtime <- stan_lmer(log(disper_unweighted_unifrac) ~ -1 + species*time + (1|indivID)
                                                  , data=mf_all_noinfect
                                                  , prior_intercept = normal(location = 0,scale = 5, autoscale = TRUE)
                                                  , prior = normal(location=0, scale=5, autoscale=TRUE)
                                                  , seed= 29473
  )
  
  save(lmer_disper_unweighted_unifrac_all, file="./4_Bayesian_models/lmer_disper_unweighted_unifrac_all.RData")
  save(lmer_disper_unweighted_unifrac_all_wtime, file="./4_Bayesian_models/lmer_disper_unweighted_unifrac_all_wtime.RData")
} else {
    load("./4_Bayesian_models/lmer_disper_unweighted_unifrac_all.RData")
    load("./4_Bayesian_models/lmer_disper_unweighted_unifrac_all_wtime.RData")
}
prior_summary(lmer_disper_unweighted_unifrac_all)
prior_summary(lmer_disper_unweighted_unifrac_all_wtime)

disper_unweighted_unifrac_processed_all <- process_glmer_all(g_lmer = lmer_disper_unweighted_unifrac_all
                                                     , dep= "disper_unweighted_unifrac"
                                                     , name_dep="log_disper_unweighted_unifrac"
                                                     , transform_raw_func="log"
                                                     , intercept_present=F
                                                     , time_factor = T
                                                     , fit_distr="normal"
                                                     , mf_con=mf_con
                                                     , mf_treat=mf_treat
)
disper_unweighted_unifrac_processed_all$gg_model_Distribution_of_all_values
disper_unweighted_unifrac_processed_all$gg_p
disper_unweighted_unifrac_processed_all$gg_ExpectedDistribution_and_Bd_exposed
disper_unweighted_unifrac_processed_all$gg_ExpectedDistribution_controls


disper_unweighted_unifrac_processed_all_wtime <- process_glmer_all_wtime(g_lmer = lmer_disper_unweighted_unifrac_all_wtime
                                                             , dep= "disper_unweighted_unifrac"
                                                             , name_dep="log_disper_unweighted_unifrac"
                                                             , transform_raw_func="log"
                                                             , intercept_present=F
                                                             , fit_distr="normal"
                                                             , mf_con=mf_con
                                                             , mf_treat=mf_treat
                                                             , time_factor = T
                                                             , time_inter = T
)
disper_unweighted_unifrac_processed_all_wtime$gg_model_Distribution_of_all_values
disper_unweighted_unifrac_processed_all_wtime$gg_p
disper_unweighted_unifrac_processed_all_wtime$gg_ExpectedDistribution_and_Bd_exposed
disper_unweighted_unifrac_processed_all_wtime$gg_ExpectedDistribution_controls



##### DISPERSION (weighted Unifrac) (II) #######
if ( RERUN_DISP ) {
  
  lmer_disper_weighted_unifrac_all <- stan_lmer(log(disper_weighted_unifrac) ~ -1 + species + (1|indivID) + time
                                            , data=mf_all_noinfect
                                            , prior_intercept = normal(location = 0,scale = 5, autoscale = TRUE)
                                            , prior = normal(location=0, scale=5, autoscale=TRUE)
                                            , seed= 29473
  )
  
  lmer_disper_weighted_unifrac_all_wtime <- stan_lmer(log(disper_weighted_unifrac) ~ -1 + species*time + (1|indivID) 
                                                , data=mf_all_noinfect
                                                , prior_intercept = normal(location = 0,scale = 5, autoscale = TRUE)
                                                , prior = normal(location=0, scale=5, autoscale=TRUE)
                                                , seed= 29473
  )
  
  save(lmer_disper_weighted_unifrac_all, file="./4_Bayesian_models/lmer_disper_weighted_unifrac_all.RData")
  save(lmer_disper_weighted_unifrac_all_wtime, file="./4_Bayesian_models/lmer_disper_weighted_unifrac_all_wtime.RData")
} else {
    load("./4_Bayesian_models/lmer_disper_weighted_unifrac_all.RData")
    load("./4_Bayesian_models/lmer_disper_weighted_unifrac_all_wtime.RData")
}
prior_summary(lmer_disper_weighted_unifrac_all)
prior_summary(lmer_disper_weighted_unifrac_all_wtime)

disper_weighted_unifrac_processed_all <- process_glmer_all(g_lmer = lmer_disper_weighted_unifrac_all
                                                   , dep= "disper_weighted_unifrac"
                                                   , name_dep="log_disper_weighted_unifrac"
                                                   , transform_raw_func="log"
                                                   , intercept_present=F
                                                   , time_factor = T
                                                   , fit_distr="normal"
                                                   , mf_con=mf_con
                                                   , mf_treat=mf_treat
)
disper_weighted_unifrac_processed_all$gg_model_Distribution_of_all_values
disper_weighted_unifrac_processed_all$gg_p
disper_weighted_unifrac_processed_all$gg_ExpectedDistribution_and_Bd_exposed
disper_weighted_unifrac_processed_all$gg_ExpectedDistribution_controls

disper_weighted_unifrac_processed_all_wtime <- process_glmer_all_wtime(g_lmer = lmer_disper_weighted_unifrac_all_wtime
                                                           , dep= "disper_weighted_unifrac"
                                                           , name_dep="log_disper_weighted_unifrac"
                                                           , transform_raw_func="log"
                                                           , intercept_present=F
                                                           , fit_distr="normal"
                                                           , mf_con=mf_con
                                                           , mf_treat=mf_treat
                                                           , time_factor = T
                                                           , time_inter = T
)
disper_weighted_unifrac_processed_all_wtime$gg_model_Distribution_of_all_values
disper_weighted_unifrac_processed_all_wtime$gg_p
disper_weighted_unifrac_processed_all_wtime$gg_ExpectedDistribution_and_Bd_exposed
disper_weighted_unifrac_processed_all_wtime$gg_ExpectedDistribution_controls



##### DISTANCE (Bray-curtis) (II) #######
if ( RERUN_DIST ) {
  
  glmer_dist_braycurtis_all <- stan_glmer(dist_braycurtis ~ -1 + species + (1|indivID)
                                      , data=mf_all_noinfect
                                      , family =mgcv::betar
                                      , prior_intercept = normal(location = 0.5,scale = 2.5, autoscale = TRUE)
                                      , prior = normal(location=0.5, scale=2.5, autoscale=TRUE)
                                      , seed= 623445
  )
  
  glmer_dist_braycurtis_all_wtime <- stan_glmer(dist_braycurtis ~ -1 + species*time + (1|indivID)
                                          , data=mf_all_noinfect
                                          , family =mgcv::betar
                                          , prior_intercept = normal(location = 0.5,scale = 2.5, autoscale = TRUE)
                                          , prior = normal(location=0.5, scale=2.5, autoscale=TRUE)
                                          , seed= 623445
  )
  save(glmer_dist_braycurtis_all, file="./4_Bayesian_models/glmer_dist_braycurtis_all.RData")
  save(glmer_dist_braycurtis_all_wtime, file="./4_Bayesian_models/glmer_dist_braycurtis_all_wtime.RData")
} else {
    load("./4_Bayesian_models/glmer_dist_braycurtis_all.RData")
    load("./4_Bayesian_models/glmer_dist_braycurtis_all_wtime.RData")
}
prior_summary(glmer_dist_braycurtis_all)
prior_summary(glmer_dist_braycurtis_all_wtime)


dist_braycurtis_processed_all <- process_glmer_all(g_lmer = glmer_dist_braycurtis_all
                                           , dep= "dist_braycurtis"
                                           , name_dep="dist_braycurtis"
                                           , transform_raw_func="beta"
                                           , intercept_present=F
                                           , fit_distr="Beta"
                                           , mf_con=mf_con
                                           , mf_treat=mf_treat
)

dist_braycurtis_processed_all$gg_model_Distribution_of_all_values
dist_braycurtis_processed_all$gg_p
dist_braycurtis_processed_all$gg_ExpectedDistribution_and_Bd_exposed
dist_braycurtis_processed_all$gg_ExpectedDistribution_controls

dist_braycurtis_processed_all_wtime <- process_glmer_all_wtime(g_lmer = glmer_dist_braycurtis_all_wtime
                                                   , dep= "dist_braycurtis"
                                                   , name_dep="dist_braycurtis"
                                                   , transform_raw_func="beta"
                                                   , intercept_present=F
                                                   , fit_distr="Beta"
                                                   , mf_con=mf_con
                                                   , mf_treat=mf_treat
                                                   , time_factor = T
                                                   , time_inter = T
)

dist_braycurtis_processed_all_wtime$gg_model_Distribution_of_all_values
dist_braycurtis_processed_all_wtime$gg_p
dist_braycurtis_processed_all_wtime$gg_ExpectedDistribution_and_Bd_exposed
dist_braycurtis_processed_all_wtime$gg_ExpectedDistribution_controls


##### DISTANCE (unweighted Unifrac) (II) #######
if ( RERUN_DIST ) {
  
  glmer_dist_unweighted_unifrac_all <- stan_glmer(dist_unweighted_unifrac ~ -1 + species + (1|indivID)
                                              , data=mf_all_noinfect
                                              , family =mgcv::betar
                                              , prior_intercept = normal(location = 0.5,scale = 2.5, autoscale = TRUE)
                                              , prior = normal(location=0.5, scale=2.5, autoscale=TRUE)
                                              , seed= 623445
  )
  
  
  glmer_dist_unweighted_unifrac_all_wtime <- stan_glmer(dist_unweighted_unifrac ~ -1 + species*time + (1|indivID)
                                                  , data=mf_all_noinfect
                                                  , family =mgcv::betar
                                                  , prior_intercept = normal(location = 0.5,scale = 2.5, autoscale = TRUE)
                                                  , prior = normal(location=0.5, scale=2.5, autoscale=TRUE)
                                                  , seed= 623445
  )
  save(glmer_dist_unweighted_unifrac_all, file="./4_Bayesian_models/glmer_dist_unweighted_unifrac_all.RData")
  save(glmer_dist_unweighted_unifrac_all_wtime, file="./4_Bayesian_models/glmer_dist_unweighted_unifrac_all_wtime.RData")
} else {
  load("./4_Bayesian_models/glmer_dist_unweighted_unifrac_all.RData")
    load("./4_Bayesian_models/glmer_dist_unweighted_unifrac_all_wtime.RData")
    
}
prior_summary(glmer_dist_unweighted_unifrac_all)
prior_summary(glmer_dist_unweighted_unifrac_all_wtime)

dist_unweighted_unifrac_processed_all <- process_glmer_all(g_lmer = glmer_dist_unweighted_unifrac_all
                                                   , dep= "dist_unweighted_unifrac"
                                                   , name_dep="dist_unweighted_unifrac"
                                                   , transform_raw_func="beta"
                                                   , intercept_present=F
                                                   , fit_distr="Beta"
                                                   , mf_con=mf_con
                                                   , mf_treat=mf_treat
)

dist_unweighted_unifrac_processed_all$gg_model_Distribution_of_all_values
dist_unweighted_unifrac_processed_all$gg_p
dist_unweighted_unifrac_processed_all$gg_ExpectedDistribution_and_Bd_exposed
dist_unweighted_unifrac_processed_all$gg_ExpectedDistribution_controls

dist_unweighted_unifrac_processed_all_wtime <- process_glmer_all_wtime(g_lmer = glmer_dist_unweighted_unifrac_all_wtime
                                                           , dep= "dist_unweighted_unifrac"
                                                           , name_dep="dist_unweighted_unifrac"
                                                           , transform_raw_func="beta"
                                                           , intercept_present=F
                                                           , fit_distr="Beta"
                                                           , mf_con=mf_con
                                                           , mf_treat=mf_treat
                                                           , time_factor = T
                                                           , time_inter = T
)

dist_unweighted_unifrac_processed_all_wtime$gg_model_Distribution_of_all_values
dist_unweighted_unifrac_processed_all_wtime$gg_p
dist_unweighted_unifrac_processed_all_wtime$gg_ExpectedDistribution_and_Bd_exposed
dist_unweighted_unifrac_processed_all_wtime$gg_ExpectedDistribution_controls

##### DISTANCE (weighted Unifrac) (II) #######
if ( RERUN_DIST ) {
  
  glmer_dist_weighted_unifrac_all <- stan_glmer(dist_weighted_unifrac ~ -1 + species + (1|indivID)
                                            , data=mf_all_noinfect
                                            , family =mgcv::betar
                                            , prior_intercept = normal(location = 0.5,scale = 2.5, autoscale = TRUE)
                                            , prior = normal(location=0.5, scale=2.5, autoscale=TRUE)
                                            , seed= 623445
  )
  
  
  glmer_dist_weighted_unifrac_all_wtime <- stan_glmer(dist_weighted_unifrac ~ -1 + species*time + (1|indivID)
                                                , data=mf_all_noinfect
                                                , family =mgcv::betar
                                                , prior_intercept = normal(location = 0.5,scale = 2.5, autoscale = TRUE)
                                                , prior = normal(location=0.5, scale=2.5, autoscale=TRUE)
                                                , seed= 623445
  )
  save(glmer_dist_weighted_unifrac_all, file="./4_Bayesian_models/glmer_dist_weighted_unifrac_all.RData")
  save(glmer_dist_weighted_unifrac_all_wtime, file="./4_Bayesian_models/glmer_dist_weighted_unifrac_all_wtime.RData")
} else {
    load("./4_Bayesian_models/glmer_dist_weighted_unifrac_all.RData")
    load("./4_Bayesian_models/glmer_dist_weighted_unifrac_all_wtime.RData")
}
prior_summary(glmer_dist_weighted_unifrac_all)
prior_summary(glmer_dist_weighted_unifrac_all_wtime)

dist_weighted_unifrac_processed_all <- process_glmer_all(g_lmer = glmer_dist_weighted_unifrac_all
                                                 , dep= "dist_weighted_unifrac"
                                                 , name_dep="dist_weighted_unifrac"
                                                 , transform_raw_func="beta"
                                                 , intercept_present=F
                                                 , fit_distr="Beta"
                                                 , mf_con=mf_con
                                                 , mf_treat=mf_treat
)

dist_weighted_unifrac_processed_all$gg_model_Distribution_of_all_values
dist_weighted_unifrac_processed_all$gg_p
dist_weighted_unifrac_processed_all$gg_ExpectedDistribution_and_Bd_exposed
dist_weighted_unifrac_processed_all$gg_ExpectedDistribution_controls


dist_weighted_unifrac_processed_all_wtime <- process_glmer_all_wtime(g_lmer = glmer_dist_weighted_unifrac_all_wtime
                                                         , dep= "dist_weighted_unifrac"
                                                         , name_dep="dist_weighted_unifrac"
                                                         , transform_raw_func="beta"
                                                         , intercept_present=F
                                                         , fit_distr="Beta"
                                                         , mf_con=mf_con
                                                         , mf_treat=mf_treat
                                                         , time_factor = T
                                                         , time_inter = T
)

dist_weighted_unifrac_processed_all_wtime$gg_model_Distribution_of_all_values
dist_weighted_unifrac_processed_all_wtime$gg_p
dist_weighted_unifrac_processed_all_wtime$gg_ExpectedDistribution_and_Bd_exposed
dist_weighted_unifrac_processed_all_wtime$gg_ExpectedDistribution_controls

##### PERCENT INHIBITORY (II) #######
if ( RERUN_PERCINHIB ) {
  
  glmer_percInhib_all <- stan_glmer(percInhib ~ -1 + species + (1|indivID)
                                , data=mf_all_noinfect
                                , family =mgcv::betar
                                , prior_intercept = normal(location = 0.5,scale = 2.5, autoscale = TRUE)
                                , prior = normal(location=0.5, scale=2.5, autoscale=TRUE)
                                , seed= 59283
  )
  
  
  glmer_percInhib_all_wtime <- stan_glmer(percInhib ~ -1 + species*time + (1|indivID)
                                    , data=mf_all_noinfect
                                    , family =mgcv::betar
                                    , prior_intercept = normal(location = 0.5,scale = 2.5, autoscale = TRUE)
                                    , prior = normal(location=0.5, scale=2.5, autoscale=TRUE)
                                    , seed= 59283
  )
  save(glmer_percInhib_all, file="./4_Bayesian_models/glmer_percInhib_all.RData")
  save(glmer_percInhib_all_wtime, file="./4_Bayesian_models/glmer_percInhib_all_wtime.RData")
} else {
    load("./4_Bayesian_models/glmer_percInhib_all.RData")
    load("./4_Bayesian_models/glmer_percInhib_all_wtime.RData")
    
}
prior_summary(glmer_percInhib_all)
prior_summary(glmer_percInhib_all_wtime)

percInhib_processed_all <- process_glmer_all(g_lmer = glmer_percInhib_all
                                     , dep= "percInhib"
                                     , name_dep="percInhib"
                                     , transform_raw_func="beta"
                                     , intercept_present=F
                                     , fit_distr="Beta"
                                     , mf_con=mf_con
                                     , mf_treat=mf_treat
)

percInhib_processed_all$gg_model_Distribution_of_all_values
percInhib_processed_all$gg_p
percInhib_processed_all$gg_ExpectedDistribution_and_Bd_exposed
percInhib_processed_all$gg_ExpectedDistribution_controls


percInhib_processed_all_wtime <- process_glmer_all_wtime(g_lmer = glmer_percInhib_all_wtime
                                             , dep= "percInhib"
                                             , name_dep="percInhib"
                                             , transform_raw_func="beta"
                                             , intercept_present=F
                                             , fit_distr="Beta"
                                             , mf_con=mf_con
                                             , mf_treat=mf_treat
                                             , time_factor = T
                                             , time_inter = T
)

percInhib_processed_all_wtime$gg_model_Distribution_of_all_values
percInhib_processed_all_wtime$gg_p
percInhib_processed_all_wtime$gg_ExpectedDistribution_and_Bd_exposed
percInhib_processed_all_wtime$gg_ExpectedDistribution_controls

##### INHIBITORY RICHNESS (II) #######
if ( RERUN_INHIBRICH ) {
  # Here, I include a random variable of sample because of over-dispersion from poisson
  glmer_inhibRich_all <- stan_glmer(inhibRich ~ species + (1|indivID) + (1|SampleID) + time, data=mf_all_noinfect
                                , prior = normal(0, 10, autoscale = TRUE)
                                , family= poisson(link="identity")
                                , seed = 5423409)
  
  glmer_inhibRich_all_wtime <- stan_glmer(inhibRich ~ species*time + (1|indivID) + (1|SampleID), data=mf_all_noinfect
                                    , prior = normal(0, 10, autoscale = TRUE)
                                    , family= poisson(link="identity")
                                    , seed = 5423409)
  
  save(glmer_inhibRich_all, file="./4_Bayesian_models/glmer_inhibRich_all.RData")
  save(glmer_inhibRich_all_wtime, file="./4_Bayesian_models/glmer_inhibRich_all_wtime.RData")
} else {
    load("./4_Bayesian_models/glmer_inhibRich_all.RData")
    load("./4_Bayesian_models/glmer_inhibRich_all_wtime.RData")
}
prior_summary(glmer_inhibRich_all)
prior_summary(glmer_inhibRich_all_wtime)

inhibRich_processed_all <- process_glmer_all(g_lmer = glmer_inhibRich_all
                                     , dep= "inhibRich"
                                     , name_dep="inhibRich"
                                     , transform_raw_func="None"
                                     , intercept_present=T
                                     , fit_distr="Poisson"
                                     , mf_con=mf_con
                                     , mf_treat=mf_treat
                                     , time_factor=T
)

inhibRich_processed_all$gg_model_Distribution_of_all_values
inhibRich_processed_all$gg_p
inhibRich_processed_all$gg_ExpectedDistribution_and_Bd_exposed
inhibRich_processed_all$gg_ExpectedDistribution_controls


inhibRich_processed_all_wtime <- process_glmer_all_wtime(g_lmer = glmer_inhibRich_all_wtime
                                             , dep= "inhibRich"
                                             , name_dep="inhibRich"
                                             , transform_raw_func="None"
                                             , intercept_present=T
                                             , fit_distr="Poisson"
                                             , mf_con=mf_con
                                             , mf_treat=mf_treat
                                             , time_factor=T
                                             , time_inter = T
)

inhibRich_processed_all_wtime$gg_model_Distribution_of_all_values
inhibRich_processed_all_wtime$gg_p
inhibRich_processed_all_wtime$gg_ExpectedDistribution_and_Bd_exposed
inhibRich_processed_all_wtime$gg_ExpectedDistribution_controls


# Combine all p values and save work
all_tabs <- c("observed_otus"
              , "chao1"
              , "shannon"
              , "faith_pd"
              , "disper_braycurtis"
              , "disper_unweighted_unifrac"
              , "disper_weighted_unifrac"
              , "dist_braycurtis"
              , "dist_unweighted_unifrac"
              , "dist_weighted_unifrac"
              , "percInhib"
              , "inhibRich")
# Change names of columns so they're different
for ( tab in all_tabs) {
  temp_tab <- get(paste0(tab,"_processed_all_wtime"))$all_p 
  name1 <- paste0("exp_",tab)
  name2 <- paste0("p_",tab)
  colnames(temp_tab) <- c("indivID",name1, name2,"max_Bd_load")
  assign(paste0("all_p_pred_",tab), temp_tab)
}


# combine all_p's
all_p_pred <- data.frame(indivID=observed_otus_processed_all_wtime$all_p$indivID, infect=observed_otus_processed_all_wtime$all_p$infect) 
for ( tab in all_tabs ) {
  all_p_pred <- get(paste0("all_p_pred_",tab)) %>%
    dplyr::select(indivID,paste0("exp_",tab), paste0("p_",tab)) %>%
    left_join(all_p_pred)
}

save(all_p_pred, file="./4_Bayesian_models/all_p_pred.RData")

write.table(all_p, file="./4_Bayesian_models/all_p.txt", quote=FALSE, row.names = FALSE
            , sep="\t")



