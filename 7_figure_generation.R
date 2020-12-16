library(ggplot2)
library(MASS) # for fit distr
library(tidyverse)
library(gridExtra)
library(RColorBrewer) # colors for barplots
# setwd("~/Documents/PhD/Project_5Species/5Sp_analysis_Nov2019")
### Creating figures for 5 species manuscript ####

# Make directory for figures
dir.create("./7_figures")

# Mapping file
load("./3_5sp_mapping_otu_downstream/mf_alt_filt_final.RData")
mf_con <- mf_alt_filt_final %>%
  filter(Bd_exposure == "Control")
mf_treat <- mf_alt_filt_final %>%
  filter(Bd_exposure == "Bd-exposed")
load("./3_5sp_mapping_otu_downstream/otu_filt.RData")
load("./3_5sp_mapping_otu_downstream/otu_filt_inhibOnly.RData")
load("./3_5sp_mapping_otu_downstream/taxonomy.RData")

load("./4_Bayesian_models/all_p_con.RData")
load("./4_Bayesian_models/all_p.RData")
load("./4_Bayesian_models/all_p_combined.RData")

load("./4_Bayesian_models/all_p_pred.RData")
# load("./5_random_forest/all_RF_predictBD.RData")

# Updated random forest LOO results
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
  separate(indivID, into=c("species","indiv"), remove=FALSE)
all_p_con <- all_p_con %>%
  mutate(PABD=NA, Bd_load=NA) %>%
  separate(indivID, into=c("species","indiv"), remove=FALSE)
all_p_pred <- all_p_pred %>%
  mutate(PABD=ifelse(infect>0,1,0)) %>%
  rename(Bd_load=infect) %>%
  separate(indivID, into=c("species","indiv"), remove=FALSE)

#### experimental design ####
plot_experimental_design <-mf_alt_filt_final %>%
  mutate(Contaminated = factor(ifelse(orig_contam ==1, "!Contaminated",NA), levels=c("!Contaminated"))
         , Bd_logload = (Bd_load)) %>%
  ggplot(aes(x=time, y=indiv)) +
  geom_line(aes(group=indivID, col=Bd_exposure)) +
  geom_point(aes(group=indivID, bg=Bd_logload, col=Bd_exposure), cex=4, pch=21)+
  scale_color_manual(values=c("black","blue","orange")) +
  scale_fill_gradient(low = "white", high = "red") +
  geom_vline(aes(xintercept=5.5), col="orange")+
  geom_point(aes(group=indivID, col=Contaminated), cex=1, pch=19)+ ## NEW LINE
  facet_wrap(~species, nrow=5) +
  xlab("Time (Day)") +
  ylab("Individual Amphibian") +
  theme(strip.text = element_text(size=14), 
        axis.title = element_text(size=14))+
  labs(fill=expression(paste(italic("Bd")," load (log+1)"))
       , col = expression(paste(italic("Bd"), " exposure")))

# plot_experimental_design <- mf_alt_filt_final %>%
#   mutate(Contaminated = factor(ifelse(orig_contam ==1, "!Contaminated",NA), levels=c("!Contaminated"))
#          , Bd_logload = (Bd_load)) %>%
#   ggplot(aes(x=time, y=indiv)) +
#   geom_line(aes(group=indivID, col=Bd_exposure)) +
#   geom_point(aes(group=indivID, bg=Bd_logload), cex=4, pch=21)+
#   scale_color_manual(values=c("black","blue","orange")) +
#   scale_fill_gradient(low = "white", high = "red") +
#   geom_vline(aes(xintercept=5.5), col="orange")+
#   geom_point(aes(group=indivID, col=Contaminated), cex=1, pch=19)+ ## NEW LINE
#   facet_wrap(~species, nrow=5) +
#   xlab("Time (Day)") +
#   ylab("Individual Amphibian") +
#   theme(strip.text = element_text(size=14), 
#         axis.title = element_text(size=14))
#plot_experimental_design

pdf(file = "./7_figures/experimental_design.pdf", width = 4.5,height = 10.5)
plot_experimental_design
dev.off()


#### Control data: how does it vary? ####
gg_NMDS <- mf_con %>%
  ggplot(aes(x=NMDS1_braycurtis, y=NMDS2_braycurtis)) +
  geom_point(aes(col=species, alpha=time), cex=3, show.legend = FALSE) +
  ylab("NMDS 2 (Bray-curtis)") +
  xlab("NMDS 1 (Bray-curtis)") +
  geom_text(aes(x=-0.75, y=0.75, label=paste0("Stress: ",round(unique(NMDS_stress_braycurtis))/100)))

# gg_NMDS <- mf_con %>%
#   ggplot(aes(x=NMDS1_unweighted_unifrac, y=NMDS2_unweighted_unifrac)) +
#   geom_point(aes(col=species, alpha=time), cex=3, show.legend = FALSE) +
#   ylab("NMDS 2 (Unweighted Unifrac)") +
#   xlab("NMDS 1 (Unweighted Unifrac)") +
#   geom_text(aes(x=-0.75, y=0.75, label=paste0("Stress: ",round(unique(NMDS_stress_unweighted_unifrac))/100)))

gg_infect <- mf_treat  %>%
  group_by(species, indivID) %>%
  summarize("Max Bd load (log)"=max(Bd_load)) %>%
  ggplot(aes(x=species, y=`Max Bd load (log)`)) +
  geom_point(aes(col=species), cex=3, position = position_jitter(width=0.1, height=0.05), show.legend = FALSE) +
  xlab("Species") + ylab(expression(paste("Max ",italic("Bd")," load (log+1)")))

temp1a <-  mf_con %>%
  dplyr::select(species, shannon) %>%
  mutate(metric="Shannon_diversity") %>%
  rename(value=shannon)
temp1b <-  mf_con %>%
  dplyr::select(species, observed_otus) %>%
  mutate(metric="Bacterial_richness") %>%
  rename(value=observed_otus)
temp2 <- mf_con %>%
  dplyr::select(species, inhibRich) %>%
  mutate(metric="Richness_of_putative_inhibitory_bacteria")%>%
  rename(value=inhibRich)
temp3 <- mf_con %>%
  dplyr::select(species, percInhib) %>%
  mutate(metric="Proportion_of_putative_inhibitory_bacteria")%>%
  rename(value=percInhib)
temp4 <- mf_con %>%
  dplyr::select(species, disper_unweighted_unifrac) %>%
  mutate(metric="Dispersion (log_unweighted_Unifrac)")%>%
  rename(value=disper_unweighted_unifrac)
temp5 <- mf_con %>%
  dplyr::select(species, dist_unweighted_unifrac) %>%
  mutate(metric="Instability (unweighted_Unifrac)")%>%
  rename(value=dist_unweighted_unifrac)


gg_all <- rbind(temp1a,temp1b,temp2,temp3,temp4, temp5) %>%
  rename(Species=species) %>%
  mutate(Metric = gsub("_"," ",metric, fixed=TRUE)) %>%
  mutate(Metric = factor(Metric, levels=c("Bacterial richness","Shannon diversity","Richness of putative inhibitory bacteria","Proportion of putative inhibitory bacteria","Dispersion (log unweighted Unifrac)", "Instability (unweighted Unifrac)"))) %>%
  ggplot(aes(x=Species, y=value)) +
  geom_boxplot() +
  geom_point(aes(col=Species), position = position_jitter(width=0.1, height=0), alpha=1/3)+
  facet_wrap(Metric~., scales = "free",nrow=3) +
  ylab("")+
  xlab("Species") 
lay <- rbind(c(1,2,2),
             c(3,2,2))

# options(repr.plot.height=8, repr.plot.width=12)
plot_microbial_overview <- grid.arrange(gg_NMDS, gg_all, gg_infect, layout_matrix = lay)
plot_microbial_overview
ggsave(filename = "./7_figures/data_summary_controls.pdf",plot_microbial_overview
       , height=8, width=12)

#### Inihibitory Richness affects chances of infection ####
glm_PABD_inhibRich <- glm(PABD ~ p_inhibRich, family = binomial(link="logit"), data=all_p)

# x.predict <- data.frame(p_inhibRich=rep(seq(0,1, length.out = 100),times=5), species=rep(c("Anbo","Rhma","Lica","Lipi","Osse"), each=100))
x.predict <- data.frame(p_inhibRich=seq(0,1, length.out = 100))
y.predict <- predict(glm_PABD_inhibRich, newdata = x.predict, type = "link", se.fit = TRUE)

# Combine the hypothetical data and predicted values
new.data <- cbind(x.predict, y.predict)

# Calculate confidence intervals
std <- qnorm(0.95 / 2 + 0.5)
new.data$ymin <- glm_PABD_inhibRich$family$linkinv(new.data$fit - std * new.data$se)
new.data$ymax <- glm_PABD_inhibRich$family$linkinv(new.data$fit + std * new.data$se)
new.data$fit <- glm_PABD_inhibRich$family$linkinv(new.data$fit)  # Rescale to 0-1
plot_PABD_inhibRich <- all_p %>%
  rename(Species=species) %>%
  ggplot() +
  geom_point(aes(x=p_inhibRich, y=PABD,col=Species), position=position_jitter(height=0.05, width=0), cex=3) +
  geom_ribbon(data=new.data, aes(x=p_inhibRich, ymin=ymin, ymax=ymax), alpha=0.2, col="lightgrey") +
  geom_line(data=new.data, aes(x=p_inhibRich, y=fit)) +
  ylab("Probability of infection") + 
  xlab("Richness of inhibitory bacteria (Percentile within amph. species)") +
  theme_classic()

pdf("./7_figures/EFFECT_inhibRich_PABD.pdf",height=4, width=6)
plot_PABD_inhibRich
dev.off()

#### Percent inhibitory predicts Bd load ####
plot_PABD_percInhib <- all_p %>%
  rename(Species=species) %>%
  ggplot(aes(x=p_percInhib, y=Bd_load))+
  geom_point(aes(col=Species), position=position_jitter(height=0, width=0), cex=3) +
  # geom_ribbon(data=new.data, aes(x=p_inhibRich, ymin=ymin, ymax=ymax), alpha=0.2, col="lightgrey") +
  # geom_line(data=new.data, aes(x=p_inhibRich, y=fit)) +
  geom_smooth(method="lm", col="black")+
  ylab("Infection intensity (log Bd load)") + 
  xlab("Proportion of inhibitory bacteria (Percentile within amph. species)") +
  theme_classic()
# plot_PABD_percInhib
pdf("./7_figures/EFFECT_percRich_Bdload.pdf",height=4, width=6)
plot_PABD_percInhib
dev.off()



#### Dispersion affects chances of infection ####

#### UWU ####
glm_PABD_disper_unweighted_unifrac <- glm(PABD ~ p_disper_unweighted_unifrac, family = binomial(link="logit"), data=all_p)

# x.predict <- data.frame(p_inhibRich=rep(seq(0,1, length.out = 100),times=5), species=rep(c("Anbo","Rhma","Lica","Lipi","Osse"), each=100))
x.predict <- data.frame(p_disper_unweighted_unifrac=seq(0,1, length.out = 100))
y.predict <- predict(glm_PABD_disper_unweighted_unifrac, newdata = x.predict, type = "link", se.fit = TRUE)

# Combine the hypothetical data and predicted values
new.data <- cbind(x.predict, y.predict)

# Calculate confidence intervals
std <- qnorm(0.95 / 2 + 0.5)
new.data$ymin <- glm_PABD_disper_unweighted_unifrac$family$linkinv(new.data$fit - std * new.data$se)
new.data$ymax <- glm_PABD_disper_unweighted_unifrac$family$linkinv(new.data$fit + std * new.data$se)
new.data$fit <- glm_PABD_disper_unweighted_unifrac$family$linkinv(new.data$fit)  # Rescale to 0-1
plot_PABD_disper <- all_p %>%
  rename(Species=species) %>%
  ggplot() +
  geom_point(aes(x=p_disper_unweighted_unifrac, y=PABD,col=Species), position=position_jitter(height=0.05, width=0), cex=3) +
  geom_ribbon(data=new.data, aes(x=p_disper_unweighted_unifrac, ymin=ymin, ymax=ymax), alpha=0.2, col="lightgrey") +
  geom_line(data=new.data, aes(x=p_disper_unweighted_unifrac, y=fit)) +
  ylab("Probability of infection") + 
  xlab("Dispersion (Unweighted Unifrac,(Percentile within amph. species)") +
  theme_classic()
# plot_PABD_disper

pdf("./7_figures/EFFECT_disper_PABD_UWU.pdf",height=4, width=6)
plot_PABD_disper
dev.off()

#### BC ####
glm_PABD_disper_braycurtis <- glm(PABD ~ p_disper_braycurtis, family = binomial(link="logit"), data=all_p)

# x.predict <- data.frame(p_inhibRich=rep(seq(0,1, length.out = 100),times=5), species=rep(c("Anbo","Rhma","Lica","Lipi","Osse"), each=100))
x.predict <- data.frame(p_disper_braycurtis=seq(0,1, length.out = 100))
y.predict <- predict(glm_PABD_disper_braycurtis, newdata = x.predict, type = "link", se.fit = TRUE)

# Combine the hypothetical data and predicted values
new.data <- cbind(x.predict, y.predict)

# Calculate confidence intervals
std <- qnorm(0.95 / 2 + 0.5)
new.data$ymin <- glm_PABD_disper_braycurtis$family$linkinv(new.data$fit - std * new.data$se)
new.data$ymax <- glm_PABD_disper_braycurtis$family$linkinv(new.data$fit + std * new.data$se)
new.data$fit <- glm_PABD_disper_braycurtis$family$linkinv(new.data$fit)  # Rescale to 0-1
plot_PABD_disper <- all_p %>%
  rename(Species=species) %>%
  ggplot() +
  geom_point(aes(x=p_disper_braycurtis, y=PABD,col=Species), position=position_jitter(height=0.05, width=0), cex=3) +
  geom_ribbon(data=new.data, aes(x=p_disper_braycurtis, ymin=ymin, ymax=ymax), alpha=0.2, col="lightgrey") +
  geom_line(data=new.data, aes(x=p_disper_braycurtis, y=fit)) +
  ylab("Probability of infection") + 
  xlab("Dispersion (Bray-curtis) (Percentile within amph. species)") +
  theme_classic()

pdf("./7_figures/EFFECT_disper_PABD_BC.pdf",height=4, width=6)
plot_PABD_disper
dev.off()


#### WU ####
glm_PABD_disper_weighted_unifrac <- glm(PABD ~ p_disper_weighted_unifrac, family = binomial(link="logit"), data=all_p)

# x.predict <- data.frame(p_inhibRich=rep(seq(0,1, length.out = 100),times=5), species=rep(c("Anbo","Rhma","Lica","Lipi","Osse"), each=100))
x.predict <- data.frame(p_disper_weighted_unifrac=seq(0,1, length.out = 100))
y.predict <- predict(glm_PABD_disper_weighted_unifrac, newdata = x.predict, type = "link", se.fit = TRUE)

# Combine the hypothetical data and predicted values
new.data <- cbind(x.predict, y.predict)

# Calculate confidence intervals
std <- qnorm(0.95 / 2 + 0.5)
new.data$ymin <- glm_PABD_disper_weighted_unifrac$family$linkinv(new.data$fit - std * new.data$se)
new.data$ymax <- glm_PABD_disper_weighted_unifrac$family$linkinv(new.data$fit + std * new.data$se)
new.data$fit <- glm_PABD_disper_weighted_unifrac$family$linkinv(new.data$fit)  # Rescale to 0-1
plot_PABD_disper <- all_p %>%
  rename(Species=species) %>%
  ggplot() +
  geom_point(aes(x=p_disper_weighted_unifrac, y=PABD,col=Species), position=position_jitter(height=0.05, width=0), cex=3) +
  geom_ribbon(data=new.data, aes(x=p_disper_weighted_unifrac, ymin=ymin, ymax=ymax), alpha=0.2, col="lightgrey") +
  geom_line(data=new.data, aes(x=p_disper_weighted_unifrac, y=fit)) +
  ylab("Probability of infection") + 
  xlab("Dispersion (Weighted Unifrac) (Percentile within amph. species)") +
  theme_classic()

pdf("./7_figures/EFFECT_disper_PABD_WU.pdf",height=4, width=6)
plot_PABD_disper
dev.off()

#### Distance affects infection intensity ####
#### WU ####
plot_PABD_dist <- all_p %>%
  rename(Species=species) %>%
  ggplot(aes(x=p_dist_weighted_unifrac, y=Bd_load))+
  geom_point(aes(col=Species), position=position_jitter(height=0, width=0), cex=3) +
  # geom_ribbon(data=new.data, aes(x=p_inhibRich, ymin=ymin, ymax=ymax), alpha=0.2, col="lightgrey") +
  # geom_line(data=new.data, aes(x=p_inhibRich, y=fit)) +
  geom_smooth(method="lm", col="black")+
  ylab("Infection intensity (log BD load)") + 
  xlab("Instability (Weighted Unifrac) (Percentile within amph. species)") +
  theme_classic()
# plot_PABD_dist
pdf("./7_figures/EFFECT_dist_BDload_WU.pdf",height=4, width=6)
plot_PABD_dist
dev.off()

#### UWU ####
plot_PABD_dist <- all_p %>%
  rename(Species=species) %>%
  ggplot(aes(x=p_dist_unweighted_unifrac, y=Bd_load))+
  geom_point(aes(col=Species), position=position_jitter(height=0, width=0), cex=3) +
  # geom_ribbon(data=new.data, aes(x=p_inhibRich, ymin=ymin, ymax=ymax), alpha=0.2, col="lightgrey") +
  # geom_line(data=new.data, aes(x=p_inhibRich, y=fit)) +
  geom_smooth(method="lm", col="black")+
  ylab("Infection intensity (log BD load)") + 
  xlab("Instability (Unweighted Unifrac) (Percentile within amph. species)") +
  theme_classic()
# plot_PABD_dist
pdf("./7_figures/EFFECT_dist_BDload_UWU.pdf",height=4, width=6)
plot_PABD_dist
dev.off()

#### WU ####
plot_PABD_dist <- all_p %>%
  rename(Species=species) %>%
  ggplot(aes(x=p_dist_braycurtis, y=Bd_load))+
  geom_point(aes(col=Species), position=position_jitter(height=0, width=0), cex=3) +
  # geom_ribbon(data=new.data, aes(x=p_inhibRich, ymin=ymin, ymax=ymax), alpha=0.2, col="lightgrey") +
  # geom_line(data=new.data, aes(x=p_inhibRich, y=fit)) +
  geom_smooth(method="lm", col="black")+
  ylab("Infection intensity (log BD load)") + 
  xlab("Instability (Bray-curtis) (Percentile within amph. species)") +
  theme_classic()
# plot_PABD_dist
pdf("./7_figures/EFFECT_dist_BDload_BC.pdf",height=4, width=6)
plot_PABD_dist
dev.off()


#### Correlating ####

## INHIBRICH vs OVERALLRICH
all_p_combined <- rbind(all_p, all_p_con)
plot_corr_inhibRich_ASVrich <- all_p_combined %>%
  separate(indivID, into=c("Species","indiv"), remove=FALSE) %>%
  ggplot() +
  geom_point(aes(x=p_chao1, y=p_inhibRich, col=Species), cex=4) +
  # geom_smooth(aes(x=p_chao1, y=p_inhibRich), lty=2, col="grey", method="lm", se = FALSE)+
  xlab("Total Bacterial Richness (Percentile within species)")+
  ylab("Richness of inhib. bact. (Percentile within species)")+
  theme_classic()
# plot_corr_inhibRich_ASVrich
pdf("./7_figures/corr_observed_otus_inhibRich.pdf",height=4, width=6)
plot_corr_inhibRich_ASVrich
dev.off()

plot_corr_percInhib_inhibRich <- all_p_combined %>%
  separate(indivID, into=c("Species","indiv"), remove=FALSE) %>%
  ggplot() +
  geom_point(aes(x=p_percInhib, y=p_inhibRich, col=Species), cex=4) +
  # geom_smooth(aes(x=p_percInhib, y=p_inhibRich), lty=2, col="grey", method="lm", se = FALSE)+
  xlab("Proportion of inhib. bact. (Percentile within species)")+
  ylab("Richness of inhib. bact. (Percentile within species)") +
  theme_classic()
# plot_corr_ASVrich_inhibRich
pdf("./7_figures/corr_inhibRich_percInhib.pdf",height=4, width=6)
plot_corr_percInhib_inhibRich
dev.off()

plot_corr_percInhib_obsotus <- all_p_combined %>%
  separate(indivID, into=c("Species","indiv"), remove=FALSE) %>%
  ggplot() +
  geom_point(aes(x=p_chao1, y=p_percInhib, col=Species), cex=4) +
  # geom_smooth(aes(x=p_chao1, y=p_percInhib), lty=2, col="grey", method="lm", se=FALSE) +
  xlab("Total Bacterial Richness (Percentile within species)")+
  ylab("Proportion of inhib. bact. (Percentile within species)")+
  theme_classic()
# plot_corr_percInhib_obsotus
pdf("./7_figures/corr_observedotus_percInhib.pdf",height=4, width=6)
plot_corr_percInhib_obsotus
dev.off()

plot_corr_dist_disper <- all_p_combined %>%
  separate(indivID, into=c("Species","indiv"), remove=FALSE) %>%
  ggplot() +
  geom_point(aes(x=p_dist_unweighted_unifrac, y=p_disper_unweighted_unifrac, col=Species), cex=4) +
  xlab("Instability (Unweighted Unifrac) (Percentile within species)")+
  ylab("Dispersion (Unweighted Unifrac) (Percentile within species)")+
  theme_classic()
# plot_corr_percInhib_obsotus
pdf("./7_figures/corr_dist_disper.pdf",height=4, width=6)
plot_corr_dist_disper
dev.off()



#### Taxa summaries plot #####
# Change OTU sequences to names
# Get full taxonomy
getTaxa <- function(t, uppr_taxa_naming="c", naming_sep="__") {
  library(tidyverse)
  taxa_full <- t %>%
    separate(taxonomy, sep="; ", into = c("k","p","c","o","f","g","s"), remove=FALSE) %>%
    mutate(Taxa=NA)
  # manually "fill left" since it's glitching out
  for ( r in 1:nrow(taxa_full)) {
    lowest<- NA
    factored_levels<-factor( rev(c("k","p","c","o","f","g","s")), levels=rev( c("k","p","c","o","f","g","s")))
    current<-1
    while( is.na(lowest) ) {
      lowest <- taxa_full[r,as.character(factored_levels[current])]
      if ( is.na(lowest) ) {
        current <- current+1
      }
    }
    taxa_full[r,as.vector(factored_levels[1:current])] <- lowest
  }
  
  # Go through and fill in taxonomy as required
  for ( r in 1:nrow(taxa_full)) {
    utemp <- taxa_full[r,uppr_taxa_naming] 
    if (is.na(utemp)) {
      final <- "Unassigned"
    } else {
      # Get species
      stemp <- taxa_full[r,'s'] 
      sep_stemp <- unlist(strsplit(stemp, split=naming_sep ))
      if ( sep_stemp[1] =="s") {
        gtemp <- taxa_full[r,"g"] 
        sep_gtemp <- unlist(strsplit(gtemp, split=naming_sep ))
        # What if genus doesn't exist?
        lvl <- 2
        while (is.na(sep_gtemp[2])) {
          current_lvl <- c("s","g","f","o","c","p","k")[lvl]
          gtemp <- taxa_full[r,current_lvl] 
          sep_gtemp <- unlist(strsplit(gtemp, split=naming_sep ))
          lvl <- lvl+1
        }
        
        if (sep_gtemp[1]=="g") {
          if ( is.na(sep_stemp[2]) ) {
            sfinal <- paste0(sep_gtemp[2],"_sp")
          } else {
            sfinal <- paste0(sep_gtemp[2],"_",sep_stemp[2])
          }
        } else {
          if ( is.na(sep_stemp[2]) ) {
            sfinal <- paste0(gtemp,"_sp")
          } else {
            sfinal <- paste0(gtemp,"_",sep_stemp[2])
          }
        }
        
      } else {
        sfinal <- stemp
      }
      
      # Get upper level
      sep_utemp <- unlist(strsplit(utemp, split=naming_sep ))
      if ( sep_utemp[1] ==uppr_taxa_naming ) {
        final <- paste0(sep_utemp[2],":",sfinal)
      } else {
        final <- utemp
      }
      
    }
    # Save final name
    taxa_full[r,"Taxa"] <- make.unique(final, sep="-")
  }
  taxa_full[,"Taxa"] <- make.unique(taxa_full[,"Taxa"], sep="-")
  
  return(as.data.frame(taxa_full))
}
taxa <- getTaxa(t=taxonomy, uppr_taxa_naming = "c")

# Create OTU table of ONLY non-inhibitory OTUs
otu_filt_noinhib <- otu_filt %>%
  filter(!(`#OTU ID` %in% otu_filt_inhibOnly[,"#OTU ID"]))
otu_filt_inhib <- otu_filt_inhibOnly

# Change OTU names to taxonomy and create two mapping files; one for inhibitory and one for non-inhibitory OTUs
rownames(otu_filt_noinhib) <- taxa[match(otu_filt_noinhib$`#OTU ID`, taxa$Sequence),"Taxa"]
otu_filt_noinhib$`#OTU ID` <- NULL
rownames(otu_filt_inhib) <- taxa[match(otu_filt_inhib$`#OTU ID`, taxa$Sequence), "Taxa"]
otu_filt_inhib$`#OTU ID` <- NULL

mf_temp <- mf_alt_filt_final %>%
  rename(Time_Point=time) %>%
  dplyr::select(SampleID, species, indivID, Time_Point, Bd_exposure, Bd_load)

# Collapse by genus because there are too many isolates
mf_withoutinhib  <- otu_filt_noinhib %>%
  t() %>%
  as.data.frame() %>%
  mutate(SampleID=colnames(otu_filt_noinhib)) %>%
  right_join(mf_temp) %>%
  gather(-c(SampleID, species, indivID, Time_Point, Bd_exposure, Bd_load), key="OTU", value="Reads") 
# %>%
#   separate(OTU, sep="-", into=c("OTU_only","v"), remove=TRUE) %>%
#   group_by(SampleID, species, indivID, time, Bd_exposure, Bd_load, OTU_only) %>%
#   summarize(Reads=sum(Reads)) %>%
#   ungroup() %>%
#   rename(OTU=OTU_only)

mf_withinhib <- otu_filt_inhib %>%
  t() %>%
  as.data.frame() %>%
  mutate(SampleID=colnames(otu_filt_inhib)) %>%
  right_join(mf_temp)%>%
  gather(-c(SampleID, species, indivID, Time_Point, Bd_exposure, Bd_load), key="OTU", value="Reads") 
# %>%
#   separate(OTU, sep="-", into=c("OTU_only","v"), remove=TRUE) %>%
#   group_by(SampleID, species, indivID, time, Bd_exposure, Bd_load, OTU_only) %>%
#   summarize(Reads=sum(Reads))%>%
#   ungroup() %>%
#   rename(OTU=OTU_only)

# Get most abundant OTUs for inhibitory and non-inhibitory-- if over 5% overall abundance, then keep (calculate by species)
thresh <- 0.02
percReads_noninhib <- mf_withoutinhib %>%
  group_by(OTU, SampleID, species) %>%
  summarize(propReads=sum(Reads)/10000) %>%
  ungroup()
OVERALL_order_abundance_noninhib <- mf_withoutinhib %>%
  group_by(OTU) %>%
  summarize(sumReads=sum(Reads)) %>%
  ungroup() %>%
  mutate(percReads=sumReads/sum(sumReads)) %>%
  arrange(-percReads) %>%
  pull(OTU)
# NON inhib abundance
Anbo_noninhib <- percReads_noninhib %>%
  filter(species=="Anbo") %>%
  group_by(OTU) %>%
  summarize(percReads=mean(propReads)) %>%
  arrange(-percReads) %>%
  filter(percReads>thresh) %>%
  pull(OTU)
Rhma_noninhib <- percReads_noninhib %>%
  filter(species=="Rhma") %>%
  group_by(OTU) %>%
  summarize(percReads=mean(propReads)) %>%
  arrange(-percReads) %>%
  filter(percReads>thresh)%>%
  pull(OTU)
Osse_noninhib <- percReads_noninhib %>%
  filter(species=="Osse") %>%
  group_by(OTU) %>%
  summarize(percReads=mean(propReads)) %>%
  arrange(-percReads) %>%
  filter(percReads>thresh)%>%
  pull(OTU)
Raca_noninhib <- percReads_noninhib %>%
  filter(species=="Raca") %>%
  group_by(OTU) %>%
  summarize(percReads=mean(propReads)) %>%
  arrange(-percReads) %>%
  filter(percReads>thresh)%>%
  pull(OTU)
Rapi_noninhib <- percReads_noninhib %>%
  filter(species=="Rapi") %>%
  group_by(OTU) %>%
  summarize(percReads=mean(propReads)) %>%
  arrange(-percReads) %>%
  filter(percReads>thresh)%>%
  pull(OTU)
most_abund_noninhib <- unique(c(Anbo_noninhib, Rhma_noninhib, Osse_noninhib, Raca_noninhib, Rapi_noninhib))

percReads_inhib <- mf_withinhib %>%
  group_by(OTU, SampleID, species) %>%
  summarize(propReads=sum(Reads)/10000) %>%
  ungroup()
OVERALL_order_abundance_inhib <- mf_withinhib %>%
  group_by(OTU) %>%
  summarize(sumReads=sum(Reads)) %>%
  ungroup() %>%
  mutate(percReads=sumReads/sum(sumReads)) %>%
  arrange(-percReads) %>%
  pull(OTU)
# INHIB abundance
Anbo_inhib <- percReads_inhib %>%
  filter(species=="Anbo") %>%
  group_by(OTU) %>%
  summarize(percReads=mean(propReads)) %>%
  arrange(-percReads) %>%
  filter(percReads>thresh) %>%
  pull(OTU)
Rhma_inhib <- percReads_inhib %>%
  filter(species=="Rhma") %>%
  group_by(OTU) %>%
  summarize(percReads=mean(propReads)) %>%
  arrange(-percReads) %>%
  filter(percReads>thresh)%>%
  pull(OTU)
Osse_inhib <- percReads_inhib %>%
  filter(species=="Osse") %>%
  group_by(OTU) %>%
  summarize(percReads=mean(propReads)) %>%
  arrange(-percReads) %>%
  filter(percReads>thresh)%>%
  pull(OTU)
Raca_inhib <- percReads_inhib %>%
  filter(species=="Raca") %>%
  group_by(OTU) %>%
  summarize(percReads=mean(propReads)) %>%
  arrange(-percReads) %>%
  filter(percReads>thresh)%>%
  pull(OTU)
Rapi_inhib <- percReads_inhib %>%
  filter(species=="Rapi") %>%
  group_by(OTU) %>%
  summarize(percReads=mean(propReads)) %>%
  arrange(-percReads) %>%
  filter(percReads>thresh)%>%
  pull(OTU)
most_abund_inhib <- unique(c(Anbo_inhib, Rhma_inhib, Osse_inhib, Raca_inhib, Rapi_inhib))
# 
# # Add colors (old 5% threshold)
# non_inhib_OTUs_abundant_col <- data.frame(OTU=most_abund_noninhib, colours=c(
#   "goldenrod4" # [1] "[Saprospirae]:f__Chitinophagaceae_sp.6"    
#   , "darkred" # [2] "Gammaproteobacteria:Acinetobacter_sp.2"    
#   , "green"# [3] "Flavobacteriia:Chryseobacterium_sp"        
#   , "lightgoldenrod"# [4] "Flavobacteriia:o__Flavobacteriales"        
#   , "deepskyblue3"# [5] "Betaproteobacteria:f__Rhodocyclaceae"      
#   , "plum1"# [6] "Gammaproteobacteria:Perlucidibaca_sp"      
#   , "darkviolet"# [7] "Verrucomicrobiae:f__Verrucomicrobiaceae_sp"
#   , "yellow3"# [8] "Clostridia:o__Clostridiales"               
#   , "yellow"# [9] "Clostridia:g__Clostridium.1"               
#   , "brown2"# [10] "Gammaproteobacteria:Acinetobacter_sp.3"    
#   , "olivedrab"# [11] "Betaproteobacteria:f__Alcaligenaceae.1" 
#  ))
# 
# inhib_OTUs_abundant_col <- data.frame(OTU=most_abund_inhib, colours=c(
#   "magenta" # [1] "Gammaproteobacteria:g__Pseudomonas.3"        
#   , "lightslateblue" # [2] "Betaproteobacteria:f__Comamonadaceae.1"   
#   , "firebrick" # [7] "Alphaproteobacteria:f__Rhizobiaceae.1"       
#   , "mediumspringgreen" # [3] "Flavobacteriia:Chryseobacterium_sp.1"        
#   , "yellow2" # [4] "Gammaproteobacteria:Acinetobacter_guillouiae"
#   , "cyan1" # [10] "Betaproteobacteria:Comamonas_sp"        
#   , "hotpink3" # [6] "Gammaproteobacteria:g__Pseudomonas"          
#   , "coral3"# [8] "Alphaproteobacteria:f__Rhizobiaceae.2"       
#   , "goldenrod" # [9] "Gammaproteobacteria:Azorhizophilus_sp"       
#   , "orange" # [5] "Betaproteobacteria:o__Burkholderiales.1"     
# ))

# Add colors (2% threshold)
non_inhib_OTUs_abundant_col <- data.frame(OTU=most_abund_noninhib, colours=c(
  "goldenrod4"# "[Saprospirae]:f__Chitinophagaceae_sp-6"     
  , "darkred"# "Gammaproteobacteria:Acinetobacter_sp-2"    
  , "green"# "Flavobacteriia:Chryseobacterium_sp"         
  , "chartreuse" # "Actinobacteria:o__Actinomycetales-3"       
  , "darkgreen" # "Cytophagia:Emticicia_sp"                    
  , "plum1" # "Gammaproteobacteria:Perlucidibaca_sp"      
  , "lightgoldenrod" # "Flavobacteriia:o__Flavobacteriales"         
  , "deepskyblue3" # "Betaproteobacteria:f__Rhodocyclaceae"      
  , "darkorchid" # "Sphingobacteriia:f__Sphingobacteriaceae"    
  , "darkorchid1" # "Sphingobacteriia:f__Sphingobacteriaceae-1" 
  , "darkmagenta" # "Alphaproteobacteria:f__Sphingomonadaceae-2" 
  , "darkslategray4" # "Betaproteobacteria:Comamonas_terrigena"    
  , "maroon" # "Verrucomicrobiae:f__Verrucomicrobiaceae_sp" 
  , "yellow3"# "Clostridia:o__Clostridiales"               
  , "yellow"# "Clostridia:g__Clostridium-1"                
  , "brown2"# "Gammaproteobacteria:Acinetobacter_sp-3"    
  , "mediumpurple1" # "Sphingobacteriia:g__Pedobacter"             
  , "olivedrab"# "Betaproteobacteria:f__Alcaligenaceae-1"    
  , "aquamarine2"# "Betaproteobacteria:f__Comamonadaceae-26"    
  , "lightsalmon" # "Flavobacteriia:Fluviicola_sp-1"      
  #   
  
))

inhib_OTUs_abundant_col <- data.frame(OTU=most_abund_inhib, colours=c(
  "magenta" # "Gammaproteobacteria:g__Pseudomonas-3"         
  , "lightslateblue" # "Betaproteobacteria:f__Comamonadaceae-1"      
  , "firebrick" # "Alphaproteobacteria:f__Rhizobiaceae-1"        
  , "mediumspringgreen" # "Flavobacteriia:Chryseobacterium_sp-1"        
  , "yellow2" # "Gammaproteobacteria:Acinetobacter_guillouiae" 
  , "hotpink3" # "Gammaproteobacteria:g__Pseudomonas"          
  , "coral3" # "Alphaproteobacteria:f__Rhizobiaceae-2"        
  , "goldenrod" # "Gammaproteobacteria:Azorhizophilus_sp"       
  , "orange" # "Betaproteobacteria:o__Burkholderiales-1"
))



# Reorder overall ranks
OVERALL_order_abundance_inhib <- c(OVERALL_order_abundance_inhib[match(inhib_OTUs_abundant_col$OTU,OVERALL_order_abundance_inhib)], OVERALL_order_abundance_inhib[-match(inhib_OTUs_abundant_col$OTU,OVERALL_order_abundance_inhib)])
OVERALL_order_abundance_noninhib <- c(OVERALL_order_abundance_noninhib[match(non_inhib_OTUs_abundant_col$OTU,OVERALL_order_abundance_noninhib)], OVERALL_order_abundance_noninhib[-match(non_inhib_OTUs_abundant_col$OTU,OVERALL_order_abundance_noninhib)])

mf_withinhib <- mf_withinhib %>%
  left_join(inhib_OTUs_abundant_col) %>%
  mutate(OTU = factor(OTU, levels=OVERALL_order_abundance_inhib)) %>%
  filter(Reads!=0)

mf_withoutinhib <- mf_withoutinhib %>%
  left_join(non_inhib_OTUs_abundant_col) %>%
  mutate(OTU = factor(OTU, levels=OVERALL_order_abundance_noninhib))%>%
  filter(Reads!=0)

# Get named vector
set.seed(1248)
grey_random <- sample(grey.colors(n=length(OVERALL_order_abundance_noninhib)-nrow(non_inhib_OTUs_abundant_col)))
non_inhib_col_vector <- c(as.vector(non_inhib_OTUs_abundant_col$colours), grey_random[1:(length(OVERALL_order_abundance_noninhib)-nrow(non_inhib_OTUs_abundant_col))])
names(non_inhib_col_vector) <- OVERALL_order_abundance_noninhib

inhib_col_vector <- c(as.vector(inhib_OTUs_abundant_col$colours), grey_random[1:(length(OVERALL_order_abundance_inhib)-nrow(inhib_OTUs_abundant_col))] )
names(inhib_col_vector) <- OVERALL_order_abundance_inhib


# Create a treatment and control
mf_withoutinhib_con <- mf_withoutinhib %>%
  filter(Bd_exposure=="Control") %>%
  mutate(OTU = factor(OTU, levels=OVERALL_order_abundance_noninhib))
mf_withoutinhib_treat <- mf_withoutinhib %>%
  filter(Bd_exposure=="Bd-exposed")%>%
  mutate(OTU = factor(OTU, levels=OVERALL_order_abundance_noninhib))
mf_withinhib_con <- mf_withinhib %>%
  filter(Bd_exposure=="Control")%>%
  mutate(OTU = factor(OTU, levels=OVERALL_order_abundance_inhib))
mf_withinhib_treat <- mf_withinhib %>%
  filter(Bd_exposure=="Bd-exposed")%>%
  mutate(OTU = factor(OTU, levels=OVERALL_order_abundance_inhib))

#### Plotting taxa summaries-- separated ####
dir.create("./7_figures/TaxaSummaries")

# Get legend
# Non-inhibitory legend
g_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  legend
} 

pdf("./7_figures/TaxaSummaries/NonInhibitory_legend.pdf")
grid.arrange(
  g_legend(
    non_inhib_OTUs_abundant_col%>%
      rename("Non-Inhibitory OTUs"=OTU) %>%
      ggplot() +
      geom_bar(aes(x=`Non-Inhibitory OTUs`, y=1, fill=`Non-Inhibitory OTUs`), stat="identity") +
      scale_fill_manual(values=non_inhib_col_vector[1:nrow(non_inhib_OTUs_abundant_col)]) 
  )
)
dev.off()
pdf("./7_figures/TaxaSummaries/Inhibitory_legend.pdf")
grid.arrange(
  g_legend(
    inhib_OTUs_abundant_col%>%
      rename("Inhibitory OTUs"=OTU) %>%
      ggplot() +
      geom_bar(aes(x=`Inhibitory OTUs`, y=1, fill=`Inhibitory OTUs`), stat="identity") +
      scale_fill_manual(values=inhib_col_vector[1:nrow(inhib_OTUs_abundant_col)]) 
  )
)
dev.off()

####### CONTROL
# ANBO
anbo_con_noninhib <- mf_withoutinhib_con %>%
  # separate(indivID, into=c(NA,"indiv"), remove=FALSE) %>%
  filter(species=="Anbo") %>%
  ggplot() +
  geom_bar(aes(x=Time_Point, y=Reads/10000, fill=factor(OTU, levels=rev(OVERALL_order_abundance_noninhib))), stat="identity", show.legend = FALSE) +
  scale_fill_manual(values=non_inhib_col_vector) +
  facet_wrap(indivID~., nrow=1) + 
  scale_y_reverse() +
  xlab(NULL) +
  ylab("Proportion of reads") +
  theme(axis.text.x = element_blank()
        , axis.ticks.x = element_blank()
        )
anbo_con_inhib <- mf_withinhib_con %>%
  # separate(indivID, into=c(NA,"indiv"), remove=FALSE) %>%
  filter(species=="Anbo") %>%
  ggplot() +
  geom_bar(aes(x=Time_Point, y=Reads/10000, fill=factor(OTU, levels=rev(OVERALL_order_abundance_inhib))), stat="identity", show.legend = FALSE) +
  scale_fill_manual(values=inhib_col_vector) +
  facet_wrap(indivID~., nrow=1) +
  ylab("Proportion of reads") +
  theme(strip.background=element_blank()
        , strip.text.x = element_blank())
# Get width of plot
panels <- mf_withinhib_con %>%
  filter(species=="Anbo") %>%
  pull(indivID) %>%
  unique() %>%
  length()
pdf("./7_figures/TaxaSummaries/Anbo_con.pdf", height=4, width=panels*2+1)
grid.arrange(anbo_con_noninhib, anbo_con_inhib, nrow=2)
dev.off()

# RHMA
rhma_con_noninhib <- mf_withoutinhib_con %>%
  # separate(indivID, into=c(NA,"indiv"), remove=FALSE) %>%
  filter(species=="Rhma") %>%
  ggplot() +
  geom_bar(aes(x=Time_Point, y=Reads/10000, fill=factor(OTU, levels=rev(OVERALL_order_abundance_noninhib))), stat="identity", show.legend = FALSE) +
  scale_fill_manual(values=non_inhib_col_vector) +
  facet_wrap(indivID~., nrow=1) + 
  scale_y_reverse() +
  xlab(NULL) +
  ylab("Proportion of reads") +
  theme(axis.text.x = element_blank()
        , axis.ticks.x = element_blank())
rhma_con_inhib <- mf_withinhib_con %>%
  # separate(indivID, into=c(NA,"indiv"), remove=FALSE) %>%
  filter(species=="Rhma") %>%
  ggplot() +
  geom_bar(aes(x=Time_Point, y=Reads/10000, fill=factor(OTU, levels=rev(OVERALL_order_abundance_inhib))), stat="identity", show.legend = FALSE) +
  scale_fill_manual(values=inhib_col_vector) +
  facet_wrap(indivID~., nrow=1) +
  ylab("Proportion of reads") +
  theme(strip.background=element_blank()
        , strip.text.x = element_blank())

panels <- mf_withinhib_con %>%
  filter(species=="Rhma") %>%
  pull(indivID) %>%
  unique() %>%
  length()
pdf("./7_figures/TaxaSummaries/Rhma_con.pdf", height=4, width=panels*2+1)
 grid.arrange(rhma_con_noninhib, rhma_con_inhib, nrow=2)
dev.off()

# OSSE
osse_con_noninhib <- mf_withoutinhib_con %>%
  # separate(indivID, into=c(NA,"indiv"), remove=FALSE) %>%
  filter(species=="Osse") %>%
  ggplot() +
  geom_bar(aes(x=Time_Point, y=Reads/10000, fill=factor(OTU, levels=rev(OVERALL_order_abundance_noninhib))), stat="identity", show.legend = FALSE) +
  scale_fill_manual(values=non_inhib_col_vector) +
  facet_wrap(indivID~., nrow=1) + 
  scale_y_reverse() +
  xlab(NULL) +
  ylab("Proportion of reads") +
  theme(axis.text.x = element_blank()
        , axis.ticks.x = element_blank())
osse_con_inhib <- mf_withinhib_con %>%
  # separate(indivID, into=c(NA,"indiv"), remove=FALSE) %>%
  filter(species=="Osse") %>%
  ggplot() +
  geom_bar(aes(x=Time_Point, y=Reads/10000, fill=factor(OTU, levels=rev(OVERALL_order_abundance_inhib))), stat="identity", show.legend = FALSE) +
  scale_fill_manual(values=inhib_col_vector) +
  facet_wrap(indivID~., nrow=1) +
  ylab("Proportion of reads") +
  theme(strip.background=element_blank()
        , strip.text.x = element_blank())

panels <- mf_withinhib_con %>%
  filter(species=="Osse") %>%
  pull(indivID) %>%
  unique() %>%
  length()
pdf("./7_figures/TaxaSummaries/Osse_con.pdf", height=4, width=panels*2+1)
grid.arrange(osse_con_noninhib, osse_con_inhib, nrow=2)
dev.off()


# Raca
raca_con_noninhib <- mf_withoutinhib_con %>%
  # separate(indivID, into=c(NA,"indiv"), remove=FALSE) %>%
  filter(species=="Raca") %>%
  ggplot() +
  geom_bar(aes(x=Time_Point, y=Reads/10000, fill=factor(OTU, levels=rev(OVERALL_order_abundance_noninhib))), stat="identity", show.legend = FALSE) +
  scale_fill_manual(values=non_inhib_col_vector) +
  facet_wrap(indivID~., nrow=1) + 
  scale_y_reverse() +
  xlab(NULL) +
  ylab("Proportion of reads") +
  theme(axis.text.x = element_blank()
        , axis.ticks.x = element_blank())
raca_con_inhib <- mf_withinhib_con %>%
  # separate(indivID, into=c(NA,"indiv"), remove=FALSE) %>%
  filter(species=="Raca") %>%
  ggplot() +
  geom_bar(aes(x=Time_Point, y=Reads/10000, fill=factor(OTU, levels=rev(OVERALL_order_abundance_inhib))), stat="identity", show.legend = FALSE) +
  scale_fill_manual(values=inhib_col_vector) +
  facet_wrap(indivID~., nrow=1) +
  ylab("Proportion of reads") +
  theme(strip.background=element_blank()
        , strip.text.x = element_blank())

panels <- mf_withinhib_con %>%
  filter(species=="Raca") %>%
  pull(indivID) %>%
  unique() %>%
  length()
pdf("./7_figures/TaxaSummaries/Raca_con.pdf", height=4, width=panels*2+1)
grid.arrange(raca_con_noninhib, raca_con_inhib, nrow=2)
dev.off()


# RAPI
rapi_con_noninhib <- mf_withoutinhib_con %>%
  # separate(indivID, into=c(NA,"indiv"), remove=FALSE) %>%
  filter(species=="Rapi") %>%
  ggplot() +
  geom_bar(aes(x=Time_Point, y=Reads/10000, fill=factor(OTU, levels=rev(OVERALL_order_abundance_noninhib))), stat="identity", show.legend = FALSE) +
  scale_fill_manual(values=non_inhib_col_vector) +
  facet_wrap(indivID~., nrow=1) + 
  scale_y_reverse() +
  xlab(NULL) +
  ylab("Proportion of reads") +
  theme(axis.text.x = element_blank()
        , axis.ticks.x = element_blank())
rapi_con_inhib <- mf_withinhib_con %>%
  # separate(indivID, into=c(NA,"indiv"), remove=FALSE) %>%
  filter(species=="Rapi") %>%
  ggplot() +
  geom_bar(aes(x=Time_Point, y=Reads/10000, fill=factor(OTU, levels=rev(OVERALL_order_abundance_inhib))), stat="identity", show.legend = FALSE) +
  scale_fill_manual(values=inhib_col_vector) +
  facet_wrap(indivID~., nrow=1) +
  ylab("Proportion of reads") +
  theme(strip.background=element_blank()
        , strip.text.x = element_blank())

panels <- mf_withinhib_con %>%
  filter(species=="Rapi") %>%
  pull(indivID) %>%
  unique() %>%
  length()
pdf("./7_figures/TaxaSummaries/Rapi_con.pdf", height=4, width=panels*2+1)
grid.arrange(rapi_con_noninhib, rapi_con_inhib, nrow=2)
dev.off()



##### TREATMENT
# ANBO
anbo_treat_noninhib <- mf_withoutinhib_treat %>%
  # separate(indivID, into=c(NA,"indiv"), remove=FALSE) %>%
  filter(species=="Anbo") %>%
  ggplot() +
  geom_bar(aes(x=Time_Point, y=Reads/10000, fill=factor(OTU, levels=rev(OVERALL_order_abundance_noninhib))), stat="identity", show.legend = FALSE) +
  scale_fill_manual(values=non_inhib_col_vector) +
  facet_wrap(indivID~., nrow=1) + 
  scale_y_reverse() +
  xlab(NULL) +
  ylab("Proportion of reads") +
  theme(axis.text.x = element_blank()
        , axis.ticks.x = element_blank())
anbo_treat_inhib <- mf_withinhib_treat %>%
  # separate(indivID, into=c(NA,"indiv"), remove=FALSE) %>%
  filter(species=="Anbo") %>%
  ggplot() +
  geom_bar(aes(x=Time_Point, y=Reads/10000, fill=factor(OTU, levels=rev(OVERALL_order_abundance_inhib))), stat="identity", show.legend = FALSE) +
  scale_fill_manual(values=inhib_col_vector) +
  facet_wrap(indivID~., nrow=1) +
  ylab("Proportion of reads") +
  theme(strip.background=element_blank()
        , strip.text.x = element_blank())

panels <- mf_withinhib_treat %>%
  filter(species=="Anbo") %>%
  pull(indivID) %>%
  unique() %>%
  length()
pdf("./7_figures/TaxaSummaries/Anbo_treat.pdf", height=4, width=panels*2+1)
grid.arrange(anbo_treat_noninhib, anbo_treat_inhib, nrow=2)
dev.off()


# RHMA
rhma_treat_noninhib <- mf_withoutinhib_treat %>%
  # separate(indivID, into=c(NA,"indiv"), remove=FALSE) %>%
  filter(species=="Rhma") %>%
  ggplot() +
  geom_bar(aes(x=Time_Point, y=Reads/10000, fill=factor(OTU, levels=rev(OVERALL_order_abundance_noninhib))), stat="identity", show.legend = FALSE) +
  scale_fill_manual(values=non_inhib_col_vector) +
  facet_wrap(indivID~., nrow=1) + 
  scale_y_reverse() +
  xlab(NULL) +
  ylab("Proportion of reads") +
  theme(axis.text.x = element_blank()
        , axis.ticks.x = element_blank())
rhma_treat_inhib <- mf_withinhib_treat %>%
  # separate(indivID, into=c(NA,"indiv"), remove=FALSE) %>%
  filter(species=="Rhma") %>%
  ggplot() +
  geom_bar(aes(x=Time_Point, y=Reads/10000, fill=factor(OTU, levels=rev(OVERALL_order_abundance_inhib))), stat="identity", show.legend = FALSE) +
  scale_fill_manual(values=inhib_col_vector) +
  facet_wrap(indivID~., nrow=1) +
  ylab("Proportion of reads") +
  theme(strip.background=element_blank()
        , strip.text.x = element_blank())

panels <- mf_withinhib_treat %>%
  filter(species=="Rhma") %>%
  pull(indivID) %>%
  unique() %>%
  length()
pdf("./7_figures/TaxaSummaries/Rhma_treat.pdf", height=4, width=panels*2+1)
grid.arrange(rhma_treat_noninhib, rhma_treat_inhib, nrow=2)
dev.off()

# OSSE
osse_treat_noninhib <- mf_withoutinhib_treat %>%
  # separate(indivID, into=c(NA,"indiv"), remove=FALSE) %>%
  filter(species=="Osse") %>%
  ggplot() +
  geom_bar(aes(x=Time_Point, y=Reads/10000, fill=factor(OTU, levels=rev(OVERALL_order_abundance_noninhib))), stat="identity", show.legend = FALSE) +
  scale_fill_manual(values=non_inhib_col_vector) +
  facet_wrap(indivID~., nrow=1) + 
  scale_y_reverse() +
  xlab(NULL) +
  ylab("Proportion of reads") +
  theme(axis.text.x = element_blank()
        , axis.ticks.x = element_blank())
osse_treat_inhib <- mf_withinhib_treat %>%
  # separate(indivID, into=c(NA,"indiv"), remove=FALSE) %>%
  filter(species=="Osse") %>%
  ggplot() +
  geom_bar(aes(x=Time_Point, y=Reads/10000, fill=factor(OTU, levels=rev(OVERALL_order_abundance_inhib))), stat="identity", show.legend = FALSE) +
  scale_fill_manual(values=inhib_col_vector) +
  facet_wrap(indivID~., nrow=1) +
  ylab("Proportion of reads") +
  theme(strip.background=element_blank()
        , strip.text.x = element_blank())


panels <- mf_withinhib_treat %>%
  filter(species=="Osse") %>%
  pull(indivID) %>%
  unique() %>%
  length()
pdf("./7_figures/TaxaSummaries/Osse_treat.pdf", height=4, width=panels*2+1)
grid.arrange(osse_treat_noninhib, osse_treat_inhib, nrow=2)
dev.off()


# Raca
raca_treat_noninhib <- mf_withoutinhib_treat %>%
  # separate(indivID, into=c(NA,"indiv"), remove=FALSE) %>%
  filter(species=="Raca") %>%
  ggplot() +
  geom_bar(aes(x=Time_Point, y=Reads/10000, fill=factor(OTU, levels=rev(OVERALL_order_abundance_noninhib))), stat="identity", show.legend = FALSE) +
  scale_fill_manual(values=non_inhib_col_vector) +
  facet_wrap(indivID~., nrow=1) + 
  scale_y_reverse() +
  xlab(NULL) +
  ylab("Proportion of reads") +
  theme(axis.text.x = element_blank()
        , axis.ticks.x = element_blank())
raca_treat_inhib <- mf_withinhib_treat %>%
  # separate(indivID, into=c(NA,"indiv"), remove=FALSE) %>%
  filter(species=="Raca") %>%
  ggplot() +
  geom_bar(aes(x=Time_Point, y=Reads/10000, fill=factor(OTU, levels=rev(OVERALL_order_abundance_inhib))), stat="identity", show.legend = FALSE) +
  scale_fill_manual(values=inhib_col_vector) +
  facet_wrap(indivID~., nrow=1) +
  ylab("Proportion of reads") +
  theme(strip.background=element_blank()
        , strip.text.x = element_blank())


panels <- mf_withinhib_treat %>%
  filter(species=="Raca") %>%
  pull(indivID) %>%
  unique() %>%
  length()
pdf("./7_figures/TaxaSummaries/Raca_treat.pdf", height=4, width=panels*2+1)
grid.arrange(raca_treat_noninhib, raca_treat_inhib, nrow=2)
dev.off()


# RAPI
rapi_treat_noninhib <- mf_withoutinhib_treat %>%
  # separate(indivID, into=c(NA,"indiv"), remove=FALSE) %>%
  filter(species=="Rapi") %>%
  ggplot() +
  geom_bar(aes(x=Time_Point, y=Reads/10000, fill=factor(OTU, levels=rev(OVERALL_order_abundance_noninhib))), stat="identity", show.legend = FALSE) +
  scale_fill_manual(values=non_inhib_col_vector) +
  facet_wrap(indivID~., nrow=1) + 
  scale_y_reverse() +
  xlab(NULL) +
  ylab("Proportion of reads") +
  theme(axis.text.x = element_blank()
        , axis.ticks.x = element_blank())
rapi_treat_inhib <- mf_withinhib_treat %>%
  # separate(indivID, into=c(NA,"indiv"), remove=FALSE) %>%
  filter(species=="Rapi") %>%
  ggplot() +
  geom_bar(aes(x=Time_Point, y=Reads/10000, fill=factor(OTU, levels=rev(OVERALL_order_abundance_inhib))), stat="identity", show.legend = FALSE) +
  scale_fill_manual(values=inhib_col_vector) +
  facet_wrap(indivID~., nrow=1) +
  ylab("Proportion of reads") +
  theme(strip.background=element_blank()
        , strip.text.x = element_blank())


panels <- mf_withinhib_treat %>%
  filter(species=="Rapi") %>%
  pull(indivID) %>%
  unique() %>%
  length()
pdf("./7_figures/TaxaSummaries/Rapi_treat.pdf", height=4, width=panels*2+1)
grid.arrange(rapi_treat_noninhib, rapi_treat_inhib, nrow=2)
dev.off()

#### Updated LOO Random forest ####
## Accuracy
errorRate_test_all <- rbind(cbind(compare_PABD_onlyp_nosp_LOO, species=FALSE), cbind(compare_PABD_onlyp_wsp_LOO, species=TRUE)
                            , cbind(compare_PABD_count_nosp_LOO, species=FALSE), cbind(compare_PABD_count_wsp_LOO, species=TRUE)
                            ,cbind(compare_PABD_PA_nosp_LOO, species=FALSE), cbind(compare_PABD_PA_wsp_LOO, species=TRUE)) %>%
  mutate(Training_data = ifelse(otu_type=="NA", "Community traits",ifelse(otu_type=="count","Bacterial counts","Bacterial prevalence")))

errorRate_test_all %>%
  group_by(species, Training_data) %>%
  summarize(correct=mean(Test_pred==Test_obs))

observed_only <- errorRate_test_all %>%
  select(indivID, Test_obs, species) %>%
  group_by(indivID, Test_obs, species) %>%
  summarize(Infection_status=unique(Test_obs)) %>%
  ungroup() %>% select(-Test_obs) %>% mutate(Training_data="Observed outcome")

ggsave("./7_figures/RF_PABD_predictions_withspecies.pdf", height=3.5, width=5
       ,errorRate_test_all %>%
         filter(species==TRUE) %>%
         select(indivID, species, Test_pred, Training_data) %>%
         rename(Infection_status=Test_pred) %>%
         rbind(observed_only) %>%
         mutate(Training_data= factor(Training_data, levels=c("Observed outcome", "Community traits","Bacterial counts", "Bacterial prevalence")))%>%
         ggplot() + geom_tile(aes(x=Training_data, y=indivID, fill=Infection_status)
                              , col="black", width=0.9) +
         scale_fill_manual(values=c("darkred","lightblue"), name="Infection status") +
         theme_bw() +
         theme(axis.text.x = element_text(angle=90)) +
         scale_x_discrete(breaks=c("Observed outcome","Community traits","Bacterial counts","Bacterial prevalence")
                          , labels=(c("Observed\noutcome","Community\ntraits","Bacterial\ncounts","Bacterial\nprevalence")))+
         ylab("Amphibian ID") +xlab("Training Dataset") +
         geom_vline(aes(xintercept=1.5), lwd=2)
)
ggsave("./7_figures/RF_PABD_predictions_nospecies.pdf", height=3.5, width=5
       ,errorRate_test_all %>%
         filter(species==FALSE) %>%
         select(indivID, species, Test_pred, Training_data) %>%
         rename(Infection_status=Test_pred) %>%
         rbind(observed_only) %>%
         mutate(Training_data= factor(Training_data, levels=c("Observed outcome", "Community traits","Bacterial counts", "Bacterial prevalence")))%>%
         ggplot() + geom_tile(aes(x=Training_data, y=indivID, fill=Infection_status)
                              , col="black", width=0.9) +
         scale_fill_manual(values=c("darkred","lightblue"), name="Infection status") +
         theme_bw() +
         theme(axis.text.x = element_text(angle=90)) +
         scale_x_discrete(breaks=c("Observed outcome","Community traits","Bacterial counts","Bacterial prevalence")
                          , labels=(c("Observed\noutcome","Community\ntraits","Bacterial\ncounts","Bacterial\nprevalence")))+
         ylab("Amphibian ID") +xlab("Training Dataset") +
         geom_vline(aes(xintercept=1.5), lwd=2)
)


## Plotting MSE

MSE_test_all <- rbind(cbind(compare_infect_onlyp_nosp_LOO, species="No species predictor"), cbind(compare_infect_onlyp_wsp_LOO, species="With species as predictor")
                      , cbind(compare_infect_count_nosp_LOO, species="No species predictor"), cbind(compare_infect_count_wsp_LOO, species="With species as predictor")
                      ,cbind(compare_infect_PA_nosp_LOO, species="No species predictor"), cbind(compare_infect_PA_wsp_LOO, species="With species as predictor")) %>%
  group_by(species, otu_type) %>%
  summarize(MSE=mean((Test_obs-Test_pred)^2)) %>%
  unite(otu_type, species, col=group, remove=FALSE) %>%
  mutate(Training_data = ifelse(otu_type=="NA", "Community traits",ifelse(otu_type=="count","Bacterial counts","Bacterial prevalence")))
MSE_train_all_long <- rbind(cbind(MSE_infect_onlyp_wsp_LOO, species="With species as predictor")
                            , cbind(MSE_infect_onlyp_nosp_LOO, species="No species predictor")
                            , cbind(MSE_infect_count_wsp_LOO, species="With species as predictor")
                            , cbind(MSE_infect_count_nosp_LOO, species="No species predictor")
                            , cbind(MSE_infect_PA_wsp_LOO, species="With species as predictor")
                            , cbind(MSE_infect_PA_nosp_LOO, species="No species predictor")) %>%
  unite(otu_type, species, col=group, remove=FALSE) %>%
  mutate(Training_data = ifelse(otu_type=="NA", "Community traits",ifelse(otu_type=="count","Bacterial counts","Bacterial prevalence"))) %>%
  rename(MSE=error)

ggsave("./7_figures/MSE_diagnostics_all.pdf", height=3, width=5
       ,MSE_train_all_long %>%
         ggplot() + geom_violin(aes(x=Training_data, y=MSE)) + geom_jitter(aes(x=Training_data, y=MSE), height=0, width=0.2, alpha=0.2)+
         geom_point(data=MSE_test_all, aes(x=Training_data, y=MSE), col="red") +
         facet_wrap(.~species)+
         theme_bw()+
         theme(axis.text.x = element_text(angle=90)) +
         scale_x_discrete(labels=c("Bacterial\ncounts","Bacterial\nprevalence","Community\ntraits"))
)

infect_predict_all <- rbind(cbind(compare_infect_onlyp_wsp_LOO, species="With species as predictor")
                            , cbind(compare_infect_onlyp_nosp_LOO, species="No species predictor")
                            , cbind(compare_infect_count_wsp_LOO, species="With species as predictor")
                            , cbind(compare_infect_count_nosp_LOO, species="No species predictor")
                            , cbind(compare_infect_PA_wsp_LOO, species="With species as predictor")
                            , cbind(compare_infect_PA_nosp_LOO, species="No species predictor")) %>%
  unite(otu_type, species, col=group, remove=FALSE) %>%
  mutate(Training_data = ifelse(otu_type=="NA", "Community traits",ifelse(otu_type=="count","Bacterial counts","Bacterial prevalence"))) 

ggsave("./7_figures/RF_infect_predictions_all.pdf", height=3, width=6
       ,infect_predict_all %>%
         ggplot(aes(x=Test_obs, y=Test_pred, col=Training_data)) +geom_point() +
         geom_abline(aes(slope=1, intercept=0), col="black", lty=2) +
         geom_smooth(method="lm") +
         facet_wrap(.~species) +
         xlim(0,8)+ylim(0,8)+xlab(expression(paste("Observed ",italic("Bd")," load")))+ ylab(expression(paste("Random Forest Predicted ",italic("Bd")," load")))+
         theme_bw() +
         scale_color_manual(values=c("green","darkgreen","blue"), name="Training data")
)

ggsave("./7_figures/RF_infect_predictions_withspecies.pdf", height=3, width=5
       ,infect_predict_all %>%
         filter(species=="With species as predictor") %>%
         ggplot(aes(x=Test_obs, y=Test_pred, col=Training_data)) +geom_point() +
         geom_abline(aes(slope=1, intercept=0), col="black", lty=2) +
         geom_smooth(method="lm") +
         xlim(0,8)+ylim(0,8)+xlab(expression(paste("Observed ",italic("Bd")," load")))+ ylab(expression(paste("Random Forest Predicted ",italic("Bd")," load")))+
         theme_bw() +
         scale_color_manual(values=c("green","darkgreen","blue"), name="Training data")
)

ggsave("./7_figures/RF_infect_predictions_nospecies.pdf", height=3, width=5
       ,infect_predict_all %>%
         filter(species=="No species predictor") %>%
         ggplot(aes(x=Test_obs, y=Test_pred, col=Training_data)) +geom_point() +
         geom_abline(aes(slope=1, intercept=0), col="black", lty=2) +
         geom_smooth(method="lm") +
         xlim(0,8)+ylim(0,8)+xlab(expression(paste("Observed ",italic("Bd")," load")))+ ylab(expression(paste("Random Forest Predicted ",italic("Bd")," load")))+
         theme_bw() +
         scale_color_manual(values=c("green","darkgreen","blue"), name="Training data")
)

#### Importance from RF LOO ####
rename_p <- data.frame(old=c("p_disper_unweighted_unifrac","p_inhibRich","species","p_faith_pd","p_disper_weighted_unifrac"
        ,"p_disper_braycurtis","p_shannon","p_observed_otus","p_percInhib","p_dist_unweighted_unifrac"
        ,"p_dist_weighted_unifrac","p_chao1","p_dist_braycurtis"  )
      ,new= c("Dispersion (Unweighted Unifrac)" , "Richness of Bd-inhib. bact.", "Amphibian species", "Alpha diversity (Faith's PD)", "Dispersion (Weighted Unifrac)"
          , "Dispersion (Bray-curtis)", "Alpha diversity (Shannon)" , "Alpha diversity (Observed OTUs)", "Proportion of Bd-inhib. bact." , "Instability (Unweighted Unifrac)"
          , "Instability (Weighted Unifrac)", "Alpha diversity (Chao1)", "Instability (Bray-curtis)")
      ,metric=c("Dispersion","Inhibitory","Host","Alpha Diversity","Dispersion"
                 ,"Dispersion","Alpha Diversity","Alpha Diversity","Inhibitory","Instability"
                 ,"Instability","Alpha Diversity","Instability")) 
color_metrics <- c(Host="black",`Alpha Diversity`="red", Dispersion="blue",Instability="darkblue",Inhibitory="purple")
factored_p <- c("Amphibian species"
                , "Alpha diversity (Observed OTUs)"
                , "Alpha diversity (Chao1)"
                , "Alpha diversity (Shannon)"
                , "Alpha diversity (Faith's PD)"
                , "Richness of Bd-inhib. bact."
                , "Proportion of Bd-inhib. bact."
                , "Dispersion (Bray-curtis)" 
                ,"Dispersion (Unweighted Unifrac)" 
                , "Dispersion (Weighted Unifrac)"
                , "Instability (Bray-curtis)"
                , "Instability (Unweighted Unifrac)"
                , "Instability (Weighted Unifrac)"
                
)
#### Community level combined 
ggsave("./7_figures/importance_PABD_onlyp.pdf", height=4, width=3
       ,importance_PABD_onlyp_wsp_LOO %>%
         mutate(taxonomy= rename_p[match(taxonomy, rename_p$old),"new"])%>%
         arrange(-MeanDecreaseAccuracy) %>% mutate(taxonomy=factor(taxonomy, levels=unique(taxonomy))) %>%
         ggplot() + geom_jitter(aes(x=taxonomy, y=MeanDecreaseAccuracy), height=0, width=0.1, col="darkgrey")+
         geom_hline(aes(yintercept=0), col="red", alpha=0.2)+
         theme_bw()+
         theme(axis.text.x=element_text(angle=90, hjust=1))+
         ylab("Decr. in accuracy when absent")+xlab("Predictor")
)

ggsave("./7_figures/importance_infect_onlyp.pdf", height=4, width=3
       ,importance_infect_onlyp_wsp_LOO %>%
         mutate(taxonomy= rename_p[match(taxonomy, rename_p$old),"new"] )%>%
         arrange(-X.IncMSE) %>% mutate(taxonomy=factor(taxonomy, levels=unique(taxonomy))) %>%
         ggplot() + geom_jitter(aes(x=taxonomy, y=X.IncMSE), height=0, width=0.1, col="darkgrey")+
         geom_hline(aes(yintercept=0), col="red", alpha=0.2)+
         theme_bw()+
         theme(axis.text.x=element_text(angle=90, hjust=1))+
         ylab("% Increase MSE when absent")+xlab("Predictor")
)

gg_PABD <- importance_PABD_onlyp_wsp_LOO %>%
  mutate(taxonomy= rename_p[match(taxonomy, rename_p$old),"new"] , Metric_type= rename_p[match(taxonomy, rename_p$new),"metric"])%>%
  mutate(taxonomy=factor(taxonomy, levels=factored_p), Metric_type=factor(Metric_type, levels=c("Host","Alpha Diversity","Inhibitory","Dispersion","Instability"))) %>%
  ggplot() + geom_jitter(aes(x=taxonomy, y=MeanDecreaseAccuracy,col=Metric_type), height=0, width=0.1)+
  geom_hline(aes(yintercept=0), col="red", alpha=0.2)+
  theme_bw()+
  theme(axis.text.x=element_blank())+
  scale_color_manual(values=color_metrics, name="Predictor type")+
  ylab(expression(paste("Decr. in accuracy (Pres/Abs ",italic("Bd"),")")))+xlab(NULL)
gg_infect <- importance_infect_onlyp_wsp_LOO %>%
  mutate(taxonomy= rename_p[match(taxonomy, rename_p$old),"new"] , Metric_type= rename_p[match(taxonomy, rename_p$new),"metric"])%>%
  mutate(taxonomy=factor(taxonomy, levels=factored_p), Metric_type=factor(Metric_type, levels=c("Host","Alpha Diversity","Inhibitory","Dispersion","Instability"))) %>%
  ggplot() + geom_jitter(aes(x=taxonomy, y=X.IncMSE, col=Metric_type), height=0, width=0.1)+
  geom_hline(aes(yintercept=0), col="red", alpha=0.2)+
  theme_bw()+
  theme(axis.text.x=element_text(angle=270, hjust=0))+
  scale_color_manual(values=color_metrics, name="Predictor type")+
  ylab(expression(paste("% Increase MSE (",italic("Bd")," load)")))+xlab("Predictor")

pdf("./7_figures/importance_onlyp_combined.pdf")
grid.arrange(gg_PABD, gg_infect, ncol=1, layout_matrix=rbind(1,1,1,1,2,2,2,2,2,2,2))
dev.off()


## PABD

# Combined count and PA
importance_PABD_combined <- rbind(importance_PABD_count_wsp_LOO
                                        , importance_PABD_PA_wsp_LOO)
# Filter to exclude zeros
importance_PABD_combined_temp <- importance_PABD_combined %>%
  group_by(taxonomy, inhibitory, otu_type) %>%
  summarize(meanMDA=mean(MeanDecreaseAccuracy), sdMDA=sd(MeanDecreaseAccuracy)) %>%
  arrange(-meanMDA) %>% filter(meanMDA>0) 
# Keep only top 95% of all
toKeep_PABD <- importance_PABD_combined_temp %>% filter(meanMDA>quantile(importance_PABD_combined_temp$meanMDA, prob=0.95)) %>% pull(taxonomy) %>% unique()
# q95_count <- quantile(importance_PABD_count_combined%>%filter(X.IncMSE>0)%>%pull(X.IncMSE), probs = c(0.95))
ggsave("./7_figures/importance_PABD_combinedASV.pdf",height=7, width=5
       ,importance_PABD_combined%>%
         filter(taxonomy%in%toKeep_PABD) %>%
         mutate(Inhibitory=ifelse(inhibitory==1,"Yes","No"), otu_type=ifelse(otu_type=="count","Bacterial counts","Bacterial prevalence")) %>%
         arrange(-MeanDecreaseAccuracy) %>% mutate(taxonomy=factor(taxonomy, levels=unique(taxonomy))) %>%
         ggplot() + geom_jitter(aes(x=taxonomy, y=MeanDecreaseAccuracy, col=Inhibitory), height=0, width=0.1)+
         scale_color_manual(values=c("black","purple"), na.value="darkgrey")+
         theme_bw()+
         theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))+
         facet_grid(otu_type~.)+ylab("Decr. in accuracy of infection presence/absence")+xlab("Predictor")
)


## Infect

# Combined count and PA
importance_infect_combined <- rbind(importance_infect_count_wsp_LOO
                                    ,importance_infect_PA_wsp_LOO)
# Filter to exclude zeros
importance_infect_combined_temp <- importance_infect_combined %>%
  group_by(taxonomy, inhibitory, otu_type) %>%
  summarize(meanIncMSE=mean(X.IncMSE), sdMDA=sd(X.IncMSE)) %>%
  arrange(-meanIncMSE) %>% filter(meanIncMSE>0) 
# Keep only top 95% of all
toKeep_infect <- importance_infect_combined_temp %>% filter(meanIncMSE>quantile(importance_infect_combined_temp$meanIncMSE, prob=0.95)) %>% pull(taxonomy) %>% unique()
# q95_count <- quantile(importance_PABD_count_combined%>%filter(X.IncMSE>0)%>%pull(X.IncMSE), probs = c(0.95))
ggsave("./7_figures/importance_infect_combinedASV.pdf",height=7, width=5
       ,importance_infect_combined%>%
         filter(taxonomy%in%toKeep_infect) %>%
         mutate(Inhibitory=ifelse(inhibitory==1,"Yes","No"), otu_type=ifelse(otu_type=="count","Bacterial counts","Bacterial prevalence")) %>%
         arrange(-X.IncMSE) %>% mutate(taxonomy=factor(taxonomy, levels=unique(taxonomy))) %>%
         ggplot() + geom_jitter(aes(x=taxonomy, y=X.IncMSE, col=Inhibitory), height=0, width=0.1)+
         scale_color_manual(values=c("black","purple"), na.value="darkgrey")+
         theme_bw()+
         theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))+
         facet_grid(otu_type~.)+ylab(expression(paste("Incr. % MSE of predicted ",italic("Bd")," load")))+xlab("Predictor")
)




###### Random Forest diagnostics from parameter sweep ##########

# Diagnostics

load("./5_random_forest_diagnostics/infect_test_train_comparison_diag.RData")
load("./5_random_forest_diagnostics/PABD_test_train_comparison_diag.RData")
load("./5_random_forest_diagnostics/PABD_train_diag.RData")
load("./5_random_forest_diagnostics/infect_train_MSE_diag.RData")
load("./5_random_forest_diagnostics/infect_train_R2_diag.RData")

# Output
load("./5a_random_forest_output/all_importance_infect.RData")
load("./5a_random_forest_output/all_importance_PABD.RData")
load("./5a_random_forest_output/all_test_train_comparisons_infect.RData")
load("./5a_random_forest_output/all_test_train_comparisons_PABD.RData")
load("./5a_random_forest_output/summary_error_infect_MSER2.RData")
load("./5a_random_forest_output/summary_error_infect.RData")
load("./5a_random_forest_output/summary_error_PABD.RData")
load("./5a_random_forest_output/summary_indiv_infect.RData")
load("./5a_random_forest_output/all_importance_PABD.RData")
load("./5a_random_forest_output/summary_indiv_PABD.RData")


summary_error_infect_diag <- infect_test_train_comparison_diag %>%
  group_by(p_type, prop, otu_type, repl) %>%
  mutate(meanObs = mean(Test_obs)) %>%
  mutate(ressq=(Test_pred-Test_obs)^2, ressq_total = (Test_obs-meanObs)^2) %>%
  summarize(MSE=mean(ressq), R2=1-(sum(ressq)/sum(ressq_total))) %>%
  mutate(set="Test") %>%
  ungroup()

summary_error_PABD_diag <- PABD_test_train_comparison_diag %>%
  group_by(p_type, prop, otu_type, repl) %>%
  summarize(nCorrect=mean(Test_pred==Test_obs)) %>%
  mutate(set="Test") %>%
  ungroup()

infect_train_MSER2_diag <- infect_train_MSE_diag %>%
  rename(MSE=error) %>%
  dplyr::select(-error_type) %>%
  left_join(infect_train_R2_diag) %>%
  mutate(R2=error, set="Train") %>%
  dplyr::select(p_type, prop, otu_type, repl, MSE, R2, set)

PABD_train_diag <- PABD_train_diag %>%
  mutate(nCorrect=1-error) %>%
  dplyr::select(p_type, prop, otu_type, repl, nCorrect) %>%
  mutate(set = "Train")

## Plotting

ggsave(file="./7_figures/RF_diag_allInfect.pdf"
     , infect_test_train_comparison_diag %>%
       mutate(Proportion_used_for_training=factor(prop)) %>%
       mutate(p_type = ifelse(p_type=="onlyp","Community trait predictors",ifelse(p_type=="nop", "Bacterial count predictors", "Community & Bacterial predictors"))
              , otu_type = ifelse(otu_type=="PA","Bacterial prevalence","Summed Bacterial counts")) %>%
       mutate(p_type = factor(p_type, levels=c("Community trait predictors","Bacterial count predictors","Community & Bacterial predictors" ))) %>%
       ggplot(aes(x=Test_obs, y=Test_pred)) + geom_point(aes(col=Proportion_used_for_training)) +
       geom_abline(aes(intercept=0, slope=1)) +
       geom_smooth(aes(col=Proportion_used_for_training), method="lm") +
       ylim(0,9)+xlim(0,9) +
       facet_grid(otu_type~p_type) +
       xlab("Observed Bd load") + ylab("Random Forst Predicted Bd load")+
       theme_bw()
)

ggsave(file="./7_figures/RF_diag_MSE.pdf", height=3, width=6
       , rbind(summary_error_infect_diag, infect_train_MSER2_diag) %>%
         mutate(p_type = ifelse(p_type=="onlyp","Community traits",ifelse(p_type=="nop", "Bacterial counts", "Community & Bacterial"))
                , otu_type = ifelse(otu_type=="PA","Bacterial prevalence","Summed Bacterial counts")) %>%
         mutate(p_type = factor(p_type, levels=c("Community traits","Bacterial counts","Community & Bacterial" ))) %>%
         rename(Predictors_used = p_type, Set=set) %>%
         ggplot(aes(x=prop, y=MSE)) + geom_jitter(aes(col=Predictors_used), height=0, width=0.01) +
         geom_smooth(aes(col=Predictors_used, lty=Set), se=FALSE, method='lm') +
         facet_grid(.~otu_type) +
         xlab("Proportion used for training")+
         theme_bw()
)


ggsave(file="./7_figures/RF_diag_R2.pdf", height=3, width=6
       ,rbind(summary_error_infect_diag, infect_train_MSER2_diag) %>%
         mutate(p_type = ifelse(p_type=="onlyp","Community traits",ifelse(p_type=="nop", "Bacterial counts", "Community & Bacterial"))
                , otu_type = ifelse(otu_type=="PA","Bacterial prevalence","Summed Bacterial counts")) %>%
         mutate(p_type = factor(p_type, levels=c("Community traits","Bacterial counts","Community & Bacterial" ))) %>%
         rename(Predictors_used = p_type, Set=set) %>%
         mutate(R2 = ifelse(R2<0, 0,R2)) %>%
         ggplot(aes(x=prop, y=R2)) + geom_jitter(aes(col=Predictors_used), height=0, width=0.01) +
         geom_smooth(aes(col=Predictors_used, lty=Set), se=FALSE, method="lm") +
         facet_grid(.~otu_type) +
         ylim(0,1) +
         xlab("Proportion used for training")+
         theme_bw()
)

ggsave(file="7_figures/RF_diag_propCorrect.pdf", height=3, width=6
       , rbind(summary_error_PABD_diag, PABD_train_diag) %>%
         mutate(p_type = ifelse(p_type=="onlyp","Community traits",ifelse(p_type=="nop", "Bacterial counts", "Community & Bacterial"))
                , otu_type = ifelse(otu_type=="PA","Bacterial prevalence","Summed Bacterial counts")) %>%
         mutate(p_type = factor(p_type, levels=c("Community traits","Bacterial counts","Community & Bacterial" ))) %>%
         rename(Predictors_used = p_type, Set=set) %>%
         ggplot(aes(x=prop, y=nCorrect)) + geom_jitter(aes(col=Predictors_used), height=0, width=0.025) +
         geom_smooth(aes(col=Predictors_used, lty=Set), se=FALSE, method="lm") +
         facet_grid(.~otu_type) +
         xlab("Proportion used for training") + ylab("Proportion predictions correct")+
         theme_bw()

)

# 
# 
# #### RF old simulations ####
# #
# ## Simple diagnostics
# # Mean squared Error
# summary_error_infect_MSER2_filt <- summary_error_infect_MSER2 %>%
#   filter(model_otu_type %in% c("ponly_NA", "nop_count")) %>%
#   mutate(model_otu_type = ifelse(model_otu_type=="ponly_NA","Community traits","Bacterial counts"))
# summary_error_infect_MSER2 %>%
#   mutate(sdMSE=sqrt(varMSE)) %>% dplyr::select(meanMSE, sdMSE, everything())
# ggsave("./7_figures/RF_MSE_sim.pdf", height=3, width=5
#        , summary_error_infect %>%
#          filter(model_otu_type %in% c("ponly_NA", "nop_count")) %>%
#          mutate(model_otu_type = factor(ifelse(model_otu_type=="ponly_NA","Community traits","Bacterial counts"),levels=c("Community traits","Bacterial counts"))) %>%
#          ggplot() + geom_violin(aes(x=model_otu_type, y=MSE)) +
#          geom_jitter(aes(x=model_otu_type, y=MSE), width = 0.25, height=0) +
#          geom_point(data=summary_error_infect_MSER2_filt, aes(x=model_otu_type, y=meanMSE), col="red", cex=3) +
#          geom_segment(data=summary_error_infect_MSER2_filt, aes(x=model_otu_type, xend=model_otu_type, y=meanMSE-sqrt(varMSE), yend=meanMSE+sqrt(varMSE)), col="red") +
#          xlab("Predictors used") + ylab("Mean Squared Error")+
#          theme_classic()
# )
# 
# # Proportion correct
# summary_error_PABD_mean_filt <- summary_error_PABD %>%
#   group_by(model_otu_type) %>%
#   summarize(meanPropCorrect=mean(propCorrect), sd=sd(propCorrect))%>%
#   filter(model_otu_type %in% c("ponly_NA", "nop_count")) %>%
#   mutate(model_otu_type = factor(ifelse(model_otu_type=="ponly_NA","Community traits","Bacterial counts"),levels=c("Community traits","Bacterial counts")))
# summary_error_PABD %>%
#   group_by(model_otu_type) %>%
#   summarize(meanPropCorrect=mean(propCorrect), sd=sd(propCorrect))
# ggsave("./7_figures/RF_propCorrect_sim.pdf", height=3, width=5
#        , summary_error_PABD %>%
#          filter(model_otu_type %in% c("ponly_NA", "nop_count")) %>%
#          mutate(model_otu_type = factor(ifelse(model_otu_type=="ponly_NA","Community traits","Bacterial counts"),levels=c("Community traits","Bacterial counts"))) %>%
#          ggplot() + geom_violin(aes(x=model_otu_type, y=propCorrect)) +
#          geom_jitter(aes(x=model_otu_type, y=propCorrect), height=0.05, width=0.25) +
#          geom_point(data=summary_error_PABD_mean_filt, aes(x=model_otu_type, y=meanPropCorrect), col="red") +
#          geom_segment(data=summary_error_PABD_mean_filt, aes(x=model_otu_type, xend=model_otu_type, y=meanPropCorrect-sd, yend=meanPropCorrect+sd), col="red") +
#          theme_classic() + xlab("Predictors used") +ylab("Proportion correct")
# )
# 
# # Plotting only certain models
# # Infect
# all_test_train_comparisons_infect_filt <- all_test_train_comparisons_infect %>%
#   filter(model %in% c("ponly","nop"), otu_type %in% c("count",NA)) %>%
#   unite(model, otu_type, col = model_otu_type, remove=FALSE) %>%
#   rename(Predictors_used = model_otu_type) %>%
#   mutate(Predictors_used = factor(ifelse(Predictors_used=="ponly_NA","Community traits","Bacterial counts"),levels=c("Community traits","Bacterial counts")))
# summary_indiv_infect_filt <- summary_indiv_infect %>%
#   filter(model %in% c("ponly","nop"), otu_type %in% c("count",NA)) %>%
#   rename(Predictors_used = model_otu_type) %>%
#   mutate(Predictors_used = factor(ifelse(Predictors_used=="ponly_NA","Community traits","Bacterial counts"),levels=c("Community traits","Bacterial counts")))
# ggsave(filename="./7_figures/RF_full_bdload_filt.pdf", height=3, width=5
#        , all_test_train_comparisons_infect_filt %>%
#          ggplot(aes(x=Test_obs, y=Test_pred)) +
#          geom_jitter(aes(col=Predictors_used ),alpha=0.2, height=0, width=0.2) +
#          # geom_smooth(aes(col=Predictors_used), method="lm") +
#          geom_segment(data=summary_indiv_infect_filt, aes(x=Test_obs, xend=Test_obs, y=meanTestPred-sqrt(varTestPred), yend=meanTestPred+sqrt(varTestPred), col=Predictors_used)) +
#          geom_point(data=summary_indiv_infect_filt,aes(x=Test_obs, y=meanTestPred, col=Predictors_used)) +
#          geom_smooth(data=summary_indiv_infect_filt,aes(x=Test_obs, y=meanTestPred, col=Predictors_used), method="lm") +
#          geom_abline(aes(intercept=0, slope=1)) +
#          ylim(0,8) + xlim(0,8)+
#          theme_bw() +
#          ylab("Random Forest Predicted Bd load") + xlab("Observed Bd load")
# )
# 
# 
# summary_indiv_PABD_filt <- summary_indiv_PABD %>%
#   filter(model %in% c("ponly","nop"), otu_type %in% c("count",NA)) %>%
#   rename(Predictors_used = model_otu_type) %>%
#   mutate(Predictors_used = factor(ifelse(Predictors_used == "ponly_NA","Community traits","Bacterial counts"),levels=c("Community traits","Bacterial counts"))) %>%
#   separate(indivID, into=c("species","indiv"), remove=FALSE) %>%
#   arrange(species, indiv) %>%
#   mutate(indivID=factor(indivID, levels=unique(indivID)))
# # summary_error_PABD_mean_filt
# 
# observed_PABD <- summary_indiv_PABD_filt %>%
#   dplyr::select(Predictors_used,indivID, Observed ) %>%
#   filter(Predictors_used =="Bacterial counts") %>%
#   mutate(Predictors_used = "Observed outcome") %>%
#   mutate(Observed=ifelse(Observed=="INFECTED",1,0)) %>%
#   rename(Proportion_infected=Observed)
# addline_format <- function(x,...){
#   gsub('\\s','\n',x)
# }
# 
# ggsave("./7_figures/RF_full_PABD_filt.pdf", height=3, width=5
#   , summary_indiv_PABD_filt %>%
#     rename(Proportion_infected=infectRatio) %>%
#     dplyr::select(Predictors_used, indivID, Proportion_infected) %>%
#     rbind(observed_PABD) %>%
#     mutate(Predictors_used= factor(Predictors_used, levels=c("Observed outcome","Community traits","Bacterial counts"))) %>%
#     ggplot() +
#     geom_tile(aes(x=Predictors_used, y=indivID, fill=Proportion_infected), width=0.9, height=0.9) +
#     geom_vline(aes(xintercept=1.5), lwd=2)+
#     scale_fill_gradient2(low="lightblue",mid="white",high="darkred", midpoint=0.5,name="Proportion \npredicted \nto be infected (RF)") +
#     theme_classic() +
#     # theme(axis.text.x = element_text(angle=90)) +
#     xlab("Predictors used") +ylab("Amphibian ID") +
#     scale_x_discrete(breaks=c("Observed outcome","Community traits","Bacterial counts")
#                      , labels=addline_format(c("Observed outcome","Community traits","Bacterial counts")))
#   )
# 
# 
# ##### Supp
# 
# ggsave("./7_figures//RF_residuals_by_individual.pdf"
#        ,all_test_train_comparisons_infect %>%
#          group_by(indivID, model, otu_type, rep) %>%
#          summarize(ressq = (Test_pred-Test_obs)^2) %>%
#          ungroup() %>%
#          separate(indivID, into=c("species","indiv"), remove=FALSE) %>%
#          arrange(species, indiv) %>%
#          mutate(indivID =factor(indivID, levels=unique(indivID))) %>%
#          rename(Predictors_used=model) %>%
#          mutate(Predictors_used = factor(ifelse(Predictors_used=="ponly","Community traits", ifelse(Predictors_used=="nop", "Bacteria","both community and bacteria"))
#                                          ,levels=c("Community traits","ASVs","both community and ASV"))
#                 , otu_type = ifelse(otu_type=="count", "ASV count", ifelse(otu_type=="prev","Bacterial prevalence","No Bacteria"))) %>%
#          ggplot() +
#          geom_jitter(aes(x=indivID, y=sqrt(ressq), col=Predictors_used), height=0, width=0.2) +
#          facet_grid(otu_type~.)+
#          theme(axis.text.x = element_text(angle=90))+
#          ylab("Residuals")+ xlab("Amphibian ID")
# )
# 
# 
# ggsave("./7_figures/RF_error_by_individual.pdf", height=3, width=10
#        ,summary_indiv_PABD_filt %>%
#          separate(indivID, into=c("species","Individual")) %>%
#          rename(Observed_outcome=Observed) %>%
#          ggplot()  +
#          geom_jitter(aes(x=species, y=propCorrect,col=Individual, pch=Observed_outcome), height=0.0, width=0.25, cex=2) +
#          facet_grid(.~Predictors_used) +
#          theme(axis.text.x=element_text(angle=90)) +
#          xlab("Amphibian species") + ylab("Proportion correct")
# 
# )
# 
# ggsave("./7_figures/RF_infectfit_by_individual.pdf", height=3, width=6
#        , summary_indiv_infect %>%
#          filter(model_otu_type %in% c("ponly_NA","nop_count")) %>%
#          mutate(model_otu_type= factor(ifelse(model_otu_type=="ponly_NA","Community traits","Bacterial counts"), levels=c("Community traits","Bacterial counts"))) %>%
#          separate(indivID, into=c("species","indiv"), remove=FALSE) %>%
#          arrange(species, indiv) %>%
#          mutate(indivID = factor(indivID, levels=unique(indivID)), species=factor(species, levels=c("Anbo","Rhma","Osse","Raca","Rapi"))) %>%
#          ggplot(aes(x=Test_obs, y=meanTestPred)) +geom_point(aes(col=species)) +
#          # geom_jitter(data=all_test_train_comparisons_infect,aes(x=Test_obs, y=Test_pred), alpha=0.2, height=0, width=0.2) +
#          # geom_smooth(aes(col=model_otu_type), method="lm") +
#          geom_segment(aes(x=Test_obs, xend=Test_obs, y=meanTestPred-sqrt(varTestPred), yend=meanTestPred+sqrt(varTestPred), col=species)) +
#          geom_abline(aes(intercept=0, slope=1)) +
#          ylim(0,8) + xlim(0,8)  +
#          facet_wrap(.~model_otu_type) +
#          ylab("Mean predicted Bd load (by individual)") + xlab("Observed Bd load")+
#          theme_bw()
# 
#        )
# #



