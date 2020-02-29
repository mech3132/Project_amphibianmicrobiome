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
load("./5_random_forest/all_RF_predictBD.RData")

# load("./5_logratio_tests/mf_alt_filt_with_inhibOTUs.RData")

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
plot_experimental_design <- mf_alt_filt_final %>%
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
  geom_text(aes(x=-0.75, y=0.75, label=paste0("Stress: ",round(unique(mf_con$NMDS_stress_braycurtis))/100)))

gg_infect <- mf_treat  %>%
  group_by(species, indivID) %>%
  summarize("Max Bd load (log)"=max(Bd_load)) %>%
  ggplot(aes(x=species, y=`Max Bd load (log)`)) +
  geom_point(aes(col=species), cex=3, position = position_jitter(width=0.1, height=0.05), show.legend = FALSE) +
  xlab("Species")

temp1a <-  mf_con %>%
  dplyr::select(species, shannon) %>%
  mutate(metric="Shannon_diversity") %>%
  rename(value=shannon)
temp1b <-  mf_con %>%
  dplyr::select(species, observed_otus) %>%
  mutate(metric="Observed_OTUs_richness") %>%
  rename(value=observed_otus)
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
  mutate(metric="Dispersion_from_centroid_(log_unweighted_Unifrac)")%>%
  rename(value=disper_unweighted_unifrac)
temp5 <- mf_con %>%
  dplyr::select(species, dist_weighted_unifrac) %>%
  mutate(metric="Distance_from_previous_timepoint_(weighted_Unifrac)")%>%
  rename(value=dist_weighted_unifrac)


gg_all <- rbind(temp1a,temp1b,temp2,temp3,temp4, temp5) %>%
  rename(Species=species) %>%
  mutate(Metric = gsub("_"," ",metric, fixed=TRUE)) %>%
  mutate(Metric = factor(Metric, levels=c("Observed OTUs richness","Shannon diversity","Inhibitory OTU Richness","Percent Inhibitory","Dispersion from centroid (log unweighted Unifrac)", "Distance from previous timepoint (weighted Unifrac)"))) %>%
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
  xlab("Inhibitory ASV Richness (Percentile)") +
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
  ylab("Infection intensity (log BD load)") + 
  xlab("Proportion inhibitory (Percentile)") +
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
  xlab("Unweighted Unifrac dispersion from centroid (Percentile)") +
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
  xlab("Bray Curtis dispersion from centroid (Percentile)") +
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
  xlab("Weighted Unifrac dispersion from centroid (Percentile)") +
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
  xlab("Weighted Unifrac distance between time points (Percentile)") +
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
  xlab("Unweighted Unifrac distance between time points (Percentile)") +
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
  xlab("Bray-curtis distance between time points (Percentile)") +
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
  xlab("Overall ASV Richness (Percentile)")+
  ylab("Inhibitory Richness (Percentile)")+
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
  xlab("Proportion Inhibitory (Percentile)")+
  ylab("Inhibitory Richness (Percentile)") +
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
  xlab("Overall ASV Richness (Percentile)")+
  ylab("Proportion Inhibitory (Percentile)")+
  theme_classic()
# plot_corr_percInhib_obsotus
pdf("./7_figures/corr_observedotus_percInhib.pdf",height=4, width=6)
plot_corr_percInhib_obsotus
dev.off()

plot_corr_dist_disper <- all_p_combined %>%
  separate(indivID, into=c("Species","indiv"), remove=FALSE) %>%
  ggplot() +
  geom_point(aes(x=p_dist_unweighted_unifrac, y=p_disper_unweighted_unifrac, col=Species), cex=4) +
  xlab("Unweighted Unifrac distance (Percentile)")+
  ylab("Unweighted Unifrac dispersion (Percentile)")+
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
  select(SampleID, species, indivID, Time_Point, Bd_exposure, Bd_load)

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


###### Random Forest ##########

p <- all_p %>%
  dplyr::select(-c(grep("exp_*", colnames(all_p)))) 
# Get all p's
p_colnames <- colnames(p)
comm_metrics <- factor(p_colnames[grep("p_", p_colnames)], levels=c("p_observed_otus","p_chao1"
                                                                    ,"p_faith_pd","p_shannon"
                                                                    ,"p_dist_braycurtis"
                                                                    ,"p_dist_unweighted_unifrac"
                                                                    ,"p_dist_weighted_unifrac"
                                                                    ,"p_disper_braycurtis"
                                                                    ,"p_disper_unweighted_unifrac"
                                                                    ,"p_disper_weighted_unifrac"
                                                                    ,"p_inhibRich"
                                                                    ,"p_percInhib"
))
# colors for metrics
comm_metrics_col <- c(
  "red" # "p_inhibRich"                 
  , "darkred" # "p_percInhib"                 
  ,"darkblue"# "p_dist_weighted_unifrac"    
  , "mediumblue"# "p_dist_unweighted_unifrac"   
  , "lightblue"# "p_dist_braycurtis"           
  ,"mediumpurple4"# "p_disper_weighted_unifrac"  
  ,"purple" # "p_disper_unweighted_unifrac" 
  , "darkorchid1" # "p_disper_braycurtis"         
  , "yellow"# "p_faith_pd"                 
  , "green"# "p_shannon"                   
  , "orange" # "p_chao1"                     
  , "darkorange4"# "p_observed_otus"      
)
names(comm_metrics_col) <- comm_metrics

### Plotting 
ggsave(filename = "7_figures/randomForest_rank_infect.pdf", height=4, width=7,
       all_RF_predictBD$rank_withp_infect$importance_RF_rank_withp_infect %>%
         filter(!is.na(inhibitory)) %>%
         mutate(Inhibitory=ifelse(inhibitory==1, "Inhibitory","Not inhibitory")) %>%
         ggplot() +
         geom_point(aes(x=X.IncMSE, y=Inhibitory), position = position_jitter(width=0, height=0.1)) +
         geom_vline(data=data.frame(metric=comm_metrics, value=all_RF_predictBD$rank_withp_infect$importance_RF_rank_withp_infect[match(comm_metrics,all_RF_predictBD$rank_withp_infect$importance_RF_rank_withp_infect$taxa),"X.IncMSE" ] )
                    , mapping=aes(xintercept=value, col=comm_metrics)
                    , lwd=2, alpha=0.5) +
         scale_color_manual(values=comm_metrics_col)
)


ggsave(filename = "7_figures/randomForest_rank_PABD.pdf", height=4, width=7,
       all_RF_predictBD$rank_withp_PABD$importance_RF_rank_withp_PABD %>%
         filter(!is.na(inhibitory)) %>%
         mutate(Inhibitory=ifelse(inhibitory==1, "Inhibitory","Not inhibitory")) %>%
         ggplot() +
         geom_point(aes(x=MeanDecreaseAccuracy, y=Inhibitory), position = position_jitter(width=0, height=0.1)) +
         geom_vline(data=data.frame(metric=comm_metrics, value=all_RF_predictBD$rank_withp_PABD$importance_RF_rank_withp_PABD[match(comm_metrics,all_RF_predictBD$rank_withp_PABD$importance_RF_rank_withp_PABD$taxa),"MeanDecreaseAccuracy" ] )
                    , mapping=aes(xintercept=value, col=comm_metrics)
                    , lwd=2, alpha=0.5) +
         scale_color_manual(values=comm_metrics_col)
)

# Training set
ggsave(filename = "7_figures/randomForest_validation_PABD.pdf", height=3, width=5,
       data.frame(Observed=all_RF_predictBD$rank_onlyp_PABD$test_train_comparisonrank_onlyp_PABD[,2]
                  ,Community_level_only=all_RF_predictBD$rank_onlyp_PABD$test_train_comparisonrank_onlyp_PABD[,1]
                  ,Community_and_OTU_level=all_RF_predictBD$rank_withp_PABD$test_train_comparisonrank_withp_PABD[,1]
                  ,OTU_level_only=all_RF_predictBD$rank_nop_PABD$test_train_comparisonrank_nop_PABD[,1]
       ) %>%
         mutate(Individual=factor(seq(1,7))) %>%
         gather(-Individual, key=TrainingSet, value=Bd_State) %>%
         mutate(TrainingSet=factor(TrainingSet, levels=c("Observed"
                                                         , "Community_level_only"
                                                         , "OTU_level_only"
                                                         , "Community_and_OTU_level"))) %>%
         ggplot() + geom_tile(aes(x=TrainingSet, y=Individual, fill=Bd_State),col="black" )+
         theme(axis.text.x = element_text(angle=45, vjust = 1, hjust = 1)
               , rect = element_blank()) +
         geom_vline(aes(xintercept=1.5), lwd=4))


ggsave(filename = "./7_figures/randomForest_validation_infect.pdf", height=3, width=5,
       data.frame(Observed=all_RF_predictBD$rank_onlyp_infect$test_train_comparisonrank_onlyp_infect[,2]
                  ,Community_level_only=all_RF_predictBD$rank_onlyp_infect$test_train_comparisonrank_onlyp_infect[,1]
                  ,Community_and_OTU_level=all_RF_predictBD$rank_withp_infect$test_train_comparisonrank_withp_infect[,1]
                  ,OTU_level_only=all_RF_predictBD$rank_nop_infect$test_train_comparisonrank_nop_infect[,1]
       ) %>%
         gather(Community_level_only,Community_and_OTU_level,OTU_level_only, key=Model, value=Predicted) %>%
         ggplot() +
         geom_point(aes(x=Observed, y=Predicted, col=Model)) +
         geom_abline(aes(intercept=0, slope=1))+
         xlim(0,8)+ylim(0,8)+
         xlab("Observed Bd load") + ylab("Predicted Bd load")
)



