#!bin/bash

######## Individual ASV analysis ######
library(rstanarm)
library(bayestestR)
library(vegan) # for vegdist, adonis
library(MASS) # For isoMDS
library(ape) # for pcoa
library(tidyverse) # for data manipulation
library(gridExtra)
library(ggpubr)
##### Load data ######
load("./3_5sp_mapping_otu_downstream/mf_alt_filt_final.RData")
load("./3_5sp_mapping_otu_downstream/otu_filt.RData")
load("./3_5sp_mapping_otu_downstream/braycurtis_filt.RData")
# bc_dist <- read.delim("./2_diversity_outputs/exported_braycurtis/distance-matrix.tsv", sep="\t", header=TRUE, row.names = 1)
wu_dist <- read.delim("./2_diversity_outputs/exported_weighted_unifrac/distance-matrix.tsv", sep="\t", header=TRUE, row.names = 1)
# colnames(bc_dist) <- gsub("^X","",colnames(bc_dist))
colnames(wu_dist) <- gsub("^X","",colnames(wu_dist))
taxonomy <- read.delim("3_5sp_mapping_otu_downstream/taxonomy.txt")
source("./msc_randomforest_codebits.R")
allTaxa <- getTaxa(taxonomy)
listTaxa <- allTaxa$Taxa
names(listTaxa) <- allTaxa$Sequence
# load("./3_5sp_mapping_otu_downstream/otu_filt_inhibOnly.RData")
###### Set up mf ######
mf_con <- mf_alt_filt_final %>%
  filter(Bd_exposure == "Control")
# otu_con <- otu_filt %>% select(one_of(mf_con$SampleID))
# rownames(otu_con) <- otu_filt$`#OTU ID`
# bc_dist_con <- bc_dist[mf_con$SampleID,mf_con$SampleID]
wu_dist_con <- wu_dist[mf_con$SampleID,mf_con$SampleID]
wu_dist_treat <- wu_dist[mf_treat$SampleID,mf_treat$SampleID]

# bc_dist_con <- 

mf_treat <- mf_alt_filt_final %>%
  filter(Bd_exposure == "Bd-exposed")

mf_nonbd <- mf_treat %>% filter(prepost == "Pre") %>% full_join(mf_con)


mf_collapsed_bd <- mf_alt_filt_final %>%
  filter(Bd_exposure == "Bd-exposed") %>% 
  group_by(species, indivID, prepost) %>%
  summarize(maxBd = max(Bd_load, na.rm = TRUE), medBd = median(Bd_load, na.rm = TRUE), fracBd = sum(PABD, na.rm=FALSE)/n(), sdBdMax = Bd_load_sd[which.max(Bd_load)]
            , maxBd2 = max(Bd_load2, na.rm=TRUE), medBd2 = median(Bd_load2, na.rm=TRUE), sdBdMax2 = Bd_load_sd[which.max(Bd_load2)]) %>%
  group_by(species, indivID) %>% mutate(maxBd_outcome = max(maxBd, na.rm=TRUE), fracBd_outcome = max(fracBd, na.rm=TRUE), sdBdMax = sdBdMax[which.max(maxBd)]
                                        , maxBd_outcome2 = max(maxBd2, na.rm=TRUE), sdBdMax2 = sdBdMax2[which.max(maxBd2)]) %>%
  filter(prepost=="Pre") %>%
  mutate(PABD = ifelse(fracBd_outcome>0,1,0)) %>%
  select(species, indivID, maxBd_outcome, maxBd_outcome2, fracBd_outcome, PABD, sdBdMax, sdBdMax2) %>% #rename(Species = species) %>%
  ungroup() 


###### CHOSEN MODELS
allVar <- c("observed_otus","shannon"
            ,"inhibRich","percInhib"
            ,"dist_braycurtis","dist_unweighted_unifrac","dist_weighted_unifrac"
            ,"disper_braycurtis","disper_unweighted_unifrac","disper_weighted_unifrac")
chosen_models <- c(observed_otus = "lognormal", chao1 = "lognormal", shannon = "normal", faith_pd = "gamma"
                   , disper_braycurtis = "lognormal", disper_unweighted_unifrac = "lognormal", disper_weighted_unifrac = "lognormal"
                   , dist_braycurtis = "beta", dist_unweighted_unifrac = "beta", dist_weighted_unifrac = "beta"
                   , percInhib = "beta", inhibRich = "poisson")
allVarNames <- c("ASV richness", "Shannon diversity"
                 , "Richness of\nputative inhibitory bacteria", "Proportion\nputative inhibitory bacteria"
                 , "Instability through time\n(Bray-curtis)", "Instability through time\n(unweighted Unifrac)", "Instability through time\n(weighted Unifrac)"
                 , "Dispersion from centroid\n(Bray-curtis)", "Dispersion from centroid\n(unweighted Unifrac)", "Dispersion from centroid\n(weighted Unifrac)")
names(allVarNames) <- allVar
allVarFilt <- c("observed_otus","shannon","inhibRich","percInhib","dist_weighted_unifrac","disper_weighted_unifrac")

dir.create("5_Msc_Stats")
dir.create("5_Msc_Stats/plots")
dir.create("5_Msc_Stats/stats")

######### Basic figures ###########

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
plot_experimental_design


pdf(file = "./5_Msc_Stats/plots/exp_design.pdf", width = 4.5,height = 10.5)
plot_experimental_design
dev.off()
png(file = "./5_Msc_Stats/plots/exp_design.png", width = 4.5,height = 10.5)
plot_experimental_design
dev.off()

##### General outcomes plotting #####
gg_infectionOnly <- mf_collapsed_bd %>%
  ggplot() + geom_point(aes(x=fracBd_outcome, y=maxBd_outcome, col=species), position=position_dodge2(width=0.1), cex = 2) +
  xlab("Fraction of timepoints infected") + ylab("Maximum Bd load (log)") + labs(col="Amphibian\nSpecies")
ggsave(filename = "5_Msc_Stats/plots/infection_outcomes.png", height=3, width=4
       ,gg_infectionOnly)

mf_con %>%
  select(species, indivID, time, one_of(allVarFilt)) %>%
  pivot_longer(-c(species, indivID, time), names_to = "Metric", values_to = "value") %>%
  ggplot() + geom_point(aes(x=time, y=value, col=species), position=position_dodge2(width=0.2, )) +
  facet_grid(Metric ~ species, scales = "free", switch="y")

##### Adonis #######
if (!file.exists("5_Msc_Stats/stats/permanova_exp_treat.txt")){
  permanova_sptime_con <- adonis2(as.dist(wu_dist_con) ~ species + time + species:time, data=mf_con)
  permanova_exp_treat <- adonis2(as.dist(wu_dist[mf_treat$SampleID,mf_treat$SampleID]) ~ species*time + prepost + PABD + Bd_load, data=mf_treat, by="terms")

  capture.output(permanova_sptime_con, file="5_Msc_Stats/stats/permanova_sptime_con.txt")
  capture.output(permanova_exp_treat, file="5_Msc_Stats/stats/permanova_exp_treat.txt")
} 

##### @~~~~~~~~~~~~ PCoA ~~~~~~~~~~~~@ #####
pcoa_con <- pcoa(as.dist(wu_dist_con))
# nmds_con <- isoMDS(as.dist(wu_dist_con), k = 2)
axisVar <- pcoa_con$values$Relative_eig
ax_plot <- c(1,2)
gg_pcoa <- as.data.frame(pcoa_con$vectors) %>%
  rownames_to_column(var="SampleID") %>% left_join(mf_con) %>%
  rename_at(vars(all_of(colnames(as.data.frame(pcoa_con$vectors))[ax_plot])), ~c(paste0("PCoA ", ax_plot[1]), paste0("PCoA ", ax_plot[2]))) %>%
  ggplot() + geom_point(aes(x=get(paste0("PCoA ", ax_plot[1])), y=get(paste0("PCoA ", ax_plot[2])), col=species), cex=2) +
  xlab(paste0("PCoA ", ax_plot[1], " (", round(axisVar[ax_plot[1]]*100,1), "% of variation)"))+ ylab(paste0("PCoA ", ax_plot[2], " (", round(axisVar[ax_plot[2]]*100,1), "% of variation)")) +
  labs(col="Species")
ggsave(filename = "5_Msc_Stats/plots/pcoa_con.png", height=4, width=5
       ,gg_pcoa)

pcoa_treat <- pcoa(as.dist(wu_dist_treat))
# nmds_con <- isoMDS(as.dist(wu_dist_con), k = 2)
axisVar <- pcoa_treat$values$Relative_eig
ax_plot <- c(1,2)
gg_pcoa_treat <- as.data.frame(pcoa_treat$vectors) %>%
  rownames_to_column(var="SampleID") %>% left_join(mf_treat) %>%
  mutate(Bd_load=ifelse(Bd_load==0, NA, Bd_load)) %>%
  filter(prepost=="Post") %>%
  rename_at(vars(all_of(colnames(as.data.frame(pcoa_treat$vectors))[ax_plot])), ~c(paste0("PCoA ", ax_plot[1]), paste0("PCoA ", ax_plot[2]))) %>%
  ggplot() + geom_point(aes(x=get(paste0("PCoA ", ax_plot[1])), y=get(paste0("PCoA ", ax_plot[2])), col=Bd_load),cex=2) +
  xlab(paste0("PCoA ", ax_plot[1], " (", round(axisVar[ax_plot[1]]*100,1), "% of variation)"))+ ylab(paste0("PCoA ", ax_plot[2], " (", round(axisVar[ax_plot[2]]*100,1), "% of variation)")) +
  labs(col=expression(paste("Max ",italic("Bd")," load"))) +scale_color_gradient(low="grey", high="red", na.value = "black") +
  facet_wrap(.~species)
gg_pcoa_treat
ggsave(filename = "5_Msc_Stats/plots/pcoa_treat.png", height=4, width=6
       ,gg_pcoa_treat)

##### Multipanel figure summary #####
# 
# gg_betaplot <- nmds_con$points%>% as.data.frame() %>% rownames_to_column(var="SampleID") %>%
#   left_join(mf_con) %>%
#   ggplot() +geom_point(aes(x=V1, y=V2, col=species, alpha=time), cex=3, show.legend=FALSE) +
#   ylab("NMDS 2 (Bray-curtis)") +
#   xlab("NMDS 1 (Bray-curtis)") +
#   geom_text(aes(x=-0.25, y=-0.6, label=paste0("Stress: ",round(unique(nmds_con$stress))/100)))

gg_expdesign <-mf_alt_filt_final %>%
  mutate(Contaminated = factor(ifelse(orig_contam ==1, "!Contaminated",NA), levels=c("!Contaminated"))
         , Bd_logload = (Bd_load)) %>%
  ggplot(aes(x=time, y=indiv)) +
  geom_line(aes(group=indivID, col=Bd_exposure)) +
  geom_point(aes(group=indivID, bg=Bd_logload, col=Bd_exposure), cex=2, pch=21)+
  scale_color_manual(values=c("black","blue","orange")) +
  scale_fill_gradient(low = "white", high = "red") +
  geom_vline(aes(xintercept=5.5), col="orange")+
  # geom_point(aes(group=indivID, col=Contaminated), cex=1, pch=19)+ ## NEW LINE
  facet_wrap(~species, nrow=5) +
  xlab("Time point") +
  ylab("Individual Amphibian") +
  theme(strip.text = element_text(size=14), 
        axis.title = element_text(size=14))+
  labs(fill=expression(paste(italic("Bd")," load (log+1)"))
       , col = expression(paste(italic("Bd"), " exposure"))) + theme_grey() + theme(strip.background = element_rect(colour="black"))
gg_expdesign

gg_betaplot <- as.data.frame(pcoa_con$vectors) %>%
  rownames_to_column(var="SampleID") %>% left_join(mf_con) %>%
  rename_at(vars(all_of(colnames(as.data.frame(pcoa_con$vectors))[ax_plot])), ~c(paste0("PCoA ", ax_plot[1]), paste0("PCoA ", ax_plot[2]))) %>%
  ggplot() + geom_point(aes(x=get(paste0("PCoA ", ax_plot[1])), y=get(paste0("PCoA ", ax_plot[2])), col=species), cex=2, show.legend = FALSE) +
  xlab(paste0("PCoA ", ax_plot[1], " (", round(axisVar[ax_plot[1]]*100,1), "% of variation)"))+ ylab(paste0("PCoA ", ax_plot[2], " (", round(axisVar[ax_plot[2]]*100,1), "% of variation)")) 
gg_betaplot

gg_infect <- mf_treat  %>%
  group_by(species, indivID) %>%
  summarize("Max Bd load (log)"=max(Bd_load)) %>%
  ggplot(aes(x=species, y=`Max Bd load (log)`)) +
  geom_point(aes(col=species), cex=3, position = position_dodge2(width=0.5), show.legend = FALSE) +
  xlab("Amphibian Species") + ylab(expression(paste("Max ",italic("Bd")," load (log+1)"))) + theme_bw()
gg_infect

temp1a <-  mf_con %>%
  dplyr::select(species, shannon) %>%
  mutate(Metric=allVarNames["shannon"]) %>%
  rename(value=shannon)
temp1b <-  mf_con %>%
  dplyr::select(species, observed_otus) %>%
  mutate(Metric=allVarNames["observed_otus"]) %>%
  rename(value=observed_otus)
temp2 <- mf_con %>%
  dplyr::select(species, inhibRich) %>%
  mutate(Metric=allVarNames["inhibRich"])%>%
  rename(value=inhibRich)
temp3 <- mf_con %>%
  dplyr::select(species, percInhib) %>%
  mutate(Metric=allVarNames["percInhib"])%>%
  rename(value=percInhib)
temp4 <- mf_con %>%
  dplyr::select(species, disper_weighted_unifrac) %>%
  mutate(Metric=allVarNames["disper_weighted_unifrac"])%>%
  rename(value=disper_weighted_unifrac)
temp5 <- mf_con %>%
  dplyr::select(species, dist_weighted_unifrac) %>%
  mutate(Metric=allVarNames["dist_weighted_unifrac"])%>%
  rename(value=dist_weighted_unifrac)


gg_all <- rbind(temp1a,temp1b,temp2,temp3,temp4, temp5) %>%
  rename(Species=species) %>%
  # mutate(Metric = gsub("_"," ",metric, fixed=TRUE)) %>%
  mutate(Metric = factor(Metric, levels=allVarNames[allVarFilt])) %>%
  ggplot(aes(x=Species, y=value)) +
  geom_boxplot() +
  geom_point(aes(col=Species), position = position_jitter(width=0.1, height=0), alpha=1/3, show.legend = FALSE)+
  facet_wrap(Metric~., scales = "free",nrow=3, strip.position = "left") +
  ylab("")+
  xlab("Amphibian Species") +theme_bw() + labs(col="Amphibian\nSpecies")
gg_all

# Legend
gg_legendspecies <- mf_con %>% ggplot(aes(x=NA, y=NA, col=species)) +
  geom_point() + labs(col="Amphibian\nSpecies")
g_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  legend
} 
g_speciesLeg <- g_legend(gg_legendspecies)

lay <- rbind(c(4,4,4,1,1,1,3,3,3,5),
             c(4,4,4,2,2,2,2,2,2,5),
             c(4,4,4,2,2,2,2,2,2,5)
)



# options(repr.plot.height=8, repr.plot.width=12)
plot_microbial_overview <- grid.arrange(gg_betaplot, gg_all, gg_infect,gg_expdesign,g_speciesLeg, layout_matrix = lay)
plot_microbial_overview
ggsave(filename = "./5_Msc_Stats/plots/data_summary_controls_multipanel.pdf",plot_microbial_overview
       , height=8, width=11)
ggsave(filename = "./5_Msc_Stats/plots/data_summary_controls_multipanel.png",plot_microbial_overview
       , height=8, width=11)

##### @~~~~~~~~~~~~ Centroid ~~~~~~~~~~~~@ #####

#### Calculate centroid by speciesxtime ####
# groups by species, time
cenGroups_spt <- mf_con %>% select(SampleID, species, time) %>% unite(species, time, col="speciestime") 
centroids_control_sptime <- matrix(nrow=length(unique(cenGroups_spt$speciestime)), ncol=ncol(pcoa_con$vectors), dimnames = list(rownames=unique(cenGroups_spt$speciestime), colnames=colnames(pcoa_con$vectors)))
for ( c in unique(cenGroups_spt$speciestime) ) {
  tempSampleID <- cenGroups_spt %>% filter(speciestime==c) %>% pull(SampleID)
  centroids_control_sptime[c,] <- colMeans(rbind(pcoa_con$vectors[tempSampleID,],pcoa_con$vectors[tempSampleID,]))
}

#### Distance between centroids ####
# Get distance between centroids
centroidDistances <- matrix(nrow=nrow(centroids_control_sptime), ncol=nrow(centroids_control_sptime), dimnames=list(rownames(centroids_control_sptime), rownames(centroids_control_sptime)))
for ( i in 1:(nrow(centroids_control_sptime))) {
  for (j in 1:(nrow(centroids_control_sptime))) {
    iname <- rownames(centroids_control_sptime)[i]
    jname <- rownames(centroids_control_sptime)[j]
    centroidDistances[iname,jname] <- sqrt(sum((centroids_control_sptime[i,] - centroids_control_sptime[j,])^2))
  }
}

centroidDistances_long <- centroidDistances %>% as.data.frame() %>%rownames_to_column(var="speciestime") %>%
  pivot_longer(-c(speciestime), names_to = "Centroid2", values_to = "distance") %>%
  rename(Centroid1 = speciestime) %>%
  separate(Centroid1, into=c("species1", "time1"), remove=FALSE, sep="_") %>%
  separate(Centroid2, into=c("species2", "time2"), remove=FALSE, sep="_")  %>%
  filter(!is.na(distance)) %>% mutate(time1=as.numeric(time1), time2=as.numeric(time2)) %>%
  mutate(species1 = factor(species1, levels=c("Anbo","Rhma","Osse","Raca","Rapi")), species2 = factor(species2, levels=c("Anbo","Rhma","Osse","Raca","Rapi")))

#### Difference between individual and speciesxtime centroid ####
dist_from_centroid_con <- mf_con %>% select(SampleID, species, indivID, time) %>% right_join(cenGroups_spt) %>% mutate(dist_from_centroid=NA) %>% as.data.frame()
for ( s in 1:nrow(dist_from_centroid_con) ){
  sid_temp <- (dist_from_centroid_con[s,"SampleID"])
  spt_temp <- (dist_from_centroid_con[s,"speciestime"])
  dist_from_centroid_con[which(dist_from_centroid_con$SampleID==sid_temp),"dist_from_centroid"] <- sqrt(sum((centroids_control_sptime[spt_temp,] - pcoa_con$vectors[sid_temp,])^2))
}


##### Difference between individual and STARTING species centroid ####
start_centroid <- centroids_control_sptime %>% as.data.frame() %>%
  rownames_to_column(var="speciestime") %>% separate(speciestime, into=c("species","time"), sep="_", remove=FALSE) %>%
  group_by(species) %>% mutate(mintime = min(as.numeric(time))) %>% ungroup() %>%
  filter(mintime==as.numeric(time)) %>% select(-c(speciestime,time,mintime)) %>%
  data.frame(row.names=1)
dist_from_startcentroid_con <- mf_con %>% select(SampleID, species, indivID, time) %>% right_join(cenGroups_spt) %>% mutate(dist_from_start=NA) %>% as.data.frame()
for ( s in 1:nrow(dist_from_startcentroid_con) ){
  sid_temp <- (dist_from_startcentroid_con[s,"SampleID"])
  sp_temp <- paste0(dist_from_startcentroid_con[s,"species"])
  dist_from_startcentroid_con[which(dist_from_startcentroid_con$SampleID==sid_temp),"dist_from_start"] <- sqrt(sum((start_centroid[sp_temp,] - pcoa_con$vectors[sid_temp,])^2))
}



##### Plots for centroid ######
### Distance from starting community 
centroidDistances_long %>%
  group_by(species1, species2) %>% mutate(smallesttime=min(time1)) %>% ungroup() %>%
  filter(species1==species2, time1==smallesttime) %>% filter(time2!=1) %>%
  ggplot() + geom_line(aes(x=time2, y=distance, col=species1))
test_dat_temp <- centroidDistances_long %>%
  group_by(species1, species2) %>% mutate(smallesttime=min(time1)) %>% ungroup() %>%
  filter(species1==species2, time1==smallesttime) %>% filter(time2!=1) 
stan_distance_from_starting <- stan_glm(distance ~ -1 +species1 +time2 , data=test_dat_temp)
sink("5_Msc_Stats/stats/centroid_stats.txt", append = FALSE)
# plot(stan_distance_from_starting)
"=================Effect of time on distance from starting community (species-level):================="
citemp <- ci(rstan::extract(stan_distance_from_starting$stanfit)$beta[,6], method="HDI")
"Median:"
median(rstan::extract(stan_distance_from_starting$stanfit)$beta[,6])
"PD: "
sum(rstan::extract(stan_distance_from_starting$stanfit)$beta[,6]>0)/length(rstan::extract(stan_distance_from_starting$stanfit)$beta[,6])
"95%HDI: "
c(citemp$CI_low, citemp$CI_high)
sink()

### Distance from previous time point
centroidDistances_long %>%
  group_by(species1, species2) %>% mutate(smallesttime=min(time1)) %>% ungroup() %>%
  filter(species1==species2, (time2-time1)==1) %>% 
  ggplot() + geom_line(aes(x=time1, y=distance, col=species1))
test_dat_temp <- centroidDistances_long %>%
  group_by(species1, species2) %>% mutate(smallesttime=min(time1)) %>% ungroup() %>%
  filter(species1==species2, (time2-time1)==1)
stan_distance_from_previous <- stan_glm(distance ~ -1 +species1 +time1 , data=test_dat_temp)
plot(stan_distance_from_previous)
citemp <- ci(rstan::extract(stan_distance_from_previous$stanfit)$beta[,6], method="HDI")

sink("5_Msc_Stats/stats/centroid_stats.txt", append = TRUE)
# plot(stan_distance_from_starting)
"=================Effect of time on distance from previous adjacent community (species-level):================="
"Median:"
median(rstan::extract(stan_distance_from_previous$stanfit)$beta[,6])
"PD: "
sum(rstan::extract(stan_distance_from_previous$stanfit)$beta[,6]>0)/length(rstan::extract(stan_distance_from_previous$stanfit)$beta[,6])
"95%HDI: "
c(citemp$CI_low, citemp$CI_high)
sink()


### Differences between species
centroidDistances_long %>%
  group_by(species1, species2) %>% mutate(smallesttime=min(time1)) %>% ungroup() %>%
  filter(species1!=species2, time1==time2) %>%
  unite(species1, col="sp_compare",species2, sep="_", remove=FALSE) %>%
  ggplot() + geom_line(aes(x=(time1), y=distance, col=species2)) +
  facet_grid(.~species1)
centroidDistances_long %>%
  group_by(species1, species2) %>% mutate(smallesttime=min(time1)) %>% ungroup() %>%
  filter(species1!=species2, time1==time2) %>%
  unite(species1, col="sp_compare",species2, sep="_", remove=FALSE) %>%
  ggplot() + geom_boxplot(aes(x=factor(time1), y=distance)) 

test_dat_temp <- centroidDistances_long %>%
  group_by(species1, species2) %>% mutate(smallesttime=min(time1)) %>% ungroup() %>%
  filter(species1!=species2, time1==time2) %>%
  unite(species1, col="sp_compare",species2, sep="_", remove=FALSE)
stan_distance_betweenspecies <- stan_glm(distance ~ -1 +species2*time1 , data=test_dat_temp)
plot(stan_distance_betweenspecies)
namesTemp <- names(fixef(stan_distance_betweenspecies))
citemp <- apply(rstan::extract(stan_distance_betweenspecies$stanfit)$beta[,c(6,7,8,9,10)], 2, function(x) c(ci(x,method="HDI")$CI_low,ci(x,method="HDI")$CI_high ))
colnames(citemp ) <- namesTemp[c(6,7,8,9,10)]
sink("5_Msc_Stats/stats/centroid_stats.txt", append = TRUE)
"=================Difference between species at each time point (species-level):================="
"Median:"
medtemp <- apply(rstan::extract(stan_distance_betweenspecies$stanfit)$beta[,c(6,7,8,9,10)], 2, function(x) median(x))
names(medtemp) <- namesTemp[c(6,7,8,9,10)]
medtemp
"PD: "
pdtemp <- apply(rstan::extract(stan_distance_betweenspecies$stanfit)$beta[,c(6,7,8,9,10)],2,function(x) sum(x>0)/length(x))
names(pdtemp) <- namesTemp[c(6,7,8,9,10)]
pdtemp
"95%HDI: "
citemp
sink() 

### Individual distance from centroid by tp
dist_from_centroid_con %>%
  filter(dist_from_centroid!=0) %>%
  ggplot() + geom_line(aes(x=time, y=dist_from_centroid, group=indivID, col=species)) +
  facet_grid(.~species)

### Individual distance from START centroid
dist_from_startcentroid_con %>%
  filter(dist_from_start!=0) %>%
  ggplot() + geom_line(aes(x=time, y=dist_from_start, group=indivID, col=species)) +
  facet_grid(.~species)

############### @~~~~~~~~~~~~ Distance within infect/noninfect ~~~~~~~~~~~~@ ############
mf_treatpost <- mf_alt_filt_final %>% filter(prepost=="Post", Bd_exposure == "Bd-exposed")
# wu_dist_treatpost <- wu_dist[mf_treatpost$SampleID,mf_treatpost$SampleID ]
if (!file.exists("5_Msc_Stats/paired_dist_byinfect.RData")) {
  paired_dist_byinfect <- data.frame()
  for ( sp in unique(mf_treatpost$species) ) {
    mf_temp <- mf_treatpost %>% filter(species==sp)
    wu_dist_temp <- wu_dist[mf_temp$SampleID, mf_temp$SampleID]
    wu_dist_temp[lower.tri(wu_dist_temp)] <- NA
    dist_long_temp <- wu_dist_temp %>% as.data.frame() %>% rownames_to_column(var="samp1") %>%
      pivot_longer(-samp1, names_to = "samp2", values_to="dist_paired") %>%
      filter(!is.na(dist_paired), dist_paired>0) %>%
      left_join(mf_temp %>% rename(samp1=SampleID, sp1=species, indiv1=indivID, t1=time, PABD1=PABD) %>% select(samp1, sp1, indiv1, t1, PABD1)) %>%
      left_join(mf_temp %>% rename(samp2=SampleID, sp2=species, indiv2=indivID, t2=time, PABD2=PABD) %>% select(samp2, sp2, indiv2, t2, PABD2))
    
    paired_dist_byinfect <- rbind(paired_dist_byinfect,dist_long_temp)
    #   
    # print(sp)
    # total <- nrow(mf_temp)*(nrow(mf_temp)-1)
    # pb <- txtProgressBar(min = 0, max = total, style = 3)
    # i <- 1
    # for ( r1 in 1:(nrow(mf_temp)-1) ) {
    #   indiv1 <-  unlist(mf_temp[r1,"indivID"])
    #   t1 <-  unlist(mf_temp[r1,"time"])
    #   PABD1 <- unlist(mf_temp[r1,"PABD"])
    #   for ( r2 in 2:nrow(mf_temp) ) {
    #     indiv2 <-  unlist(mf_temp[r2,"indivID"])
    #     t2 <-  unlist(mf_temp[r2,"time"])     
    #     PABD2 <- unlist(mf_temp[r2,"PABD"])
    #     
    #     dist_paired <- wu_dist[unlist(mf_temp[r1,"SampleID"]), unlist(mf_temp[r2,"SampleID"])]
    #     
    #     paired_dist_byinfect <- rbind(paired_dist_byinfect, data.frame(species=sp, PABD1 = PABD1, PABD2=PABD2, dist_paired = dist_paired, indiv1 = indiv1, indiv2=indiv2, t1=t1, t2=t2))
    #     
    #     #------- PROGRESS BAR ------- #
    #     Sys.sleep(0.1)
    #     i <- i+1
    #     # update progress bar
    #     setTxtProgressBar(pb, i)
    #     #-----------------------------#
    #   }
    # }
    # 
  }
  save(paired_dist_byinfect, file = "5_Msc_Stats/paired_dist_byinfect.RData")
} else {
  load("5_Msc_Stats/paired_dist_byinfect.RData")
}

paired_dist_byinfect %>%
  filter(t1==t2, indiv1!=indiv2, PABD1==PABD2) %>%
  ggplot(aes(x=t1, y=dist_paired, col=factor(PABD1))) + geom_point() +
  geom_smooth() +
  facet_grid(sp1~.)

##### @~~~~~~~~~~~~ ASV Turnover ~~~~~~~~~~~~@ #####

# Mf for combining
mf_toCombine <- mf_alt_filt_final %>% select(species, indivID, time, prepost, Bd_exposure, Bd_load, PABD) %>%
  rename(PABD_single = PABD) %>%
  left_join(mf_collapsed_bd)


##### Individuals; adjacent timepoints and total timepoints ####
indivList <- unique(mf_alt_filt_final$indivID)
summary_turnover_adjacent <- data.frame()
ASV_turnover_summary_adjacent <- data.frame()
summary_turnover_allTime <- data.frame()
ASV_turnover_summary_allTime <- data.frame()
for ( id in indivList) {
  subSamp_temp <- mf_alt_filt_final %>% filter(indivID ==id ) %>% arrange(time) %>%
    select(SampleID, time)
  # By species
  sp_temp <- strsplit(id, "_")[[1]][1]
  for ( r in 1:(nrow(subSamp_temp)-1)) {
    # r=1
    timeTemp <- pull(subSamp_temp[c(r),"time"])
    sampTemp <- pull(subSamp_temp[c(r,r+1),"SampleID"])
    datTemp <- otu_filt[,sampTemp]/colSums(otu_filt[,sampTemp]) 
    summary_temp <- datTemp %>%as_tibble() %>%
      rename_at(vars(all_of(sampTemp)), ~c("samp1","samp2")) %>%
      filter(!(samp1==0 & samp2==0)) %>%
      mutate(ASVpresent = samp1>0,lostASV = (samp2==0 & samp1>0), gainASV = (samp1==0 & samp2>0),keptASV = (samp1>0 & samp2>0)) %>%
      rowwise() %>%mutate(relAbundsame = min(c(samp2, samp1))) %>% 
      select(-c(samp1, samp2)) %>%
      colSums()
    
    byASV_temp <- datTemp %>% as_tibble() %>%
      cbind(data.frame(ASVID=otu_filt$`#OTU ID`)) %>%
      rename_at(vars(all_of(sampTemp)), ~c("samp1","samp2")) %>%
      filter(!(samp1==0 & samp2==0)) %>%
      mutate(change = abs(samp2-samp1), sign = sign(samp2-samp1)) %>%
      select(ASVID, samp1, change, sign)
    
    summary_turnover_adjacent <- rbind(summary_turnover_adjacent, data.frame(species = sp_temp, indivID = id, time = timeTemp, t(summary_temp)))
    ASV_turnover_summary_adjacent <- rbind(ASV_turnover_summary_adjacent, data.frame(species = sp_temp, indivID = id, time = timeTemp, byASV_temp))

  }
  ### All timepoints
  allSamps_temp <- pull(subSamp_temp[,"SampleID"])
  dat_tempall <- otu_filt[,allSamps_temp]/colSums(otu_filt[,allSamps_temp])
  allTimePoints_temp <- dat_tempall %>% mutate(fracConstant = apply(dat_tempall,1,min), ASVID=otu_filt$`#OTU ID`, total = rowSums(dat_tempall)) %>%
    filter(total > 0) %>%mutate(indivID = id, species = sp_temp)
  
  summaryAll_temp <- allTimePoints_temp %>% select(indivID, species, ASVID, fracConstant) %>% group_by(indivID, species) %>%
    dplyr::summarise(relAbundsame = sum(fracConstant), ASVpresent=n(), nASVssam = sum(fracConstant>0)) %>% ungroup()
  
  asvsummary_temp <- allTimePoints_temp %>% select(species, indivID, ASVID, fracConstant)%>% filter(fracConstant>0)
  
  ASV_turnover_summary_allTime <- rbind(ASV_turnover_summary_allTime, asvsummary_temp)
  summary_turnover_allTime<- rbind(summary_turnover_allTime,summaryAll_temp)
}
summary_turnover_adjacent <- summary_turnover_adjacent %>% mutate(species = factor(species, levels=c("Anbo","Rhma","Osse","Raca","Rapi")))
ASV_turnover_summary_adjacent <- ASV_turnover_summary_adjacent %>% mutate(species = factor(species, levels=c("Anbo","Rhma","Osse","Raca","Rapi")))
summary_turnover_allTime <- summary_turnover_allTime %>% mutate(species = factor(species, levels=c("Anbo","Rhma","Osse","Raca","Rapi")))
ASV_turnover_summary_allTime <- ASV_turnover_summary_allTime %>% mutate(species = factor(species, levels=c("Anbo","Rhma","Osse","Raca","Rapi")))

# head(ASV_turnover_summary_adjacent)

### Relationship between persistence and abundance

#### Do ASVs that stay have higher relative abundance?
ASV_turnover_summary_adjacent %>% group_by(species,indivID, ASVID) %>%
  dplyr::summarise(aveRelAbund = mean(samp1), sdRelAbund = sd(samp1), nPrev = sum(samp1>0)/n()) %>%
  ungroup() %>% filter(nPrev>0) %>%
  ggplot(aes(x=aveRelAbund, y=nPrev))+
  # geom_line() +
  geom_point() + facet_grid(.~species)

gg_asvbyindiv <- summary_turnover_allTime %>% left_join(mf_collapsed_bd) %>%
  ggplot() + geom_boxplot(aes(x=species, y=ASVpresent)) +
  geom_jitter(aes(x=species, y=nASVssam*(max(summary_turnover_allTime$ASVpresent)/max(summary_turnover_allTime$nASVssam))+200), col="blue", width=0.2, height=0) +
  labs(col="Maximum\nBd load") +xlab("Species") +
  scale_y_continuous(name="Total number of ASVs\nobserved on an individual (log)", sec.axis = sec_axis(trans=~./(max(summary_turnover_allTime$ASVpresent)/max(summary_turnover_allTime$nASVssam)), name="Number ASVs persistent\nthrough all time points on an individual")) +
  # scale_color_gradient(low="grey",high="red", na.value = "black") + 
  theme(axis.title.y.right = element_text(color="blue")
           , axis.text.y.right = element_text(color="blue"))
gg_asvbyindiv
ggsave(filename = "5_Msc_Stats/plots/ASV_on_individual_over_time.png", height=4, width=6
       ,gg_asvbyindiv)

summary_turnover_allTime %>% left_join(mf_collapsed_bd) %>% mutate(Bd_exposure = ifelse(is.na(PABD), "Control","Bd-exposed")) %>%
  ggplot() + geom_jitter(aes(x=species, y=relAbundsame), width=0.2, height=0) +
  labs(col="Maximum\nBd load") +xlab("Species") +ylab("Fraction of reads that\nremained consistent through time") +
  # scale_color_gradient(low="grey",high="red", na.value = "black") + 
  theme_bw()

summary_turnover_allTime %>% left_join(mf_collapsed_bd) %>% mutate(Bd_exposure = ifelse(is.na(PABD), "Control","Bd-exposed")) %>%
  ggplot() + geom_jitter(aes(x=species, y=nASVssam, col=Bd_exposure), width=0.2, height=0) +
  labs(col="Treatment") +xlab("Species") +ylab("Number of ASVs that\nremained consistent through time") +
  # scale_color_gradient(low="grey",high="red", na.value = "black") + 
  theme_bw()
gg_asvthroughtime_treatment <- summary_turnover_allTime %>% left_join(mf_collapsed_bd) %>% mutate(Bd_exposure = ifelse(is.na(PABD), "Control","Bd-exposed")) %>%
  mutate(percASV = nASVssam/ASVpresent) %>%
  ggplot() + geom_jitter(aes(x=species, y=percASV, col=Bd_exposure), width=0.2, height=0) +
  labs(col="Treatment") +xlab("Species") +ylab("Fraction of ASVs that\nremained consistent through time") +
  # scale_color_gradient(low="grey",high="red", na.value = "black") + 
  theme_bw()
gg_asvthroughtime_treatment
ggsave(filename = "5_Msc_Stats/plots/ASVthroughtime_bytreatment.png",height=4, width=6,
       gg_asvthroughtime_treatment)
##### ASVs kept over time by individual #####
summary_asvs_through_time <- summary_turnover_allTime %>% group_by(species) %>%
  dplyr::summarise(ASVs_kept_over_time = mean(nASVssam), , sd_ASVs_kept_over_time = sd(nASVssam), percASV_kept_over_time = mean(nASVssam/ASVpresent), sd_percASV_kept_over_time = sd(nASVssam/ASVpresent))
summary_asvs_through_time
write.table(summary_asvs_through_time, file = "5_Msc_Stats/stats/summary_asvs_through_time_by_indiv.txt", quote=FALSE, row.names=FALSE, sep="\t")

summary_turnover_allTime_adj <- summary_turnover_allTime %>% left_join(mf_collapsed_bd) %>% mutate(Bd_exposure = ifelse(is.na(PABD), 0,1))

##### STAN: number ASVs across all time; treatment and infection outcome ####
stan_coreASVsTreatControl <- stan_glmer(nASVssam ~ Bd_exposure + (1|species) + (1|indivID),data=summary_turnover_allTime_adj
                                        , family = poisson(link="identity")
                                        , iter=5000, adapt_delta = )

plot(stan_coreASVsTreatControl, pars=c("Bd_exposure"))
plot(stan_coreASVsTreatControl)
samps <- rstan::extract(stan_coreASVsTreatControl$stanfit)
# Average ASVs across all species
sink("5_Msc_Stats/stats/ASV_stats.txt", append = FALSE)
"==============Effect of Bd exposure on number of ASVs over all time points============="
"Intercept; average ASVs overall"
"median: "
median(samps$alpha)
"median Per time point"
median(mf_con$observed_otus)
"PD"
min(c(sum(samps$alpha>0)/length(samps$alpha),sum(samps$alpha<0)/length(samps$alpha)))
"95% CI"
c(ci(samps$alpha, method="HDI")$CI_low, ci(samps$alpha, method="HDI")$CI_high)
"PD"
"Effect of Bd on ASVs across all timepoints"
effects_samps <- samps$beta
"median"
median(effects_samps)
"PD"
min(c(sum(effects_samps>0)/length(effects_samps),sum(effects_samps<0)/length(effects_samps)))
"95% CI"
c(ci(effects_samps, method="HDI")$CI_low,ci(effects_samps, method="HDI")$CI_high)
sink()

set.seed(9823)
stan_coreASVsTreatControl_infect <- stan_glmer(nASVssam ~ Bd_exposure:PABD + Bd_exposure:PABD:maxBd_outcome + (1|species) + (1|indivID),data=summary_turnover_allTime_adj
                                        , family = poisson(link="identity"))
plot(stan_coreASVsTreatControl_infect, pars=c("Bd_exposure:PABD","Bd_exposure:PABD:maxBd_outcome"))
plot(stan_coreASVsTreatControl_infect)
samps <- rstan::extract(stan_coreASVsTreatControl_infect$stanfit)
effects_samps <- samps$beta
effect_inter <- samps$alpha
sink("5_Msc_Stats/stats/ASV_stats.txt", append = TRUE)
"================Effect of bd on number of ASVs ============="
## INTERCEPT
"Effect of PABD"
"median:"
median(effect_inter[,1])
"PD:"
min(c(sum(effect_inter[,1]>0)/length(effect_inter[,1]),sum(effect_inter[,1]<0)/length(effect_inter[,1])))
"95%CI"
c(ci(effect_inter[,1], method="HDI")$CI_low, ci(effect_inter[,1], method="HDI")$CI_high)

##PABD
"Effect of PABD"
"median:"
median(effects_samps[,1])
"PD:"
min(c(sum(effects_samps[,1]>0)/length(effects_samps[,1]),sum(effects_samps[,1]<0)/length(effects_samps[,1])))
"95%CI"
c(ci(effects_samps[,1], method="HDI")$CI_low, ci(effects_samps[,1], method="HDI")$CI_high)

##load
"Effect of Bd load"
"median:"
median(effects_samps[,2])
"PD:"
c(min(sum(effects_samps[,2]>0)/length(effects_samps[,2]),sum(effects_samps[,2]<0)/length(effects_samps[,2])))
"95%CI"
c(ci(effects_samps[,2], method="HDI")$CI_low, ci(effects_samps[,2], method="HDI")$CI_high)
## Number of ASVs kept between time points
"Number ASVs kept between time points"
median(summary_turnover_adjacent$keptASV)
median(effects_samps$alpha)
sink()

ASVlabels <- c(percASVLost = "ASVs lost\nsince last time point", percASVGain = "ASVs gained\nsince last time point", percASVKept = "ASVs retained\nsince last time point", relAbundsame="Fraction of reads\nthat remained consistent\nbetween time points")
speciesLabels <- c(Anbo="Anbo", Rhma = "Rhma", Osse = "Osse", Raca= "Raca", Rapi = "Rapi")
mf_for_indiv <- mf_alt_filt_final %>% select(species, indivID, time, Bd_exposure, Bd_load, PABD) 
summary_turnover_adjacent %>% as_tibble() %>%
  left_join(mf_for_indiv) %>%
  mutate(percASVLost = lostASV/ASVpresent, percASVGain = gainASV/ASVpresent, percASVKept = keptASV/ASVpresent) %>%
  select(species, indivID, time, Bd_exposure, percASVLost, percASVGain, percASVKept, relAbundsame, Bd_load, PABD) %>%
  pivot_longer(-c(species, indivID, time,Bd_exposure, Bd_load, PABD), names_to = "Metric", values_to = "Percent") %>%
  mutate(group = paste0(indivID, Metric)) %>%
  rowwise() %>% mutate(Metric = ASVlabels[Metric]) %>%
  # filter(Metric=="percASVGain") %>%
  ggplot() + geom_line(aes(x=time, y=Percent, col=Metric, group=group)) +
  geom_point(aes(x=time, y=Percent, fill=Bd_load), pch=21) +
  facet_grid(species~Bd_exposure, scales = "free") +
  scale_color_manual(values=c("darkgreen","orange","blue","magenta"))+
  scale_fill_gradient(low="grey", high="red") +
  labs(col = "Metric", fill="Bd load (log)")
gg_ASVsretained <- summary_turnover_adjacent %>% as_tibble() %>%
  left_join(mf_for_indiv) %>%
  mutate(percASVLost = lostASV/ASVpresent, percASVGain = gainASV/ASVpresent, percASVKept = keptASV/ASVpresent) %>%
  select(species, indivID, time, keptASV, Bd_load, PABD, Bd_exposure) %>%
  # pivot_longer(-c(species, indivID, time,Bd_exposure, Bd_load, PABD), names_to = "Metric", values_to = "Percent") %>%
  # mutate(group = paste0(indivID, Metric)) %>%
  # filter(Metric == "percASVKept") %>%
  # rowwise() %>% mutate(Metric = ASVlabels[Metric]) %>%
  # filter(Metric=="percASVGain") %>%
  ggplot() + geom_line(aes(x=time, y=keptASV, group=indivID)) +
  geom_point(aes(x=time, y=keptASV, fill=Bd_load), pch=21) +
  facet_grid(species~Bd_exposure, scales = "free") +
  # scale_color_manual(values=c("darkgreen","orange","blue","magenta"))+
  scale_fill_gradient(low="grey", high="red") +
  labs(fill=expression(paste(italic("Bd")," load (log)"))) + xlab("Time") + ylab("Percent of ASVs retained since last time point")
gg_ASVsretained
ggsave(filename = "5_Msc_Stats/plots/ASVs_retained_through_time_treatment.png", height=6, width=5,
       gg_ASVsretained)

ASVlabels <- c(percASVLost = "ASVs lost\nsince last time point", percASVGain = "ASVs gained\nsince last time point", percASVKept = "ASVs retained\nsince last time point", relAbundsame="Fraction of reads\nthat remained consistent\nbetween time points")
speciesLabels <- c(Anbo="Anbo", Rhma = "Rhma", Osse = "Osse", Raca= "Raca", Rapi = "Rapi")
summary_turnover_adjacent %>% as_tibble() %>%
  mutate(percASVLost = lostASV/ASVpresent, percASVGain = gainASV/ASVpresent, percASVKept = keptASV/ASVpresent) %>%
  select(species, indivID, time, percASVLost, percASVGain, percASVKept, relAbundsame) %>%
  pivot_longer(-c(species, indivID, time), names_to = "Metric", values_to = "Percent") %>%
  mutate(group = paste0(indivID, Metric)) %>%
  # filter(Metric=="percASVGain") %>%
  ggplot(aes(x=time, y=Percent, col=Metric, group=group)) + geom_line() +
  facet_grid(Metric~species, scales = "free", labeller = as_labeller(c(ASVlabels, speciesLabels)), switch="y")

##### ASVs in common across individuals and through time #####
timeList <- sort(unique(mf_alt_filt_final$time))
otu_filt_mat <- otu_filt
rownames(otu_filt_mat) <- otu_filt$`#OTU ID`
otu_filt_mat <- otu_filt_mat[,-match("#OTU ID", colnames(otu_filt_mat))]

if (!file.exists("5_Msc_Stats/allSharedASVs_byindiv.RData")) {
  allSharedASVs_byindiv <- data.frame()
  allSharedASVs_bysp <- data.frame()
  for ( t in timeList ) {
    # t=1
    # t = timeList[1]
    dat_t_temp <- mf_alt_filt_final %>% filter(time == t) 
    indivList_temp <- dat_t_temp$indivID
    print(paste0("TIME: ",t))
    print("individual by individual")
    
    total <- length(indivList_temp)^2
    pb <- txtProgressBar(min = 0, max = total, style = 3)
    i <- 1
    for ( i1 in indivList_temp ) {
      sid_1 <- dat_t_temp %>% filter(indivID==i1) %>% pull(SampleID)
      for ( i2 in indivList_temp ) {
        sid_2 <- dat_t_temp %>% filter(indivID == i2) %>% pull(SampleID)
        otu_temp <- otu_filt_mat[,c( sid_1, sid_2)]
        # i2 = indivList_temp[2]
        sharedFrac <- sum(apply(otu_temp,1,min))/colSums(otu_temp)[1]
        sharedASVs <- sum(rowSums(otu_temp>0)==2)
        totalASVs <- sum(otu_temp[,sid_1]>0)
        sharedFracASVs <- sharedASVs/totalASVs
        
        allSharedASVs_byindiv <- rbind(allSharedASVs_byindiv, data.frame(indivID_1 =i1, indivID_2 = i2, time = t, sharedFrac = sharedFrac, sharedASVs = sharedASVs, sharedFracASVs = sharedFracASVs, totalASVs= totalASVs))
        #------- PROGRESS BAR ------- #
        Sys.sleep(0.1)
        i <- i+1
        # update progress bar
        setTxtProgressBar(pb, i)
        #-----------------------------#
      }
    }
    
    print("Species by species")
    spList <- unique(dat_t_temp %>% pull(species))
    total <- length(spList)
    pb <- txtProgressBar(min = 0, max = total, style = 3)
    i <- 1
    for ( i1 in spList ) {
      sid_1 <- dat_t_temp %>% filter(species==i1) %>% pull(SampleID)
      otu_temp <- otu_filt_mat[,c( sid_1)]
      sharedFrac <- sum(apply(otu_temp,1,min))/colSums(otu_temp)[1]
      sharedASVs <- sum(rowSums(otu_temp>0)==length(sid_1))
      totalASVs <- sum(otu_temp[,sid_1]>0)
      sharedFracASVs <- sharedASVs/totalASVs
      
      allSharedASVs_bysp <- rbind(allSharedASVs_bysp, data.frame(species =i1, time = t, sharedFrac = sharedFrac, sharedASVs = sharedASVs, sharedFracASVs = sharedFracASVs, totalASVs = totalASVs))
      #------- PROGRESS BAR ------- #
      Sys.sleep(0.1)
      i <- i+1
      # update progress bar
      setTxtProgressBar(pb, i)
      #-----------------------------#
      
    }
  }
  save(allSharedASVs_byindiv, file="5_Msc_Stats/allSharedASVs_byindiv.RData")
  save(allSharedASVs_bysp, file="5_Msc_Stats/allSharedASVs_bysp.RData")
  
} else {
  load("5_Msc_Stats/allSharedASVs_byindiv.RData")
  load("5_Msc_Stats/allSharedASVs_bysp.RData")
}

allSharedASVs_byindiv_adj <- allSharedASVs_byindiv %>%
  separate(indivID_1, into = c("sp_1","indiv"), remove=FALSE, sep="_") %>%
  separate(indivID_2, into = c("sp_2","indiv"), remove=FALSE, sep="_") %>%
  select(-indiv)
  # rowwise() %>% mutate(sp_1 = pull(mf_alt_filt_final[match(indivID_1, mf_alt_filt_final$SampleID),"species"]), sp_2 = pull(mf_alt_filt_final[match(indivID_2, mf_alt_filt_final$SampleID),"species"])) %>%
  # ungroup()

allSharedASVs_byindiv_adj %>% filter(sp_1==sp_2, indivID_1 !=indivID_2) %>%
  ggplot() + geom_point(aes(x=time, y=sharedFrac, col=sp_2)) +
  facet_grid(.~sp_1)
allSharedASVs_byindiv_adj %>% filter(sp_1==sp_2, indivID_1 !=indivID_2) %>%
  ggplot() + geom_point(aes(x=time, y=sharedFracASVs, col=sp_2)) +
  facet_grid(.~sp_1)
allSharedASVs_byindiv_adj %>% filter(sp_1==sp_2, indivID_1 !=indivID_2) %>%
  ggplot() + geom_point(aes(x=time, y=sharedASVs, col=sp_2)) +
  facet_grid(.~sp_1)
allSharedASVs_byindiv_adj %>% filter(sp_1==sp_2, indivID_1 !=indivID_2) %>%
    ggplot() + geom_point(aes(x=time, y=totalASVs))  +
  facet_grid(.~sp_1)

allSharedASVs_bysp %>%
  ggplot() + geom_point(aes(x=time, y=sharedFrac, col=species))

allSharedASVs_bysp %>%
  ggplot() + geom_point(aes(x=time, y=sharedASVs, col=species)) +
  geom_line(aes(x=time, y=sharedASVs, col=species))


allSharedASVs_bysp %>%
  ggplot() + geom_point(aes(x=time, y=sharedFracASVs, col=species))

#### Shared fraction of ASVS within species through time #####
# lm_shared_btwn_species_through_time <- stan_glm(sharedFracASVs ~ species*time, data=allSharedASVs_bysp
#                                                 , family=gaussian)
# plot(lm_shared_btwn_species_through_time)
# plot(lm_shared_btwn_species_through_time, pars=c("time"))

 
# lm_shared_btwn_indiv_through_time <- stan_glmer(sharedFracASVs ~ sp_1 + time + (1|indivID_1), data=allSharedASVs_byindiv_filt
#                                                 , family=gaussian)
# plot(lm_shared_btwn_indiv_through_time)
# ### No effect of time for shared ASVs within species

#### STAN: Shared fraction of ASVS between individual through time #####

allSharedASVs_byindiv_filt <-allSharedASVs_byindiv_adj %>% filter(indivID_1!=indivID_2, sp_1==sp_2)
if (!file.exists("5_Msc_Stats/lm_shared_btwn_indiv_through_time.RData") ){
  set.seed(2934)
  lm_shared_btwn_indiv_through_time <- stan_glmer(sharedFracASVs ~ time + (1|indivID_1) + (1|sp_1), data=allSharedASVs_byindiv_filt
                                                  ,family=gaussian
                                                  , iter = 4000, adapt_delta=0.999)
  save(lm_shared_btwn_indiv_through_time, file="5_Msc_Stats/lm_shared_btwn_indiv_through_time.RData")
} else {
  load("5_Msc_Stats/lm_shared_btwn_indiv_through_time.RData")
}

plot(lm_shared_btwn_indiv_through_time, pars=c("time"))
samps<- rstan::extract(lm_shared_btwn_indiv_through_time$stanfit)
sink("5_Msc_Stats/stats/ASV_stats.txt", append = TRUE)
"================== Fraction of shared ASVs between individuals (of same species) through time ==============="
# IN general, overall decrease in shared ASVs between individuals over time, but this changes depending on individual.
"Effect of time"
"median"
median(samps$beta[,1])
"PD"
min(c(sum(samps$beta[,1]>0)/length(samps$beta[,1]),sum(samps$beta[,1]<0)/length(samps$beta[,1])))
"95% CI"
c(ci(samps$beta[,1], method="HDI")$CI_low, ci(samps$beta[,1], method="HDI")$CI_high)

"Intercept: shared between individuals"
"median"
median(samps$alpha[,1])
"PD"
min(c(sum(samps$alpha[,1]>0)/length(samps$alpha[,1]),sum(samps$alpha[,1]<0)/length(samps$alpha[,1])))
"95% CI"
c(ci(samps$alpha[,1], method="HDI")$CI_low, ci(samps$alpha[,1], method="HDI")$CI_high)

sink()

gg_sharedindiv_throughtime <- allSharedASVs_byindiv_filt %>%
  unite(indivID_1, indivID_2, col="bothIndiv", remove=FALSE) %>%
  # filter(indivID_1 == "Anbo_2") %>%
  ggplot() + geom_point(aes(x=time, y=sharedFracASVs, col=sp_1)) +
  geom_line(aes(x=time, y=sharedFracASVs, group=bothIndiv), alpha=0.1)+
  facet_grid(.~sp_1)+ labs(col="Amphibian\nSpecies") + ylab("Fraction of shared ASVs\nbetween individuals of same species") +
  xlab("Time point")
gg_sharedindiv_throughtime
ggsave(filename = "5_Msc_Stats/plots/FracASVs_shared_by_indiv.png", width=10, height=4
       ,gg_sharedindiv_throughtime)

# ## by ASVs
# if (!file.exists("5_Msc_Stats/lm_shared_btwn_indiv_through_time_ASV.RData")) {
#   lm_shared_btwn_indiv_through_time_ASV <- stan_glmer(sharedASVs ~ time + (1|indivID_1) + (1|sp_1), data=allSharedASVs_byindiv_filt
#                                                       ,family=poisson(link="log")
#                                                       , iter=4000
#                                                       , adapt_delta=  0.999)
#   save(lm_shared_btwn_indiv_through_time_ASV, file = "5_Msc_Stats/lm_shared_btwn_indiv_through_time_ASV.RData")
# } else {
#   load("5_Msc_Stats/lm_shared_btwn_indiv_through_time_ASV.RData")
# }
# 
# plot(lm_shared_btwn_indiv_through_time_ASV, pars=c("time"))
# plot(lm_shared_btwn_indiv_through_time_ASV)
# samps<- rstan::extract(lm_shared_btwn_indiv_through_time_ASV$stanfit)
# # IN general, overall decrease in shared ASVs between individuals over time, but this changes depending on individual.
# sink("5_Msc_Stats/stats/ASV_stats.txt", append = TRUE)
# "================== NUMBER of shared ASVs between individuals (of same species) through time ==============="
# # IN general, overall decrease in shared ASVs between individuals over time, but this changes depending on individual.
# "Effect of time"
# "median"
# median(samps$beta[,1])
# "PD"
# min(c(sum(samps$beta[,1]>0)/length(samps$beta[,1]),sum(samps$beta[,1]<0)/length(samps$beta[,1])))
# "95% CI"
# c(ci(samps$beta[,1], method="HDI")$CI_low, ci(samps$beta[,1], method="HDI")$CI_high)
# 
# "Intercept: shared between individuals"
# "median"
# median(samps$alpha[,1])
# "PD"
# min(c(sum(samps$alpha[,1]>0)/length(samps$alpha[,1]),sum(samps$alpha[,1]<0)/length(samps$alpha[,1])))
# "95% CI"
# c(ci(samps$alpha[,1], method="HDI")$CI_low, ci(samps$alpha[,1], method="HDI")$CI_high)
# 
# sink()

allSharedASVs_byindiv_filt %>%
  unite(indivID_1, indivID_2, col="bothIndiv", remove=FALSE) %>%
  # filter(indivID_1 == "Anbo_2") %>%
  ggplot() + geom_point(aes(x=time, y=sharedASVs, col=sp_1)) +
  geom_line(aes(x=time, y=sharedFracASVs, group=bothIndiv))+
  facet_grid(.~sp_1)


indivList <- unique(mf_alt_filt_final$indivID)
if (!file.exists("5_Msc_Stats/allSharedASVs_throughtime.RData")) {
  allSharedASVs_throughtime <- data.frame()
  allSharedASVs_throughtimeAll <- data.frame()
  for ( id in indivList ) {

    dat_id_temp <- mf_alt_filt_final %>% filter(indivID == id) 
    timeList_temp <- sort(dat_id_temp$time)
    print(id)
    
    print("Timepoints within individual")
    total <- length(timeList_temp)^2
    pb <- txtProgressBar(min = 0, max = total, style = 3)
    i <- 1
    for ( t1 in timeList_temp ) {
      sid_1 <- dat_id_temp %>% filter(time==t1) %>% pull(SampleID)
      for ( t2 in timeList_temp ) {
        sid_2 <- dat_id_temp %>% filter(time == t2) %>% pull(SampleID)
        otu_temp <- otu_filt_mat[,c( sid_1, sid_2)]
        sharedFrac <- sum(apply(otu_temp,1,min))/colSums(otu_temp)[1]
        sharedASVs <- sum(rowSums(otu_temp>0)==2)
        totalASVs <- sum(otu_temp[,sid_1]>0)
        sharedFracASVs <- sharedASVs/totalASVs
        
        allSharedASVs_throughtime <- rbind(allSharedASVs_throughtime, data.frame(time_1 =t1, time_2 = t2, indivID = id, sharedFrac = sharedFrac, sharedASVs = sharedASVs, sharedFracASVs = sharedFracASVs, totalASVs= totalASVs))
        #------- PROGRESS BAR ------- #
        Sys.sleep(0.1)
        i <- i+1
        # update progress bar
        setTxtProgressBar(pb, i)
        #-----------------------------#
      }
    }
    
    print("Across all time")

      ind_all <- dat_id_temp %>% pull(SampleID)
      otu_temp <- otu_filt_mat[,c( ind_all)]
      sharedFrac <- sum(apply(otu_temp,1,min))/colSums(otu_temp)[1]
      sharedASVs <- sum(rowSums(otu_temp>0)==length(ind_all))
      totalASVs <- sum(otu_temp[,ind_all]>0)
      sharedFracASVs <- sharedASVs/totalASVs
      
      allSharedASVs_throughtimeAll <- rbind(allSharedASVs_throughtimeAll, data.frame(indivID =id,  sharedFrac = sharedFrac, sharedASVs = sharedASVs, sharedFracASVs = sharedFracASVs, totalASVs = totalASVs))

    
  }
  save(allSharedASVs_throughtime, file="5_Msc_Stats/allSharedASVs_throughtime.RData")
  save(allSharedASVs_throughtimeAll, file="5_Msc_Stats/allSharedASVs_throughtimeAll.RData")
  
} else {
  load("5_Msc_Stats/allSharedASVs_throughtime.RData")
  load("5_Msc_Stats/allSharedASVs_throughtimeAll.RData")
  
}

allSharedASVs_throughtimeAll %>% left_join(mf_collapsed_bd) %>% mutate(Treatment = ifelse(is.na(PABD),0,1)) %>%
  separate(indivID, into=c("species","indiv"), remove=FALSE, sep = "_") %>%
 ggplot() + geom_point(aes(x=species, y=sharedASVs, col=Treatment))
allSharedASVs_throughtimeAll %>%left_join(mf_collapsed_bd) %>% mutate(Treatment = ifelse(is.na(PABD),0,1)) %>%
  separate(indivID, into=c("species","indiv"), remove=FALSE, sep = "_") %>%
  ggplot() + geom_point(aes(x=species, y=sharedFracASVs, col=Treatment))
allSharedASVs_throughtimeAll %>%left_join(mf_collapsed_bd) %>% mutate(Treatment = ifelse(is.na(PABD),0,1)) %>%
  separate(indivID, into=c("species","indiv"), remove=FALSE, sep = "_") %>%
  ggplot() + geom_point(aes(x=species, y=sharedFrac, col=Treatment))


##### STAN: Individual compared to first time point ####
allSharedASVs_throughtime_adj <- allSharedASVs_throughtime %>% group_by(indivID) %>% mutate(minTime = min(time_2)) %>% ungroup() %>%
  filter(time_2 == minTime) %>%
  separate(indivID, into=c("species","indiv"), sep="_", remove=FALSE) %>%
  filter(time_1!=1) %>% left_join(mf_collapsed_bd) %>% mutate(Treatment = ifelse(is.na(PABD), 0,1))
allSharedASVs_throughtime_adj %>%
  ggplot(aes(x=time_1, y=sharedFracASVs, group=indivID, col=species)) + geom_point() +
  geom_line(aes( lty=factor(Treatment)))
gg_sharedFromStart <- allSharedASVs_throughtime_adj %>%
  ggplot(aes(x=time_1, y=sharedASVs, group=indivID, col=species)) + geom_point() +
  geom_line(aes( lty=factor(Treatment)))+ylab("Number of shared ASVs\nwith first time point") +xlab("Time point") +
  labs(col="Amphibian\nSpecies", lty="Exposed to \nBd or control")
gg_sharedFromStart
ggsave(filename="5_Msc_Stats/plots/ASV_sharedFromStart.png", height=4, width=6
       , gg_sharedFromStart)


stan_sharedASVSindivthroughtime <- stan_glmer(sharedASVs ~ time_1 + Treatment + (1|species) + (1|indivID), data = allSharedASVs_throughtime_adj)
plot(stan_sharedASVSindivthroughtime, pars = c("time_1","Treatment"))

samps <- rstan::extract(stan_sharedASVSindivthroughtime$stanfit)
#EFfect time
sink("5_Msc_Stats/stats/ASV_stats.txt", append = TRUE)
"================== Number of shared ASVs with first time point ==============="
# IN general, overall decrease in shared ASVs between individuals over time, but this changes depending on individual.
"Effect of time"
"median"
median(samps$beta[,1])
"PD"
min(c(sum(samps$beta[,1]>0)/length(samps$beta[,1]),sum(samps$beta[,1]<0)/length(samps$beta[,1])))
"95% CI"
c(ci(samps$beta[,1], method="HDI")$CI_low, ci(samps$beta[,1], method="HDI")$CI_high)

"Effect of treatment group (should be zero)"
"median"
median(samps$beta[,2])
"PD"
min(c(sum(samps$beta[,2]>0)/length(samps$beta[,2]),sum(samps$beta[,2]<0)/length(samps$beta[,2])))
"95% CI"
c(ci(samps$beta[,2], method="HDI")$CI_low, ci(samps$beta[,2], method="HDI")$CI_high)

"Intercept: shared from first time point"
"median"
median(samps$alpha[,1])
"PD"
min(c(sum(samps$alpha[,1]>0)/length(samps$alpha[,1]),sum(samps$alpha[,1]<0)/length(samps$alpha[,1])))
"95% CI"
c(ci(samps$alpha[,1], method="HDI")$CI_low, ci(samps$alpha[,1], method="HDI")$CI_high)

sink()

##### STAN: Individual compared to previous time point ####
allSharedASVs_throughtime_sumBychange <- allSharedASVs_throughtime%>% group_by(indivID) %>% mutate(minTime = min(time_2)) %>% ungroup() %>%
  filter(time_2 == time_1+1) %>%
  separate(indivID, into=c("species","indiv"), sep="_", remove=FALSE) %>% rename(time = time_2) %>%
  left_join(mf_toCombine) %>%
  mutate(prepost = factor(prepost, levels=c("Pre","Post")), Bd_exposure = factor(Bd_exposure, levels=c("Control","Bd-exposed")))

gg_betweentimepoints_ASV <- allSharedASVs_throughtime_sumBychange %>%
  # filter(time_1!=1) %>%
  ggplot(aes(x=time_1, y=sharedFracASVs, group=indivID, col=Bd_load)) + geom_point() +
  geom_line()+
  facet_wrap(.~species, nrow=1) +
  scale_color_gradient(low="grey",high="red") +  ylab("Fraction of ASVs\nthat persist between time points") +xlab("Time point") +
  labs(col=expression(paste(italic("Bd")," load (log+1)")))+
  theme_bw()
gg_betweentimepoints_ASV

ggsave(filename = "5_Msc_Stats/plots/ASV_fractionpersist_through_time.png",height=4, width=10
       , gg_betweentimepoints_ASV)

allSharedASVs_throughtime_sumBychange %>%
  # filter(time_1!=1) %>%
  ggplot(aes(x=time_1, y=sharedASVs, group=indivID, col=Bd_load)) + geom_point() +
  geom_line()+
  facet_wrap(.~species, nrow=1) +
  scale_color_gradient(low="grey",high="red") +
  theme_bw()

#### Exposure only

if ( !file.exists("5_Msc_Stats/stan_sharedFracASVs_byindiv_bd.RData")) {
  stan_sharedFracASVs_byindiv_bd <- stan_glmer(sharedFracASVs ~  time + Bd_exposure + prepost:Bd_exposure + prepost:Bd_exposure:PABD_single + prepost:Bd_exposure:PABD_single:Bd_load + (1|species) + (1|indivID), data=allSharedASVs_throughtime_sumBychange,
                                               adapt_delta = 0.999, iter=5000)
  save(stan_sharedFracASVs_byindiv_bd, file = "5_Msc_Stats/stan_sharedFracASVs_byindiv_bd.RData")
} else {
  load("5_Msc_Stats/stan_sharedFracASVs_byindiv_bd.RData")
}

plot(stan_sharedFracASVs_byindiv_bd, pars=c("time", "Bd_exposureBd-exposed","Bd_exposureControl:prepostPost", "Bd_exposureBd-exposed:prepostPost", "Bd_exposureBd-exposed:prepostPost:PABD_single", "Bd_exposureBd-exposed:prepostPost:PABD_single:Bd_load"))
samps <- rstan::extract(stan_sharedFracASVs_byindiv_bd$stanfit)
predsamps <- samps$beta
colnames(predsamps) <- names(fixef(stan_sharedFracASVs_byindiv_bd))[-1]

sink("5_Msc_Stats/stats/ASV_stats.txt", append = TRUE)
"================== Fraction of shared ASVs with adjacent time point ==============="
apply(predsamps, 2, function(x) c(median=median(x), lower95=ci(x,method="HDI")$CI_low, upper95=ci(x,method="HDI")$CI_high, pcomp=min(c(sum(x>0)/length(x),sum(x<0)/length(x)))))

sink()



if ( !file.exists("5_Msc_Stats/stan_sharedASVs_byindiv_bd.RData")) {
  stan_sharedASVs_byindiv_bd <- stan_glmer(sharedASVs ~  time + Bd_exposure + prepost:Bd_exposure + prepost:Bd_exposure:PABD_single + prepost:Bd_exposure:PABD_single:Bd_load + (1|species) + (1|indivID), data=allSharedASVs_throughtime_sumBychange
                                               , family=poisson(link="log")
                                               ,adapt_delta = 0.999, iter=5000)
  save(stan_sharedASVs_byindiv_bd, file = "5_Msc_Stats/stan_sharedASVs_byindiv_bd.RData")
} else {
  load("5_Msc_Stats/stan_sharedASVs_byindiv_bd.RData")
}

plot(stan_sharedASVs_byindiv_bd, pars=c("time", "Bd_exposureBd-exposed","Bd_exposureControl:prepostPost", "Bd_exposureBd-exposed:prepostPost", "Bd_exposureBd-exposed:prepostPost:PABD_single", "Bd_exposureBd-exposed:prepostPost:PABD_single:Bd_load"))


if ( !file.exists("5_Msc_Stats/stan_sharedFrac_byindiv_bd.RData")) {
  stan_sharedFrac_byindiv_bd <- stan_glmer(sharedFrac ~  time + Bd_exposure + prepost:Bd_exposure + prepost:Bd_exposure:PABD_single + prepost:Bd_exposure:PABD_single:Bd_load + (1|species) + (1|indivID), data=allSharedASVs_throughtime_sumBychange,
                                               adapt_delta = 0.999, iter=5000)
  save(stan_sharedFrac_byindiv_bd, file = "5_Msc_Stats/stan_sharedFrac_byindiv_bd.RData")
} else {
  load("5_Msc_Stats/stan_sharedFrac_byindiv_bd.RData")
}

plot(stan_sharedFrac_byindiv_bd, pars=c("time", "Bd_exposureBd-exposed","Bd_exposureControl:prepostPost", "Bd_exposureBd-exposed:prepostPost", "Bd_exposureBd-exposed:prepostPost:PABD_single", "Bd_exposureBd-exposed:prepostPost:PABD_single:Bd_load"))
samps <- rstan::extract(stan_sharedFrac_byindiv_bd$stanfit)
predsamps <- samps$beta
colnames(predsamps) <- names(fixef(stan_sharedFracASVs_byindiv_bd))[-1]

sink("5_Msc_Stats/stats/ASV_stats.txt", append = TRUE)
"================== Fraction of shared TOTAL FRACTION with adjacent time point ==============="
apply(predsamps, 2, function(x) c(median=median(x), lower95=ci(x,method="HDI")$CI_low, upper95=ci(x,method="HDI")$CI_high, pcomp=min(c(sum(x>0)/length(x),sum(x<0)/length(x)))))

sink()


# 
# samps <- rstan::extract(stan_sharedASVs_byindiv_bd$stanfit)
# predsamps <- samps$beta
# colnames(predsamps) <- names(fixef(stan_sharedASVs_byindiv_bd))[-1]

# sink("5_Msc_Stats/stats/ASV_stats.txt", append = TRUE)
# "================== Number of shared ASVs between adjacent time points ==============="
# apply(predsamps, 2, function(x) c(median=median(x), lower95=ci(x,method="HDI")$CI_low, upper95=ci(x,method="HDI")$CI_high, pcomp=min(c(sum(x>0)/length(x),sum(x<0)/length(x)))))
# 
# sink()

# stan_sharedASVs_byindiv_bd <- stan_glmer(sharedASVs ~ Bd_exposure + prepost:Bd_exposure + prepost:Bd_exposure:PABD_single + prepost:Bd_exposure:PABD_single:Bd_load + (1|species) + (1|indivID), data=allSharedASVs_throughtime_sumBychange
#                                          , family=poisson(link="identity")
#                                          ,adapt_delta = 0.999, iter=5000)
# plot(stan_sharedASVs_byindiv_bd, pars=c("Bd_exposureBd-exposed","Bd_exposureControl:prepostPost", "Bd_exposureBd-exposed:prepostPost", "Bd_exposureBd-exposed:prepostPost:PABD_single", "Bd_exposureBd-exposed:prepostPost:PABD_single:Bd_load"))
# 
# 
# stan_sharedFrac_byindiv_bd <- stan_glmer(sharedFrac ~ Bd_exposure + prepost:Bd_exposure + prepost:Bd_exposure:PABD_single + prepost:Bd_exposure:PABD_single:Bd_load + (1|species) + (1|indivID), data=allSharedASVs_throughtime_sumBychange
#                                          ,adapt_delta = 0.999, iter=5000)
# plot(stan_sharedFrac_byindiv_bd, pars=c("Bd_exposureBd-exposed","Bd_exposureControl:prepostPost", "Bd_exposureBd-exposed:prepostPost", "Bd_exposureBd-exposed:prepostPost:PABD_single", "Bd_exposureBd-exposed:prepostPost:PABD_single:Bd_load"))

#### Number of ASVs retained throughout all time points by individual ####
allSharedASVs_throughtimeAll %>%
  separate(indivID, into=c("species","indiv"), sep="_", remove=FALSE) %>%
  ggplot() + geom_point(aes(x=species, y=sharedFracASVs))
allSharedASVs_throughtimeAll %>%
  separate(indivID, into=c("species","indiv"), sep="_", remove=FALSE) %>%
  ggplot() + geom_point(aes(x=species, y=sharedASVs))


############### @~~~~~~~~~~~~ ASVs associated with infection ~~~~~~~~~~~~@ ############

#### Exposure vs non-exposure ####

total <- 5*2*2
pb <- txtProgressBar(min = 0, max = total, style = 3)
i <- 1
allPrevbyASV <- data.frame()
allReadsbyASV <- data.frame()
allVarbyASV <- data.frame()
# allPrevbyASV_infection <- data.frame()
# allReadsbyASV_infection <- data.frame()
# allVarbyASV_infection <- data.frame()
for (sp in unique(mf_alt_filt_final$species) ) {
  for ( pp in c("Pre","Post")) {
    for ( be in c("Bd-exposed","Control")) {
      mf_temp <- mf_alt_filt_final %>% filter(prepost==pp, Bd_exposure==be, species == sp)
      otu_temp <- otu_filt_mat[,mf_temp$SampleID]
      prevByASV <- apply(otu_temp, 1, function(x) (sum(x!=0)/length(x))) %>% as.data.frame() %>% rownames_to_column(var="ASVID") %>%
        mutate(species=sp, prepost=pp, Bd_exposure=be, Grouping="ppbe", PABD = NA) %>% rename_at(vars("."), ~"prevalence") %>% as_tibble()
      readsByASV <- apply(otu_temp, 1, function(x) sum(x)) %>% as.data.frame() %>% rownames_to_column(var="ASVID") %>%
        mutate(species=sp, prepost=pp, Bd_exposure=be, Grouping="ppbe", PABD = NA)%>% rename_at(vars("."), ~"reads") %>% as_tibble()
      varByASV <- apply(otu_temp, 1, function(x) sd(x)) %>% as.data.frame() %>% rownames_to_column(var="ASVID") %>%
        mutate(species=sp, prepost=pp, Bd_exposure=be, Grouping="ppbe", PABD = NA)%>% rename_at(vars("."), ~"sd") %>% as_tibble()
      # plot(prevByASV~log(abundByASV))
      # Save prevalence and abundance info
      allPrevbyASV <- rbind(allPrevbyASV, prevByASV)
      allReadsbyASV <- rbind(allReadsbyASV, readsByASV)
      allVarbyASV <- rbind(allVarbyASV, varByASV)
      if (be == "Bd-exposed" & pp == "Post") {
        for ( bd in c(1,0)) {
          mf_temp2 <- mf_temp %>% filter(prepost==pp, PABD==bd, species == sp, Bd_exposure==be)
          otu_temp2 <- otu_filt_mat[,mf_temp2$SampleID]
          prevByASV2 <- apply(otu_temp2, 1, function(x) (sum(x!=0)/length(x))) %>% as.data.frame() %>% rownames_to_column(var="ASVID") %>%
            mutate(species=sp, prepost=pp, Bd_exposure=be, Grouping="bd", PABD=bd) %>% rename_at(vars("."), ~"prevalence") %>% as_tibble()
          readsByASV2 <- apply(otu_temp2, 1, function(x) sum(x)) %>% as.data.frame() %>% rownames_to_column(var="ASVID") %>%
            mutate(species=sp, prepost=pp, Bd_exposure=be, Grouping="bd",PABD=bd)%>% rename_at(vars("."), ~"reads") %>% as_tibble()
          varByASV2 <- apply(otu_temp2, 1, function(x) sd(x)) %>% as.data.frame() %>% rownames_to_column(var="ASVID") %>%
            mutate(species=sp, prepost=pp, Bd_exposure=be, Grouping="bd",PABD=bd)%>% rename_at(vars("."), ~"sd") %>% as_tibble()
          # plot(prevByASV~log(abundByASV))
          # Save prevalence and abundance info
          allPrevbyASV <- rbind(allPrevbyASV, prevByASV2)
          allReadsbyASV <- rbind(allReadsbyASV, readsByASV2)
          allVarbyASV <- rbind(allVarbyASV, varByASV2)

        }
        #------- PROGRESS BAR ------- #
        Sys.sleep(0.1)
        i <- i+1
        # update progress bar
        setTxtProgressBar(pb, i)
        #-----------------------------#
      }
      
    }

  }
}

## Add taxonomy
toMergeTaxa <-allTaxa %>% rename(ASVID=Sequence) %>% select(ASVID, inhibitory, Taxa, g, f, o, c)
allReadsbyASV_adj <- allReadsbyASV %>%
  left_join(toMergeTaxa) %>%
  left_join(allVarbyASV) %>%
  left_join(allPrevbyASV)

allReadsbyASV_adj %>%
  group_by(prepost, Bd_exposure, PABD, species) %>% summarise(med_prev = median(prevalence[prevalence!=0], na.rm=TRUE), med_reads = median(reads[reads!=0], na.rm=TRUE), med_sd = median(sd[sd!=0], na.rm=TRUE))

allReadsbyASV_adj %>%
  filter(prevalence!=0) %>%
  unite(prepost, Bd_exposure, PABD, col="Group", sep="_", remove=FALSE) %>%
  ggplot() + geom_boxplot(aes(x=Group, y=log(reads))) +
  facet_wrap(.~species) + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

### Get top 5 ASVs for each group
set.seed(5734)
allReadsbyASV_adj %>%
  group_by(prepost, Bd_exposure, PABD, species) %>% 
  mutate(ASVrank = rank(-(reads), na.last = TRUE), relAbundReads = reads/sum(reads)) %>%
  ungroup() %>%
  filter(ASVrank<=5) %>%
  unite(prepost, Bd_exposure, PABD, col="Group", sep="_", remove=FALSE) %>%
  mutate(Group = factor(Group, levels=c("Pre_Control_NA", "Pre_Bd-exposed_NA", "Post_Control_NA", "Post_Bd-exposed_NA", "Post_Bd-exposed_0", "Post_Bd-exposed_1" ))) %>%
  mutate(species = factor(species, levels=c("Anbo","Rhma","Osse","Raca","Rapi"))) %>%
  ggplot() + geom_bar(aes(x=Group, y=relAbundReads, fill=Taxa), stat = "identity") +
  facet_grid(species~., drop=TRUE) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  scale_fill_manual(values=sample(colors()[-grep("grey|gray|white|cream",colors())]))

### Plot nicely
allReadsbyASV_adj_filt <- allReadsbyASV_adj %>%
  group_by(prepost, Bd_exposure, PABD, species) %>% 
  mutate(ASVrank = rank(-(reads), na.last = TRUE), relAbundReads = reads/sum(reads)) %>%
  ungroup() %>%
  filter(ASVrank<=5) %>%
  unite(prepost, Bd_exposure, PABD, col="Group", sep="_", remove=FALSE) %>%
  mutate(Group = factor(Group, levels=c("Pre_Control_NA", "Pre_Bd-exposed_NA", "Post_Control_NA", "Post_Bd-exposed_NA", "Post_Bd-exposed_0", "Post_Bd-exposed_1" ))) %>%
  mutate(species = factor(species, levels=c("Anbo","Rhma","Osse","Raca","Rapi"))) %>%
  mutate(Treatment = Bd_exposure, prepost = ifelse(prepost=="Pre", "Pre-Bd-exposure", "Post-Bd-exposure"))
set.seed(5734)
col_for_taxa <- sample(colors()[-grep("grey|gray|white|cream",colors())],size=length(unique(allReadsbyASV_adj_filt$Taxa)))
names(col_for_taxa) <- unique(allReadsbyASV_adj_filt$Taxa)
gg_forLegend <- g_legend(allReadsbyASV_adj_filt %>%
                            ggplot() + geom_bar(aes(x=Bd_exposure, y=relAbundReads, fill=Taxa), stat = "identity") +
                            scale_fill_manual(values=col_for_taxa) 

)

gg_pre <- allReadsbyASV_adj_filt %>%
  filter(prepost=="Pre-Bd-exposure") %>%
  ggplot() + geom_bar(aes(x=Treatment, y=relAbundReads, fill=Taxa), stat = "identity", show.legend = FALSE) +
  facet_grid(species~prepost) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  scale_fill_manual(values=col_for_taxa) + ylab("Total relative abundance of reads")
gg_post <- allReadsbyASV_adj_filt %>%
  filter(prepost=="Post-Bd-exposure", is.na(PABD)) %>%
  ggplot() + geom_bar(aes(x=Treatment, y=relAbundReads, fill=Taxa), stat = "identity", show.legend = FALSE) +
  facet_grid(species~prepost) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  scale_fill_manual(values=col_for_taxa) + ylab("Total relative abundance of reads")
gg_PABD <- allReadsbyASV_adj_filt %>%
  filter(prepost=="Post-Bd-exposure", !is.na(PABD)) %>%
  mutate(Infection = ifelse(PABD==1, "Bd-exposed\n(Infected)","Bd-exposed\n(Not infected)")) %>%
  ggplot() + geom_bar(aes(x=Infection, y=relAbundReads, fill=Taxa), stat = "identity", show.legend = FALSE) +
  facet_grid(species~Treatment, drop=TRUE) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  scale_fill_manual(values=col_for_taxa) + ylab("Total relative abundance of reads")
lay <- rbind(c(1,2,3,4,4,4)
             ,c(1,2,3,4,4,4)
             ,c(1,2,3,4,4,4)
             ,c(1,2,3,4,4,4)
)

gg_taxasummary <- grid.arrange(gg_pre, gg_post, gg_PABD, gg_forLegend, layout_matrix = lay)
ggsave(filename = "5_Msc_Stats/plots/taxasummary.png",height=6, width=13
       ,gg_taxasummary)

allReadsbyASV_adj %>%
  group_by(prepost, Bd_exposure, PABD, species) %>% 
  mutate(ASVrank = rank(-(reads), na.last = TRUE), relAbundReads = reads/sum(reads)) %>%
  ungroup() %>%
  filter(ASVrank<=5) %>%
  unite(prepost, Bd_exposure, PABD, col="Group", sep="_", remove=FALSE) %>%
  mutate(Group = factor(Group, levels=c("Pre_Control_NA", "Pre_Bd-exposed_NA", "Post_Control_NA", "Post_Bd-exposed_NA", "Post_Bd-exposed_0", "Post_Bd-exposed_1" ))) %>%
  mutate(species = factor(species, levels=c("Anbo","Rhma","Osse","Raca","Rapi"))) %>%
  ggplot() + geom_bar(aes(x=Group, y=relAbundReads, fill=factor(inhibitory)), stat = "identity") +
  facet_grid(species~., scales="free", drop=TRUE) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))# +
  # scale_fill_manual(values=sample(colors()[-grep("grey|gray|white|cream",colors())]))



allPrevbyASV %>% left_join(allReadsbyASV) %>%
  filter(Grouping=="bd") %>%
  ggplot() + geom_point(aes(x=prevalence, y=log(reads), col=species))  +
  facet_grid(PABD~.)

allReadsbyASV %>% 
  ggplot()

##### ASV composition by individual ######
mf_toCombined_withSampleID <- mf_alt_filt_final%>%select(SampleID, species, indivID, time) %>% left_join(mf_toCombine)

# Top 3
top3ASV <- otu_filt %>% rename(ASVID="#OTU ID") %>% pivot_longer(-ASVID, names_to = "SampleID", values_to = "Reads") %>%
  group_by(SampleID) %>% mutate(rank = rank(-Reads)) %>% ungroup() %>% filter(rank<=3) %>%
  pull(ASVID) %>% unique()

mf_and_otu_dominant <- otu_filt %>% rename(ASVID="#OTU ID") %>% pivot_longer(-ASVID, names_to = "SampleID", values_to = "Reads") %>%
  filter(ASVID %in% top3ASV) %>%
  left_join(allTaxa %>% rename(ASVID=Sequence)) %>%
  left_join(mf_toCombined_withSampleID) %>%
  mutate(RelAbund = Reads/10000)
mf_and_otu_dominant %>%
  ggplot() + geom_bar(aes(x=factor(time), y=RelAbund, fill=g, group=indivID), stat="identity", position = position_dodge()) +
  facet_grid(species~Bd_exposure) +
  geom_vline(aes(xintercept=5.5), col="orange")

mf_and_otu_dominant %>%
  filter(species=="Anbo") %>%
  ggplot() + geom_bar(aes(x=factor(time), y=RelAbund, fill=g), stat="identity", show.legend = FALSE) +
  facet_grid(indivID~.) +
  geom_point(aes(x=factor(time), y=1, col = Bd_load)) +
  scale_color_gradient(low="white", high="red")
              
mf_and_otu_dominant %>%
  filter(species=="Rhma") %>%
  ggplot() + geom_bar(aes(x=factor(time), y=RelAbund, fill=g), stat="identity", show.legend = FALSE) +
  facet_grid(indivID~.) +
  geom_point(aes(x=factor(time), y=1, col = Bd_load)) +
  scale_color_gradient(low="white", high="red")

mf_and_otu_dominant %>%
  filter(species=="Osse") %>%
  ggplot() + geom_bar(aes(x=factor(time), y=RelAbund, fill=g), stat="identity", show.legend = FALSE) +
  facet_grid(indivID~.) +
  geom_point(aes(x=factor(time), y=1, col = Bd_load)) +
  scale_color_gradient(low="white", high="red")
