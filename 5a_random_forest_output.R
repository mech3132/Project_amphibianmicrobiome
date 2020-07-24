#!/bin/bash Rscript

#### Loading data #####
library(tidyverse)
dir.create("5a_random_forest_output")

# Prevalence
load("./5_random_forest_prevalence/infect_test_train_comparison_ponly.RData")
load("./5_random_forest_prevalence/infect_test_train_comparison_ponly2.RData")
load("./5_random_forest_prevalence/infect_test_train_comparison_nop.RData")
load("./5_random_forest_prevalence/infect_test_train_comparison_nop2.RData")
load("./5_random_forest_prevalence/infect_test_train_comparison_withp.RData")
load("./5_random_forest_prevalence/infect_test_train_comparison_withp2.RData")
load("./5_random_forest_prevalence/PABD_test_train_comparison_ponly.RData")
load("./5_random_forest_prevalence/PABD_test_train_comparison_ponly2.RData")
load("./5_random_forest_prevalence/PABD_test_train_comparison_nop.RData")
load("./5_random_forest_prevalence/PABD_test_train_comparison_nop2.RData")
load("./5_random_forest_prevalence/PABD_test_train_comparison_withp.RData")
load("./5_random_forest_prevalence/PABD_test_train_comparison_withp2.RData")

infect_test_train_comparison_nop_prev <- infect_test_train_comparison_nop
infect_test_train_comparison_nop_prev2 <- infect_test_train_comparison_nop2
infect_test_train_comparison_withp_prev <- infect_test_train_comparison_withp
infect_test_train_comparison_withp_prev2 <- infect_test_train_comparison_withp2

PABD_test_train_comparison_nop_prev <- PABD_test_train_comparison_nop
PABD_test_train_comparison_nop_prev2 <- PABD_test_train_comparison_nop2
PABD_test_train_comparison_withp_prev <- PABD_test_train_comparison_withp
PABD_test_train_comparison_withp_prev2 <- PABD_test_train_comparison_withp2

# Count
load("./5_random_forest_count/infect_test_train_comparison_nop.RData")
load("./5_random_forest_count/infect_test_train_comparison_nop2.RData")
load("./5_random_forest_count/infect_test_train_comparison_withp.RData")
load("./5_random_forest_count/infect_test_train_comparison_withp2.RData")
load("./5_random_forest_count/PABD_test_train_comparison_nop.RData")
load("./5_random_forest_count/PABD_test_train_comparison_nop2.RData")
load("./5_random_forest_count/PABD_test_train_comparison_withp.RData")
load("./5_random_forest_count/PABD_test_train_comparison_withp2.RData")

infect_test_train_comparison_nop_count <- infect_test_train_comparison_nop
infect_test_train_comparison_nop_count2 <- infect_test_train_comparison_nop2
infect_test_train_comparison_withp_count <- infect_test_train_comparison_withp
infect_test_train_comparison_withp_count2 <- infect_test_train_comparison_withp2

PABD_test_train_comparison_nop_count <- PABD_test_train_comparison_nop
PABD_test_train_comparison_nop_count2 <- PABD_test_train_comparison_nop2
PABD_test_train_comparison_withp_count <- PABD_test_train_comparison_withp
PABD_test_train_comparison_withp_count2 <- PABD_test_train_comparison_withp2

load("./5_random_forest_prevalence/importance_ponly_PABD.RData")
load("./5_random_forest_prevalence/importance_ponly_PABD2.RData")
load("./5_random_forest_prevalence/importance_ponly.RData")
load("./5_random_forest_prevalence/importance_ponly2.RData")

remove(infect_test_train_comparison_nop)
remove(infect_test_train_comparison_nop2)
remove(infect_test_train_comparison_withp)
remove(infect_test_train_comparison_withp2)
remove(PABD_test_train_comparison_nop)
remove(PABD_test_train_comparison_nop2)
remove(PABD_test_train_comparison_withp)
remove(PABD_test_train_comparison_withp2)


load("./5_random_forest_prevalence/importance_ponly_PABD.RData")
load("./5_random_forest_prevalence/importance_ponly_PABD2.RData")
load("./5_random_forest_prevalence/importance_ponly.RData")
load("./5_random_forest_prevalence/importance_ponly2.RData")

#### Combine all datasets and save ####
all_test_train_comparisons_infect <- rbind(
  cbind(infect_test_train_comparison_ponly, model="ponly",otu_type=NA)
  ,cbind(infect_test_train_comparison_ponly2, model="ponly",otu_type=NA)
  ,cbind(infect_test_train_comparison_nop_prev, model="nop", otu_type="prev")
  ,cbind(infect_test_train_comparison_nop_prev2, model="nop", otu_type="prev")
  ,cbind(infect_test_train_comparison_nop_count, model="nop", otu_type="count")
  ,cbind(infect_test_train_comparison_nop_count2, model="nop", otu_type="count")
  ,cbind(infect_test_train_comparison_withp_prev, model="withp", otu_type="prev")
  ,cbind(infect_test_train_comparison_withp_prev2, model="withp", otu_type="prev")
  ,cbind(infect_test_train_comparison_withp_count, model="withp", otu_type="count")
  ,cbind(infect_test_train_comparison_withp_count2, model="withp", otu_type="count")
)
save(all_test_train_comparisons_infect, file="./5a_random_forest_output/all_test_train_comparisons_infect.RData")

all_test_train_comparisons_PABD <- rbind(
  cbind(PABD_test_train_comparison_ponly, model="ponly",otu_type=NA)
  ,cbind(PABD_test_train_comparison_ponly2, model="ponly",otu_type=NA)
  ,cbind(PABD_test_train_comparison_nop_prev, model="nop", otu_type="prev")
  ,cbind(PABD_test_train_comparison_nop_prev2, model="nop", otu_type="prev")
  ,cbind(PABD_test_train_comparison_nop_count, model="nop", otu_type="count")
  ,cbind(PABD_test_train_comparison_nop_count2, model="nop", otu_type="count")
  ,cbind(PABD_test_train_comparison_withp_prev, model="withp", otu_type="prev")
  ,cbind(PABD_test_train_comparison_withp_prev2, model="withp", otu_type="prev")
  ,cbind(PABD_test_train_comparison_withp_count, model="withp", otu_type="count")
  ,cbind(PABD_test_train_comparison_withp_count2, model="withp", otu_type="count")
)
save(all_test_train_comparisons_PABD, file="./5a_random_forest_output/all_test_train_comparisons_PABD.RData")

all_importance_infect <- rbind(importance_ponly, importance_ponly2)
all_importance_PABD <- rbind(importance_ponly_PABD, importance_ponly_PABD2)
save(all_importance_infect, file="./5a_random_forest_output/all_importance_infect.RData")
save(all_importance_PABD, file="./5a_random_forest_output/all_importance_PABD.RData")

#### Calculating errors infect ####
summary_error_infect <- all_test_train_comparisons_infect %>%
  unite(model, otu_type, col = model_otu_type, remove=FALSE) %>%
  group_by(model, otu_type,model_otu_type, rep) %>%
  mutate(ave_obs=mean(Test_obs)) %>%
  ungroup() %>%
  mutate(ressq = (Test_pred-Test_obs)^2, varsq = (Test_obs-ave_obs)^2) %>%
  group_by(model, otu_type, model_otu_type,rep) %>%
  summarize(MSE=mean(ressq), MSE_res=sum(ressq), MSE_total=sum(varsq), R2 = 1-(sum(ressq)/sum(varsq))) %>%
  ungroup() %>%
  mutate(model_otu_type = factor(model_otu_type, levels=c("ponly_NA","nop_prev","nop_count","withp_prev","withp_count")))
save(summary_error_infect, file="./5a_random_forest_output/summary_error_infect.RData")

summary_error_infect_MSER2 <- summary_error_infect %>%
  group_by(model_otu_type) %>%
  summarize(meanMSE=mean(MSE), medianMSE=median(MSE), varMSE=var(MSE),meanR2=mean(R2), medianR2=median(R2), varR2=var(R2)) 
save(summary_error_infect_MSER2, file="./5a_random_forest_output/summary_error_infect_MSER2.RData")

summary_indiv_infect <- all_test_train_comparisons_infect %>%
  unite(model, otu_type, col = model_otu_type, remove=FALSE) %>%
  group_by(indivID, model, otu_type, model_otu_type, Test_obs) %>%
  mutate(meanTestPred=mean(Test_pred), varTestPred=var(Test_pred))  %>%
  ungroup() %>%
  mutate(model_otu_type = factor(model_otu_type, levels=c("ponly_NA","nop_prev","nop_count","withp_prev","withp_count")))
save(summary_indiv_infect, file="./5a_random_forest_output/summary_indiv_infect.RData")

#### Calcuating errors PABD ####

summary_error_PABD <- all_test_train_comparisons_PABD %>% 
  unite(model, otu_type, col = model_otu_type, remove=FALSE) %>%
  group_by(model, otu_type,model_otu_type, rep) %>%
  summarize(propCorrect = mean(Test_pred==Test_obs)) %>%
  ungroup() %>%
  mutate(model_otu_type = factor(model_otu_type, levels=c("ponly_NA","nop_prev","nop_count","withp_prev","withp_count")))
save(summary_error_PABD, file="./5a_random_forest_output/summary_error_PABD.RData")


summary_indiv_PABD <- all_test_train_comparisons_PABD %>% 
  unite(model, otu_type, col = model_otu_type, remove=FALSE) %>%
  group_by(model, otu_type,model_otu_type, indivID) %>%
  summarize(propCorrect = mean(Test_pred==Test_obs), Observed=unique(Test_obs)) %>%
  ungroup() %>%
  mutate(model_otu_type = factor(model_otu_type, levels=c("ponly_NA","nop_prev","nop_count","withp_prev","withp_count")))
save(summary_indiv_PABD, file="./5a_random_forest_output/summary_indiv_PABD.RData")



summary_indiv_PABD <- all_test_train_comparisons_PABD %>% 
  unite(model, otu_type, col = model_otu_type, remove=FALSE) %>%
  group_by(model, otu_type,model_otu_type, indivID) %>%
  summarize(propCorrect = mean(Test_pred==Test_obs), Observed=unique(Test_obs), infectRatio=(sum(Test_pred=="INFECTED")/n())) %>%
  ungroup() %>%
  mutate(model_otu_type = factor(model_otu_type, levels=c("ponly_NA","nop_prev","nop_count","withp_prev","withp_count")))
save(summary_indiv_PABD, file="./5a_random_forest_output/summary_indiv_PABD.RData")



#### Plotting infect ####

all_test_train_comparisons_infect %>%
  unite(model, otu_type, col = model_otu_type, remove=FALSE) %>%
  ggplot(aes(x=Test_obs, y=Test_pred)) +geom_jitter(aes(col=model_otu_type ), height=0, width=0.2) +
  geom_smooth(aes(col=model_otu_type), method="lm") +
  geom_abline(aes(intercept=0, slope=1)) +
  ylim(0,8) + xlim(0,8)

ggsave("./5a_random_forest_output/RF_MSE_sim.pdf"
  ,summary_error_infect %>%
    ggplot() + geom_violin(aes(x=model_otu_type, y=MSE)) +
    geom_jitter(aes(x=model_otu_type, y=MSE), width = 0.25, height=0) +
    geom_point(data=summary_error_infect_MSER2, aes(x=model_otu_type, y=meanMSE), col="red", cex=3) +
    geom_segment(data=summary_error_infect_MSER2, aes(x=model_otu_type, xend=model_otu_type, y=meanMSE-sqrt(varMSE), yend=meanMSE+sqrt(varMSE)), col="red")
  
)
summary_error_infect %>%
  ggplot() + geom_violin(aes(x=model_otu_type, y=R2)) +
  geom_jitter(aes(x=model_otu_type, y=R2), width = 0.25, height=0) +
  geom_point(data=summary_error_infect_MSER2, aes(x=model_otu_type, y=meanR2), col="red", cex=3) +
  geom_segment(data=summary_error_infect_MSER2, aes(x=model_otu_type, xend=model_otu_type, y=meanR2-sqrt(varMSE), yend=meanR2+sqrt(varMSE)), col="red")

summary_indiv_infect %>%
  ggplot(aes(x=Test_obs, y=meanTestPred)) +geom_point(aes(col=model_otu_type )) +
  # geom_jitter(data=all_test_train_comparisons_infect,aes(x=Test_obs, y=Test_pred), alpha=0.2, height=0, width=0.2) +
  geom_smooth(aes(col=model_otu_type), method="lm") +
  geom_segment(aes(x=Test_obs, xend=Test_obs, y=meanTestPred-sqrt(varTestPred), yend=meanTestPred+sqrt(varTestPred), col=model_otu_type)) +
  geom_abline(aes(intercept=0, slope=1)) +
  ylim(0,8) + xlim(0,8)

ggsave(filename="./5a_random_forest_output/RF_allmodels_bdload.pdf"
       , all_test_train_comparisons_infect %>%
         unite(model, otu_type, col = model_otu_type, remove=FALSE) %>%
         ggplot(aes(x=Test_obs, y=Test_pred)) +
         geom_jitter(aes(col=model_otu_type ),alpha=0.2, height=0, width=0.2) +
         geom_smooth(aes(col=model_otu_type), method="lm") +
         geom_segment(data=summary_indiv_infect, aes(x=Test_obs, xend=Test_obs, y=meanTestPred-sqrt(varTestPred), yend=meanTestPred+sqrt(varTestPred), col=model_otu_type)) +
         geom_point(data=summary_indiv_infect,aes(x=Test_obs, y=meanTestPred, col=model_otu_type)) +
         geom_abline(aes(intercept=0, slope=1)) +
         ylim(0,8) + xlim(0,8)
       )

## Plotting only certain models
all_test_train_comparisons_infect_filt <- all_test_train_comparisons_infect %>%
  filter(model %in% c("ponly","nop"), otu_type %in% c("prev",NA))
summary_indiv_infect_filt <- summary_indiv_infect %>%
  filter(model %in% c("ponly","nop"), otu_type %in% c("prev",NA)) 
ggsave(filename="./5a_random_forest_output/RF_allmodels_bdload_filt.pdf"
       , all_test_train_comparisons_infect_filt %>%
         unite(model, otu_type, col = model_otu_type, remove=FALSE) %>%
         ggplot(aes(x=Test_obs, y=Test_pred)) +
         geom_jitter(aes(col=model_otu_type ),alpha=0.2, height=0, width=0.2) +
         geom_smooth(aes(col=model_otu_type), method="lm") +
         geom_segment(data=summary_indiv_infect_filt, aes(x=Test_obs, xend=Test_obs, y=meanTestPred-sqrt(varTestPred), yend=meanTestPred+sqrt(varTestPred), col=model_otu_type)) +
         geom_point(data=summary_indiv_infect_filt,aes(x=Test_obs, y=meanTestPred, col=model_otu_type)) +
         geom_abline(aes(intercept=0, slope=1)) +
         ylim(0,8) + xlim(0,8)
)

# By species?

ggsave("./5a_random_forest_output/RF_residuals_by_individual.pdf"
      ,all_test_train_comparisons_infect %>%
        group_by(indivID, model, otu_type, rep) %>%
        summarize(ressq = (Test_pred-Test_obs)^2) %>%
        ungroup() %>%
        separate(indivID, into=c("species","indiv"), remove=FALSE) %>%
        arrange(species, indiv) %>%
        mutate(indivID =factor(indivID, levels=unique(indivID))) %>%
        ggplot() + 
        geom_jitter(aes(x=indivID, y=sqrt(ressq), col=model), height=0, width=0.2) +
        facet_grid(otu_type~.)+
        theme(axis.text.x = element_text(angle=90))
      )

summary_indiv_infect %>%
  separate(indivID, into=c("species","indiv"), remove=FALSE) %>%
  arrange(species, indiv) %>%
  mutate(indivID = factor(indivID, levels=unique(indivID)), species=factor(species, levels=c("Anbo","Rhma","Osse","Raca","Rapi"))) %>%
  ggplot(aes(x=Test_obs, y=meanTestPred)) +geom_point(aes(pch=indiv,col=species)) +
  # geom_jitter(data=all_test_train_comparisons_infect,aes(x=Test_obs, y=Test_pred), alpha=0.2, height=0, width=0.2) +
  # geom_smooth(aes(col=model_otu_type), method="lm") +
  geom_segment(aes(x=Test_obs, xend=Test_obs, y=meanTestPred-sqrt(varTestPred), yend=meanTestPred+sqrt(varTestPred), col=species)) +
  geom_abline(aes(intercept=0, slope=1)) +
  ylim(0,8) + xlim(0,8)  +
  facet_wrap(.~model_otu_type)



#### Plotting PABD ####
summary_error_PABD_mean <- summary_error_PABD %>%
  group_by(model_otu_type) %>%
  summarize(meanPropCorrect=mean(propCorrect), sd=sd(propCorrect))
ggsave("./5a_random_forest_output/RF_propCorrect_sim.pdf"
       , summary_error_PABD %>%
         ggplot() + geom_violin(aes(x=model_otu_type, y=propCorrect)) +
         geom_jitter(aes(x=model_otu_type, y=propCorrect), height=0.05, width=0.25) +
         geom_point(data=summary_error_PABD_mean, aes(x=model_otu_type, y=meanPropCorrect), col="red") +
         geom_segment(data=summary_error_PABD_mean, aes(x=model_otu_type, xend=model_otu_type, y=meanPropCorrect-sd, yend=meanPropCorrect+sd), col="red")
       
       )
ggsave("./5a_random_forest_output/RF_error_by_individual.pdf", height=3, width=10
       ,summary_indiv_PABD %>%
         separate(indivID, into=c("species","indiv")) %>%
         ggplot()  +
         geom_jitter(aes(x=species, y=propCorrect,col=indiv, pch=Observed), height=0.0, width=0.25, cex=2) +
         facet_grid(.~model_otu_type) +
         theme(axis.text.x=element_text(angle=90))
       
       )

summary_indiv_PABD %>%
  filter(model=="ponly") %>%
  separate(indivID, into=c("species","indiv")) %>%
  dplyr::select(species) %>%table()


#### Look at importance ####


rbind(importance_ponly_PABD, importance_ponly_PABD2) %>%
  group_by(taxa) %>%
  summarize(infectedMean=mean(INFECTED), infectedSD=sd(INFECTED)
            , notinfectMean=mean(NOT_INFECTED), notinfectSD=sd(NOT_INFECTED)
            , meanMean=mean(MeanDecreaseAccuracy), meanSD=sd(MeanDecreaseAccuracy)) %>%
  arrange(-meanMean)

rbind(importance_ponly, importance_ponly2) %>%
  group_by(taxa) %>%
  summarize( meanIncMSE=mean(X.IncMSE), sdIncMSE=sd(X.IncMSE)) %>%
  arrange(-meanIncMSE)

