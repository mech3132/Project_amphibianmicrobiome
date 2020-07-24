#!/bin/bash Rscript

#### Loading data #####
library(randomForest)
library(tidyverse)
load("3_5sp_mapping_otu_downstream/mf_alt_filt_final.RData")
load("3_5sp_mapping_otu_downstream/otu_filt.RData")
load("3_5sp_mapping_otu_downstream/taxonomy.RData")
load("4_Bayesian_models/all_p.RData")
source("msc_randomforest_codebits.R")
dir.create("5_random_forest_2")

###### Create master file of predictors #########
otu <- otu_filt
mf <- mf_alt_filt_final %>%
  select(SampleID, prepost, indivID, Bd_exposure) 
p <- all_p %>%
  select(-c(grep("exp_*", colnames(all_p)))) 
p_infectonly <- p %>%
  select(-c(grep("p_*", colnames(p))))

# Change OTU names to nice ones
Taxa <- getTaxa(t=taxonomy)
# otu$`#OTU ID`<- Taxa[match(otu_filt$`#OTU ID`, Taxa$Sequence), "Taxa"]
otu <- otu_filt

# Make OTU into ranks
otu_rank <- count_to_percentileranking(OTU=otu)
# Make OTU into presence absence
otu_PA <- count_to_presenceabsence(OTU=otu)
# Just abundance
otu_count <- otu[,-which(colnames(otu)=="#OTU ID")]
rownames(otu_count) <- otu$`#OTU ID`

# Add OTU table to mf 
mf_rank <- mf_otu_combine(otu_rank, mf)
mf_PA <- mf_otu_combine(otu_PA, mf)
mf_count <- mf_otu_combine(otu_count, mf)


#--------- Random forest data setup -------------#
all_RF_predictBD <- list()
RERUN_RF <- FALSE
if ( RERUN_RF ) {
  total <- 3*3*2*7*1
  # create progress bar
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  i <- 1
  for ( otu_type in c("rank","PA","count") ) {
    for ( p_type in c("onlyp","withp","nop") ) {
      for ( bd_type in c("infect","PABD")) {
        for ( prop in c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)) {
          for ( repl in c(1) ) {
            all_RF_predictBD[[paste(otu_type, p_type, bd_type, prop, repl, sep="_")]] <- list()
            # Get a formatted mapping file
            mf_pred_temp_bd <- get_mf_pred_temp(p_type=p_type, otu_type = otu_type, bd_type=bd_type, p=p, mf_temp=mf_temp)
            
            # Mutate infection depending on question
            Error_metric <- ifelse(bd_type=="PABD","MeanDecreaseAccuracy","X.IncMSE" )
            
            # Pick training and validation set
            set.seed(423*repl)
            r <- sample(nrow(mf_pred_temp_bd), prop*nrow(mf_pred_temp_bd), replace=FALSE)
            trainSet <- mf_pred_temp_bd[r,]
            testSet <- mf_pred_temp_bd[-r,]
            # Run random forest model
            RF <- randomForest(response ~ ., data=trainSet,importance=TRUE)
            # Run prediction
            testSet_predictions <- predict(RF, testSet, type="class")
            # Get error of predicted set
            error <- ifelse(bd_type == "PABD", mean(testSet_predictions==testSet$response), mean((testSet_predictions-testSet$response)^2))
            # importance of factors
            importance_RF <- data.frame(importance(RF)) %>%
              mutate(taxa=rownames(importance(RF))) %>%
              select(taxa, everything()) %>%
              arrange(-get(paste0(Error_metric)))
            
            # Add inhibitory to importance
            # Match inhibitory to not
            importance_RF[,c("taxonomy","inhibitory")] <- Taxa[match(importance_RF$taxa, Taxa$Sequence), c("Taxa","inhibitory")]
            importance_RF <- importance_RF %>%
              mutate(taxonomy=ifelse(is.na(taxonomy), taxa, taxonomy))
            
            # Save data
            all_RF_predictBD[[paste(otu_type, p_type, bd_type, prop,repl, sep="_")]][["mf_used"]] <- mf_pred_temp_bd
            all_RF_predictBD[[paste(otu_type, p_type, bd_type, prop,repl, sep="_")]][["RF"]] <- RF
            all_RF_predictBD[[paste(otu_type, p_type, bd_type, prop,repl, sep="_")]][["importance"]] <- importance_RF
            all_RF_predictBD[[paste(otu_type, p_type, bd_type, prop,repl, sep="_")]][["testSetError"]] <- error
            all_RF_predictBD[[paste(otu_type, p_type, bd_type, prop,repl, sep="_")]][["test_train_comparison"]] <- data.frame(Test_pred=testSet_predictions, Test_obs=testSet$response)
            all_RF_predictBD[[paste(otu_type, p_type, bd_type, prop,repl, sep="_")]][["R2"]] <- getR2(Test_pred=testSet_predictions, Test_obs=testSet$response)
            
            # Insert into data summary plot
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
  }
  close(pb)
  save(all_RF_predictBD, file="./5_random_forest_2/all_RF_predictBD.RData")
}
load("./5_random_forest_2/all_RF_predictBD.RData")

# What is the best model?
## Extract validation R2 from infect models
R2_summary <-rbind(extractError(all_RF_predictBD, bd_type="infect", otu_type = "rank", error_type = "R2")
      , extractError(all_RF_predictBD, bd_type="infect", otu_type = "count", error_type = "R2")
      , extractError(all_RF_predictBD, bd_type="infect", otu_type = "PA", error_type = "R2")
      
)  %>%
  mutate(error_type="% Unexplained variance (load)"
         , onlyp=1-onlyp
         , nop=1-nop
         , withp=1-withp)
MSE_summary <-rbind(extractError(all_RF_predictBD, bd_type="infect", otu_type = "rank", error_type = "MSE")
                   , extractError(all_RF_predictBD, bd_type="infect", otu_type = "count", error_type = "MSE")
                   , extractError(all_RF_predictBD, bd_type="infect", otu_type = "PA", error_type = "MSE")
                   
)%>%
  mutate(error_type="MSE (load)")
errorRate_summary <- rbind(extractError(all_RF_predictBD, bd_type="PABD", otu_type = "rank", error_type = "errorRate")
                           , extractError(all_RF_predictBD, bd_type="PABD", otu_type = "count", error_type = "errorRate")
                           , extractError(all_RF_predictBD, bd_type="PABD", otu_type = "PA", error_type = "errorRate")
                           
)%>%
  mutate(error_type="errorRate (pres/absence)")

ggsave("5_random_forest_2/model_diagnostics.pdf", height=8, width=10,
       rbind(R2_summary, MSE_summary, errorRate_summary) %>%
         gather(onlyp, nop, withp, key=TrainingParameters, value=value) %>%
         mutate(TrainingParameters=ifelse(TrainingParameters=="onlyp","Community metrics only",ifelse(TrainingParameters=="nop","ASV metrics only","Both community and ASV metrics")))%>%
         mutate(TrainingParameters=factor(TrainingParameters, levels=c("Community metrics only","ASV metrics only","Both community and ASV metrics"))) %>%
         rename(Dataset=type) %>%
         mutate(otu_type = ifelse(otu_type=="rank","ASV rank",ifelse(otu_type=="count","ASV frequency (rarefied)", "ASV presence/absence"))) %>%
         mutate(otu_type = factor(otu_type, levels=c("ASV presence/absence","ASV rank","ASV frequency (rarefied)"))) %>%
         mutate(error_type=factor(error_type, levels=c("errorRate (pres/absence)","MSE (load)","% Unexplained variance (load)"))) %>%
         ggplot() + geom_line(aes(x=prop_train, y=value, col=TrainingParameters, lty=Dataset)) + facet_grid(error_type~otu_type, scales = "free_y", switch = "y")+
         xlab("Proportion used in training set")+ylab("")
       
       )


rbind(extractError(all_RF_predictBD, bd_type="infect", otu_type = "rank", error_type = "R2")
      , extractError(all_RF_predictBD, bd_type="infect", otu_type = "count", error_type = "R2")
      , extractError(all_RF_predictBD, bd_type="infect", otu_type = "PA", error_type = "R2")
) %>%
  gather(onlyp, nop, withp, key=Model, value=R2) %>%
  mutate(Model=factor(Model, levels=c("onlyp","nop","withp"))) %>%
  ggplot() + geom_line(aes(x=prop_train, y=1-R2, col=Model, lty=type)) + facet_wrap(.~otu_type, ncol=1)

rbind(extractError(all_RF_predictBD, bd_type="infect", otu_type = "rank", error_type = "MSE")
      , extractError(all_RF_predictBD, bd_type="infect", otu_type = "count", error_type = "MSE")
      , extractError(all_RF_predictBD, bd_type="infect", otu_type = "PA", error_type = "MSE")
) %>%
  gather(onlyp, nop, withp, key=Model, value=MSE) %>%
  mutate(Model=factor(Model, levels=c("onlyp","nop","withp"))) %>%
  ggplot() + geom_line(aes(x=prop_train, y=MSE, col=Model, lty=type)) + facet_wrap(.~otu_type, ncol=1)


rbind(extractError(all_RF_predictBD, bd_type="PABD", otu_type = "rank", error_type = "errorRate")
      , extractError(all_RF_predictBD, bd_type="PABD", otu_type = "count", error_type = "errorRate")
      , extractError(all_RF_predictBD, bd_type="PABD", otu_type = "PA", error_type = "errorRate")
) %>%
  gather(onlyp, nop, withp, key=Model, value=errorRate) %>%
  mutate(Model=factor(Model, levels=c("onlyp","nop","withp"))) %>%
  ggplot() + geom_line(aes(x=prop_train, y=errorRate, col=Model, lty=type)) + facet_wrap(.~otu_type, ncol=1)


### Are parameter importance constant across?

extractImportance(all_RF_predictBD, bd_type="infect", otu_type="rank", top=3) %>%
  ggplot() +geom_line(aes(x=prop, y=importance, col=taxonomy), show.legend = FALSE) +facet_wrap(.~model)

extractImportance(all_RF_predictBD, bd_type="infect", otu_type="count", top=3) %>%
  ggplot() +geom_line(aes(x=prop, y=importance, col=taxonomy), show.legend = FALSE) +facet_wrap(.~model)

extractImportance(all_RF_predictBD, bd_type="infect", otu_type="PA", top=3) %>%
  ggplot() +geom_line(aes(x=prop, y=importance, col=taxonomy), show.legend = FALSE) +facet_wrap(.~model)

extractImportance(all_RF_predictBD, bd_type="PABD", otu_type="rank", top=3) %>%
  ggplot() +geom_line(aes(x=prop, y=importance, col=taxonomy), show.legend = FALSE) +facet_wrap(.~model)

extractImportance(all_RF_predictBD, bd_type="PABD", otu_type="count", top=3) %>%
  ggplot() +geom_line(aes(x=prop, y=importance, col=taxonomy), show.legend = FALSE) +facet_wrap(.~model)

extractImportance(all_RF_predictBD, bd_type="PABD", otu_type="PA", top=3) %>%
  ggplot() +geom_line(aes(x=prop, y=importance, col=taxonomy), show.legend = FALSE) +facet_wrap(.~model)

## Extract accuracy from PABD models

write.table(make_summary_RF(all_output = all_RF_predictBD, otu_type = "rank", bd_type = "infect", prop=0.8, repl = 1)
            , quote = FALSE
            , sep = "\t"
            , row.names = FALSE
            , col.names=TRUE
            , file="./5_random_forest_2/RF_rank_infect.txt")
write.table(make_summary_RF(all_output = all_RF_predictBD, otu_type = "rank", bd_type = "PABD", prop=0.8, repl = 1)
            , quote = FALSE
            , sep = "\t"
            , row.names = FALSE
            , col.names=TRUE
            , file="./5_random_forest_2/RF_rank_PABD.txt")
write.table(make_summary_RF(all_output = all_RF_predictBD, otu_type = "PA", bd_type = "infect", prop=0.8, repl = 1)
            , quote = FALSE
            , sep = "\t"
            , row.names = FALSE
            , col.names=TRUE
            , file="./5_random_forest_2/RF_PA_infect.txt")
write.table(make_summary_RF(all_output = all_RF_predictBD, otu_type = "PA", bd_type = "PABD", prop=0.8, repl = 1)
            , quote = FALSE
            , sep = "\t"
            , row.names = FALSE
            , col.names=TRUE
            , file="./5_random_forest_2/RF_PA_PABD.txt")
write.table(make_summary_RF(all_output = all_RF_predictBD, otu_type = "count", bd_type = "infect", prop=0.8, repl = 1)
            , quote = FALSE
            , sep = "\t"
            , row.names = FALSE
            , col.names=TRUE
            , file="./5_random_forest_2/RF_count_infect.txt")
write.table(make_summary_RF(all_output = all_RF_predictBD, otu_type = "count", bd_type = "PABD", prop=0.8, repl = 1)
            , quote = FALSE
            , sep = "\t"
            , row.names = FALSE
            , col.names=TRUE
            , file="./5_random_forest_2/RF_count_PABD.txt")


#### Are the OTUs predicting infection just OTUs that correlate with that particular species?
all_RF_predictSpecies <- list()

RERUN_species=FALSE
if (RERUN_species) {
  total <- 2*2
  # create progress bar
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  i <- 1
  for ( otu_type in c("rank","PA","count") ){
    # for ( bd_type in c("PABD","infect")) {
    mf_temp <- get(paste0("mf_",otu_type)) %>% 
      separate(indivID, into=c("species","indiv")) %>%
      select(-c(SampleID, indiv, prepost, Bd_exposure)) %>%
      mutate(species=factor(species))
    
    Error_metric <- "MeanDecreaseAccuracy"
    
    # Presence absence
    set.seed(423)
    r <- sample(nrow(mf_temp), 0.7*nrow(mf_temp), replace=FALSE)
    trainSet <- mf_temp[r,]
    testSet <- mf_temp[-r,]
    
    RF <- randomForest(species ~ ., data=trainSet,importance=TRUE)
    testSet_predictions <- predict(RF, testSet, type="class")
    # Get error of predicted set
    # error <- ifelse(bd_type == "PABD", mean(testSet_predictions==testSet$response), mean((testSet_predictions-testSet$response)^2))
    error <- 1-mean(testSet_predictions==testSet$species)
    # importance of factors
    importance_RF <- data.frame(importance(RF)) %>%
      mutate(taxa=rownames(importance(RF))) %>%
      select(taxa, everything()) %>%
      arrange(-get(paste0(Error_metric)))
    
    # Add inhibitory to importance
    # Match inhibitory to not
    importance_RF[,c("taxonomy","inhibitory")] <- Taxa[match(importance_RF$taxa, Taxa$Sequence), c("Taxa","inhibitory")]
    importance_RF <- importance_RF %>%
      mutate(taxonomy=ifelse(is.na(taxonomy), taxa, taxonomy))
    
    # Save data
    all_RF_predictSpecies[[paste(otu_type)]][[paste0("mf_pred_",otu_type)]] <- mf_temp
    all_RF_predictSpecies[[paste(otu_type)]][[paste0("RF_",otu_type)]] <- RF
    all_RF_predictSpecies[[paste(otu_type)]][[paste0("importance_RF_",otu_type)]] <- importance_RF
    all_RF_predictSpecies[[paste(otu_type)]][[paste0("testSet_error",otu_type)]] <- error
    all_RF_predictSpecies[[paste(otu_type)]][[paste0("test_train_comparison")]] <- data.frame(Test_pred=testSet_predictions, Test_obs=testSet$species)
    
    # Insert into data summary plot
    #------- PROGRESS BAR ------- #
    Sys.sleep(0.1)
    i <- i+1
    # update progress bar
    setTxtProgressBar(pb, i)
    #-----------------------------#
    
  }
  
  save(all_RF_predictSpecies, file="./5_random_forest_2/all_RF_predictSpecies.RData")
}
load("./5_random_forest_2/all_RF_predictSpecies.RData")
all_RF_predictSpecies$PA$importance_RF_PA %>%
  select(c("taxonomy","inhibitory", "MeanDecreaseAccuracy")) %>%
  head(20)
ggsave("./5_random_forest_2/corr_spModel_bdModel_PA.pdf"
       ,all_RF_predictBD$PA_nop_PABD_0.7_1$importance %>%
         select(c("taxonomy", "MeanDecreaseAccuracy")) %>%
         rename(bdModel_MeanDecreaseAccuracy=MeanDecreaseAccuracy) %>%
         left_join(all_RF_predictSpecies$PA$importance_RF_PA) %>%
         rename(spModel_MeanDecreaseAccuracy=MeanDecreaseAccuracy) %>%
         select("taxonomy","bdModel_MeanDecreaseAccuracy","spModel_MeanDecreaseAccuracy") %>%
         ggplot() + geom_point(aes(x=spModel_MeanDecreaseAccuracy, y=bdModel_MeanDecreaseAccuracy))
)

all_RF_predictSpecies$rank$importance_RF_rank %>%
  select(c("taxonomy","inhibitory", "MeanDecreaseAccuracy")) %>%
  head(20)


all_RF_predictSpecies$count$importance_RF_count %>%
  select(c("taxonomy","inhibitory", "MeanDecreaseAccuracy")) %>%
  head(20)
ggsave("./5_random_forest_2/corr_spModel_bdModel_count.pdf"
       , all_RF_predictBD$count_nop_PABD_0.7_1$importance %>%
         select(c("taxonomy", "MeanDecreaseAccuracy")) %>%
         rename(bdModel_MeanDecreaseAccuracy=MeanDecreaseAccuracy) %>%
         left_join(all_RF_predictSpecies$count$importance_RF_count) %>%
         rename(spModel_MeanDecreaseAccuracy=MeanDecreaseAccuracy) %>%
         select("taxonomy","bdModel_MeanDecreaseAccuracy","spModel_MeanDecreaseAccuracy") %>%
         ggplot() + geom_point(aes(x=spModel_MeanDecreaseAccuracy, y=bdModel_MeanDecreaseAccuracy))
       
       )


### Let's inspect the OTUs that got pulled out

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


### Checking whether rank and PA correlate
all_RF_predictBD$rank_withp_infect_0.7_1$importance %>%
  mutate(PercIncrMSE_rank_infect = X.IncMSE) %>%
  select(taxonomy, PercIncrMSE_rank_infect) %>%
  full_join(all_RF_predictBD$PA_withp_infect_0.7_1$importance) %>%
  mutate(PercIncrMSE_PA_infect = X.IncMSE) %>%
  ggplot() +
  geom_point(aes(x=PercIncrMSE_PA_infect, y=PercIncrMSE_rank_infect))
all_RF_predictBD$rank_withp_PABD_0.7_1$importance %>%
  mutate(MeanDecrAccuracy_rank_PABD = MeanDecreaseAccuracy) %>%
  select(taxonomy, MeanDecrAccuracy_rank_PABD) %>%
  full_join(all_RF_predictBD$PA_withp_PABD_0.7_1$importance) %>%
  mutate(MeanDecrAccuracy_PA_PABD = MeanDecreaseAccuracy) %>%
  ggplot() +
  geom_point(aes(x=MeanDecrAccuracy_PA_PABD, y=MeanDecrAccuracy_rank_PABD))

all_RF_predictBD$count_withp_infect_0.7_1$importance %>%
  mutate(PercIncrMSE_count_infect = X.IncMSE) %>%
  select(taxonomy, PercIncrMSE_count_infect) %>%
  full_join(all_RF_predictBD$PA_withp_infect_0.7_1$importance) %>%
  mutate(PercIncrMSE_PA_infect = X.IncMSE) %>%
  ggplot() +
  geom_point(aes(x=PercIncrMSE_PA_infect, y=PercIncrMSE_count_infect))
all_RF_predictBD$count_withp_PABD_0.7_1$importance %>%
  mutate(MeanDecrAccuracy_count_PABD = MeanDecreaseAccuracy) %>%
  select(taxonomy, MeanDecrAccuracy_count_PABD) %>%
  full_join(all_RF_predictBD$PA_withp_PABD_0.7_1$importance) %>%
  mutate(MeanDecrAccuracy_PA_PABD = MeanDecreaseAccuracy) %>%
  ggplot() +
  geom_point(aes(x=MeanDecrAccuracy_PA_PABD, y=MeanDecrAccuracy_count_PABD))

all_RF_predictBD$count_withp_infect_0.7_1$importance %>%
  mutate(PercIncrMSE_count_infect = X.IncMSE) %>%
  select(taxonomy, PercIncrMSE_count_infect) %>%
  full_join(all_RF_predictBD$rank_withp_infect_0.7_1$importance) %>%
  mutate(PercIncrMSE_rank_infect = X.IncMSE) %>%
  ggplot() +
  geom_point(aes(x=PercIncrMSE_rank_infect, y=PercIncrMSE_count_infect))
all_RF_predictBD$count_withp_PABD_0.7_1$importance %>%
  mutate(MeanDecrAccuracy_count_PABD = MeanDecreaseAccuracy) %>%
  select(taxonomy, MeanDecrAccuracy_count_PABD) %>%
  full_join(all_RF_predictBD$rank_withp_PABD_0.7_1$importance) %>%
  mutate(MeanDecrAccuracy_rank_PABD = MeanDecreaseAccuracy) %>%
  ggplot() +
  geom_point(aes(x=MeanDecrAccuracy_rank_PABD, y=MeanDecrAccuracy_count_PABD))

## Everything

all_RF_predictBD$count_withp_infect_0.7_1$importance %>%
  filter(!is.na(inhibitory)) %>%
  ggplot() +
  geom_point(aes(x=X.IncMSE, y=inhibitory), position = position_jitter(width=0, height=0.1)) +
  geom_vline(data=data.frame(metric=comm_metrics, value=all_RF_predictBD$count_withp_infect_0.7_1$importance[match(comm_metrics,all_RF_predictBD$count_withp_infect_0.7_1$importance$taxa),"X.IncMSE" ] )
             , mapping=aes(xintercept=value, col=comm_metrics)
             , lwd=2, alpha=0.5) +
  scale_color_manual(values=comm_metrics_col)
# 
# all_RF_predictBD$rank_withp_infect$importance_RF_rank_withp_infect %>%
#   filter(!is.na(inhibitory)) %>%
#   ggplot() +
#   geom_point(aes(x=X.IncMSE, y=inhibitory), position = position_jitter(width=0, height=0.1)) +
#   geom_vline(data=data.frame(metric=comm_metrics, value=all_RF_predictBD$rank_withp_infect$importance_RF_rank_withp_infect[match(comm_metrics,all_RF_predictBD$rank_withp_infect$importance_RF_rank_withp_infect$taxa),"X.IncMSE" ] )
#             , mapping=aes(xintercept=value, col=comm_metrics)
#             , lwd=2, alpha=0.5) +
#   scale_color_manual(values=comm_metrics_col)

# all_RF_predictBD$PA_withp_infect$importance_RF_PA_withp_infect %>%
#   filter(!is.na(inhibitory)) %>%
#   ggplot() +
#   geom_point(aes(x=X.IncMSE, y=inhibitory), position = position_jitter(width=0, height=0.1)) +
#   geom_vline(data=data.frame(metric=comm_metrics, value=all_RF_predictBD$PA_withp_infect$importance_RF_PA_withp_infect[match(comm_metrics,all_RF_predictBD$PA_withp_infect$importance_RF_PA_withp_infect$taxa),"X.IncMSE" ] )
#              , mapping=aes(xintercept=value, col=comm_metrics)
#              , lwd=2, alpha=0.5) +
#   scale_color_manual(values=comm_metrics_col)

all_RF_predictBD$count_withp_PABD_0.7_1$importance %>%
  filter(!is.na(inhibitory)) %>%
  ggplot() +
  geom_point(aes(x=MeanDecreaseAccuracy, y=inhibitory), position = position_jitter(width=0, height=0.1)) +
  geom_vline(data=data.frame(metric=comm_metrics, value=all_RF_predictBD$count_withp_PABD_0.7_1$importance[match(comm_metrics,all_RF_predictBD$count_withp_PABD_0.7_1$importance$taxa),"MeanDecreaseAccuracy" ] )
             , mapping=aes(xintercept=value, col=comm_metrics)
             , lwd=2, alpha=0.5) +
  scale_color_manual(values=comm_metrics_col)
# 
# all_RF_predictBD$rank_withp_PABD$importance_RF_rank_withp_PABD %>%
#   filter(!is.na(inhibitory)) %>%
#   ggplot() +
#   geom_point(aes(x=MeanDecreaseAccuracy, y=inhibitory), position = position_jitter(width=0, height=0.1)) +
#   geom_vline(data=data.frame(metric=comm_metrics, value=all_RF_predictBD$rank_withp_PABD$importance_RF_rank_withp_PABD[match(comm_metrics,all_RF_predictBD$rank_withp_PABD$importance_RF_rank_withp_PABD$taxa),"MeanDecreaseAccuracy" ] )
#              , mapping=aes(xintercept=value, col=comm_metrics)
#              , lwd=2, alpha=0.5) +
#   scale_color_manual(values=comm_metrics_col)
# 
# 
# all_RF_predictBD$PA_withp_PABD$importance_RF_PA_withp_PABD %>%
#   filter(!is.na(inhibitory)) %>%
#   ggplot() +
#   geom_point(aes(x=MeanDecreaseAccuracy, y=inhibitory), position = position_jitter(width=0, height=0.1)) +
#   geom_vline(data=data.frame(metric=comm_metrics, value=all_RF_predictBD$PA_withp_PABD$importance_RF_PA_withp_PABD[match(comm_metrics,all_RF_predictBD$PA_withp_PABD$importance_RF_PA_withp_PABD$taxa),"MeanDecreaseAccuracy" ] )
#              , mapping=aes(xintercept=value, col=comm_metrics)
#              , lwd=2, alpha=0.5) +
#   scale_color_manual(values=comm_metrics_col)



#### Looking at importance #######


#### RANK ONLY #####

# Rank: PABD
testSet <- as.numeric(rownames(all_RF_predictBD$rank_onlyp_PABD_0.7_1$test_train_comparison))
species_order <- all_RF_predictBD$rank_onlyp_PABD_0.7_1$mf_used[testSet,"species"]
## Plot
ggsave("5_random_forest_2/predictions_PABD_rank.pdf",height=5, width=7,
data.frame(species=species_order
           ,Observed=all_RF_predictBD$rank_onlyp_PABD_0.7_1$test_train_comparison[,2]
           ,Community_metrics_only=all_RF_predictBD$rank_onlyp_PABD_0.7_1$test_train_comparison[,1]
           ,Community_and_ASV_metrics=all_RF_predictBD$rank_withp_PABD_0.7_1$test_train_comparison[,1]
           ,ASV_metrics_only=all_RF_predictBD$rank_nop_PABD_0.7_1$test_train_comparison[,1]
) %>%
  mutate(Individual=factor(paste0(seq(1,7), " (",species,")"))) %>%
  select(-species) %>%
  gather(-Individual, key=TrainingSet, value=Predicted_Bd_State) %>%
  mutate(TrainingSet=factor(TrainingSet, levels=c("Observed"
                                                  , "Community_metrics_only"
                                                  , "ASV_metrics_only"
                                                  , "Community_and_ASV_metrics"))) %>%
  ggplot() + geom_tile(aes(x=TrainingSet, y=Individual, fill=Predicted_Bd_State),col="black" )+
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust = 1)
        , rect = element_blank()) +
  geom_vline(aes(xintercept=1.5), lwd=4)
)

# Rank: infect

ggsave("5_random_forest_2/predictions_infect_rank.pdf",height=3, width=5,
       
data.frame(Observed=all_RF_predictBD$rank_onlyp_infect_0.7_1$test_train_comparison[,2]
           ,Community_metrics_only=all_RF_predictBD$rank_onlyp_infect_0.7_1$test_train_comparison[,1]
           ,Community_and_ASV_metrics=all_RF_predictBD$rank_withp_infect_0.7_1$test_train_comparison[,1]
           ,ASV_metrics_only=all_RF_predictBD$rank_nop_infect_0.7_1$test_train_comparison[,1]
) %>%
  gather(Community_metrics_only,Community_and_ASV_metrics,ASV_metrics_only, key=Model, value=Predicted) %>%
  ggplot() +
  geom_point(aes(x=Observed, y=Predicted, col=Model)) +
  geom_abline(aes(intercept=0, slope=1))+
  xlim(0,8)+ylim(0,8)+
  xlab("Observed Bd load") + ylab("Predicted Bd load")
)
data.frame(Observed=all_RF_predictBD$rank_onlyp_infect_0.7_1$test_train_comparison[,2]
           ,Community_level_only=all_RF_predictBD$rank_onlyp_infect_0.7_1$test_train_comparison[,1]
           ,Community_and_OTU_level=all_RF_predictBD$rank_withp_infect_0.7_1$test_train_comparison[,1]
           ,OTU_level_only=all_RF_predictBD$rank_nop_infect_0.7_1$test_train_comparison[,1]
) %>%
  gather(Community_level_only,Community_and_OTU_level,OTU_level_only, key=Model, value=Predicted) %>%
  mutate(Residuals=Predicted-Observed) %>%
  ggplot() +
  geom_point(aes(x=Model, y=Residuals)) 

#### PRESENCE ABSENCE ONLY #####
## PA: PABD
## Plot
ggsave("5_random_forest_2/predictions_PABD_PA.pdf",height=5, width=7,
       data.frame(species=species_order
                  ,Observed=all_RF_predictBD$PA_onlyp_PABD_0.7_1$test_train_comparison[,2]
                  ,Community_metrics_only=all_RF_predictBD$PA_onlyp_PABD_0.7_1$test_train_comparison[,1]
                  ,Community_and_ASV_metrics=all_RF_predictBD$PA_withp_PABD_0.7_1$test_train_comparison[,1]
                  ,ASV_metrics_only=all_RF_predictBD$PA_nop_PABD_0.7_1$test_train_comparison[,1]
       ) %>%
         mutate(Individual=factor(paste0(seq(1,7), " (",species,")"))) %>%
         select(-species) %>%
         gather(-Individual, key=TrainingSet, value=Predicted_Bd_State) %>%
         mutate(TrainingSet=factor(TrainingSet, levels=c("Observed"
                                                         , "Community_metrics_only"
                                                         , "ASV_metrics_only"
                                                         , "Community_and_ASV_metrics"))) %>%
         ggplot() + geom_tile(aes(x=TrainingSet, y=Individual, fill=Predicted_Bd_State),col="black" )+
         theme(axis.text.x = element_text(angle=45, vjust = 1, hjust = 1)
               , rect = element_blank()) +
         geom_vline(aes(xintercept=1.5), lwd=4)
       
       )


## PA: infect

ggsave("5_random_forest_2/predictions_infect_PA.pdf", height=3, width=5
       , data.frame(Observed=all_RF_predictBD$PA_onlyp_infect_0.7_1$test_train_comparison[,2]
                    ,Community_metrics_only=all_RF_predictBD$PA_onlyp_infect_0.7_1$test_train_comparison[,1]
                    ,Community_and_ASV_metrics=all_RF_predictBD$PA_withp_infect_0.7_1$test_train_comparison[,1]
                    ,ASV_metrics_only=all_RF_predictBD$PA_nop_infect_0.7_1$test_train_comparison[,1]
       ) %>%
         gather(Community_metrics_only,Community_and_ASV_metrics,ASV_metrics_only, key=Model, value=Predicted) %>%
         ggplot() +
         geom_point(aes(x=Observed, y=Predicted, col=Model)) +
         geom_abline(aes(intercept=0, slope=1))+
         xlim(0,8)+ylim(0,8)+
         xlab("Observed Bd load") +ylab("Predicted Bd load")
       )



#### COUNT ONLY #####

# count: PABD
## Plot
ggsave("5_random_forest_2/preditions_PABD_count.pdf", height=5, width=7,
       data.frame(species=species_order 
                  ,Observed=all_RF_predictBD$count_onlyp_PABD_0.7_1$test_train_comparison[,2]
                  ,Community_metrics_only=all_RF_predictBD$count_onlyp_PABD_0.7_1$test_train_comparison[,1]
                  ,Community_and_ASV_metrics=all_RF_predictBD$count_withp_PABD_0.7_1$test_train_comparison[,1]
                  ,ASV_metrics_only=all_RF_predictBD$count_nop_PABD_0.7_1$test_train_comparison[,1]) %>%
         mutate(Individual=factor(paste0(seq(1,7), " (",species,")"))) %>%
         select(-species) %>%
         gather(-Individual, key=TrainingSet, value=Predicted_Bd_State) %>%
         mutate(TrainingSet=factor(TrainingSet, levels=c("Observed"
                                                         , "Community_metrics_only"
                                                         , "ASV_metrics_only"
                                                         , "Community_and_ASV_metrics"))) %>%
         ggplot() +geom_tile(aes(x=TrainingSet, y=Individual, fill=Predicted_Bd_State),col="black" ) +
         geom_vline(aes(xintercept=1.5), lwd=4) +
         theme(axis.text.x = element_text(angle=45, vjust = 1, hjust = 1), rect = element_blank())
       )


  
  
# count: infect
ggsave("5_random_forest_2/predictions_infect_count.pdf", height=3, width=5,
       data.frame(Observed=all_RF_predictBD$count_onlyp_infect_0.7_1$test_train_comparison[,2]
                  ,Community_metrics_only=all_RF_predictBD$count_onlyp_infect_0.7_1$test_train_comparison[,1]
                  ,Community_and_ASV_metrics=all_RF_predictBD$count_withp_infect_0.7_1$test_train_comparison[,1]
                  ,ASV_metrics_only=all_RF_predictBD$count_nop_infect_0.7_1$test_train_comparison[,1]
       ) %>%
         gather(Community_metrics_only,Community_and_ASV_metrics,ASV_metrics_only, key=Model, value=Predicted) %>%
         ggplot() +
         geom_point(aes(x=Observed, y=Predicted, col=Model)) +
         geom_abline(aes(intercept=0, slope=1))+
         xlim(0,8)+ylim(0,8)+
         xlab("Observed Bd load") + ylab("Predicted Bd load")
       )


data.frame(Observed=all_RF_predictBD$count_onlyp_infect_0.7_1$test_train_comparison[,2]
           ,Community_level_only=all_RF_predictBD$count_onlyp_infect_0.7_1$test_train_comparison[,1]
           ,Community_and_OTU_level=all_RF_predictBD$count_withp_infect_0.7_1$test_train_comparison[,1]
           ,OTU_level_only=all_RF_predictBD$count_nop_infect_0.7_1$test_train_comparison[,1]
) %>%
  gather(Community_level_only,Community_and_OTU_level,OTU_level_only, key=Model, value=Predicted) %>%
  mutate(Residuals=Predicted-Observed) %>%
  ggplot() +
  geom_point(aes(x=Model, y=Residuals)) 

# ### Are inhibitory bacteria more likely to be predictive?
# rank_infect_filtered <- all_RF_predictBD$rank_nop_infect$importance_RF_rank_nop_infect %>%
#   filter(X.IncMSE>0, !is.na(inhibitory)) %>%
#   mutate(inhibitory=factor(inhibitory))
# anova(lm(X.IncMSE ~ inhibitory, data=rank_infect_filtered))
# summary(glm(inhibitory ~ X.IncMSE, data=rank_infect_filtered, family = binomial))
# rank_infect_filtered %>%
#   ggplot() +geom_point(aes(x=inhibitory, y=X.IncMSE), position = position_jitter(height=0, width=0.1))
# 
# rank_PABD_filtered <- all_RF_predictBD$rank_nop_PABD$importance_RF_rank_nop_PABD %>%
#   filter(MeanDecreaseAccuracy>0, !is.na(inhibitory)) %>%
#   mutate(inhibitory=factor(inhibitory))
# anova(lm(MeanDecreaseAccuracy ~ inhibitory, data=rank_PABD_filtered))
# summary(glm(inhibitory ~ MeanDecreaseAccuracy, data=rank_PABD_filtered, family = binomial))
# rank_PABD_filtered %>%
#   ggplot() +geom_point(aes(x=inhibitory, y=MeanDecreaseAccuracy), position=position_jitter(height=0, width=0.1))
# 
# 
# 
# 


####### Save PA table for otu ##########

mf_tosave <- mf_PA %>%
  filter(prepost=="Pre", Bd_exposure=="Bd-exposed") %>%
  select(-c("SampleID","prepost","Bd_exposure")) %>%
  gather(-c(indivID),key=OTU, value=PA) %>%
  mutate(PA=as.numeric(PA)) %>%
  # filter(indivID==c("Anbo_4","Anbo_2"), OTU=="TACGTAGGGTACGAGCGTTGTCCGGAATTATTGGGCGTAAAGAGCTCGTAGGTGGTTGGTCACGTCTGCTGTGGAAACGCAACGCTTAAC") %>%
  group_by(indivID, OTU) %>%
  summarize(avePA=mean(PA)) %>%
  ungroup() %>%
  separate(indivID, into=c("species","indiv"), remove=FALSE) %>%
  select(-indiv) %>%
  # filter(!is.na(avePA)) %>%
  spread(key=OTU, value=avePA)  %>% 
  select(-species) %>%
  as.data.frame()

rownames(mf_tosave) <- mf_tosave$indivID
mf_tosave <- mf_tosave %>%
  select(-indivID) %>%
  t() %>%as.data.frame()

otu_PA_tosave <- cbind("#OTU ID"=rownames(mf_tosave),mf_tosave) %>%
  select("#OTU ID",everything())

# Save rank and prevelance 
write.table(otu_PA_tosave, file="./5_random_forest_2/OTUTable_prevalence.txt", quote=FALSE, sep="\t"
            , row.names = FALSE, col.names = TRUE)



######## LOAD QIIME2 DATA #######
for ( type in c("rarefied_mean","rarefied_median","rarefied_sum","stratified-rare")) {
  # type="rarefied_mean"
  all_props <- data.frame()
  for ( prop in c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)) {
    if (file.exists(paste0("./5a_QIIME2_RF/ALL_",type,"/predictions_",prop,"/predictions.tsv"))) {
      temp <- read.delim(paste0("./5a_QIIME2_RF/ALL_",type,"/predictions_",prop,"/predictions.tsv"))
      temp$prop <- prop
      all_props <- rbind(all_props, temp)
    }
  }
  if ( type == "stratified-rare") {
   indivList <- pull(mf_alt_filt_final[match(all_props$SampleID,mf_alt_filt_final$SampleID),"indivID"])
   all_props <- all_props %>%
     mutate(SampleID = indivList)
  }
  combined <- p %>%
    select(indivID, infect) %>%
    rename(SampleID=indivID, Observed=infect) %>%
    right_join(all_props) %>%
    arrange(SampleID, prop) %>%
    mutate(prop_test=1-prop) %>%
    rename(prop_train=prop)
  
  name <- gsub("-","_",type)
  assign(paste0("predobs_",name), combined)
  
  R2 <- combined %>%
    group_by(prop_train) %>%
    mutate(meanObs=mean(Observed)) %>%
    ungroup() %>%
    mutate(res_sq = (Observed-prediction)^2, tot_sq=(Observed-meanObs)^2) %>%
    group_by(prop_train) %>%
    summarize(var_res=sum(res_sq)/sum(tot_sq), mse = mean(res_sq)) %>%
    rename(var_unexplained = var_res) %>%
    mutate(prop_test=1-prop_train)%>%
    mutate(model=name)
  assign(paste0("R2_",name), R2)
  
}

rbind(R2_rarefied_mean, R2_rarefied_median, R2_rarefied_sum, R2_stratified_rare) %>%
  rename(Aggregate_method=model) %>%
  ggplot() + geom_line(aes(x=prop_train, y=var_unexplained, col=Aggregate_method)) +
  ylab("% Unexplained variation (load)") + xlab("Proportion used for training set")

rbind(R2_rarefied_mean, R2_rarefied_median, R2_rarefied_sum, R2_stratified_rare) %>%
  rename(Aggregate_method=model) %>%
  ggplot() + geom_line(aes(x=prop_train, y=mse, col=Aggregate_method)) +
  ylab("Mean Squared Error (load)") + xlab("Proportion used for training set")


######### LOAD QIIME2 data binary ########
for ( type in c("rarefied_mean_bin","rarefied_median_bin","rarefied_sum_bin","stratified-rare_bin")) {
  all_props <- data.frame()
  for ( prop in c(0.1,0.2,0.3,0.4,0.5,0.6)) {
    if ( file.exists(paste0("./5a_QIIME2_RF/ALL_",type,"/predictions_",prop,"/predictions.tsv"))) {
      temp <- read.delim(paste0("./5a_QIIME2_RF/ALL_",type,"/predictions_",prop,"/predictions.tsv"))
      temp$prop <- prop
      if ( type == "stratified-rare_bin") {
        indivList <- pull(mf_alt_filt_final[match(temp$SampleID,mf_alt_filt_final$SampleID),"indivID"])
        temp <- temp %>%
          mutate(SampleID = indivList)
      }
      all_props <- rbind(all_props, temp)
      
    } else {
      next
    }
  }
  combined <- p %>%
    select(indivID, infect) %>%
    rename(SampleID=indivID, Observed=infect) %>%
    mutate(Observed=ifelse(Observed>0, TRUE, FALSE)) %>%
    right_join(all_props) %>%
    arrange(SampleID, prop) %>%
    mutate(prop_test=1-prop) %>%
    rename(prop_train=prop)
  
  name <- gsub("-","_",type)
  combined$type <- name
  
  assign(paste0("predobs_bin_",name), combined)
  
  errorRate <- combined %>%
    group_by(prop_train, type) %>%
    summarise(correct=mean(Observed==prediction))%>%
    ungroup()
  assign(paste0("errRate_",name), errorRate)
  
  
}

rbind(predobs_bin_rarefied_mean_bin, predobs_bin_rarefied_median_bin, predobs_bin_rarefied_sum_bin, predobs_bin_stratified_rare_bin) %>%
  group_by(prop, type) %>%
  summarise(correct=mean(Observed==prediction)) %>%
  ggplot() + geom_line(aes(x=prop, y=correct, col=type)) +
  ylab("Proportion Correct (Pres/Abs)") + xlab("Proportion used for training set")


#### Comparison between R and QIIME2

## Combine R2 for all QIIME
all_R2 <- rbind(R2_rarefied_sum, R2_stratified_rare, R2_rarefied_median, R2_rarefied_mean) %>%
  mutate(R2=1-var_unexplained) %>%
  select(prop_train, R2, model)

# Plot R2 against prop_train for all models
ggsave("5_random_forest_2/qiimevsR_R2.pdf", height=3, width=5
       , extractError(all_RF_predictBD, bd_type = "infect", otu_type = "count", error_type="R2") %>%
         filter(type=="Validation") %>%
         select(prop_train, nop)%>%
         rename(R2=nop) %>%
         mutate(model="R (ASV count)") %>%
         rbind(all_R2) %>%
         mutate(model=ifelse(model=="rarefied_mean","QIIME2 (rarefied_mean)"
                             , ifelse(model=="rarefied_median", "QIIME2 (rarefied_median)",
                                      ifelse(model=="rarefied_sum","QIIME2 (rarefied_sum)"
                                             , ifelse(model=="stratified_rare","QIIME2 (rarefied_stratified)", "R (ASV count summed)"))))) %>%
         ggplot() +
         # + geom_bar(aes(x=prop_train, y=R2, fill=model)
         # , stat = "identity", position = "dodge")
         geom_line(aes(x=prop_train, y=R2, col=model))+
         coord_cartesian(ylim=c(0, 1))+
         ylab("R^2 of model")+xlab("Proportion of data used for training")
)


predobs_stratified_rare %>%
  filter(prop_train==0.4) %>%
  select(Observed, prediction) %>%
  rename(Test_pred=prediction, Test_obs=Observed) %>%
  mutate(Type="QIIME (rarefied stratified)") %>%
  rbind(data.frame(all_RF_predictBD$count_nop_infect_0.5_1$test_train_comparison, Type="R")) %>%
  ggplot() +geom_point(aes(x=Test_obs, y=Test_pred, col=Type)) + geom_abline(aes(slope=1, intercept=0))

ggsave("5_random_forest_2/qiimevsR_proptrain0.4.pdf"
       ,predobs_rarefied_median %>%
         filter(prop_train==0.4) %>%
         select(Observed, prediction) %>%
         rename(Test_pred=prediction, Test_obs=Observed) %>%
         mutate(Type="QIIME (rarefied median)") %>%
         rbind(data.frame(all_RF_predictBD$count_nop_infect_0.4_1$test_train_comparison, Type="R")) %>%
         ggplot() +geom_point(aes(x=Test_obs, y=Test_pred, col=Type)) + geom_abline(aes(slope=1, intercept=0))+
         ylim(0,8)+xlim(0,8)
)
ggsave("5_random_forest_2/qiimevsR_proptrain0.7.pdf"
       ,predobs_rarefied_median %>%
         filter(prop_train==0.7) %>%
         select(Observed, prediction) %>%
         rename(Test_pred=prediction, Test_obs=Observed) %>%
         mutate(Type="QIIME (rarefied median)") %>%
         rbind(data.frame(all_RF_predictBD$count_nop_infect_0.7_1$test_train_comparison, Type="R")) %>%
         ggplot() +geom_point(aes(x=Test_obs, y=Test_pred, col=Type)) + geom_abline(aes(slope=1, intercept=0))+
         ylim(0,8)+xlim(0,8)
)


## Combine R2 for all QIIME
all_errRate <- rbind(errRate_stratified_rare_bin, errRate_rarefied_mean_bin, errRate_rarefied_median_bin, errRate_rarefied_sum_bin) %>%
  rename(model=type) %>%
  select(prop_train, correct, model) 

extractError(all_RF_predictBD, bd_type = "PABD", otu_type = "count", error_type="errorRate") %>%
  filter(type=="Validation") %>%
  mutate(correct=1-nop) %>%
  select(prop_train, correct)%>%
  mutate(model="R (ASV count)") %>%
  rbind(all_errRate) %>%
  mutate(model=ifelse(model=="rarefied_mean_bin","QIIME2 (rarefied_mean)"
                      , ifelse(model=="rarefied_median_bin", "QIIME2 (rarefied_median)",
                               ifelse(model=="rarefied_sum_bin","QIIME2 (rarefied_sum)"
                                      , ifelse(model=="stratified_rare_bin","QIIME2 (rarefied_stratified)", "R (ASV count summed)"))))) %>%
  ggplot() +
  # + geom_bar(aes(x=prop_train, y=R2, fill=model)
  # , stat = "identity", position = "dodge")
  geom_line(aes(x=prop_train, y=correct, col=model))+
  coord_cartesian(ylim=c(0, 1))

ggsave("5_random_forest_2/qiimevsR_binary.pdf", height=3, width=5
       , extractError(all_RF_predictBD, bd_type = "PABD", otu_type = "count", error_type="errorRate") %>%
         filter(type=="Validation") %>%
         mutate(correct=1-nop) %>%
         select(prop_train, correct)%>%
         mutate(model="R (ASV count)") %>%
         rbind(all_errRate) %>%
         mutate(model=ifelse(model=="rarefied_mean_bin","QIIME2 (rarefied_mean)"
                             , ifelse(model=="rarefied_median_bin", "QIIME2 (rarefied_median)",
                                      ifelse(model=="rarefied_sum_bin","QIIME2 (rarefied_sum)"
                                             , ifelse(model=="stratified_rare_bin","QIIME2 (rarefied_stratified)", "R (ASV count summed)"))))) %>%
         ggplot() +
         # + geom_bar(aes(x=prop_train, y=R2, fill=model)
         # , stat = "identity", position = "dodge")
         geom_line(aes(x=prop_train, y=correct, col=model))+
         coord_cartesian(ylim=c(0, 1)) +
         ylab("Proportion correct") + xlab("Proportion used for training set")
)


