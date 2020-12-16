#!/bin/bash Rscript

#### Loading data #####
library(randomForest)
library(tidyverse)
load("3_5sp_mapping_otu_downstream/mf_alt_filt_final.RData")
load("3_5sp_mapping_otu_downstream/otu_filt.RData")
load("3_5sp_mapping_otu_downstream/taxonomy.RData")
load("4_Bayesian_models/all_p.RData")
source("msc_randomforest_codebits.R")
dir.create("5_random_forest_diagnostics")

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

#### Run PA, all p, infect ####
RERUN_RF <- FALSE
if ( RERUN_RF ) {
  total <- 1*2*3*1*8*3
  
  # create progress bar
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  i <- 1
  for ( v in c(1)) {
    # v=7
    all_RF_predictBD <- list()
    for ( otu_type in c("count","PA") ) {
      for ( p_type in c("onlyp","withp","nop") ) {
        # for ( p_type in c("onlyp") ) {
        for ( bd_type in c("infect")) {
          for ( prop in c(0.2,0.3,0.4,0.5,0.6,0.7,0.8)) {
            for ( repl in c(1:3) ) {
              
              all_RF_predictBD[[paste(otu_type, p_type, bd_type, prop, repl, sep="_")]] <- list()
              # Get a formatted mapping file
              get_mf_pred_temp_bd <- get_mf_pred_temp(p_type=p_type, otu_type = otu_type, bd_type=bd_type, p=p)
              mf_pred_temp_bd <- get_mf_pred_temp_bd[[1]]
              indiv_list <- get_mf_pred_temp_bd[[2]]
              # Mutate infection depending on question
              Error_metric <- ifelse(bd_type=="PABD","MeanDecreaseAccuracy","X.IncMSE" )
              
              # Pick training and validation set
              # set.seed(423*repl)
              r <- sample(nrow(mf_pred_temp_bd), prop*nrow(mf_pred_temp_bd), replace=FALSE)
              trainSet <- mf_pred_temp_bd[r,]
              testSet <- mf_pred_temp_bd[-r,]
              testSet_indiv <- as.character(unlist(indiv_list[-r,]))
              
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
              # all_RF_predictBD[[paste(otu_type, p_type, bd_type, prop,repl, sep="_")]][["mf_used"]] <- mf_pred_temp_bd
              all_RF_predictBD[[paste(otu_type, p_type, bd_type, prop,repl, sep="_")]][["RF"]] <- RF
              all_RF_predictBD[[paste(otu_type, p_type, bd_type, prop,repl, sep="_")]][["importance"]] <- importance_RF
              all_RF_predictBD[[paste(otu_type, p_type, bd_type, prop,repl, sep="_")]][["testSetError"]] <- error
              all_RF_predictBD[[paste(otu_type, p_type, bd_type, prop,repl, sep="_")]][["test_train_comparison"]] <- data.frame(indivID=testSet_indiv, Test_pred=testSet_predictions, Test_obs=testSet$response)
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
    assign(paste0("RF_infect_paramsweep_",v), all_RF_predictBD)
    # get(paste0("RF_infect_10_",v)) %>%
    # save(file=paste0("./5_random_forest_prevalence/","RF_infect_10_",v,".RData"))
  }
  close(pb)
  # 
  # load("./5_random_forest_diagnostics/RF_infect_paramsweep_1.RData" )
  
  save(RF_infect_paramsweep_1,file="./5_random_forest_diagnostics/RF_infect_paramsweep_1.RData" )
  infect_test_train_comparison_diag <- extractComparison(allOutput = RF_infect_paramsweep_1, replicates=c(1:3), bd_type = "infect")
  infect_train_MSE_diag <- extractTrainError(allOutput = RF_infect_paramsweep_1, replicates=c(1:3), bd_type = "infect", error_type = "MSE")
  infect_train_R2_diag <- extractTrainError(allOutput = RF_infect_paramsweep_1, replicates=c(1:3), bd_type = "infect", error_type = "R2")
  save(infect_test_train_comparison_diag,file="./5_random_forest_diagnostics/infect_test_train_comparison_diag.RData" )
  save(infect_train_MSE_diag,file="./5_random_forest_diagnostics/infect_train_MSE_diag.RData" )
  save(infect_train_R2_diag,file="./5_random_forest_diagnostics/infect_train_R2_diag.RData" )
  
}

remove(RF_infect_paramsweep_1)
#### Run PA, all p, PABD. ####
# RERUN_RF <- TRUE
if ( RERUN_RF ) {
  total <- 1*2*3*1*8*3
  # create progress bar
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  i <- 1
  for ( v in c(1)) {
    # v=1
    all_RF_predictBD <- list()
    
    for ( otu_type in c("count","PA") ) {
      for ( p_type in c("onlyp","withp","nop") ) {
        # for ( p_type in c("onlyp","nop") ) {
        for ( bd_type in c("PABD")) {
          for ( prop in c(0.2,0.3,0.4,0.5,0.6,0.7,0.8)) {
            for ( repl in c(1:3) ) {
              
              all_RF_predictBD[[paste(otu_type, p_type, bd_type, prop, repl, sep="_")]] <- list()
              # Get a formatted mapping file
              get_mf_pred_temp_bd <- get_mf_pred_temp(p_type=p_type, otu_type = otu_type, bd_type=bd_type, p=p)
              mf_pred_temp_bd <- get_mf_pred_temp_bd[[1]]
              indiv_list <- get_mf_pred_temp_bd[[2]]
              # Mutate infection depending on question
              Error_metric <- ifelse(bd_type=="PABD","MeanDecreaseAccuracy","X.IncMSE" )
              
              # Pick training and validation set
              # set.seed(423*repl)
              r <- sample(nrow(mf_pred_temp_bd), prop*nrow(mf_pred_temp_bd), replace=FALSE)
              trainSet <- mf_pred_temp_bd[r,]
              # Need to make sure there are at least 2 classes
              while (length(unique(trainSet$response))<2) {
                r <- sample(nrow(mf_pred_temp_bd), prop*nrow(mf_pred_temp_bd), replace=FALSE)
                trainSet <- mf_pred_temp_bd[r,]
              }
              testSet <- mf_pred_temp_bd[-r,]
              testSet_indiv <- as.character(unlist(indiv_list[-r,]))
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
              # all_RF_predictBD[[paste(otu_type, p_type, bd_type, prop,repl, sep="_")]][["mf_used"]] <- mf_pred_temp_bd
              all_RF_predictBD[[paste(otu_type, p_type, bd_type, prop,repl, sep="_")]][["RF"]] <- RF
              all_RF_predictBD[[paste(otu_type, p_type, bd_type, prop,repl, sep="_")]][["importance"]] <- importance_RF
              all_RF_predictBD[[paste(otu_type, p_type, bd_type, prop,repl, sep="_")]][["testSetError"]] <- error
              all_RF_predictBD[[paste(otu_type, p_type, bd_type, prop,repl, sep="_")]][["test_train_comparison"]] <- data.frame(indivID = testSet_indiv, Test_pred=testSet_predictions, Test_obs=testSet$response)
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
    # assign(paste0("RF_PABD_10_",v), all_RF_predictBD)
    assign(paste0("RF_PABD_parasweep_",v), all_RF_predictBD)
    
    # get(paste0("RF_PABD_10_",v)) %>%
    # save(file=paste0("./5_random_forest_prevalence/","RF_PABD_10_",v,".RData"))
  }
  close(pb)

  save(RF_PABD_parasweep_1,file="./5_random_forest_diagnostics/RF_PABD_parasweep_1.RData" )
  PABD_test_train_comparison_diag <- extractComparison(allOutput = RF_PABD_parasweep_1, replicates=c(1:3), bd_type = "PABD")
  PABD_train_diag <- extractTrainError(allOutput = RF_PABD_parasweep_1, replicates=c(1:3), bd_type="PABD", error_type = "errorRate")
  save(PABD_test_train_comparison_diag,file="./5_random_forest_diagnostics/PABD_test_train_comparison_diag.RData" )
  save(PABD_train_diag,file="./5_random_forest_diagnostics/PABD_train_diag.RData" )
  
}

# load("./5_random_forest_diagnostics/RF_PABD_parasweep_1.RData")
# load("./5_random_forest_diagnostics/RF_infect_paramsweep_1.RData")

load("./5_random_forest_diagnostics/infect_test_train_comparison_diag.RData")
load("./5_random_forest_diagnostics/PABD_test_train_comparison_diag.RData")
load("./5_random_forest_diagnostics/PABD_train_diag.RData")
# PABD_train_diag <- PABD_train_diag %>%
#   as_tibble() %>%
#   mutate(error=as.numeric(error)/100, prop=as.numeric(prop)/10, repl=as.numeric(repl))
load("./5_random_forest_diagnostics/infect_train_MSE_diag.RData")
# infect_train_MSE_diag <- infect_train_MSE_diag %>%
#   as.data.frame() %>%
#   mutate(error=as.numeric(error)/100, prop=as.numeric(prop)/10, repl=as.numeric(repl))
load("./5_random_forest_diagnostics/infect_train_R2_diag.RData")
# infect_train_R2_diag <- infect_train_R2_diag %>%
#   as.data.frame() %>%
#   mutate(error=as.numeric(error)/100, prop=as.numeric(prop)/10, repl=as.numeric(repl))

#### Error calculations ####
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
  select(-error_type) %>%
  left_join(infect_train_R2_diag) %>%
  mutate(R2=error, set="Train") %>%
  select(p_type, prop, otu_type, repl, MSE, R2, set)

PABD_train_diag <- PABD_train_diag %>%
  mutate(nCorrect=1-error) %>%
  select(p_type, prop, otu_type, repl, nCorrect) %>%
  mutate(set = "Train")




#### Plotting ####

infect_test_train_comparison_diag %>%
  ggplot(aes(x=Test_obs, y=Test_pred)) + geom_point(aes(col=factor(prop))) +
  geom_abline(aes(intercept=0, slope=1)) +
  geom_smooth(aes(col=factor(prop)), method="lm") +
  ylim(0,9)+xlim(0,9) +
  facet_grid(otu_type~p_type)

rbind(summary_error_infect_diag, infect_train_MSER2_diag) %>%
  ggplot(aes(x=prop, y=MSE)) + geom_point(aes(col=p_type)) +
  geom_smooth(aes(col=p_type, lty=set), se=FALSE) +
  facet_grid(.~otu_type)

rbind(summary_error_infect_diag, infect_train_MSER2_diag) %>%
  mutate(R2 = ifelse(R2<0, 0,R2)) %>%
  ggplot(aes(x=prop, y=R2)) + geom_jitter(aes(col=p_type), height=0, width=0.01) +
  geom_smooth(aes(col=p_type, lty=set), se=FALSE) +
  facet_grid(.~otu_type) +
  ylim(0,1)

rbind(summary_error_PABD_diag, PABD_train_diag) %>%
  ggplot(aes(x=prop, y=nCorrect)) + geom_jitter(aes(col=p_type), height=0, width=0.025) +
  geom_smooth(aes(col=p_type, lty=set), se=FALSE) +
  facet_grid(.~otu_type)
  

