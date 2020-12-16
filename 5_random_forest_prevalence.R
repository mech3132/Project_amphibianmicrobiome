#!/bin/bash Rscript

#### Loading data #####
library(randomForest)
library(tidyverse)
load("3_5sp_mapping_otu_downstream/mf_alt_filt_final.RData")
load("3_5sp_mapping_otu_downstream/otu_filt.RData")
load("3_5sp_mapping_otu_downstream/taxonomy.RData")
load("4_Bayesian_models/all_p.RData")
source("msc_randomforest_codebits.R")
dir.create("5_random_forest_prevalence")

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
  total <- 1*3*1*1*60
  
  # create progress bar
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  i <- 1
  for ( v in c(1:3)) {
    # v=7
    all_RF_predictBD <- list()
    
    for ( otu_type in c("PA") ) {
      for ( p_type in c("onlyp","withp","nop") ) {
      # for ( p_type in c("onlyp") ) {
        
        for ( bd_type in c("infect")) {
          for ( prop in c(0.8)) {
            for ( repl in c(1:10) ) {
              
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
    assign(paste0("RF_infect_10_",v), all_RF_predictBD)
    # get(paste0("RF_infect_10_",v)) %>%
    # save(file=paste0("./5_random_forest_prevalence/","RF_infect_10_",v,".RData"))
  }
  close(pb)
  # 
  save(RF_infect_10_1,file="./5_random_forest_prevalence/RF_infect_10_1.RData" )
  save(RF_infect_10_2,file="./5_random_forest_prevalence/RF_infect_10_2.RData" )
  save(RF_infect_10_3,file="./5_random_forest_prevalence/RF_infect_10_3.RData" )

  # create progress bar
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  i <- 1
  for ( v in c(4:6)) {
    # v=7
    all_RF_predictBD <- list()
    
    for ( otu_type in c("PA") ) {
      for ( p_type in c("onlyp","withp","nop") ) {
        # for ( p_type in c("onlyp") ) {
        
        for ( bd_type in c("infect")) {
          for ( prop in c(0.8)) {
            for ( repl in c(1:10) ) {
              
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
    assign(paste0("RF_infect_10_",v), all_RF_predictBD)
    # get(paste0("RF_infect_10_",v)) %>%
    # save(file=paste0("./5_random_forest_prevalence/","RF_infect_10_",v,".RData"))
  }
  close(pb)
  save(RF_infect_10_4,file="./5_random_forest_prevalence/RF_infect_10_4.RData" )
  save(RF_infect_10_5,file="./5_random_forest_prevalence/RF_infect_10_5.RData" )
  save(RF_infect_10_6,file="./5_random_forest_prevalence/RF_infect_10_6.RData" )

  # RF_infect_PA <- all_RF_predictBD
  # save(RF_infect_PA, file="./5_random_forest_prevalence/RF_infect_PA.RData")
  
}


#### Run PA, all p, PABD. ####
# RERUN_RF <- TRUE
if ( RERUN_RF ) {
  total <- 1*2*1*1*60
  # create progress bar
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  i <- 1
  for ( v in c(1:3)) {
    # v=1
    all_RF_predictBD <- list()
    
    for ( otu_type in c("PA") ) {
      for ( p_type in c("onlyp","withp","nop") ) {
      # for ( p_type in c("onlyp","nop") ) {
        for ( bd_type in c("PABD")) {
          for ( prop in c(0.8)) {
            for ( repl in c(1:10) ) {
                
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
    assign(paste0("RF_PABD_10_",v), all_RF_predictBD)
    
    # get(paste0("RF_PABD_10_",v)) %>%
    # save(file=paste0("./5_random_forest_prevalence/","RF_PABD_10_",v,".RData"))
  }
  close(pb)
  
  save(RF_PABD_10_1,file="./5_random_forest_prevalence/RF_PABD_10_1.RData" )
  save(RF_PABD_10_2,file="./5_random_forest_prevalence/RF_PABD_10_2.RData" )
  save(RF_PABD_10_3,file="./5_random_forest_prevalence/RF_PABD_10_3.RData" )

  total <- 1*2*1*1*60
  # create progress bar
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  i <- 1
  for ( v in c(4:6)) {
    # v=1
    all_RF_predictBD <- list()
    
    for ( otu_type in c("PA") ) {
      for ( p_type in c("onlyp","withp","nop") ) {
        # for ( p_type in c("onlyp","nop") ) {
        for ( bd_type in c("PABD")) {
          for ( prop in c(0.8)) {
            for ( repl in c(1:10) ) {
              
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
    assign(paste0("RF_PABD_10_",v), all_RF_predictBD)
    
    # get(paste0("RF_PABD_10_",v)) %>%
    # save(file=paste0("./5_random_forest_prevalence/","RF_PABD_10_",v,".RData"))
  }
  close(pb)
  save(RF_PABD_10_4,file="./5_random_forest_prevalence/RF_PABD_10_4.RData" )
  save(RF_PABD_10_5,file="./5_random_forest_prevalence/RF_PABD_10_5.RData" )
  save(RF_PABD_10_6,file="./5_random_forest_prevalence/RF_PABD_10_6.RData" )
  
  # RF_PABD_PA <- all_RF_predictBD
  # save(RF_PABD_PA, file="./5_random_forest_prevalence/RF_PABD_PA.RData")
}

RELOAD_RAW <- FALSE
if ( RELOAD_RAW ) {
  load("./5_random_forest_prevalence/RF_infect_10_1.RData")
  load("./5_random_forest_prevalence/RF_infect_10_2.RData")
  load("./5_random_forest_prevalence/RF_infect_10_3.RData")
  load("./5_random_forest_prevalence/RF_infect_10_4.RData")
  load("./5_random_forest_prevalence/RF_infect_10_5.RData")
  load("./5_random_forest_prevalence/RF_infect_10_6.RData")
  
  load("./5_random_forest_prevalence/RF_PABD_10_1.RData")
  load("./5_random_forest_prevalence/RF_PABD_10_2.RData")
  load("./5_random_forest_prevalence/RF_PABD_10_3.RData")
  load("./5_random_forest_prevalence/RF_PABD_10_4.RData")
  load("./5_random_forest_prevalence/RF_PABD_10_5.RData")
  load("./5_random_forest_prevalence/RF_PABD_10_6.RData")
  #### Infect extract #####
  
  infect_test_train_comparison_ponly <- extractItem(allOutput = list(RF_infect_10_1, RF_infect_10_2, RF_infect_10_3), bd_type="infect", otu_type = "PA", p_types="onlyp", reps=c(1:10)
              , extract="test_train_comparison")
  infect_test_train_comparison_ponly2 <- extractItem(allOutput = list(RF_infect_10_4, RF_infect_10_5, RF_infect_10_6), bd_type="infect", otu_type = "PA", p_types="onlyp", reps=c(1:10)
                                                    , extract="test_train_comparison", repstart=31)
  save(infect_test_train_comparison_ponly, file="./5_random_forest_prevalence/infect_test_train_comparison_ponly.RData")
  save(infect_test_train_comparison_ponly2, file="./5_random_forest_prevalence/infect_test_train_comparison_ponly2.RData")
  
  infect_test_train_comparison_nop <- extractItem(allOutput = list(RF_infect_10_1, RF_infect_10_2, RF_infect_10_3), bd_type="infect", otu_type = "PA", p_types="nop", reps=c(1:10)
                                                  , extract="test_train_comparison")
  infect_test_train_comparison_nop2 <- extractItem(allOutput = list(RF_infect_10_4, RF_infect_10_5, RF_infect_10_6), bd_type="infect", otu_type = "PA", p_types="nop", reps=c(1:10)
                                                  , extract="test_train_comparison", repstart=31)
  save(infect_test_train_comparison_nop, file="./5_random_forest_prevalence/infect_test_train_comparison_nop.RData")
  save(infect_test_train_comparison_nop2, file="./5_random_forest_prevalence/infect_test_train_comparison_nop2.RData")
  
  infect_test_train_comparison_withp <- extractItem(allOutput = list(RF_infect_10_1, RF_infect_10_2, RF_infect_10_3), bd_type="infect", otu_type = "PA", p_types="withp", reps=c(1:10)
                                                    , extract="test_train_comparison")
  infect_test_train_comparison_withp2 <- extractItem(allOutput = list(RF_infect_10_4, RF_infect_10_5, RF_infect_10_6), bd_type="infect", otu_type = "PA", p_types="withp", reps=c(1:10)
                                                    , extract="test_train_comparison", repstart=31)
  save(infect_test_train_comparison_withp, file="./5_random_forest_prevalence/infect_test_train_comparison_withp.RData")
  save(infect_test_train_comparison_withp2, file="./5_random_forest_prevalence/infect_test_train_comparison_withp2.RData")
  
  # Getting importance
  importance_ponly <- extractItem(allOutput = list(RF_infect_10_1, RF_infect_10_2, RF_infect_10_3), bd_type="infect", otu_type = "PA", p_types="onlyp", reps=c(1:10)
                                  , extract="importance")
  importance_ponly2 <- extractItem(allOutput = list(RF_infect_10_4, RF_infect_10_5, RF_infect_10_6), bd_type="infect", otu_type = "PA", p_types="onlyp", reps=c(1:10)
                                  , extract="importance", repstart=31)
  save(importance_ponly, file="./5_random_forest_prevalence/importance_ponly.RData")
  save(importance_ponly2, file="./5_random_forest_prevalence/importance_ponly2.RData")
  


  #### PABD extract####
  PABD_test_train_comparison_ponly <- extractItem(allOutput = list(RF_PABD_10_1, RF_PABD_10_2, RF_PABD_10_3), bd_type="PABD", otu_type = "PA", p_types="onlyp", reps=c(1:10)
                                                  , extract="test_train_comparison")
  PABD_test_train_comparison_ponly2 <- extractItem(allOutput = list(RF_PABD_10_4, RF_PABD_10_5, RF_PABD_10_6), bd_type="PABD", otu_type = "PA", p_types="onlyp", reps=c(1:10)
                                                  , extract="test_train_comparison", repstart=31)
  save(PABD_test_train_comparison_ponly, file="./5_random_forest_prevalence/PABD_test_train_comparison_ponly.RData")
  save(PABD_test_train_comparison_ponly2, file="./5_random_forest_prevalence/PABD_test_train_comparison_ponly2.RData")
  
  PABD_test_train_comparison_nop <- extractItem(allOutput = list(RF_PABD_10_1, RF_PABD_10_2, RF_PABD_10_3), bd_type="PABD", otu_type = "PA", p_types="nop", reps=c(1:10)
                                                , extract="test_train_comparison")
  PABD_test_train_comparison_nop2 <- extractItem(allOutput = list(RF_PABD_10_4, RF_PABD_10_5, RF_PABD_10_6), bd_type="PABD", otu_type = "PA", p_types="nop", reps=c(1:10)
                                                , extract="test_train_comparison", repstart=31)
  save(PABD_test_train_comparison_nop, file="./5_random_forest_prevalence/PABD_test_train_comparison_nop.RData")
  save(PABD_test_train_comparison_nop2, file="./5_random_forest_prevalence/PABD_test_train_comparison_nop2.RData")
  
  PABD_test_train_comparison_withp <- extractItem(allOutput = list(RF_PABD_10_1, RF_PABD_10_2, RF_PABD_10_3), bd_type="PABD", otu_type = "PA", p_types="withp", reps=c(1:10)
                                                  , extract="test_train_comparison")
  PABD_test_train_comparison_withp2 <- extractItem(allOutput = list(RF_PABD_10_4, RF_PABD_10_5, RF_PABD_10_6), bd_type="PABD", otu_type = "PA", p_types="withp", reps=c(1:10)
                                                  , extract="test_train_comparison", repstart=31)
  save(PABD_test_train_comparison_withp, file="./5_random_forest_prevalence/PABD_test_train_comparison_withp.RData")
  save(PABD_test_train_comparison_withp2, file="./5_random_forest_prevalence/PABD_test_train_comparison_withp2.RData")
  
  # Getting importance
  importance_ponly_PABD <- extractItem(allOutput = list(RF_PABD_10_1, RF_PABD_10_2, RF_PABD_10_3), bd_type="PABD", otu_type = "PA", p_types="onlyp", reps=c(1:10)
                                       , extract="importance")
  importance_ponly_PABD2 <- extractItem(allOutput = list(RF_PABD_10_4, RF_PABD_10_5, RF_PABD_10_6), bd_type="PABD", otu_type = "PA", p_types="onlyp", reps=c(1:10)
                                       , extract="importance", repstart=31)
  save(importance_ponly_PABD, file="./5_random_forest_prevalence/importance_ponly_PABD.RData")
  save(importance_ponly_PABD2, file="./5_random_forest_prevalence/importance_ponly_PABD2.RData")
  
  
  
} else {
  load("./5_random_forest_prevalence/infect_test_train_comparison_ponly.RData")
  load("./5_random_forest_prevalence/infect_test_train_comparison_nop.RData")
  load("./5_random_forest_prevalence/infect_test_train_comparison_withp.RData")
  load("./5_random_forest_prevalence/importance_ponly.RData")
  load("./5_random_forest_prevalence/infect_test_train_comparison_ponly2.RData")
  load("./5_random_forest_prevalence/infect_test_train_comparison_nop2.RData")
  load("./5_random_forest_prevalence/infect_test_train_comparison_withp2.RData")
  load("./5_random_forest_prevalence/importance_ponly2.RData")
  
  load("./5_random_forest_prevalence/PABD_test_train_comparison_ponly.RData")
  load("./5_random_forest_prevalence/PABD_test_train_comparison_nop.RData")
  load("./5_random_forest_prevalence/PABD_test_train_comparison_withp.RData")
  load("./5_random_forest_prevalence/importance_ponly_PABD.RData")
  load("./5_random_forest_prevalence/PABD_test_train_comparison_ponly2.RData")
  load("./5_random_forest_prevalence/PABD_test_train_comparison_nop2.RData")
  load("./5_random_forest_prevalence/PABD_test_train_comparison_withp2.RData")
  load("./5_random_forest_prevalence/importance_ponly_PABD2.RData")

}

# 
# #### infect plot ####
# all_test_train_comparisons <- rbind(cbind(infect_test_train_comparison_nop, model="nop")
#                                     ,cbind(infect_test_train_comparison_ponly, model="ponly")
#                                     ,cbind(infect_test_train_comparison_withp, model="withp")
#                                     ,cbind(infect_test_train_comparison_nop2, model="nop")
#                                     ,cbind(infect_test_train_comparison_ponly2, model="ponly")
#                                     ,cbind(infect_test_train_comparison_withp2, model="withp"))
# save(all_test_train_comparisons, file="./5_random_forest_prevalence/all_test_train_comparisons.RData")
# 
# all_test_train_comparisons %>%
#   ggplot() +  geom_smooth(aes(x=Test_obs, y=Test_pred, col=model), method="lm") +
#   geom_point(aes(x=Test_obs, y=Test_pred, col=model), show.legend = FALSE)+
#   geom_abline(aes(intercept=0, slope=1)) +ylim(0,8)+ xlim(0,8)
# 
# all_test_train_comparisons %>%
#   select(indivID) %>% table()
# 
# 
# infect_test_train_comparison_ponly %>%
#   ggplot() +
#   geom_point(aes(x=Test_obs, y=Test_pred, col=indivID), show.legend = FALSE)+
#   geom_abline(aes(intercept=0, slope=1)) +ylim(0,8)+ xlim(0,8)
# 
# rbind(importance_ponly, importance_ponly2) %>%
#   group_by(taxa) %>%
#   summarise(meanIncMSE= mean(X.IncMSE), sdIncMSE= sd(X.IncMSE)) %>%
#   arrange(-meanIncMSE)
# 
# #### PABD plot ####
# 
# all_test_train_comparisons_PABD <- rbind(cbind(PABD_test_train_comparison_ponly, model="ponly")
#       ,cbind(PABD_test_train_comparison_nop, model="nop")
#       , cbind(PABD_test_train_comparison_withp, model="withp")
#       ,cbind(PABD_test_train_comparison_ponly2, model="ponly")
#       ,cbind(PABD_test_train_comparison_nop2, model="nop")
#       , cbind(PABD_test_train_comparison_withp2, model="withp"))
# save(all_test_train_comparisons_PABD, file="./5_random_forest_prevalence/all_test_train_comparisons_PABD.RData")
# 
# all_test_train_comparisons_PABD%>%
#   group_by(model,rep) %>%
#   summarize(propCorrect=mean(Test_pred==Test_obs)) %>%
#   ggplot() +geom_violin(aes(x=model, y=propCorrect))+
#   geom_point(aes(x=model, y=propCorrect), position=position_jitter(height=0, width=0.25))
# 
# all_test_train_comparisons_PABD%>%
#   group_by(model,indivID) %>%
#   summarize(propCorrect=mean(Test_pred==Test_obs)) %>%
#   separate(indivID, into=c("species","indiv") , remove=FALSE) %>%
#   ggplot() +
#   # geom_violin(aes(x=indivID, y=propCorrect))+
#   geom_point(aes(x=model, y=propCorrect, col=indiv, pch=species), position=position_jitter(height=0, width=0.25))+
#   facet_wrap(.~model, scales = "free_x")
# 
# 
# rbind(importance_ponly_PABD, importance_ponly_PABD2) %>%
#   group_by(taxa) %>%
#   summarise(meanDecrAccuracy= mean(MeanDecreaseAccuracy), sdDecrAccuracy= sd(MeanDecreaseAccuracy)) %>%
#   arrange(-meanDecrAccuracy)
# 
# 


