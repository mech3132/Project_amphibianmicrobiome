#!/bin/bash Rscript

#### Loading data #####
library(randomForest)
library(tidyverse)
load("3_5sp_mapping_otu_downstream/mf_alt_filt_final.RData")
load("3_5sp_mapping_otu_downstream/otu_filt.RData")
load("3_5sp_mapping_otu_downstream/taxonomy.RData")
load("4_Bayesian_models/all_p.RData")
source("msc_randomforest_codebits.R")
dir.create("5_random_forest_LOO")

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
#### INFECT SECTION  ####
allIndiv <- mf_alt_filt_final %>% filter(PABD==1) %>% pull(indivID) %>% unique()
total <- length(allIndiv)
prop <- round(14/15,2)

#### Infect, ponly, with species, no zeros ####
RERUN_RF <- FALSE
if ( RERUN_RF ) {
  ## Set params
  p_type="onlyp"
  otu_type="NA"
  bd_type="infect"
  mf_temp=NA
  keepSp=TRUE
  # create progress bar
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  i <- 1
  RF_infect_onlyp_wsp_LOO <- list()
  # Get a formatted mapping file
  get_mf_pred_temp_bd <- get_mf_pred_temp(p_type=p_type
                                          , otu_type = otu_type
                                          , bd_type=bd_type
                                          , p=p
                                          , keepSp=keepSp)
  mf_pred_temp_bd <- get_mf_pred_temp_bd[[1]]
  indiv_list <- get_mf_pred_temp_bd[[2]]
  # Mutate infection depending on question
  Error_metric <- ifelse(bd_type=="PABD","MeanDecreaseAccuracy","X.IncMSE" )
  for ( LOO in c(1:length(allIndiv)) ) {
    RF_infect_onlyp_wsp_LOO[[paste(otu_type, p_type, bd_type, prop, LOO, sep="_")]] <- list()
    # Pick training and validation set
    trainSet <- mf_pred_temp_bd[-LOO,]
    testSet <- mf_pred_temp_bd[LOO,]
    testSet_indiv <- as.character(unlist(indiv_list[LOO,]))
    
    # Run random forest model
    RF <- randomForest(response ~ ., data=trainSet,importance=TRUE)
    # Run prediction
    testSet_predictions <- predict(RF, testSet, type="class")
    # importance of factors
    importance_RF <- data.frame(importance(RF)) %>%
      mutate(taxa=rownames(importance(RF))) %>%
      select(taxa, everything()) %>%
      arrange(-get(paste0(Error_metric)))
    
    # Add inhibitory to importance
    importance_RF[,c("taxonomy","inhibitory")] <- Taxa[match(importance_RF$taxa, Taxa$Sequence), c("Taxa","inhibitory")]
    importance_RF <- importance_RF %>%
      mutate(taxonomy=ifelse(is.na(taxonomy), taxa, taxonomy))
    
    # Save data
    RF_infect_onlyp_wsp_LOO[[paste(otu_type, p_type, bd_type, prop,LOO, sep="_")]][["RF"]] <- RF
    RF_infect_onlyp_wsp_LOO[[paste(otu_type, p_type, bd_type, prop,LOO, sep="_")]][["importance"]] <- importance_RF
    RF_infect_onlyp_wsp_LOO[[paste(otu_type, p_type, bd_type, prop,LOO, sep="_")]][["test_train_comparison"]] <- data.frame(indivID=testSet_indiv, Test_pred=testSet_predictions, Test_obs=testSet$response)
    
    # Insert into data summary plot
    #------- PROGRESS BAR ------- #
    Sys.sleep(0.1)
    i <- i+1
    # update progress bar
    setTxtProgressBar(pb, i)
    #-----------------------------#
    # }
    
  }
  
  close(pb)
  # load("./5_random_forest_LOO/RF_infect_onlyp_wsp_LOO.RData")
  save(RF_infect_onlyp_wsp_LOO,file="./5_random_forest_LOO/RF_infect_onlyp_wsp_LOO.RData" )
  compare_infect_onlyp_wsp_LOO <- extractItemLOO(RF_infect_onlyp_wsp_LOO, bd_type = bd_type, p_type=p_type, otu_type = otu_type,extract = "test_train_comparison" )
  save(compare_infect_onlyp_wsp_LOO,file="./5_random_forest_LOO/compare_infect_onlyp_wsp_LOO.RData" )
  importance_infect_onlyp_wsp_LOO <-extractItemLOO(RF_infect_onlyp_wsp_LOO, bd_type = bd_type, p_type=p_type, otu_type = otu_type,extract = "importance")
  save(importance_infect_onlyp_wsp_LOO,file="./5_random_forest_LOO/importance_infect_onlyp_wsp_LOO.RData" )
  MSE_infect_onlyp_wsp_LOO <- extractTrainErrorLOO(RF_infect_onlyp_wsp_LOO, bd_type = bd_type, p_type=p_type, otu_type = otu_type)
  save(MSE_infect_onlyp_wsp_LOO,file="./5_random_forest_LOO/MSE_infect_onlyp_wsp_LOO.RData" )
}
remove(RF_infect_onlyp_wsp_LOO)

#### Infect, ponly, no species, no zeros ####
if ( RERUN_RF ) {
  ## Set params
  p_type="onlyp"
  otu_type="NA"
  bd_type="infect"
  mf_temp=NA
  keepSp=FALSE
  # create progress bar
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  i <- 1
  RF_infect_onlyp_nosp_LOO <- list()
  # Get a formatted mapping file
  get_mf_pred_temp_bd <- get_mf_pred_temp(p_type=p_type
                                          , otu_type = otu_type
                                          , bd_type=bd_type
                                          , p=p
                                          , keepSp=keepSp)
  mf_pred_temp_bd <- get_mf_pred_temp_bd[[1]]
  indiv_list <- get_mf_pred_temp_bd[[2]]
  # Mutate infection depending on question
  Error_metric <- ifelse(bd_type=="PABD","MeanDecreaseAccuracy","X.IncMSE" )
  for ( LOO in c(1:length(allIndiv)) ) {
    RF_infect_onlyp_nosp_LOO[[paste(otu_type, p_type, bd_type, prop, LOO, sep="_")]] <- list()
    # Pick training and validation set
    trainSet <- mf_pred_temp_bd[-LOO,]
    testSet <- mf_pred_temp_bd[LOO,]
    testSet_indiv <- as.character(unlist(indiv_list[LOO,]))
    
    # Run random forest model
    RF <- randomForest(response ~ ., data=trainSet,importance=TRUE)
    # Run prediction
    testSet_predictions <- predict(RF, testSet, type="class")
    # importance of factors
    importance_RF <- data.frame(importance(RF)) %>%
      mutate(taxa=rownames(importance(RF))) %>%
      select(taxa, everything()) %>%
      arrange(-get(paste0(Error_metric)))
    
    # Add inhibitory to importance
    importance_RF[,c("taxonomy","inhibitory")] <- Taxa[match(importance_RF$taxa, Taxa$Sequence), c("Taxa","inhibitory")]
    importance_RF <- importance_RF %>%
      mutate(taxonomy=ifelse(is.na(taxonomy), taxa, taxonomy))
    
    # Save data
    RF_infect_onlyp_nosp_LOO[[paste(otu_type, p_type, bd_type, prop,LOO, sep="_")]][["RF"]] <- RF
    RF_infect_onlyp_nosp_LOO[[paste(otu_type, p_type, bd_type, prop,LOO, sep="_")]][["importance"]] <- importance_RF
    RF_infect_onlyp_nosp_LOO[[paste(otu_type, p_type, bd_type, prop,LOO, sep="_")]][["test_train_comparison"]] <- data.frame(indivID=testSet_indiv, Test_pred=testSet_predictions, Test_obs=testSet$response)
    
    # Insert into data summary plot
    #------- PROGRESS BAR ------- #
    Sys.sleep(0.1)
    i <- i+1
    # update progress bar
    setTxtProgressBar(pb, i)
    #-----------------------------#
    # }
    
  }
  
  close(pb)
  # load("./5_random_forest_LOO/RF_infect_onlyp_nosp_LOO.RData")
  save(RF_infect_onlyp_nosp_LOO,file="./5_random_forest_LOO/RF_infect_onlyp_nosp_LOO.RData" )
  compare_infect_onlyp_nosp_LOO <- extractItemLOO(RF_infect_onlyp_nosp_LOO, bd_type = bd_type, p_type=p_type, otu_type = otu_type,extract = "test_train_comparison" )
  save(compare_infect_onlyp_nosp_LOO,file="./5_random_forest_LOO/compare_infect_onlyp_nosp_LOO.RData" )
  importance_infect_onlyp_nosp_LOO <-extractItemLOO(RF_infect_onlyp_nosp_LOO, bd_type = bd_type, p_type=p_type, otu_type = otu_type,extract = "importance")
  save(importance_infect_onlyp_nosp_LOO,file="./5_random_forest_LOO/importance_infect_onlyp_nosp_LOO.RData" )
  MSE_infect_onlyp_nosp_LOO <- extractTrainErrorLOO(RF_infect_onlyp_nosp_LOO, bd_type = bd_type, p_type=p_type, otu_type = otu_type)
  save(MSE_infect_onlyp_nosp_LOO,file="./5_random_forest_LOO/MSE_infect_onlyp_nosp_LOO.RData" )
}
remove(RF_infect_onlyp_nosp_LOO)


#### Infect, count, with species, no zeros ####
if ( RERUN_RF ) {
  ## Set params
  p_type="nop"
  otu_type="count"
  bd_type="infect"
  keepSp=TRUE
  # create progress bar
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  i <- 1
  RF_infect_count_wsp_LOO <- list()
  # Get a formatted mapping file
  get_mf_pred_temp_bd <- get_mf_pred_temp(p_type=p_type
                                          , otu_type = otu_type
                                          , bd_type=bd_type
                                          , p=p
                                          , keepSp=keepSp)
  mf_pred_temp_bd <- get_mf_pred_temp_bd[[1]]
  indiv_list <- get_mf_pred_temp_bd[[2]]
  # Mutate infection depending on question
  Error_metric <- ifelse(bd_type=="PABD","MeanDecreaseAccuracy","X.IncMSE" )
  for ( LOO in c(1:length(allIndiv)) ) {
    RF_infect_count_wsp_LOO[[paste(otu_type, p_type, bd_type, prop, LOO, sep="_")]] <- list()
    # Pick training and validation set
    trainSet <- mf_pred_temp_bd[-LOO,]
    testSet <- mf_pred_temp_bd[LOO,]
    testSet_indiv <- as.character(unlist(indiv_list[LOO,]))
    
    # Run random forest model
    RF <- randomForest(response ~ ., data=trainSet,importance=TRUE)
    # Run prediction
    testSet_predictions <- predict(RF, testSet, type="class")
    # importance of factors
    importance_RF <- data.frame(importance(RF)) %>%
      mutate(taxa=rownames(importance(RF))) %>%
      select(taxa, everything()) %>%
      arrange(-get(paste0(Error_metric)))
    
    # Add inhibitory to importance
    importance_RF[,c("taxonomy","inhibitory")] <- Taxa[match(importance_RF$taxa, Taxa$Sequence), c("Taxa","inhibitory")]
    importance_RF <- importance_RF %>%
      mutate(taxonomy=ifelse(is.na(taxonomy), taxa, taxonomy))
    
    # Save data
    RF_infect_count_wsp_LOO[[paste(otu_type, p_type, bd_type, prop,LOO, sep="_")]][["RF"]] <- RF
    RF_infect_count_wsp_LOO[[paste(otu_type, p_type, bd_type, prop,LOO, sep="_")]][["importance"]] <- importance_RF
    RF_infect_count_wsp_LOO[[paste(otu_type, p_type, bd_type, prop,LOO, sep="_")]][["test_train_comparison"]] <- data.frame(indivID=testSet_indiv, Test_pred=testSet_predictions, Test_obs=testSet$response)
    
    # Insert into data summary plot
    #------- PROGRESS BAR ------- #
    Sys.sleep(0.1)
    i <- i+1
    # update progress bar
    setTxtProgressBar(pb, i)
    #-----------------------------#
    # }
    
  }
  
  close(pb)
  # load("./5_random_forest_LOO/RF_infect_count_wsp_LOO.RData")
  save(RF_infect_count_wsp_LOO,file="./5_random_forest_LOO/RF_infect_count_wsp_LOO.RData" )
  compare_infect_count_wsp_LOO <- extractItemLOO(RF_infect_count_wsp_LOO, bd_type = bd_type, p_type=p_type, otu_type = otu_type,extract = "test_train_comparison" )
  save(compare_infect_count_wsp_LOO,file="./5_random_forest_LOO/compare_infect_count_wsp_LOO.RData" )
  importance_infect_count_wsp_LOO <-extractItemLOO(RF_infect_count_wsp_LOO, bd_type = bd_type, p_type=p_type, otu_type = otu_type,extract = "importance")
  save(importance_infect_count_wsp_LOO,file="./5_random_forest_LOO/importance_infect_count_wsp_LOO.RData" )
  MSE_infect_count_wsp_LOO <- extractTrainErrorLOO(RF_infect_count_wsp_LOO, bd_type = bd_type, p_type=p_type, otu_type = otu_type)
  save(MSE_infect_count_wsp_LOO,file="./5_random_forest_LOO/MSE_infect_count_wsp_LOO.RData" )
}
remove(RF_infect_count_wsp_LOO)


#### Infect, count, no species, no zeros ####
if ( RERUN_RF ) {
  ## Set params
  p_type="nop"
  otu_type="count"
  bd_type="infect"
  keepSp=FALSE
  # create progress bar
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  i <- 1
  RF_infect_count_nosp_LOO <- list()
  # Get a formatted mapping file
  get_mf_pred_temp_bd <- get_mf_pred_temp(p_type=p_type
                                          , otu_type = otu_type
                                          , bd_type=bd_type
                                          , p=p
                                          , keepSp=keepSp)
  mf_pred_temp_bd <- get_mf_pred_temp_bd[[1]]
  indiv_list <- get_mf_pred_temp_bd[[2]]
  # Mutate infection depending on question
  Error_metric <- ifelse(bd_type=="PABD","MeanDecreaseAccuracy","X.IncMSE" )
  for ( LOO in c(1:length(allIndiv)) ) {
    RF_infect_count_nosp_LOO[[paste(otu_type, p_type, bd_type, prop, LOO, sep="_")]] <- list()
    # Pick training and validation set
    trainSet <- mf_pred_temp_bd[-LOO,]
    testSet <- mf_pred_temp_bd[LOO,]
    testSet_indiv <- as.character(unlist(indiv_list[LOO,]))
    
    # Run random forest model
    RF <- randomForest(response ~ ., data=trainSet,importance=TRUE)
    # Run prediction
    testSet_predictions <- predict(RF, testSet, type="class")
    # importance of factors
    importance_RF <- data.frame(importance(RF)) %>%
      mutate(taxa=rownames(importance(RF))) %>%
      select(taxa, everything()) %>%
      arrange(-get(paste0(Error_metric)))
    
    # Add inhibitory to importance
    importance_RF[,c("taxonomy","inhibitory")] <- Taxa[match(importance_RF$taxa, Taxa$Sequence), c("Taxa","inhibitory")]
    importance_RF <- importance_RF %>%
      mutate(taxonomy=ifelse(is.na(taxonomy), taxa, taxonomy))
    
    # Save data
    RF_infect_count_nosp_LOO[[paste(otu_type, p_type, bd_type, prop,LOO, sep="_")]][["RF"]] <- RF
    RF_infect_count_nosp_LOO[[paste(otu_type, p_type, bd_type, prop,LOO, sep="_")]][["importance"]] <- importance_RF
    RF_infect_count_nosp_LOO[[paste(otu_type, p_type, bd_type, prop,LOO, sep="_")]][["test_train_comparison"]] <- data.frame(indivID=testSet_indiv, Test_pred=testSet_predictions, Test_obs=testSet$response)
    
    # Insert into data summary plot
    #------- PROGRESS BAR ------- #
    Sys.sleep(0.1)
    i <- i+1
    # update progress bar
    setTxtProgressBar(pb, i)
    #-----------------------------#
    # }
    
  }
  
  close(pb)
  # load("./5_random_forest_LOO/RF_infect_count_nosp_LOO.RData")
  save(RF_infect_count_nosp_LOO,file="./5_random_forest_LOO/RF_infect_count_nosp_LOO.RData" )
  compare_infect_count_nosp_LOO <- extractItemLOO(RF_infect_count_nosp_LOO, bd_type = bd_type, p_type=p_type, otu_type = otu_type,extract = "test_train_comparison" )
  save(compare_infect_count_nosp_LOO,file="./5_random_forest_LOO/compare_infect_count_nosp_LOO.RData" )
  importance_infect_count_nosp_LOO <-extractItemLOO(RF_infect_count_nosp_LOO, bd_type = bd_type, p_type=p_type, otu_type = otu_type,extract = "importance")
  save(importance_infect_count_nosp_LOO,file="./5_random_forest_LOO/importance_infect_count_nosp_LOO.RData" )
  MSE_infect_count_nosp_LOO <- extractTrainErrorLOO(RF_infect_count_nosp_LOO, bd_type = bd_type, p_type=p_type, otu_type = otu_type)
  save(MSE_infect_count_nosp_LOO,file="./5_random_forest_LOO/MSE_infect_count_nosp_LOO.RData" )
}
remove(RF_infect_count_nosp_LOO)

#### Infect, PA, no species, no zeros ####
if ( RERUN_RF ) {
  ## Set params
  p_type="nop"
  otu_type="PA"
  bd_type="infect"
  keepSp=FALSE
  # create progress bar
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  i <- 1
  RF_infect_PA_nosp_LOO <- list()
  # Get a formatted mapping file
  get_mf_pred_temp_bd <- get_mf_pred_temp(p_type=p_type
                                          , otu_type = otu_type
                                          , bd_type=bd_type
                                          , p=p
                                          , keepSp=keepSp)
  mf_pred_temp_bd <- get_mf_pred_temp_bd[[1]]
  indiv_list <- get_mf_pred_temp_bd[[2]]
  # Mutate infection depending on question
  Error_metric <- ifelse(bd_type=="PABD","MeanDecreaseAccuracy","X.IncMSE" )
  for ( LOO in c(1:length(allIndiv)) ) {
    RF_infect_PA_nosp_LOO[[paste(otu_type, p_type, bd_type, prop, LOO, sep="_")]] <- list()
    # Pick training and validation set
    trainSet <- mf_pred_temp_bd[-LOO,]
    testSet <- mf_pred_temp_bd[LOO,]
    testSet_indiv <- as.character(unlist(indiv_list[LOO,]))
    
    # Run random forest model
    RF <- randomForest(response ~ ., data=trainSet,importance=TRUE)
    # Run prediction
    testSet_predictions <- predict(RF, testSet, type="class")
    # importance of factors
    importance_RF <- data.frame(importance(RF)) %>%
      mutate(taxa=rownames(importance(RF))) %>%
      select(taxa, everything()) %>%
      arrange(-get(paste0(Error_metric)))
    
    # Add inhibitory to importance
    importance_RF[,c("taxonomy","inhibitory")] <- Taxa[match(importance_RF$taxa, Taxa$Sequence), c("Taxa","inhibitory")]
    importance_RF <- importance_RF %>%
      mutate(taxonomy=ifelse(is.na(taxonomy), taxa, taxonomy))
    
    # Save data
    RF_infect_PA_nosp_LOO[[paste(otu_type, p_type, bd_type, prop,LOO, sep="_")]][["RF"]] <- RF
    RF_infect_PA_nosp_LOO[[paste(otu_type, p_type, bd_type, prop,LOO, sep="_")]][["importance"]] <- importance_RF
    RF_infect_PA_nosp_LOO[[paste(otu_type, p_type, bd_type, prop,LOO, sep="_")]][["test_train_comparison"]] <- data.frame(indivID=testSet_indiv, Test_pred=testSet_predictions, Test_obs=testSet$response)
    
    # Insert into data summary plot
    #------- PROGRESS BAR ------- #
    Sys.sleep(0.1)
    i <- i+1
    # update progress bar
    setTxtProgressBar(pb, i)
    #-----------------------------#
    # }
    
  }
  
  close(pb)
  # load("./5_random_forest_LOO/RF_infect_PA_nosp_LOO.RData")
  save(RF_infect_PA_nosp_LOO,file="./5_random_forest_LOO/RF_infect_PA_nosp_LOO.RData" )
  compare_infect_PA_nosp_LOO <- extractItemLOO(RF_infect_PA_nosp_LOO, bd_type = bd_type, p_type=p_type, otu_type = otu_type,extract = "test_train_comparison" )
  save(compare_infect_PA_nosp_LOO,file="./5_random_forest_LOO/compare_infect_PA_nosp_LOO.RData" )
  importance_infect_PA_nosp_LOO <-extractItemLOO(RF_infect_PA_nosp_LOO, bd_type = bd_type, p_type=p_type, otu_type = otu_type,extract = "importance")
  save(importance_infect_PA_nosp_LOO,file="./5_random_forest_LOO/importance_infect_PA_nosp_LOO.RData" )
  MSE_infect_PA_nosp_LOO <- extractTrainErrorLOO(RF_infect_PA_nosp_LOO, bd_type = bd_type, p_type=p_type, otu_type = otu_type)
  save(MSE_infect_PA_nosp_LOO,file="./5_random_forest_LOO/MSE_infect_PA_nosp_LOO.RData" )
}
remove(RF_infect_PA_nosp_LOO)

#### Infect, PA, with species, no zeros ####
if ( RERUN_RF ) {
  ## Set params
  p_type="nop"
  otu_type="PA"
  bd_type="infect"
  keepSp=TRUE
  # create progress bar
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  i <- 1
  RF_infect_PA_wsp_LOO <- list()
  # Get a formatted mapping file
  get_mf_pred_temp_bd <- get_mf_pred_temp(p_type=p_type
                                          , otu_type = otu_type
                                          , bd_type=bd_type
                                          , p=p
                                          , keepSp=keepSp)
  mf_pred_temp_bd <- get_mf_pred_temp_bd[[1]]
  indiv_list <- get_mf_pred_temp_bd[[2]]
  # Mutate infection depending on question
  Error_metric <- ifelse(bd_type=="PABD","MeanDecreaseAccuracy","X.IncMSE" )
  for ( LOO in c(1:length(allIndiv)) ) {
    RF_infect_PA_wsp_LOO[[paste(otu_type, p_type, bd_type, prop, LOO, sep="_")]] <- list()
    # Pick training and validation set
    trainSet <- mf_pred_temp_bd[-LOO,]
    testSet <- mf_pred_temp_bd[LOO,]
    testSet_indiv <- as.character(unlist(indiv_list[LOO,]))
    
    # Run random forest model
    RF <- randomForest(response ~ ., data=trainSet,importance=TRUE)
    # Run prediction
    testSet_predictions <- predict(RF, testSet, type="class")
    # importance of factors
    importance_RF <- data.frame(importance(RF)) %>%
      mutate(taxa=rownames(importance(RF))) %>%
      select(taxa, everything()) %>%
      arrange(-get(paste0(Error_metric)))
    
    # Add inhibitory to importance
    importance_RF[,c("taxonomy","inhibitory")] <- Taxa[match(importance_RF$taxa, Taxa$Sequence), c("Taxa","inhibitory")]
    importance_RF <- importance_RF %>%
      mutate(taxonomy=ifelse(is.na(taxonomy), taxa, taxonomy))
    
    # Save data
    RF_infect_PA_wsp_LOO[[paste(otu_type, p_type, bd_type, prop,LOO, sep="_")]][["RF"]] <- RF
    RF_infect_PA_wsp_LOO[[paste(otu_type, p_type, bd_type, prop,LOO, sep="_")]][["importance"]] <- importance_RF
    RF_infect_PA_wsp_LOO[[paste(otu_type, p_type, bd_type, prop,LOO, sep="_")]][["test_train_comparison"]] <- data.frame(indivID=testSet_indiv, Test_pred=testSet_predictions, Test_obs=testSet$response)
    
    # Insert into data summary plot
    #------- PROGRESS BAR ------- #
    Sys.sleep(0.1)
    i <- i+1
    # update progress bar
    setTxtProgressBar(pb, i)
    #-----------------------------#
    # }
    
  }
  
  close(pb)
  # load("./5_random_forest_LOO/RF_infect_PA_wsp_LOO.RData")
  save(RF_infect_PA_wsp_LOO,file="./5_random_forest_LOO/RF_infect_PA_wsp_LOO.RData" )
  compare_infect_PA_wsp_LOO <- extractItemLOO(RF_infect_PA_wsp_LOO, bd_type = bd_type, p_type=p_type, otu_type = otu_type,extract = "test_train_comparison" )
  save(compare_infect_PA_wsp_LOO,file="./5_random_forest_LOO/compare_infect_PA_wsp_LOO.RData" )
  importance_infect_PA_wsp_LOO <-extractItemLOO(RF_infect_PA_wsp_LOO, bd_type = bd_type, p_type=p_type, otu_type = otu_type,extract = "importance")
  save(importance_infect_PA_wsp_LOO,file="./5_random_forest_LOO/importance_infect_PA_wsp_LOO.RData" )
  MSE_infect_PA_wsp_LOO <- extractTrainErrorLOO(RF_infect_PA_wsp_LOO, bd_type = bd_type, p_type=p_type, otu_type = otu_type)
  save(MSE_infect_PA_wsp_LOO,file="./5_random_forest_LOO/MSE_infect_PA_wsp_LOO.RData" )
}
remove(RF_infect_PA_wsp_LOO)

#### Infect, all count, with species, no zeros ####
if ( RERUN_RF ) {
  ## Set params
  p_type="withp"
  otu_type="count"
  bd_type="infect"
  keepSp=TRUE
  # create progress bar
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  i <- 1
  RF_infect_withpcount_wsp_LOO <- list()
  # Get a formatted mapping file
  get_mf_pred_temp_bd <- get_mf_pred_temp(p_type=p_type
                                          , otu_type = otu_type
                                          , bd_type=bd_type
                                          , p=p
                                          , keepSp=keepSp)
  mf_pred_temp_bd <- get_mf_pred_temp_bd[[1]]
  indiv_list <- get_mf_pred_temp_bd[[2]]
  # Mutate infection depending on question
  Error_metric <- ifelse(bd_type=="PABD","MeanDecreaseAccuracy","X.IncMSE" )
  for ( LOO in c(1:length(allIndiv)) ) {
    RF_infect_withpcount_wsp_LOO[[paste(otu_type, p_type, bd_type, prop, LOO, sep="_")]] <- list()
    # Pick training and validation set
    trainSet <- mf_pred_temp_bd[-LOO,]
    testSet <- mf_pred_temp_bd[LOO,]
    testSet_indiv <- as.character(unlist(indiv_list[LOO,]))
    
    # Run random forest model
    RF <- randomForest(response ~ ., data=trainSet,importance=TRUE)
    # Run prediction
    testSet_predictions <- predict(RF, testSet, type="class")
    # importance of factors
    importance_RF <- data.frame(importance(RF)) %>%
      mutate(taxa=rownames(importance(RF))) %>%
      select(taxa, everything()) %>%
      arrange(-get(paste0(Error_metric)))
    
    # Add inhibitory to importance
    importance_RF[,c("taxonomy","inhibitory")] <- Taxa[match(importance_RF$taxa, Taxa$Sequence), c("Taxa","inhibitory")]
    importance_RF <- importance_RF %>%
      mutate(taxonomy=ifelse(is.na(taxonomy), taxa, taxonomy))
    
    # Save data
    RF_infect_withpcount_wsp_LOO[[paste(otu_type, p_type, bd_type, prop,LOO, sep="_")]][["RF"]] <- RF
    RF_infect_withpcount_wsp_LOO[[paste(otu_type, p_type, bd_type, prop,LOO, sep="_")]][["importance"]] <- importance_RF
    RF_infect_withpcount_wsp_LOO[[paste(otu_type, p_type, bd_type, prop,LOO, sep="_")]][["test_train_comparison"]] <- data.frame(indivID=testSet_indiv, Test_pred=testSet_predictions, Test_obs=testSet$response)
    
    # Insert into data summary plot
    #------- PROGRESS BAR ------- #
    Sys.sleep(0.1)
    i <- i+1
    # update progress bar
    setTxtProgressBar(pb, i)
    #-----------------------------#
    # }
    
  }
  
  close(pb)
  # load("./5_random_forest_LOO/RF_infect_withpcount_wsp_LOO.RData")
  save(RF_infect_withpcount_wsp_LOO,file="./5_random_forest_LOO/RF_infect_withpcount_wsp_LOO.RData" )
  compare_infect_withpcount_wsp_LOO <- extractItemLOO(RF_infect_withpcount_wsp_LOO, bd_type = bd_type, p_type=p_type, otu_type = otu_type,extract = "test_train_comparison" )
  save(compare_infect_withpcount_wsp_LOO,file="./5_random_forest_LOO/compare_infect_withpcount_wsp_LOO.RData" )
  importance_infect_withpcount_wsp_LOO <-extractItemLOO(RF_infect_withpcount_wsp_LOO, bd_type = bd_type, p_type=p_type, otu_type = otu_type,extract = "importance")
  save(importance_infect_withpcount_wsp_LOO,file="./5_random_forest_LOO/importance_infect_withpcount_wsp_LOO.RData" )
  MSE_infect_withpcount_wsp_LOO <- extractTrainErrorLOO(RF_infect_withpcount_wsp_LOO, bd_type = bd_type, p_type=p_type, otu_type = otu_type)
  save(MSE_infect_withpcount_wsp_LOO,file="./5_random_forest_LOO/MSE_infect_withpcount_wsp_LOO.RData" )
}
remove(RF_infect_withpcount_wsp_LOO)


#### Infect, all PA, with species, no zeros ####
if ( RERUN_RF ) {
  ## Set params
  p_type="withp"
  otu_type="PA"
  bd_type="infect"
  keepSp=TRUE
  # create progress bar
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  i <- 1
  RF_infect_withpPA_wsp_LOO <- list()
  # Get a formatted mapping file
  get_mf_pred_temp_bd <- get_mf_pred_temp(p_type=p_type
                                          , otu_type = otu_type
                                          , bd_type=bd_type
                                          , p=p
                                          , keepSp=keepSp)
  mf_pred_temp_bd <- get_mf_pred_temp_bd[[1]]
  indiv_list <- get_mf_pred_temp_bd[[2]]
  # Mutate infection depending on question
  Error_metric <- ifelse(bd_type=="PABD","MeanDecreaseAccuracy","X.IncMSE" )
  for ( LOO in c(1:length(allIndiv)) ) {
    RF_infect_withpPA_wsp_LOO[[paste(otu_type, p_type, bd_type, prop, LOO, sep="_")]] <- list()
    # Pick training and validation set
    trainSet <- mf_pred_temp_bd[-LOO,]
    testSet <- mf_pred_temp_bd[LOO,]
    testSet_indiv <- as.character(unlist(indiv_list[LOO,]))
    
    # Run random forest model
    RF <- randomForest(response ~ ., data=trainSet,importance=TRUE)
    # Run prediction
    testSet_predictions <- predict(RF, testSet, type="class")
    # importance of factors
    importance_RF <- data.frame(importance(RF)) %>%
      mutate(taxa=rownames(importance(RF))) %>%
      select(taxa, everything()) %>%
      arrange(-get(paste0(Error_metric)))
    
    # Add inhibitory to importance
    importance_RF[,c("taxonomy","inhibitory")] <- Taxa[match(importance_RF$taxa, Taxa$Sequence), c("Taxa","inhibitory")]
    importance_RF <- importance_RF %>%
      mutate(taxonomy=ifelse(is.na(taxonomy), taxa, taxonomy))
    
    # Save data
    RF_infect_withpPA_wsp_LOO[[paste(otu_type, p_type, bd_type, prop,LOO, sep="_")]][["RF"]] <- RF
    RF_infect_withpPA_wsp_LOO[[paste(otu_type, p_type, bd_type, prop,LOO, sep="_")]][["importance"]] <- importance_RF
    RF_infect_withpPA_wsp_LOO[[paste(otu_type, p_type, bd_type, prop,LOO, sep="_")]][["test_train_comparison"]] <- data.frame(indivID=testSet_indiv, Test_pred=testSet_predictions, Test_obs=testSet$response)
    
    # Insert into data summary plot
    #------- PROGRESS BAR ------- #
    Sys.sleep(0.1)
    i <- i+1
    # update progress bar
    setTxtProgressBar(pb, i)
    #-----------------------------#
    # }
    
  }
  
  close(pb)
  # load("./5_random_forest_LOO/RF_infect_withpPA_wsp_LOO.RData")
  save(RF_infect_withpPA_wsp_LOO,file="./5_random_forest_LOO/RF_infect_withpPA_wsp_LOO.RData" )
  compare_infect_withpPA_wsp_LOO <- extractItemLOO(RF_infect_withpPA_wsp_LOO, bd_type = bd_type, p_type=p_type, otu_type = otu_type,extract = "test_train_comparison" )
  save(compare_infect_withpPA_wsp_LOO,file="./5_random_forest_LOO/compare_infect_withpPA_wsp_LOO.RData" )
  importance_infect_withpPA_wsp_LOO <-extractItemLOO(RF_infect_withpPA_wsp_LOO, bd_type = bd_type, p_type=p_type, otu_type = otu_type,extract = "importance")
  save(importance_infect_withpPA_wsp_LOO,file="./5_random_forest_LOO/importance_infect_withpPA_wsp_LOO.RData" )
  MSE_infect_withpPA_wsp_LOO <- extractTrainErrorLOO(RF_infect_withpPA_wsp_LOO, bd_type = bd_type, p_type=p_type, otu_type = otu_type)
  save(MSE_infect_withpPA_wsp_LOO,file="./5_random_forest_LOO/MSE_infect_withpPA_wsp_LOO.RData" )
}
remove(RF_infect_withpPA_wsp_LOO)

#### PABD SECTION  ####
allIndiv <- mf_alt_filt_final %>%  filter(Bd_exposure=="Bd-exposed") %>%pull(indivID) %>% unique()
total <- length(allIndiv)
prop <- round(21/22,2)


#### PABD, ponly, with species ####
if ( RERUN_RF ) {
  ## Set params
  p_type="onlyp"
  otu_type="NA"
  bd_type="PABD"
  mf_temp=NA
  keepSp=TRUE
  # create progress bar
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  i <- 1
  RF_PABD_onlyp_wsp_LOO <- list()
  # Get a formatted mapping file
  get_mf_pred_temp_bd <- get_mf_pred_temp(p_type=p_type
                                          , otu_type = otu_type
                                          , bd_type=bd_type
                                          , p=p
                                          , keepSp=keepSp)
  mf_pred_temp_bd <- get_mf_pred_temp_bd[[1]]
  indiv_list <- get_mf_pred_temp_bd[[2]]
  # Mutate infection depending on question
  Error_metric <- ifelse(bd_type=="PABD","MeanDecreaseAccuracy","X.IncMSE" )
  for ( LOO in c(1:length(allIndiv)) ) {
    RF_PABD_onlyp_wsp_LOO[[paste(otu_type, p_type, bd_type, prop, LOO, sep="_")]] <- list()
    # Pick training and validation set
    trainSet <- mf_pred_temp_bd[-LOO,]
    testSet <- mf_pred_temp_bd[LOO,]
    testSet_indiv <- as.character(unlist(indiv_list[LOO,]))
    
    # Run random forest model
    RF <- randomForest(response ~ ., data=trainSet,importance=TRUE)
    # Run prediction
    testSet_predictions <- predict(RF, testSet, type="class")
    # importance of factors
    importance_RF <- data.frame(importance(RF)) %>%
      mutate(taxa=rownames(importance(RF))) %>%
      select(taxa, everything()) %>%
      arrange(-get(paste0(Error_metric)))
    
    # Add inhibitory to importance
    importance_RF[,c("taxonomy","inhibitory")] <- Taxa[match(importance_RF$taxa, Taxa$Sequence), c("Taxa","inhibitory")]
    importance_RF <- importance_RF %>%
      mutate(taxonomy=ifelse(is.na(taxonomy), taxa, taxonomy))
    
    # Save data
    RF_PABD_onlyp_wsp_LOO[[paste(otu_type, p_type, bd_type, prop,LOO, sep="_")]][["RF"]] <- RF
    RF_PABD_onlyp_wsp_LOO[[paste(otu_type, p_type, bd_type, prop,LOO, sep="_")]][["importance"]] <- importance_RF
    RF_PABD_onlyp_wsp_LOO[[paste(otu_type, p_type, bd_type, prop,LOO, sep="_")]][["test_train_comparison"]] <- data.frame(indivID=testSet_indiv, Test_pred=testSet_predictions, Test_obs=testSet$response)
    
    # Insert into data summary plot
    #------- PROGRESS BAR ------- #
    Sys.sleep(0.1)
    i <- i+1
    # update progress bar
    setTxtProgressBar(pb, i)
    #-----------------------------#
    # }
    
  }
  
  close(pb)
  # load("./5_random_forest_LOO/RF_PABD_onlyp_wsp_LOO.RData")
  save(RF_PABD_onlyp_wsp_LOO,file="./5_random_forest_LOO/RF_PABD_onlyp_wsp_LOO.RData" )
  compare_PABD_onlyp_wsp_LOO <- extractItemLOO(RF_PABD_onlyp_wsp_LOO, bd_type = bd_type, p_type=p_type, otu_type = otu_type,extract = "test_train_comparison" )
  save(compare_PABD_onlyp_wsp_LOO,file="./5_random_forest_LOO/compare_PABD_onlyp_wsp_LOO.RData" )
  importance_PABD_onlyp_wsp_LOO <-extractItemLOO(RF_PABD_onlyp_wsp_LOO, bd_type = bd_type, p_type=p_type, otu_type = otu_type,extract = "importance")
  save(importance_PABD_onlyp_wsp_LOO,file="./5_random_forest_LOO/importance_PABD_onlyp_wsp_LOO.RData" )
  MSE_PABD_onlyp_wsp_LOO <- extractTrainErrorLOO(RF_PABD_onlyp_wsp_LOO, bd_type = bd_type, p_type=p_type, otu_type = otu_type)
  save(MSE_PABD_onlyp_wsp_LOO,file="./5_random_forest_LOO/MSE_PABD_onlyp_wsp_LOO.RData" )
}
remove(RF_PABD_onlyp_wsp_LOO)

#### PABD, ponly, no species ####
if ( RERUN_RF ) {
  ## Set params
  p_type="onlyp"
  otu_type="NA"
  bd_type="PABD"
  mf_temp=NA
  keepSp=FALSE
  # create progress bar
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  i <- 1
  RF_PABD_onlyp_nosp_LOO <- list()
  # Get a formatted mapping file
  get_mf_pred_temp_bd <- get_mf_pred_temp(p_type=p_type
                                          , otu_type = otu_type
                                          , bd_type=bd_type
                                          , p=p
                                          , keepSp=keepSp)
  mf_pred_temp_bd <- get_mf_pred_temp_bd[[1]]
  indiv_list <- get_mf_pred_temp_bd[[2]]
  # Mutate infection depending on question
  Error_metric <- ifelse(bd_type=="PABD","MeanDecreaseAccuracy","X.IncMSE" )
  for ( LOO in c(1:length(allIndiv)) ) {
    RF_PABD_onlyp_nosp_LOO[[paste(otu_type, p_type, bd_type, prop, LOO, sep="_")]] <- list()
    # Pick training and validation set
    trainSet <- mf_pred_temp_bd[-LOO,]
    testSet <- mf_pred_temp_bd[LOO,]
    testSet_indiv <- as.character(unlist(indiv_list[LOO,]))
    
    # Run random forest model
    RF <- randomForest(response ~ ., data=trainSet,importance=TRUE)
    # Run prediction
    testSet_predictions <- predict(RF, testSet, type="class")
    # importance of factors
    importance_RF <- data.frame(importance(RF)) %>%
      mutate(taxa=rownames(importance(RF))) %>%
      select(taxa, everything()) %>%
      arrange(-get(paste0(Error_metric)))
    
    # Add inhibitory to importance
    importance_RF[,c("taxonomy","inhibitory")] <- Taxa[match(importance_RF$taxa, Taxa$Sequence), c("Taxa","inhibitory")]
    importance_RF <- importance_RF %>%
      mutate(taxonomy=ifelse(is.na(taxonomy), taxa, taxonomy))
    
    # Save data
    RF_PABD_onlyp_nosp_LOO[[paste(otu_type, p_type, bd_type, prop,LOO, sep="_")]][["RF"]] <- RF
    RF_PABD_onlyp_nosp_LOO[[paste(otu_type, p_type, bd_type, prop,LOO, sep="_")]][["importance"]] <- importance_RF
    RF_PABD_onlyp_nosp_LOO[[paste(otu_type, p_type, bd_type, prop,LOO, sep="_")]][["test_train_comparison"]] <- data.frame(indivID=testSet_indiv, Test_pred=testSet_predictions, Test_obs=testSet$response)
    
    # Insert into data summary plot
    #------- PROGRESS BAR ------- #
    Sys.sleep(0.1)
    i <- i+1
    # update progress bar
    setTxtProgressBar(pb, i)
    #-----------------------------#
    # }
    
  }
  
  close(pb)
  # load("./5_random_forest_LOO/RF_PABD_onlyp_nosp_LOO.RData" )
  save(RF_PABD_onlyp_nosp_LOO,file="./5_random_forest_LOO/RF_PABD_onlyp_nosp_LOO.RData" )
  compare_PABD_onlyp_nosp_LOO <- extractItemLOO(RF_PABD_onlyp_nosp_LOO, bd_type = bd_type, p_type=p_type, otu_type = otu_type,extract = "test_train_comparison" )
  save(compare_PABD_onlyp_nosp_LOO,file="./5_random_forest_LOO/compare_PABD_onlyp_nosp_LOO.RData" )
  importance_PABD_onlyp_nosp_LOO <-extractItemLOO(RF_PABD_onlyp_nosp_LOO, bd_type = bd_type, p_type=p_type, otu_type = otu_type,extract = "importance")
  save(importance_PABD_onlyp_nosp_LOO,file="./5_random_forest_LOO/importance_PABD_onlyp_nosp_LOO.RData" )
  MSE_PABD_onlyp_nosp_LOO <- extractTrainErrorLOO(RF_PABD_onlyp_nosp_LOO, bd_type = bd_type, p_type=p_type, otu_type = otu_type)
  save(MSE_PABD_onlyp_nosp_LOO,file="./5_random_forest_LOO/MSE_PABD_onlyp_nosp_LOO.RData" )
}
remove(RF_PABD_onlyp_nosp_LOO)


#### PABD, count, with species ####
if ( RERUN_RF ) {
  ## Set params
  p_type="nop"
  otu_type="count"
  bd_type="PABD"
  keepSp=TRUE
  # create progress bar
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  i <- 1
  RF_PABD_count_wsp_LOO <- list()
  # Get a formatted mapping file
  get_mf_pred_temp_bd <- get_mf_pred_temp(p_type=p_type
                                          , otu_type = otu_type
                                          , bd_type=bd_type
                                          , p=p
                                          , keepSp=keepSp)
  mf_pred_temp_bd <- get_mf_pred_temp_bd[[1]]
  indiv_list <- get_mf_pred_temp_bd[[2]]
  # Mutate infection depending on question
  Error_metric <- ifelse(bd_type=="PABD","MeanDecreaseAccuracy","X.IncMSE" )
  for ( LOO in c(1:length(allIndiv)) ) {
    RF_PABD_count_wsp_LOO[[paste(otu_type, p_type, bd_type, prop, LOO, sep="_")]] <- list()
    # Pick training and validation set
    trainSet <- mf_pred_temp_bd[-LOO,]
    testSet <- mf_pred_temp_bd[LOO,]
    testSet_indiv <- as.character(unlist(indiv_list[LOO,]))
    
    # Run random forest model
    RF <- randomForest(response ~ ., data=trainSet,importance=TRUE)
    # Run prediction
    testSet_predictions <- predict(RF, testSet, type="class")
    # importance of factors
    importance_RF <- data.frame(importance(RF)) %>%
      mutate(taxa=rownames(importance(RF))) %>%
      select(taxa, everything()) %>%
      arrange(-get(paste0(Error_metric)))
    
    # Add inhibitory to importance
    importance_RF[,c("taxonomy","inhibitory")] <- Taxa[match(importance_RF$taxa, Taxa$Sequence), c("Taxa","inhibitory")]
    importance_RF <- importance_RF %>%
      mutate(taxonomy=ifelse(is.na(taxonomy), taxa, taxonomy))
    
    # Save data
    RF_PABD_count_wsp_LOO[[paste(otu_type, p_type, bd_type, prop,LOO, sep="_")]][["RF"]] <- RF
    RF_PABD_count_wsp_LOO[[paste(otu_type, p_type, bd_type, prop,LOO, sep="_")]][["importance"]] <- importance_RF
    RF_PABD_count_wsp_LOO[[paste(otu_type, p_type, bd_type, prop,LOO, sep="_")]][["test_train_comparison"]] <- data.frame(indivID=testSet_indiv, Test_pred=testSet_predictions, Test_obs=testSet$response)
    
    # Insert into data summary plot
    #------- PROGRESS BAR ------- #
    Sys.sleep(0.1)
    i <- i+1
    # update progress bar
    setTxtProgressBar(pb, i)
    #-----------------------------#
    # }
    
  }
  
  close(pb)

  # load("./5_random_forest_LOO/RF_PABD_count_wsp_LOO.RData" )
  save(RF_PABD_count_wsp_LOO,file="./5_random_forest_LOO/RF_PABD_count_wsp_LOO.RData" )
  compare_PABD_count_wsp_LOO <- extractItemLOO(RF_PABD_count_wsp_LOO, bd_type = bd_type, p_type=p_type, otu_type = otu_type,extract = "test_train_comparison" )
  save(compare_PABD_count_wsp_LOO,file="./5_random_forest_LOO/compare_PABD_count_wsp_LOO.RData" )
  importance_PABD_count_wsp_LOO <-extractItemLOO(RF_PABD_count_wsp_LOO, bd_type = bd_type, p_type=p_type, otu_type = otu_type,extract = "importance")
  save(importance_PABD_count_wsp_LOO,file="./5_random_forest_LOO/importance_PABD_count_wsp_LOO.RData" )
  MSE_PABD_count_wsp_LOO <- extractTrainErrorLOO(RF_PABD_count_wsp_LOO, bd_type=bd_type, otu_type=otu_type, p_type=p_type)
  save(MSE_PABD_count_wsp_LOO,file="./5_random_forest_LOO/MSE_PABD_count_wsp_LOO.RData" )
}
remove(RF_PABD_count_wsp_LOO)


#### PABD, count, no species ####
if ( RERUN_RF ) {
  ## Set params
  p_type="nop"
  otu_type="count"
  bd_type="PABD"
  keepSp=FALSE
  # create progress bar
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  i <- 1
  RF_PABD_count_nosp_LOO <- list()
  # Get a formatted mapping file
  get_mf_pred_temp_bd <- get_mf_pred_temp(p_type=p_type
                                          , otu_type = otu_type
                                          , bd_type=bd_type
                                          , p=p
                                          , keepSp=keepSp)
  mf_pred_temp_bd <- get_mf_pred_temp_bd[[1]]
  indiv_list <- get_mf_pred_temp_bd[[2]]
  # Mutate infection depending on question
  Error_metric <- ifelse(bd_type=="PABD","MeanDecreaseAccuracy","X.IncMSE" )
  for ( LOO in c(1:length(allIndiv)) ) {
    RF_PABD_count_nosp_LOO[[paste(otu_type, p_type, bd_type, prop, LOO, sep="_")]] <- list()
    # Pick training and validation set
    trainSet <- mf_pred_temp_bd[-LOO,]
    testSet <- mf_pred_temp_bd[LOO,]
    testSet_indiv <- as.character(unlist(indiv_list[LOO,]))
    
    # Run random forest model
    RF <- randomForest(response ~ ., data=trainSet,importance=TRUE)
    # Run prediction
    testSet_predictions <- predict(RF, testSet, type="class")
    # importance of factors
    importance_RF <- data.frame(importance(RF)) %>%
      mutate(taxa=rownames(importance(RF))) %>%
      select(taxa, everything()) %>%
      arrange(-get(paste0(Error_metric)))
    
    # Add inhibitory to importance
    importance_RF[,c("taxonomy","inhibitory")] <- Taxa[match(importance_RF$taxa, Taxa$Sequence), c("Taxa","inhibitory")]
    importance_RF <- importance_RF %>%
      mutate(taxonomy=ifelse(is.na(taxonomy), taxa, taxonomy))
    
    # Save data
    RF_PABD_count_nosp_LOO[[paste(otu_type, p_type, bd_type, prop,LOO, sep="_")]][["RF"]] <- RF
    RF_PABD_count_nosp_LOO[[paste(otu_type, p_type, bd_type, prop,LOO, sep="_")]][["importance"]] <- importance_RF
    RF_PABD_count_nosp_LOO[[paste(otu_type, p_type, bd_type, prop,LOO, sep="_")]][["test_train_comparison"]] <- data.frame(indivID=testSet_indiv, Test_pred=testSet_predictions, Test_obs=testSet$response)
    
    # Insert into data summary plot
    #------- PROGRESS BAR ------- #
    Sys.sleep(0.1)
    i <- i+1
    # update progress bar
    setTxtProgressBar(pb, i)
    #-----------------------------#
    # }
    
  }
  
  close(pb)
  # load("./5_random_forest_LOO/RF_PABD_count_nosp_LOO.RData")
  save(RF_PABD_count_nosp_LOO,file="./5_random_forest_LOO/RF_PABD_count_nosp_LOO.RData" )
  compare_PABD_count_nosp_LOO <- extractItemLOO(RF_PABD_count_nosp_LOO, bd_type = bd_type, p_type=p_type, otu_type = otu_type,extract = "test_train_comparison" )
  save(compare_PABD_count_nosp_LOO,file="./5_random_forest_LOO/compare_PABD_count_nosp_LOO.RData" )
  importance_PABD_count_nosp_LOO <-extractItemLOO(RF_PABD_count_nosp_LOO, bd_type = bd_type, p_type=p_type, otu_type = otu_type,extract = "importance")
  save(importance_PABD_count_nosp_LOO,file="./5_random_forest_LOO/importance_PABD_count_nosp_LOO.RData" )
  MSE_PABD_count_nosp_LOO <- extractTrainErrorLOO(RF_PABD_count_nosp_LOO, bd_type = bd_type, p_type=p_type, otu_type = otu_type)
  save(MSE_PABD_count_nosp_LOO,file="./5_random_forest_LOO/MSE_PABD_count_nosp_LOO.RData" )
}
remove(RF_PABD_count_nosp_LOO)


#### PABD, PA, with species ####
if ( RERUN_RF ) {
  ## Set params
  p_type="nop"
  otu_type="PA"
  bd_type="PABD"
  keepSp=TRUE
  # create progress bar
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  i <- 1
  RF_PABD_PA_wsp_LOO <- list()
  # Get a formatted mapping file
  get_mf_pred_temp_bd <- get_mf_pred_temp(p_type=p_type
                                          , otu_type = otu_type
                                          , bd_type=bd_type
                                          , p=p
                                          , keepSp=keepSp)
  mf_pred_temp_bd <- get_mf_pred_temp_bd[[1]]
  indiv_list <- get_mf_pred_temp_bd[[2]]
  # Mutate infection depending on question
  Error_metric <- ifelse(bd_type=="PABD","MeanDecreaseAccuracy","X.IncMSE" )
  for ( LOO in c(1:length(allIndiv)) ) {
    RF_PABD_PA_wsp_LOO[[paste(otu_type, p_type, bd_type, prop, LOO, sep="_")]] <- list()
    # Pick training and validation set
    trainSet <- mf_pred_temp_bd[-LOO,]
    testSet <- mf_pred_temp_bd[LOO,]
    testSet_indiv <- as.character(unlist(indiv_list[LOO,]))
    
    # Run random forest model
    RF <- randomForest(response ~ ., data=trainSet,importance=TRUE)
    # Run prediction
    testSet_predictions <- predict(RF, testSet, type="class")
    # importance of factors
    importance_RF <- data.frame(importance(RF)) %>%
      mutate(taxa=rownames(importance(RF))) %>%
      select(taxa, everything()) %>%
      arrange(-get(paste0(Error_metric)))
    
    # Add inhibitory to importance
    importance_RF[,c("taxonomy","inhibitory")] <- Taxa[match(importance_RF$taxa, Taxa$Sequence), c("Taxa","inhibitory")]
    importance_RF <- importance_RF %>%
      mutate(taxonomy=ifelse(is.na(taxonomy), taxa, taxonomy))
    
    # Save data
    RF_PABD_PA_wsp_LOO[[paste(otu_type, p_type, bd_type, prop,LOO, sep="_")]][["RF"]] <- RF
    RF_PABD_PA_wsp_LOO[[paste(otu_type, p_type, bd_type, prop,LOO, sep="_")]][["importance"]] <- importance_RF
    RF_PABD_PA_wsp_LOO[[paste(otu_type, p_type, bd_type, prop,LOO, sep="_")]][["test_train_comparison"]] <- data.frame(indivID=testSet_indiv, Test_pred=testSet_predictions, Test_obs=testSet$response)
    
    # Insert into data summary plot
    #------- PROGRESS BAR ------- #
    Sys.sleep(0.1)
    i <- i+1
    # update progress bar
    setTxtProgressBar(pb, i)
    #-----------------------------#
    # }
    
  }
  
  close(pb)
  # load("./5_random_forest_LOO/RF_PABD_PA_wsp_LOO.RData")
  save(RF_PABD_PA_wsp_LOO,file="./5_random_forest_LOO/RF_PABD_PA_wsp_LOO.RData" )
  compare_PABD_PA_wsp_LOO <- extractItemLOO(RF_PABD_PA_wsp_LOO, bd_type = bd_type, p_type=p_type, otu_type = otu_type,extract = "test_train_comparison" )
  save(compare_PABD_PA_wsp_LOO,file="./5_random_forest_LOO/compare_PABD_PA_wsp_LOO.RData" )
  importance_PABD_PA_wsp_LOO <-extractItemLOO(RF_PABD_PA_wsp_LOO, bd_type = bd_type, p_type=p_type, otu_type = otu_type,extract = "importance")
  save(importance_PABD_PA_wsp_LOO,file="./5_random_forest_LOO/importance_PABD_PA_wsp_LOO.RData" )
  MSE_PABD_PA_wsp_LOO <- extractTrainErrorLOO(RF_PABD_PA_wsp_LOO, bd_type = bd_type, p_type=p_type, otu_type = otu_type, p_type="onlyp")
  save(MSE_PABD_PA_wsp_LOO,file="./5_random_forest_LOO/MSE_PABD_PA_wsp_LOO.RData" )
}
remove(RF_PABD_PA_wsp_LOO)

#### PABD, PA, no species ####
if ( RERUN_RF ) {
  ## Set params
  p_type="nop"
  otu_type="PA"
  bd_type="PABD"
  keepSp=FALSE
  # create progress bar
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  i <- 1
  RF_PABD_PA_nosp_LOO <- list()
  # Get a formatted mapping file
  get_mf_pred_temp_bd <- get_mf_pred_temp(p_type=p_type
                                          , otu_type = otu_type
                                          , bd_type=bd_type
                                          , p=p
                                          , keepSp=keepSp)
  mf_pred_temp_bd <- get_mf_pred_temp_bd[[1]]
  indiv_list <- get_mf_pred_temp_bd[[2]]
  # Mutate infection depending on question
  Error_metric <- ifelse(bd_type=="PABD","MeanDecreaseAccuracy","X.IncMSE" )
  for ( LOO in c(1:length(allIndiv)) ) {
    RF_PABD_PA_nosp_LOO[[paste(otu_type, p_type, bd_type, prop, LOO, sep="_")]] <- list()
    # Pick training and validation set
    trainSet <- mf_pred_temp_bd[-LOO,]
    testSet <- mf_pred_temp_bd[LOO,]
    testSet_indiv <- as.character(unlist(indiv_list[LOO,]))
    
    # Run random forest model
    RF <- randomForest(response ~ ., data=trainSet,importance=TRUE)
    # Run prediction
    testSet_predictions <- predict(RF, testSet, type="class")
    # importance of factors
    importance_RF <- data.frame(importance(RF)) %>%
      mutate(taxa=rownames(importance(RF))) %>%
      select(taxa, everything()) %>%
      arrange(-get(paste0(Error_metric)))
    
    # Add inhibitory to importance
    importance_RF[,c("taxonomy","inhibitory")] <- Taxa[match(importance_RF$taxa, Taxa$Sequence), c("Taxa","inhibitory")]
    importance_RF <- importance_RF %>%
      mutate(taxonomy=ifelse(is.na(taxonomy), taxa, taxonomy))
    
    # Save data
    RF_PABD_PA_nosp_LOO[[paste(otu_type, p_type, bd_type, prop,LOO, sep="_")]][["RF"]] <- RF
    RF_PABD_PA_nosp_LOO[[paste(otu_type, p_type, bd_type, prop,LOO, sep="_")]][["importance"]] <- importance_RF
    RF_PABD_PA_nosp_LOO[[paste(otu_type, p_type, bd_type, prop,LOO, sep="_")]][["test_train_comparison"]] <- data.frame(indivID=testSet_indiv, Test_pred=testSet_predictions, Test_obs=testSet$response)
    
    # Insert into data summary plot
    #------- PROGRESS BAR ------- #
    Sys.sleep(0.1)
    i <- i+1
    # update progress bar
    setTxtProgressBar(pb, i)
    #-----------------------------#
    # }
    
  }
  
  close(pb)
  # load("./5_random_forest_LOO/RF_PABD_PA_nosp_LOO.RData")
  save(RF_PABD_PA_nosp_LOO,file="./5_random_forest_LOO/RF_PABD_PA_nosp_LOO.RData" )
  compare_PABD_PA_nosp_LOO <- extractItemLOO(RF_PABD_PA_nosp_LOO, bd_type = bd_type, p_type=p_type, otu_type = otu_type,extract = "test_train_comparison" )
  save(compare_PABD_PA_nosp_LOO,file="./5_random_forest_LOO/compare_PABD_PA_nosp_LOO.RData" )
  importance_PABD_PA_nosp_LOO <-extractItemLOO(RF_PABD_PA_nosp_LOO, bd_type = bd_type, p_type=p_type, otu_type = otu_type,extract = "importance")
  save(importance_PABD_PA_nosp_LOO,file="./5_random_forest_LOO/importance_PABD_PA_nosp_LOO.RData" )
  MSE_PABD_PA_nosp_LOO <- extractTrainErrorLOO(RF_PABD_PA_nosp_LOO, bd_type = bd_type, p_type=p_type, otu_type = otu_type)
  save(MSE_PABD_PA_nosp_LOO,file="./5_random_forest_LOO/MSE_PABD_PA_nosp_LOO.RData" )
}
remove(RF_PABD_PA_nosp_LOO)

#### PABD, all count, with species ####
if ( RERUN_RF ) {
  ## Set params
  p_type="withp"
  otu_type="count"
  bd_type="PABD"
  keepSp=TRUE
  # create progress bar
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  i <- 1
  RF_PABD_withpcount_wsp_LOO <- list()
  # Get a formatted mapping file
  get_mf_pred_temp_bd <- get_mf_pred_temp(p_type=p_type
                                          , otu_type = otu_type
                                          , bd_type=bd_type
                                          , p=p
                                          , keepSp=keepSp)
  mf_pred_temp_bd <- get_mf_pred_temp_bd[[1]]
  indiv_list <- get_mf_pred_temp_bd[[2]]
  # Mutate infection depending on question
  Error_metric <- ifelse(bd_type=="PABD","MeanDecreaseAccuracy","X.IncMSE" )
  for ( LOO in c(1:length(allIndiv)) ) {
    RF_PABD_withpcount_wsp_LOO[[paste(otu_type, p_type, bd_type, prop, LOO, sep="_")]] <- list()
    # Pick training and validation set
    trainSet <- mf_pred_temp_bd[-LOO,]
    testSet <- mf_pred_temp_bd[LOO,]
    testSet_indiv <- as.character(unlist(indiv_list[LOO,]))
    
    # Run random forest model
    RF <- randomForest(response ~ ., data=trainSet,importance=TRUE)
    # Run prediction
    testSet_predictions <- predict(RF, testSet, type="class")
    # importance of factors
    importance_RF <- data.frame(importance(RF)) %>%
      mutate(taxa=rownames(importance(RF))) %>%
      select(taxa, everything()) %>%
      arrange(-get(paste0(Error_metric)))
    
    # Add inhibitory to importance
    importance_RF[,c("taxonomy","inhibitory")] <- Taxa[match(importance_RF$taxa, Taxa$Sequence), c("Taxa","inhibitory")]
    importance_RF <- importance_RF %>%
      mutate(taxonomy=ifelse(is.na(taxonomy), taxa, taxonomy))
    
    # Save data
    RF_PABD_withpcount_wsp_LOO[[paste(otu_type, p_type, bd_type, prop,LOO, sep="_")]][["RF"]] <- RF
    RF_PABD_withpcount_wsp_LOO[[paste(otu_type, p_type, bd_type, prop,LOO, sep="_")]][["importance"]] <- importance_RF
    RF_PABD_withpcount_wsp_LOO[[paste(otu_type, p_type, bd_type, prop,LOO, sep="_")]][["test_train_comparison"]] <- data.frame(indivID=testSet_indiv, Test_pred=testSet_predictions, Test_obs=testSet$response)
    
    # Insert into data summary plot
    #------- PROGRESS BAR ------- #
    Sys.sleep(0.1)
    i <- i+1
    # update progress bar
    setTxtProgressBar(pb, i)
    #-----------------------------#
    # }
    
  }
  
  close(pb)
  # load("./5_random_forest_LOO/RF_PABD_withpcount_wsp_LOO.RData")
  save(RF_PABD_withpcount_wsp_LOO,file="./5_random_forest_LOO/RF_PABD_withpcount_wsp_LOO.RData" )
  compare_PABD_withpcount_wsp_LOO <- extractItemLOO(RF_PABD_withpcount_wsp_LOO, bd_type = bd_type, p_type=p_type, otu_type = otu_type,extract = "test_train_comparison" )
  save(compare_PABD_withpcount_wsp_LOO,file="./5_random_forest_LOO/compare_PABD_withpcount_wsp_LOO.RData" )
  importance_PABD_withpcount_wsp_LOO <-extractItemLOO(RF_PABD_withpcount_wsp_LOO, bd_type = bd_type, p_type=p_type, otu_type = otu_type,extract = "importance")
  save(importance_PABD_withpcount_wsp_LOO,file="./5_random_forest_LOO/importance_PABD_withpcount_wsp_LOO.RData" )
  MSE_PABD_withpcount_wsp_LOO <- extractTrainErrorLOO(RF_PABD_withpcount_wsp_LOO, bd_type = bd_type, p_type=p_type, otu_type = otu_type)
  save(MSE_PABD_withpcount_wsp_LOO,file="./5_random_forest_LOO/MSE_PABD_withpcount_wsp_LOO.RData" )
}
remove(RF_PABD_withpcount_wsp_LOO)

#### PABD, all PA, with species ####
if ( RERUN_RF ) {
  ## Set params
  p_type="withp"
  otu_type="PA"
  bd_type="PABD"
  keepSp=TRUE
  # create progress bar
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  i <- 1
  RF_PABD_withpPA_wsp_LOO <- list()
  # Get a formatted mapping file
  get_mf_pred_temp_bd <- get_mf_pred_temp(p_type=p_type
                                          , otu_type = otu_type
                                          , bd_type=bd_type
                                          , p=p
                                          , keepSp=keepSp)
  mf_pred_temp_bd <- get_mf_pred_temp_bd[[1]]
  indiv_list <- get_mf_pred_temp_bd[[2]]
  # Mutate infection depending on question
  Error_metric <- ifelse(bd_type=="PABD","MeanDecreaseAccuracy","X.IncMSE" )
  for ( LOO in c(1:length(allIndiv)) ) {
    RF_PABD_withpPA_wsp_LOO[[paste(otu_type, p_type, bd_type, prop, LOO, sep="_")]] <- list()
    # Pick training and validation set
    trainSet <- mf_pred_temp_bd[-LOO,]
    testSet <- mf_pred_temp_bd[LOO,]
    testSet_indiv <- as.character(unlist(indiv_list[LOO,]))
    
    # Run random forest model
    RF <- randomForest(response ~ ., data=trainSet,importance=TRUE)
    # Run prediction
    testSet_predictions <- predict(RF, testSet, type="class")
    # importance of factors
    importance_RF <- data.frame(importance(RF)) %>%
      mutate(taxa=rownames(importance(RF))) %>%
      select(taxa, everything()) %>%
      arrange(-get(paste0(Error_metric)))
    
    # Add inhibitory to importance
    importance_RF[,c("taxonomy","inhibitory")] <- Taxa[match(importance_RF$taxa, Taxa$Sequence), c("Taxa","inhibitory")]
    importance_RF <- importance_RF %>%
      mutate(taxonomy=ifelse(is.na(taxonomy), taxa, taxonomy))
    
    # Save data
    RF_PABD_withpPA_wsp_LOO[[paste(otu_type, p_type, bd_type, prop,LOO, sep="_")]][["RF"]] <- RF
    RF_PABD_withpPA_wsp_LOO[[paste(otu_type, p_type, bd_type, prop,LOO, sep="_")]][["importance"]] <- importance_RF
    RF_PABD_withpPA_wsp_LOO[[paste(otu_type, p_type, bd_type, prop,LOO, sep="_")]][["test_train_comparison"]] <- data.frame(indivID=testSet_indiv, Test_pred=testSet_predictions, Test_obs=testSet$response)
    
    # Insert into data summary plot
    #------- PROGRESS BAR ------- #
    Sys.sleep(0.1)
    i <- i+1
    # update progress bar
    setTxtProgressBar(pb, i)
    #-----------------------------#
    # }
    
  }
  
  close(pb)
  # load("./5_random_forest_LOO/RF_PABD_withpPA_wsp_LOO.RData")
  save(RF_PABD_withpPA_wsp_LOO,file="./5_random_forest_LOO/RF_PABD_withpPA_wsp_LOO.RData" )
  compare_PABD_withpPA_wsp_LOO <- extractItemLOO(RF_PABD_withpPA_wsp_LOO, bd_type = bd_type, p_type=p_type, otu_type = otu_type,extract = "test_train_comparison" )
  save(compare_PABD_withpPA_wsp_LOO,file="./5_random_forest_LOO/compare_PABD_withpPA_wsp_LOO.RData" )
  importance_PABD_withpPA_wsp_LOO <-extractItemLOO(RF_PABD_withpPA_wsp_LOO, bd_type = bd_type, p_type=p_type, otu_type = otu_type,extract = "importance")
  save(importance_PABD_withpPA_wsp_LOO,file="./5_random_forest_LOO/importance_PABD_withpPA_wsp_LOO.RData" )
  MSE_PABD_withpPA_wsp_LOO <- extractTrainErrorLOO(RF_PABD_withpPA_wsp_LOO, bd_type = bd_type, p_type=p_type, otu_type = otu_type)
  save(MSE_PABD_withpPA_wsp_LOO,file="./5_random_forest_LOO/MSE_PABD_withpPA_wsp_LOO.RData" )
}
remove(RF_PABD_withpPA_wsp_LOO)

#### Re-load data ####
if ( !RERUN_RF) {
  
  # Updated random forest LOO results
  load("./5_random_forest_LOO/compare_infect_onlyp_wsp_LOO.RData")
  load("./5_random_forest_LOO/compare_infect_onlyp_nosp_LOO.RData")
  load("./5_random_forest_LOO/compare_infect_count_wsp_LOO.RData")
  load("./5_random_forest_LOO/compare_infect_count_nosp_LOO.RData")
  load("./5_random_forest_LOO/compare_infect_PA_wsp_LOO.RData")
  load("./5_random_forest_LOO/compare_infect_PA_nosp_LOO.RData")
  load("./5_random_forest_LOO/compare_infect_withpPA_wsp_LOO.RData")
  load("./5_random_forest_LOO/compare_infect_withpcount_wsp_LOO.RData")
  
  load("./5_random_forest_LOO/compare_PABD_onlyp_wsp_LOO.RData")
  load("./5_random_forest_LOO/compare_PABD_onlyp_nosp_LOO.RData")
  load("./5_random_forest_LOO/compare_PABD_count_wsp_LOO.RData")
  load("./5_random_forest_LOO/compare_PABD_count_nosp_LOO.RData")
  load("./5_random_forest_LOO/compare_PABD_PA_wsp_LOO.RData")
  load("./5_random_forest_LOO/compare_PABD_PA_nosp_LOO.RData")
  load("./5_random_forest_LOO/compare_PABD_withpPA_wsp_LOO.RData")
  load("./5_random_forest_LOO/compare_PABD_withpcount_wsp_LOO.RData")
  
  load("./5_random_forest_LOO/importance_infect_onlyp_wsp_LOO.RData")
  load("./5_random_forest_LOO/importance_infect_onlyp_nosp_LOO.RData")
  load("./5_random_forest_LOO/importance_infect_count_wsp_LOO.RData")
  load("./5_random_forest_LOO/importance_infect_count_nosp_LOO.RData")
  load("./5_random_forest_LOO/importance_infect_PA_wsp_LOO.RData")
  load("./5_random_forest_LOO/importance_infect_PA_nosp_LOO.RData")
  load("./5_random_forest_LOO/importance_infect_withpPA_wsp_LOO.RData")
  load("./5_random_forest_LOO/importance_infect_withpcount_wsp_LOO.RData")
  
  load("./5_random_forest_LOO/importance_PABD_onlyp_wsp_LOO.RData")
  load("./5_random_forest_LOO/importance_PABD_onlyp_nosp_LOO.RData")
  load("./5_random_forest_LOO/importance_PABD_count_wsp_LOO.RData")
  load("./5_random_forest_LOO/importance_PABD_count_nosp_LOO.RData")
  load("./5_random_forest_LOO/importance_PABD_PA_wsp_LOO.RData")
  load("./5_random_forest_LOO/importance_PABD_PA_nosp_LOO.RData")
  load("./5_random_forest_LOO/importance_PABD_withpPA_wsp_LOO.RData")
  load("./5_random_forest_LOO/importance_PABD_withpcount_wsp_LOO.RData")
  
  load("./5_random_forest_LOO/MSE_infect_onlyp_wsp_LOO.RData")
  load("./5_random_forest_LOO/MSE_infect_onlyp_nosp_LOO.RData")
  load("./5_random_forest_LOO/MSE_infect_count_wsp_LOO.RData")
  load("./5_random_forest_LOO/MSE_infect_count_nosp_LOO.RData")
  load("./5_random_forest_LOO/MSE_infect_PA_wsp_LOO.RData")
  load("./5_random_forest_LOO/MSE_infect_PA_nosp_LOO.RData")
  load("./5_random_forest_LOO/MSE_infect_withpPA_wsp_LOO.RData")
  load("./5_random_forest_LOO/MSE_infect_withpcount_wsp_LOO.RData")
  
  
}

#### Preliminary plots and data synthesis ####

rbind(cbind(compare_infect_onlyp_nosp_LOO, species=FALSE)
      , cbind(compare_infect_onlyp_wsp_LOO, species=TRUE)) %>%
  ggplot() + geom_point(aes(x=Test_obs, y=Test_pred, col=species)) + geom_abline(aes(intercept=0, slope=1))+
  ylim(0,8)+ xlim(0,8)

rbind(cbind(compare_infect_count_nosp_LOO, species=FALSE)
      , cbind(compare_infect_count_wsp_LOO, species=TRUE)) %>%
  ggplot() + geom_point(aes(x=Test_obs, y=Test_pred, col=species)) + geom_abline(aes(intercept=0, slope=1))+
  ylim(0,8)+ xlim(0,8)

rbind(cbind(compare_infect_PA_nosp_LOO, species=FALSE)
      , cbind(compare_infect_PA_wsp_LOO, species=TRUE)) %>%
  ggplot() + geom_point(aes(x=Test_obs, y=Test_pred, col=species)) + geom_abline(aes(intercept=0, slope=1))+
  ylim(0,8)+ xlim(0,8)



## Accuracy
errorRate_test_all <- rbind(cbind(compare_PABD_onlyp_nosp_LOO, species=FALSE), cbind(compare_PABD_onlyp_wsp_LOO, species=TRUE)
                           , cbind(compare_PABD_count_nosp_LOO, species=FALSE), cbind(compare_PABD_count_wsp_LOO, species=TRUE)
                           ,cbind(compare_PABD_PA_nosp_LOO, species=FALSE), cbind(compare_PABD_PA_wsp_LOO, species=TRUE)) %>%
  mutate(Training_data = ifelse(otu_type=="NA", "Community traits",ifelse(otu_type=="count","ASV counts","ASV prevalence")))

errorRate_test_all_summary <- errorRate_test_all %>%
  group_by(species, Training_data) %>%
  summarize(correct=mean(Test_pred==Test_obs))
errorRate_test_all_summary

observed_only <- errorRate_test_all %>%
  select(indivID, Test_obs, species) %>%
  group_by(indivID, Test_obs, species) %>%
  summarize(Infection_status=unique(Test_obs)) %>%
  ungroup() %>% select(-Test_obs) %>% mutate(Training_data="Observed outcome")

ggsave("./5_random_forest_LOO/PABD_predictions.pdf", height=3.5, width=5
       ,errorRate_test_all %>%
         select(indivID, species, Test_pred, Training_data) %>%
         rename(Infection_status=Test_pred) %>%
         rbind(observed_only) %>%
         mutate(Training_data= factor(Training_data, levels=c("Observed outcome", "Community traits","ASV counts", "ASV prevalence")))%>%
         ggplot() + geom_tile(aes(x=Training_data, y=indivID, fill=Infection_status)
                              , col="black", width=0.9) +
         scale_fill_manual(values=c("darkred","lightblue"), name="Infection status") +
         theme_bw() +
         theme(axis.text.x = element_text(angle=90)) +
         scale_x_discrete(breaks=c("Observed outcome","Community traits","ASV counts","ASV prevalence")
                          , labels=(c("Observed\noutcome","Community\ntraits","ASV\ncounts","ASV\nprevalence")))+
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
  mutate(Training_data = ifelse(otu_type=="NA", "Community traits",ifelse(otu_type=="count","ASV counts","ASV prevalence")))
MSE_test_all
MSE_train_all_long <- rbind(cbind(MSE_infect_onlyp_wsp_LOO, species="With species as predictor")
                            , cbind(MSE_infect_onlyp_nosp_LOO, species="No species predictor")
                            , cbind(MSE_infect_count_wsp_LOO, species="With species as predictor")
                            , cbind(MSE_infect_count_nosp_LOO, species="No species predictor")
                            , cbind(MSE_infect_PA_wsp_LOO, species="With species as predictor")
                            , cbind(MSE_infect_PA_nosp_LOO, species="No species predictor")) %>%
  unite(otu_type, species, col=group, remove=FALSE) %>%
  mutate(Training_data = ifelse(otu_type=="NA", "Community traits",ifelse(otu_type=="count","ASV counts","ASV prevalence"))) %>%
  rename(MSE=error)

ggsave("./5_random_forest_LOO/MSE_all.pdf", height=3, width=5
       ,MSE_train_all_long %>%
         ggplot() + geom_violin(aes(x=Training_data, y=MSE)) + geom_jitter(aes(x=Training_data, y=MSE), height=0, width=0.2, alpha=0.2)+
         geom_point(data=MSE_test_all, aes(x=Training_data, y=MSE), col="red") +
         facet_wrap(.~species)+
         theme_bw()+
         theme(axis.text.x = element_text(angle=90)) +
         scale_x_discrete(labels=c("ASV\ncounts","ASV\nprevalence","Community\ntraits"))
       )

infect_predict_all <- rbind(cbind(compare_infect_onlyp_wsp_LOO, species="With species as predictor")
                            , cbind(compare_infect_onlyp_nosp_LOO, species="No species predictor")
                            , cbind(compare_infect_count_wsp_LOO, species="With species as predictor")
                            , cbind(compare_infect_count_nosp_LOO, species="No species predictor")
                            , cbind(compare_infect_PA_wsp_LOO, species="With species as predictor")
                            , cbind(compare_infect_PA_nosp_LOO, species="No species predictor")) %>%
  unite(otu_type, species, col=group, remove=FALSE) %>%
  mutate(Training_data = ifelse(otu_type=="NA", "Community traits",ifelse(otu_type=="count","ASV counts","ASV prevalence"))) 

ggsave("./5_random_forest_LOO/infect_predictions.pdf", height=3, width=6
       ,infect_predict_all %>%
         ggplot(aes(x=Test_obs, y=Test_pred, col=Training_data)) +geom_point() +
         geom_abline(aes(slope=1, intercept=0), col="black", lty=2) +
         geom_smooth(method="lm") +
         facet_wrap(.~species) +
         xlim(0,8)+ylim(0,8)+xlab(expression(paste("Observed ",italic("Bd")," load")))+ ylab(expression(paste("Predicted ",italic("Bd")," load")))+
         theme_bw()
       )

#### Quick test with combined p and no p ####

rbind(cbind(compare_infect_onlyp_wsp_LOO, species=TRUE)
      , cbind(compare_infect_PA_wsp_LOO, species=TRUE)
      ,cbind(compare_infect_withpPA_wsp_LOO, species=TRUE)) %>%
  unite(otu_type, p_type, col=type,remove=FALSE) %>%
  ggplot() + geom_point(aes(x=Test_obs, y=Test_pred, col=type)) + geom_abline(aes(intercept=0, slope=1))+
  geom_smooth(aes(x=Test_obs, y=Test_pred, col=type), method="lm")+
  ylim(0,8)+ xlim(0,8)

rbind(cbind(compare_infect_onlyp_wsp_LOO, species=TRUE)
      , cbind(compare_infect_PA_wsp_LOO, species=TRUE)
      ,cbind(compare_infect_withpPA_wsp_LOO, species=TRUE)) %>%
  group_by(otu_type, p_type) %>%
  summarize(MSE = mean((Test_obs-Test_pred)^2))

rbind(cbind(compare_infect_onlyp_wsp_LOO, species=TRUE)
      , cbind(compare_infect_count_wsp_LOO, species=TRUE)
      ,cbind(compare_infect_withpcount_wsp_LOO, species=TRUE)) %>%
  unite(otu_type, p_type, col=type,remove=FALSE) %>%
  ggplot() + geom_point(aes(x=Test_obs, y=Test_pred, col=type)) + geom_abline(aes(intercept=0, slope=1))+
  geom_smooth(aes(x=Test_obs, y=Test_pred, col=type), method="lm")+
  ylim(0,8)+ xlim(0,8)

rbind(cbind(compare_infect_onlyp_wsp_LOO, species=TRUE)
      , cbind(compare_infect_count_wsp_LOO, species=TRUE)
      ,cbind(compare_infect_withpcount_wsp_LOO, species=TRUE)) %>%
  group_by(otu_type, p_type) %>%
  summarize(MSE = mean((Test_obs-Test_pred)^2))

rbind(cbind(compare_PABD_onlyp_wsp_LOO, species=TRUE)
      , cbind(compare_PABD_count_wsp_LOO, species=TRUE)
      ,cbind(compare_PABD_withpcount_wsp_LOO, species=TRUE)) %>%
  group_by(otu_type, p_type) %>%
  summarize(correctRate = mean(Test_obs==Test_pred))
rbind(cbind(compare_PABD_onlyp_wsp_LOO, species=TRUE)
      , cbind(compare_PABD_PA_wsp_LOO, species=TRUE)
      ,cbind(compare_PABD_withpPA_wsp_LOO, species=TRUE)) %>%
  group_by(otu_type, p_type) %>%
  summarize(correctRate = mean(Test_obs==Test_pred))



#### Looking at importance ####

## For PABD
ggsave("./5_random_forest_LOO/importance_PABD_onlyp.pdf", height=6, width=3
       ,rbind(cbind(importance_PABD_onlyp_wsp_LOO,species="With species as a predictor"),
              cbind(importance_PABD_onlyp_nosp_LOO,species="No species predictor")) %>%
         arrange(-MeanDecreaseAccuracy) %>% mutate(taxonomy=factor(taxonomy, levels=unique(taxonomy))) %>%
         ggplot() + geom_jitter(aes(x=taxonomy, y=MeanDecreaseAccuracy), height=0, width=0.1, col="darkgrey")+
         geom_hline(aes(yintercept=0), col="red", alpha=0.2)+
         theme_bw()+
         theme(axis.text.x=element_text(angle=90, hjust=1))+
         facet_grid(species~.)+ylab("Decrease in accuracy when absent")+xlab("Predictor")
)

importance_PABD_count_combined <- rbind(cbind(importance_PABD_count_wsp_LOO,species="With species as a predictor"),
                                          cbind(importance_PABD_count_nosp_LOO,species="No species predictor"))
toKeep_count <- importance_PABD_count_combined %>%
  group_by(taxonomy, inhibitory) %>%
  summarize(meanMDA=mean(MeanDecreaseAccuracy), sdMDA=sd(MeanDecreaseAccuracy)) %>%
  arrange(-meanMDA) %>% filter(meanMDA>1) %>% pull(taxonomy)
# q95_count <- quantile(importance_PABD_count_combined%>%filter(X.IncMSE>0)%>%pull(X.IncMSE), probs = c(0.95))
ggsave("./5_random_forest_LOO/importance_PABD_count.pdf",height=7, width=5
       ,importance_PABD_count_combined%>%
         filter(taxonomy%in%toKeep_count) %>%
         mutate(Inhibitory=ifelse(inhibitory==1,"Yes","No")) %>%
         arrange(-MeanDecreaseAccuracy) %>% mutate(taxonomy=factor(taxonomy, levels=unique(taxonomy))) %>%
         ggplot() + geom_jitter(aes(x=taxonomy, y=MeanDecreaseAccuracy, col=Inhibitory), height=0, width=0.1)+
         scale_color_manual(values=c("black","purple"), na.value="darkgrey")+
         theme_bw()+
         theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))+
         facet_grid(species~.)+ylab("Decrease in accuracy when absent")+xlab("Predictor")
)

importance_PABD_prevalence_combined <- rbind(cbind(importance_PABD_PA_wsp_LOO,species="With species as a predictor"),
                                               cbind(importance_PABD_PA_nosp_LOO,species="No species predictor"))
toKeep_PA <- importance_PABD_prevalence_combined %>%
  group_by(taxonomy, inhibitory) %>%
  summarize(meanMDA=mean(MeanDecreaseAccuracy), sdMDA=sd(MeanDecreaseAccuracy)) %>%
  arrange(-meanMDA) %>% filter(meanMDA>1) %>% pull(taxonomy)
# q95_PA <- quantile(importance_PABD_prevalence_combined%>%filter(X.IncMSE>0)%>%pull(X.IncMSE), probs = c(0.95))
ggsave("./5_random_forest_LOO/importance_PABD_PA.pdf",height=7, width=5
       ,importance_PABD_prevalence_combined%>%
         filter(taxonomy%in%toKeep_PA) %>%
         mutate(Inhibitory=ifelse(inhibitory==1,"Yes","No")) %>%
         arrange(-MeanDecreaseAccuracy) %>% mutate(taxonomy=factor(taxonomy, levels=unique(taxonomy))) %>%
         ggplot() + geom_jitter(aes(x=taxonomy, y=MeanDecreaseAccuracy, col=Inhibitory), height=0, width=0.1)+
         scale_color_manual(values=c("black","purple"), na.value="darkgrey")+
         theme_bw()+
         theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))+
         facet_grid(species~.)+ylab("Decrease in accuracy when absent")+xlab("Predictor")
)

## For infect
ggsave("./5_random_forest_LOO/importance_infect_onlyp.pdf", height=6, width=3
       ,rbind(cbind(importance_infect_onlyp_wsp_LOO,species="With species as a predictor"),
              cbind(importance_infect_onlyp_nosp_LOO,species="No species predictor")) %>%
         arrange(-X.IncMSE) %>% mutate(taxonomy=factor(taxonomy, levels=unique(taxonomy))) %>%
         ggplot() + geom_jitter(aes(x=taxonomy, y=X.IncMSE), height=0, width=0.1, col="darkgrey")+
         geom_hline(aes(yintercept=0), col="red", alpha=0.2)+
         theme_bw()+
         theme(axis.text.x=element_text(angle=90, hjust=1))+
         facet_grid(species~.)+ylab("Percent increase in MSE when absent")+xlab("Predictor")
       )

importance_infect_count_combined <- rbind(cbind(importance_infect_count_wsp_LOO,species="With species as a predictor"),
cbind(importance_infect_count_nosp_LOO,species="No species predictor"))
toKeep_count <- importance_infect_count_combined %>%
  group_by(taxonomy, inhibitory) %>%
  summarize(meanIncMSE=mean(X.IncMSE), sdIncMSE=sd(X.IncMSE)) %>%
  arrange(-meanIncMSE) %>% filter(meanIncMSE>1) %>% pull(taxonomy)
# q95_count <- quantile(importance_infect_count_combined%>%filter(X.IncMSE>0)%>%pull(X.IncMSE), probs = c(0.95))
ggsave("./5_random_forest_LOO/importance_infect_count.pdf",height=7, width=5
       ,importance_infect_count_combined%>%
         filter(taxonomy%in%toKeep_count) %>%
         mutate(Inhibitory=ifelse(inhibitory==1,"Yes","No")) %>%
         arrange(-X.IncMSE) %>% mutate(taxonomy=factor(taxonomy, levels=unique(taxonomy))) %>%
         ggplot() + geom_jitter(aes(x=taxonomy, y=X.IncMSE, col=Inhibitory), height=0, width=0.1)+
         scale_color_manual(values=c("black","purple"), na.value="darkgrey")+
         theme_bw()+
         theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))+
         facet_grid(species~.)+ylab("Percent increase in MSE when absent")+xlab("Predictor")
)

importance_infect_prevalence_combined <- rbind(cbind(importance_infect_PA_wsp_LOO,species="With species as a predictor"),
                                          cbind(importance_infect_PA_nosp_LOO,species="No species predictor"))
toKeep_PA <- importance_infect_prevalence_combined %>%
  group_by(taxonomy, inhibitory) %>%
  summarize(meanIncMSE=mean(X.IncMSE), sdIncMSE=sd(X.IncMSE)) %>%
  arrange(-meanIncMSE) %>% filter(meanIncMSE>1) %>% pull(taxonomy)
# q95_PA <- quantile(importance_infect_prevalence_combined%>%filter(X.IncMSE>0)%>%pull(X.IncMSE), probs = c(0.95))
ggsave("./5_random_forest_LOO/importance_infect_PA.pdf",height=7, width=5
       ,importance_infect_prevalence_combined%>%
         filter(taxonomy%in%toKeep_PA) %>%
         mutate(Inhibitory=ifelse(inhibitory==1,"Yes","No")) %>%
         arrange(-X.IncMSE) %>% mutate(taxonomy=factor(taxonomy, levels=unique(taxonomy))) %>%
         ggplot() + geom_jitter(aes(x=taxonomy, y=X.IncMSE, col=Inhibitory), height=0, width=0.1)+
         scale_color_manual(values=c("black","purple"), na.value="darkgrey")+
         theme_bw()+
         theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))+
         facet_grid(species~.)+ylab("Percent increase in MSE when absent")+xlab("Predictor")
)


