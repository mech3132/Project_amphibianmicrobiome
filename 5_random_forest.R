#!/bin/bash Rscript

#### Loading data #####
library(randomForest)
library(tidyverse)
load("3_5sp_mapping_otu_downstream/mf_alt_filt_final.RData")
load("3_5sp_mapping_otu_downstream/otu_filt.RData")
load("3_5sp_mapping_otu_downstream/taxonomy.RData")
load("4_Bayesian_models/all_p.RData")
source("msc_randomforest_codebits.R")

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

# Add OTU table to mf 
mf_rank <- mf_otu_combine(otu_rank, mf)
mf_PA <- mf_otu_combine(otu_PA, mf)

#--------- Random forest data setup -------------#

make_summary_RF <- function(all_output, otu_type=c("PA","rank"), bd_type=c("PABD","infect")) {
  # all_output <- all_RF_predictBD
  # otu_type <- "rank"
  # bd_type="infect"
  # Create empty matrix
  # Get error explained by each model into table
  base_df <- data.frame(Model=c("onlyp",rep(NA,12),"nop",rep(NA,14),"withp",rep(NA,14)))
  if (bd_type == "infect" ) {
    summary_randomforest <- cbind(base_df, data.frame(Train_set_MSE=NA, Validation_set_MSE=NA, Top_importance=NA, IncrMSE=NA))
  } else {
    summary_randomforest <- cbind(base_df, data.frame(Train_set_error_rate=NA, Validation_set_error_rate=NA, Top_importance=NA, MeanDecrAccuracy=NA))
    
  }
  
  for ( p_type in c("onlyp","nop","withp") ) {
    curr_cat <- all_output[[paste(otu_type, p_type, bd_type, sep="_")]]
    curr_model <- curr_cat[[paste0("RF_",otu_type,"_",p_type,"_",bd_type)]]
    error_rate_predict <- curr_cat[[paste0("testSet_error",otu_type,"_",p_type,"_",bd_type)]]
    if ( bd_type=="PABD" ) {
      error_rate <- curr_model$err.rate[,"OOB"][nrow(curr_model$err.rate)]*100
      var_explained <- NA
      final_error <- paste0(round(error_rate, 3), "%")
      valid_set_error <- paste0(round(error_rate_predict,3),"%")
    } else if (bd_type=="infect") {
      error_rate <- curr_model$mse[length(curr_model$mse)]
      var_explained <- curr_model$rsq[length(curr_model$rsq)]
      final_error <- paste0(round(error_rate, 3), ", ", round(var_explained*100, 2),"% var explained")
      valid_set_error <- paste0(round(error_rate_predict,3))
    }
    current_row <- ifelse(p_type=="onlyp",0,ifelse(p_type=="nop",13,28))
    
    # Finally, get top importance factors
    curr_imp <- curr_cat[[paste0("importance_RF_",otu_type,"_",p_type,"_",bd_type)]]
    if ( bd_type == "infect") {
      col_imp <- "X.IncMSE"
      col_phrase <- "%IncrMSE="
    } else if ( bd_type == "PABD" ) {
      col_imp <- "MeanDecreaseAccuracy"
      col_phrase <- "MeanDecrAccuracy="
    }
    
    # col_imp changes if it's p model
    # Only 13 factors in p model
    if ( p_type=="onlyp") {
      l <- 13
      # metric <- "Metric"
      # if ( bd=="infect" ) {
      #   col_imp <- "%IncMSE"
      # }
      metric <- "taxonomy"
    } else {
      l <- 15
      metric <- "taxonomy"
      
    }
    
    for ( tax in 1:l) {
      if (is.null(curr_imp$inhibitory[tax])){
        inhib_phrase <- ""
      } else if (is.na(curr_imp$inhibitory[tax]) ) {
        inhib_phrase <- ""
      } else if (curr_imp$inhibitory[tax]==1) {
        inhib_phrase <- " (inhibitory)"
      } else {
        inhib_phrase <- ""
      }
      
      if ( bd_type == "infect" ) {
        summary_randomforest[current_row+1, c("Train_set_MSE","Validation_set_MSE")] <- c(final_error, valid_set_error)
        summary_randomforest[current_row+tax,"Top_importance"] <- paste0(curr_imp[tax,metric], inhib_phrase)
        summary_randomforest[current_row+tax,"IncrMSE"] <- paste0(col_phrase, round(curr_imp[tax,col_imp],2))
      } else if ( bd_type == "PABD") {
        summary_randomforest[current_row+1, c("Train_set_error_rate","Validation_set_error_rate")] <- c(final_error, valid_set_error)
        summary_randomforest[current_row+tax,"Top_importance"] <- paste0(curr_imp[tax,metric], inhib_phrase)
        summary_randomforest[current_row+tax,"MeanDecrAccuracy"] <- paste0(col_phrase, round(curr_imp[tax,col_imp],2))
      }
    }
  }
  
  return(summary_randomforest)
}

# -----------------------------------------------#
all_RF_predictBD <- list()

RERUN_RF <- FALSE
if ( RERUN_RF ) {
  total <- 2*3*2
  # create progress bar
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  i <- 1
  for ( otu_type in c("rank","PA") ) {
    for ( p_type in c("onlyp","withp","nop") ) {
      for ( bd_type in c("infect","PABD")) {
        all_RF_predictBD[[paste(otu_type, p_type, bd_type, sep="_")]] <- list()
        if ( p_type == "onlyp") {
          mf_pred_temp <- p %>%
            separate(indivID, into=c("species","indiv"), remove=FALSE) %>%
            mutate(species=factor(species)) %>%
            arrange(indivID) %>%
            select(-c(indiv, indivID))
        } else {
          p_to_add <- get(ifelse(p_type == "withp", "p", "p_infectonly"))
          mf_temp <- get(ifelse(otu_type == "rank", "mf_rank", "mf_PA"))
          # Filter mf and re-arrange to be compatible with randomforest
          mf_pred_temp <- mf_temp %>%
            filter(prepost=="Pre", Bd_exposure=="Bd-exposed") %>%
            select(-c("SampleID","prepost","Bd_exposure")) %>%
            gather(-c(indivID),key=OTU, value=rank) %>%
            mutate(rank=as.numeric(rank)) %>%
            # filter(indivID==c("Anbo_4","Anbo_2"), OTU=="TACGTAGGGTACGAGCGTTGTCCGGAATTATTGGGCGTAAAGAGCTCGTAGGTGGTTGGTCACGTCTGCTGTGGAAACGCAACGCTTAAC") %>%
            group_by(indivID, OTU) %>%
            summarize(aveRank=mean(rank)) %>%
            ungroup() %>%
            separate(indivID, into=c("species","indiv"), remove=FALSE) %>%
            select(-indiv) %>%
            filter(!is.na(aveRank)) %>%
            spread(key=OTU, value=aveRank)  %>% 
            left_join(p_to_add) %>%
            arrange(indivID) %>%
            select(-indivID) %>%
            mutate(species=factor(species))
        } 
        
        # Mutate infection depending on question
        if (bd_type == "PABD" ) {
          mf_pred_temp_bd  <- mf_pred_temp %>%
            mutate(PABD = factor(ifelse(infect>0,"INFECTED","NOT_INFECTED"))) %>%
            select(-infect) %>%
            rename(response=PABD)
          Error_metric <- "MeanDecreaseAccuracy"
        } else {
          mf_pred_temp_bd <- mf_pred_temp %>%
            filter(infect >0) %>%
            rename(response=infect)
          Error_metric <- "X.IncMSE"
          
        } 
        
        # Pick training and validation set
        set.seed(423)
        r <- sample(nrow(mf_pred_temp_bd), 0.7*nrow(mf_pred_temp_bd), replace=FALSE)
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
        all_RF_predictBD[[paste(otu_type, p_type, bd_type, sep="_")]][[paste0("mf_pred_",otu_type,"_",p_type,"_",bd_type)]] <- mf_pred_temp_bd
        all_RF_predictBD[[paste(otu_type, p_type, bd_type, sep="_")]][[paste0("RF_",otu_type,"_",p_type,"_",bd_type)]] <- RF
        all_RF_predictBD[[paste(otu_type, p_type, bd_type, sep="_")]][[paste0("importance_RF_",otu_type,"_",p_type,"_",bd_type)]] <- importance_RF
        all_RF_predictBD[[paste(otu_type, p_type, bd_type, sep="_")]][[paste0("testSet_error",otu_type,"_",p_type,"_",bd_type)]] <- error
        all_RF_predictBD[[paste(otu_type, p_type, bd_type, sep="_")]][[paste0("test_train_comparison",otu_type,"_",p_type,"_",bd_type)]] <- data.frame(Test_pred=testSet_predictions, Test_obs=testSet$response)
        
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
  close(pb)
  save(all_RF_predictBD, file="./5_random_forest/all_RF_predictBD.RData")
}
load("./5_random_forest/all_RF_predictBD.RData")


write.table(make_summary_RF(all_output = all_RF_predictBD, otu_type = "rank", bd_type = "infect")
            , quote = FALSE
            , sep = "\t"
            , row.names = FALSE
            , col.names=TRUE
            , file="./5_random_forest/RF_rank_infect.txt")
write.table(make_summary_RF(all_output = all_RF_predictBD, otu_type = "rank", bd_type = "PABD")
            , quote = FALSE
            , sep = "\t"
            , row.names = FALSE
            , col.names=TRUE
            , file="./5_random_forest/RF_rank_PABD.txt")
write.table(make_summary_RF(all_output = all_RF_predictBD, otu_type = "PA", bd_type = "infect")
            , quote = FALSE
            , sep = "\t"
            , row.names = FALSE
            , col.names=TRUE
            , file="./5_random_forest/RF_PA_infect.txt")
write.table(make_summary_RF(all_output = all_RF_predictBD, otu_type = "PA", bd_type = "PABD")
            , quote = FALSE
            , sep = "\t"
            , row.names = FALSE
            , col.names=TRUE
            , file="./5_random_forest/RF_PA_PABD.txt")
#### Are the OTUs predicting infection just OTUs that correlate with that particular species?
all_RF_predictSpecies <- list()

RERUN_species=FALSE
if (RERUN_species) {
  total <- 2*2
  # create progress bar
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  i <- 1
  for ( otu_type in c("rank","PA") ){
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

  save(all_RF_predictSpecies, file="./5_random_forest/all_RF_predictSpecies.RData")
}
load("./5_random_forest/all_RF_predictSpecies.RData")
all_RF_predictSpecies$PA$importance_RF_PA %>%
  select(c("taxonomy","inhibitory", "MeanDecreaseAccuracy")) %>%
  head(20)
all_RF_predictSpecies$rank$importance_RF_rank %>%
  select(c("taxonomy","inhibitory", "MeanDecreaseAccuracy")) %>%
  head(20)

####### Do a version with "response" to BD infection
bd_load_only <- mf_alt_filt_final %>%
  select(SampleID, Bd_load)

mf_rank_post <- mf_rank %>%
  filter(prepost == "Post", Bd_exposure == "Bd-exposed") %>%
  left_join(bd_load_only) %>%
  separate(indivID, into=c("species","indiv")) %>%
  select(-c(SampleID, prepost, indiv, Bd_exposure))
mf_PA_post <- mf_PA %>%
  filter(prepost == "Post", Bd_exposure == "Bd-exposed") %>%
  left_join(bd_load_only) %>%
  separate(indivID, into=c("species","indiv")) %>%
  select(-c(SampleID, prepost, indiv, Bd_exposure))


RERUN_speciespost=FALSE
if ( RERUN_speciespost ) {
  total <- 2*2
  # create progress bar
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  i <- 1
  
  all_RF_predict_post <- list()
  for ( sp in c("Anbo","Rhma","Osse","Raca","Rapi")) {
    # Check to see if there were any infected time points
    max_load <- mf_alt_filt_final %>%
      filter(species == sp) %>%
      mutate(count = ifelse(Bd_load>0,1,0)) %>%
      select(count) %>%
      sum()
    
    if (max_load<5) {
      next
    } 
    for ( otu_type in c("rank","PA")) {
      mf_temp <- get(paste0("mf_",otu_type,"_post")) %>%
        filter(species==sp)
      for ( bd_type in c("PABD","infect")) {
        
        if ( bd_type == "PABD" ) {
          mf_temp <- mf_temp %>%
            mutate(Bd = ifelse(Bd_load>0, "INFECT","NOT_INFECT")) %>%
            mutate(Bd = factor(Bd)) %>%
            select(-Bd_load) %>%
            mutate(species = factor(species))
          Error_metric <- "MeanDecreaseAccuracy"
        } else {
          mf_temp <- mf_temp %>%
            rename(Bd = Bd_load)%>%
            mutate(species = factor(species))
          Error_metric <- "X.IncMSE"
        }
        set.seed(423)
        r <- sample(nrow(mf_temp), 0.7*nrow(mf_temp), replace=FALSE)
        trainSet <- mf_temp[r,]
        testSet <- mf_temp[-r,]
        
        RF_post <- randomForest(Bd ~ ., data=trainSet,importance=TRUE)
        testSet_predictions <- predict(RF_post, testSet, type="class")
        # table(testSet_predictions,testSet$species)
        error <- mean((testSet_predictions-testSet$Bd)^2)
        
        # Get importance
        importance_RF_post <- data.frame(importance(RF_post)) %>%
          mutate(OTU=rownames(importance(RF_post))) 
        importance_RF_post[,c("taxonomy","inhibitory")] <- Taxa[match(importance_RF_post$OTU, Taxa$Sequence), c("Taxa","inhibitory")]
        importance_RF_post <- importance_RF_post %>%
          mutate(taxonomy=ifelse(is.na(taxonomy), OTU, taxonomy)) %>%
          select(taxonomy, everything()) %>%
          arrange(-get(paste0(Error_metric))) %>%
          select(-OTU)
        
        
        # Save data
        all_RF_predict_post[[paste(otu_type,"_",bd_type)]][[sp]][[paste0("mf_pred_",otu_type,"_",bd_type)]] <- mf_temp
        all_RF_predict_post[[paste(otu_type,"_",bd_type)]][[sp]][[paste0("RF_",otu_type,"_",bd_type)]] <- RF_post
        all_RF_predict_post[[paste(otu_type,"_",bd_type)]][[sp]][[paste0("importance_RF_",otu_type,"_",bd_type)]] <- importance_RF_post
        all_RF_predict_post[[paste(otu_type,"_",bd_type)]][[sp]][[paste0("testSet_error",otu_type,"_",bd_type)]] <- error
        all_RF_predict_post[[paste(otu_type,"_",bd_type)]][[sp]][[paste0("test_train_comparison","_",bd_type)]] <- data.frame(Test_pred=testSet_predictions, Test_obs=testSet$Bd_load)
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
  save(all_RF_predict_post, file="./5_random_forest/all_RF_predict_post.RData")
}
load("./5_random_forest/all_RF_predict_post.RData")


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
all_RF_predictBD$rank_withp_infect$importance_RF_rank_withp_infect %>%
  mutate(PercIncrMSE_rank_infect = X.IncMSE) %>%
  select(taxonomy, PercIncrMSE_rank_infect) %>%
  full_join(all_RF_predictBD$PA_withp_infect$importance_RF_PA_withp_infect) %>%
  mutate(PercIncrMSE_PA_infect = X.IncMSE) %>%
  ggplot() +
  geom_point(aes(x=PercIncrMSE_PA_infect, y=PercIncrMSE_rank_infect))
all_RF_predictBD$rank_withp_PABD$importance_RF_rank_withp_PABD %>%
  mutate(MeanDecrAccuracy_rank_PABD = MeanDecreaseAccuracy) %>%
  select(taxonomy, MeanDecrAccuracy_rank_PABD) %>%
  full_join(all_RF_predictBD$PA_withp_PABD$importance_RF_PA_withp_PABD) %>%
  mutate(MeanDecrAccuracy_PA_PABD = MeanDecreaseAccuracy) %>%
  ggplot() +
  geom_point(aes(x=MeanDecrAccuracy_PA_PABD, y=MeanDecrAccuracy_rank_PABD))

## Everything
all_RF_predictBD$rank_withp_infect$importance_RF_rank_withp_infect %>%
  filter(!is.na(inhibitory)) %>%
  ggplot() +
  geom_point(aes(x=X.IncMSE, y=inhibitory), position = position_jitter(width=0, height=0.1)) +
  geom_vline(data=data.frame(metric=comm_metrics, value=all_RF_predictBD$rank_withp_infect$importance_RF_rank_withp_infect[match(comm_metrics,all_RF_predictBD$rank_withp_infect$importance_RF_rank_withp_infect$taxa),"X.IncMSE" ] )
            , mapping=aes(xintercept=value, col=comm_metrics)
            , lwd=2, alpha=0.5) +
  scale_color_manual(values=comm_metrics_col)

# all_RF_predictBD$PA_withp_infect$importance_RF_PA_withp_infect %>%
#   filter(!is.na(inhibitory)) %>%
#   ggplot() +
#   geom_point(aes(x=X.IncMSE, y=inhibitory), position = position_jitter(width=0, height=0.1)) +
#   geom_vline(data=data.frame(metric=comm_metrics, value=all_RF_predictBD$PA_withp_infect$importance_RF_PA_withp_infect[match(comm_metrics,all_RF_predictBD$PA_withp_infect$importance_RF_PA_withp_infect$taxa),"X.IncMSE" ] )
#              , mapping=aes(xintercept=value, col=comm_metrics)
#              , lwd=2, alpha=0.5) +
#   scale_color_manual(values=comm_metrics_col)


all_RF_predictBD$rank_withp_PABD$importance_RF_rank_withp_PABD %>%
  filter(!is.na(inhibitory)) %>%
  ggplot() +
  geom_point(aes(x=MeanDecreaseAccuracy, y=inhibitory), position = position_jitter(width=0, height=0.1)) +
  geom_vline(data=data.frame(metric=comm_metrics, value=all_RF_predictBD$rank_withp_PABD$importance_RF_rank_withp_PABD[match(comm_metrics,all_RF_predictBD$rank_withp_PABD$importance_RF_rank_withp_PABD$taxa),"MeanDecreaseAccuracy" ] )
             , mapping=aes(xintercept=value, col=comm_metrics)
             , lwd=2, alpha=0.5) +
  scale_color_manual(values=comm_metrics_col)


all_RF_predictBD$PA_withp_PABD$importance_RF_PA_withp_PABD %>%
  filter(!is.na(inhibitory)) %>%
  ggplot() +
  geom_point(aes(x=MeanDecreaseAccuracy, y=inhibitory), position = position_jitter(width=0, height=0.1)) +
  geom_vline(data=data.frame(metric=comm_metrics, value=all_RF_predictBD$PA_withp_PABD$importance_RF_PA_withp_PABD[match(comm_metrics,all_RF_predictBD$PA_withp_PABD$importance_RF_PA_withp_PABD$taxa),"MeanDecreaseAccuracy" ] )
             , mapping=aes(xintercept=value, col=comm_metrics)
             , lwd=2, alpha=0.5) +
  scale_color_manual(values=comm_metrics_col)



#### Looking at importance #######


# Rank: PABD
all_RF_predictBD$rank_onlyp_PABD$testSet_errorrank_onlyp_PABD
all_RF_predictBD$rank_withp_PABD$testSet_errorrank_withp_PABD
all_RF_predictBD$rank_nop_PABD$testSet_errorrank_nop_PABD

all_RF_predictBD$rank_onlyp_PABD$RF_rank_onlyp_PABD$confusion
all_RF_predictBD$rank_withp_PABD$RF_rank_withp_PABD$confusion
all_RF_predictBD$rank_nop_PABD$RF_rank_nop_PABD$confusion

all_RF_predictBD$rank_onlyp_PABD$test_train_comparisonrank_onlyp_PABD
all_RF_predictBD$rank_withp_PABD$test_train_comparisonrank_withp_PABD
all_RF_predictBD$rank_nop_PABD$test_train_comparisonrank_nop_PABD

## Plot
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
  geom_vline(aes(xintercept=1.5), lwd=4)


# Rank: infect
all_RF_predictBD$rank_onlyp_infect$testSet_errorrank_onlyp_infect
all_RF_predictBD$rank_withp_infect$testSet_errorrank_withp_infect
all_RF_predictBD$rank_nop_infect$testSet_errorrank_nop_infect

all_RF_predictBD$rank_onlyp_infect$RF_rank_onlyp_infect
all_RF_predictBD$rank_withp_infect$RF_rank_withp_infect
all_RF_predictBD$rank_nop_infect$RF_rank_nop_infect

ss_tot <- sum((all_RF_predictBD$rank_onlyp_infect$test_train_comparisonrank_onlyp_infect$Test_obs-mean(all_RF_predictBD$rank_onlyp_infect$test_train_comparisonrank_onlyp_infect$Test_obs))^2)
ss_err <- sum((all_RF_predictBD$rank_onlyp_infect$test_train_comparisonrank_onlyp_infect$Test_obs-all_RF_predictBD$rank_onlyp_infect$test_train_comparisonrank_onlyp_infect$Test_pred)^2)
1-ss_err/ss_tot

ss_tot <- sum((all_RF_predictBD$rank_withp_infect$test_train_comparisonrank_withp_infect$Test_obs-mean(all_RF_predictBD$rank_withp_infect$test_train_comparisonrank_withp_infect$Test_obs))^2)
ss_err <- sum((all_RF_predictBD$rank_withp_infect$test_train_comparisonrank_withp_infect$Test_obs-all_RF_predictBD$rank_withp_infect$test_train_comparisonrank_withp_infect$Test_pred)^2)
1-ss_err/ss_tot

ss_tot <- sum((all_RF_predictBD$rank_nop_infect$test_train_comparisonrank_nop_infect$Test_obs-mean(all_RF_predictBD$rank_nop_infect$test_train_comparisonrank_nop_infect$Test_obs))^2)
ss_err <- sum((all_RF_predictBD$rank_nop_infect$test_train_comparisonrank_nop_infect$Test_obs-all_RF_predictBD$rank_nop_infect$test_train_comparisonrank_nop_infect$Test_pred)^2)
1-ss_err/ss_tot

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

data.frame(Observed=all_RF_predictBD$rank_onlyp_infect$test_train_comparisonrank_onlyp_infect[,2]
           ,Community_level_only=all_RF_predictBD$rank_onlyp_infect$test_train_comparisonrank_onlyp_infect[,1]
           ,Community_and_OTU_level=all_RF_predictBD$rank_withp_infect$test_train_comparisonrank_withp_infect[,1]
           ,OTU_level_only=all_RF_predictBD$rank_nop_infect$test_train_comparisonrank_nop_infect[,1]
) %>%
  gather(Community_level_only,Community_and_OTU_level,OTU_level_only, key=Model, value=Predicted) %>%
  mutate(Residuals=Predicted-Observed) %>%
  ggplot() +
  geom_point(aes(x=Model, y=Residuals)) 

## PA: PABD
all_RF_predictBD$PA_onlyp_PABD$testSet_errorPA_onlyp_PABD
all_RF_predictBD$PA_withp_PABD$testSet_errorPA_withp_PABD
all_RF_predictBD$PA_nop_PABD$testSet_errorPA_nop_PABD

all_RF_predictBD$PA_onlyp_PABD$test_train_comparisonPA_onlyp_PABD
all_RF_predictBD$PA_withp_PABD$test_train_comparisonPA_withp_PABD
all_RF_predictBD$PA_nop_PABD$test_train_comparisonPA_nop_PABD

## PA: infect
all_RF_predictBD$PA_onlyp_infect$testSet_errorPA_onlyp_infect
all_RF_predictBD$PA_withp_infect$testSet_errorPA_withp_infect
all_RF_predictBD$PA_nop_infect$testSet_errorPA_nop_infect

data.frame(Observed=all_RF_predictBD$PA_onlyp_infect$test_train_comparisonPA_onlyp_infect[,2]
           ,Community_level_only=all_RF_predictBD$PA_onlyp_infect$test_train_comparisonPA_onlyp_infect[,1]
           ,Community_and_OTU_level=all_RF_predictBD$PA_withp_infect$test_train_comparisonPA_withp_infect[,1]
           ,OTU_level_only=all_RF_predictBD$PA_nop_infect$test_train_comparisonPA_nop_infect[,1]
) %>%
  gather(Community_level_only,Community_and_OTU_level,OTU_level_only, key=Model, value=Predicted) %>%
  ggplot() +
  geom_point(aes(x=Observed, y=Predicted, col=Model)) +
  geom_abline(aes(intercept=0, slope=1))+
  xlim(0,8)+ylim(0,8)

### Are inhibitory bacteria more likely to be predictive?
rank_infect_filtered <- all_RF_predictBD$rank_nop_infect$importance_RF_rank_nop_infect %>%
  filter(X.IncMSE>0, !is.na(inhibitory)) %>%
  mutate(inhibitory=factor(inhibitory))
anova(lm(X.IncMSE ~ inhibitory, data=rank_infect_filtered))
summary(glm(inhibitory ~ X.IncMSE, data=rank_infect_filtered, family = binomial))
rank_infect_filtered %>%
  ggplot() +geom_point(aes(x=inhibitory, y=X.IncMSE), position = position_jitter(height=0, width=0.1))

rank_PABD_filtered <- all_RF_predictBD$rank_nop_PABD$importance_RF_rank_nop_PABD %>%
  filter(MeanDecreaseAccuracy>0, !is.na(inhibitory)) %>%
  mutate(inhibitory=factor(inhibitory))
anova(lm(MeanDecreaseAccuracy ~ inhibitory, data=rank_PABD_filtered))
summary(glm(inhibitory ~ MeanDecreaseAccuracy, data=rank_PABD_filtered, family = binomial))
rank_PABD_filtered %>%
  ggplot() +geom_point(aes(x=inhibitory, y=MeanDecreaseAccuracy), position=position_jitter(height=0, width=0.1))



