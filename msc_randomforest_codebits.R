#!/bin/bash Rscript

###### Modular code msc ########


# Relative abundance
count_to_relabund <- function(OTU, OTUcol="#OTU ID") {
  otu_col <- which(colnames(OTU)==OTUcol)
  otus <- OTU[,otu_col]
  relabund <- OTU[,-otu_col]
  for ( c in 1:ncol(relabund)) {
    relabund[,c] <- relabund[,c]/sum(relabund[,c])
  }
  rownames(relabund) <- otus
  return(relabund)
}

# Ranking
count_to_percentileranking <- function(OTU, OTUcol="#OTU ID") {
  otu_col <- which(colnames(OTU)==OTUcol)
  otus <- OTU[,otu_col]
  ranks <- OTU[,-otu_col]
  for ( c in 1:ncol(ranks)) {
    temp_rank <- rank(-ranks[,c], ties.method = "min")
    ranks[,c] <- 1-temp_rank/max(temp_rank)
  }
  rownames(ranks) <- otus
  return(ranks)
}


count_to_presenceabsence <- function(OTU, OTUcol="#OTU ID") {
  otu_col <- which(colnames(OTU)==OTUcol)
  otus <- OTU[,otu_col]
  PA <- OTU[,-otu_col]
  PA[PA>0] <- 1
  rownames(PA) <- otus
  return(PA)
}

# Turn OTU table into variable table with mf
mf_otu_combine <- function(OTU, MF) {
  library(tidyverse)
  OTU_t <- data.frame(t(OTU))
  OTU_t$SampleID <- rownames(OTU_t)
  combined <- MF %>%
    left_join(OTU_t)
  return(combined)
}

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
    taxa_full[r,"Taxa"] <- make.unique(final, sep = "-")
  }
  taxa_full[,"Taxa"] <- make.unique(taxa_full[,"Taxa"], sep="-")
  
  return(as.data.frame(taxa_full))
}

# 
# # Random forest loop to find best mtry
# randomForest_findmtry <- function(full_set, type, response, frac_train=0.7, nTrainSets=1, nmtry=10) {
#   # First, estimate number of mtrys to start with
#   if ( type=="class") {
#     ave <- sqrt(ncol(full_set)-1)
#     lwr <- mean(c(0,ave))
#     upr <- mean(c(ave, ncol(full_set)-1 ))
#   } else if ( type == "regression") {
#     ave <- (ncol(full_set)-1)/2
#     lwr <- mean(c(0,ave))
#     upr <- mean(c(ave, ncol(full_set)-1 ))
#   }
#   # Rename response
#   full_set <- full_set %>%
#     rename(response=response)
#   for ( s in 1:nTrainSets) {
#     rTrain <- sample(nrow(full_set), size=0.7*nrow(full_set), replace=FALSE)
#     trainSet <- full_set[rTrain,]
#     testSet <- full_set[-rTrain,]
#     a=data.frame()
#     for (i in round(seq(lwr,upr, length.out=nmtry))) {
#       print(paste0("Currently doing i=",i))
#       tempModel <- randomForest(response ~ ., data = trainSet, mtry=i)
#       predValid <- predict(tempModel, testSet, type = type)
#       if ( type=="class") {
#         error <- mean(predValid != testSet$response)
#       } else if (type == "regression") {
#         error <-  sum((predValid-testSet$response)^2)
#       }
#       a <- rbind(a, data.frame(error=error, mtry=i, set=s))
#     }
#   }
#   return(a)
# }

# 
# make_summary_RF <- function(all_output, otu_type=c("PA","rank","count"), bd_type=c("PABD","infect"), prop=c(0.2,0.3,0.4,0.5,0.6,0.7,0.8), repl=1) {
#   # all_output <- all_RF_predictBD
#   # otu_type <- "rank"
#   # bd_type="infect"
#   # Create empty matrix
#   # Get error explained by each model into table
#   base_df <- data.frame(Model=c("onlyp",rep(NA,12),"nop",rep(NA,14),"withp",rep(NA,14)))
#   if (bd_type == "infect" ) {
#     summary_randomforest <- cbind(base_df, data.frame(Train_set_MSE=NA, Validation_set_MSE=NA, Top_importance=NA, IncrMSE=NA))
#   } else {
#     summary_randomforest <- cbind(base_df, data.frame(Train_set_error_rate=NA, Validation_set_error_rate=NA, Top_importance=NA, MeanDecrAccuracy=NA))
#     
#   }
#   
#   for ( p_type in c("onlyp","nop","withp") ) {
#     curr_cat <- all_output[[paste(otu_type, p_type, bd_type, prop, repl,sep="_")]]
#     curr_model <- curr_cat[["RF"]]
#     error_rate_predict <- curr_cat[["testSetError"]]
#     if ( bd_type=="PABD" ) {
#       error_rate <- curr_model$err.rate[,"OOB"][nrow(curr_model$err.rate)]*100
#       var_explained <- NA
#       final_error <- paste0(round(error_rate, 3), "%")
#       valid_set_error <- paste0(round(error_rate_predict,3),"%")
#     } else if (bd_type=="infect") {
#       error_rate <- curr_model$mse[length(curr_model$mse)]
#       var_explained <- curr_model$rsq[length(curr_model$rsq)]
#       final_error <- paste0(round(error_rate, 3), ", ", round(var_explained*100, 2),"% var explained")
#       valid_set_error <- paste0(round(error_rate_predict,3))
#     }
#     current_row <- ifelse(p_type=="onlyp",0,ifelse(p_type=="nop",13,28))
#     
#     # Finally, get top importance factors
#     curr_imp <- curr_cat[["importance"]]
#     if ( bd_type == "infect") {
#       col_imp <- "X.IncMSE"
#       col_phrase <- "%IncrMSE="
#     } else if ( bd_type == "PABD" ) {
#       col_imp <- "MeanDecreaseAccuracy"
#       col_phrase <- "MeanDecrAccuracy="
#     }
#     
#     # col_imp changes if it's p model
#     # Only 13 factors in p model
#     if ( p_type=="onlyp") {
#       l <- 13
#       # metric <- "Metric"
#       # if ( bd=="infect" ) {
#       #   col_imp <- "%IncMSE"
#       # }
#       metric <- "taxonomy"
#     } else {
#       l <- 15
#       metric <- "taxonomy"
#       
#     }
#     
#     for ( tax in 1:l) {
#       if (is.null(curr_imp$inhibitory[tax])){
#         inhib_phrase <- ""
#       } else if (is.na(curr_imp$inhibitory[tax]) ) {
#         inhib_phrase <- ""
#       } else if (curr_imp$inhibitory[tax]==1) {
#         inhib_phrase <- " (inhibitory)"
#       } else {
#         inhib_phrase <- ""
#       }
#       
#       if ( bd_type == "infect" ) {
#         summary_randomforest[current_row+1, c("Train_set_MSE","Validation_set_MSE")] <- c(final_error, valid_set_error)
#         summary_randomforest[current_row+tax,"Top_importance"] <- paste0(curr_imp[tax,metric], inhib_phrase)
#         summary_randomforest[current_row+tax,"IncrMSE"] <- paste0(col_phrase, round(curr_imp[tax,col_imp],2))
#       } else if ( bd_type == "PABD") {
#         summary_randomforest[current_row+1, c("Train_set_error_rate","Validation_set_error_rate")] <- c(final_error, valid_set_error)
#         summary_randomforest[current_row+tax,"Top_importance"] <- paste0(curr_imp[tax,metric], inhib_phrase)
#         summary_randomforest[current_row+tax,"MeanDecrAccuracy"] <- paste0(col_phrase, round(curr_imp[tax,col_imp],2))
#       }
#     }
#   }
#   
#   return(summary_randomforest)
# }

# make_summary_RF <- function(all_output, otu_type=c("PA","rank","count"), bd_type=c("PABD","infect"), repl=c(1,2,3)) {
#   # all_output <- all_RF_predictBD
#   # otu_type <- "rank"
#   # bd_type="infect"
#   # Create empty matrix
#   # Get error explained by each model into table
#   base_df <- data.frame(Model=c("onlyp",rep(NA,12),"nop",rep(NA,14),"withp",rep(NA,14)))
#   if (bd_type == "infect" ) {
#     summary_randomforest <- cbind(base_df, data.frame(Train_set_MSE=NA, Validation_set_MSE=NA, Top_importance=NA, IncrMSE=NA))
#   } else {
#     summary_randomforest <- cbind(base_df, data.frame(Train_set_error_rate=NA, Validation_set_error_rate=NA, Top_importance=NA, MeanDecrAccuracy=NA))
#     
#   }
#   
#   for ( p_type in c("onlyp","nop","withp") ) {
#     curr_cat <- all_output[[paste(otu_type, p_type, bd_type, repl,sep="_")]]
#     curr_model <- curr_cat[["RF"]]
#     error_rate_predict <- curr_cat[[testSetError]]
#     if ( bd_type=="PABD" ) {
#       error_rate <- curr_model$err.rate[,"OOB"][nrow(curr_model$err.rate)]*100
#       var_explained <- NA
#       final_error <- paste0(round(error_rate, 3), "%")
#       valid_set_error <- paste0(round(error_rate_predict,3),"%")
#     } else if (bd_type=="infect") {
#       error_rate <- curr_model$mse[length(curr_model$mse)]
#       var_explained <- curr_model$rsq[length(curr_model$rsq)]
#       final_error <- paste0(round(error_rate, 3), ", ", round(var_explained*100, 2),"% var explained")
#       valid_set_error <- paste0(round(error_rate_predict,3))
#     }
#     current_row <- ifelse(p_type=="onlyp",0,ifelse(p_type=="nop",13,28))
#     
#     # Finally, get top importance factors
#     curr_imp <- curr_cat[["importance"]]
#     if ( bd_type == "infect") {
#       col_imp <- "X.IncMSE"
#       col_phrase <- "%IncrMSE="
#     } else if ( bd_type == "PABD" ) {
#       col_imp <- "MeanDecreaseAccuracy"
#       col_phrase <- "MeanDecrAccuracy="
#     }
#     
#     # col_imp changes if it's p model
#     # Only 13 factors in p model
#     if ( p_type=="onlyp") {
#       l <- 13
#       # metric <- "Metric"
#       # if ( bd=="infect" ) {
#       #   col_imp <- "%IncMSE"
#       # }
#       metric <- "taxonomy"
#     } else {
#       l <- 15
#       metric <- "taxonomy"
#       
#     }
#     
#     for ( tax in 1:l) {
#       if (is.null(curr_imp$inhibitory[tax])){
#         inhib_phrase <- ""
#       } else if (is.na(curr_imp$inhibitory[tax]) ) {
#         inhib_phrase <- ""
#       } else if (curr_imp$inhibitory[tax]==1) {
#         inhib_phrase <- " (inhibitory)"
#       } else {
#         inhib_phrase <- ""
#       }
#       
#       if ( bd_type == "infect" ) {
#         summary_randomforest[current_row+1, c("Train_set_MSE","Validation_set_MSE")] <- c(final_error, valid_set_error)
#         summary_randomforest[current_row+tax,"Top_importance"] <- paste0(curr_imp[tax,metric], inhib_phrase)
#         summary_randomforest[current_row+tax,"IncrMSE"] <- paste0(col_phrase, round(curr_imp[tax,col_imp],2))
#       } else if ( bd_type == "PABD") {
#         summary_randomforest[current_row+1, c("Train_set_error_rate","Validation_set_error_rate")] <- c(final_error, valid_set_error)
#         summary_randomforest[current_row+tax,"Top_importance"] <- paste0(curr_imp[tax,metric], inhib_phrase)
#         summary_randomforest[current_row+tax,"MeanDecrAccuracy"] <- paste0(col_phrase, round(curr_imp[tax,col_imp],2))
#       }
#     }
#   }
#   
#   return(summary_randomforest)
# }


get_mf_pred_temp <- function(p_type, otu_type, bd_type, p, keepSp=TRUE) {
  if ( p_type == "onlyp") {
    mf_pred_temp <- p %>%
      separate(indivID, into=c("species","indiv"), remove=FALSE) %>%
      mutate(species=factor(species)) %>%
      arrange(indivID) %>%
      select(-indiv)
    
  } else {
    p_to_add <- get(ifelse(p_type == "withp", "p", "p_infectonly"))
    mf_temp <- get(paste0("mf_",otu_type))
    # Filter mf and re-arrange to be compatible with randomforest
    mf_pred_temp <- mf_temp %>%
      filter(prepost=="Pre", Bd_exposure=="Bd-exposed") %>%
      select(-c("SampleID","prepost","Bd_exposure")) %>%
      gather(-c(indivID),key=OTU, value=rank) %>%
      mutate(rank=as.numeric(rank)) %>%
      group_by(indivID, OTU) %>%
      summarize(aveRank=mean(rank)) %>%
      ungroup() %>%
      separate(indivID, into=c("species","indiv"), remove=FALSE) %>%
      select(-indiv) %>%
      filter(!is.na(aveRank)) %>%
      spread(key=OTU, value=aveRank)  %>% 
      left_join(p_to_add) %>%
      arrange(indivID) %>%
      mutate(species=factor(species))
    
  }
  
  # Mutate infection depending on question
  if (bd_type == "PABD" ) {
    mf_pred_temp_bd  <- mf_pred_temp %>%
      mutate(PABD = factor(ifelse(infect>0,"INFECTED","NOT_INFECTED"))) %>%
      select(-infect) %>%
      rename(response=PABD) %>%
      select(-c(indivID))
    indiv_list <- mf_pred_temp%>%
      mutate(PABD = factor(ifelse(infect>0,"INFECTED","NOT_INFECTED"))) %>%
      select(-infect) %>%
      rename(response=PABD) %>%
      select(indivID)
    # Error_metric <- "MeanDecreaseAccuracy"
  } else {
    mf_pred_temp_bd <- mf_pred_temp %>%
      filter(infect >0) %>%
      rename(response=infect) %>%
      select(-c(indivID))
    indiv_list <- mf_pred_temp %>%
      filter(infect >0) %>%
      rename(response=infect) %>%
      select(indivID)
  } 
  
  if ( !keepSp ) {
    mf_pred_temp_bd <- mf_pred_temp_bd %>%
      select(-species)
  }
  return(list(mf_pred_temp_bd,indiv_list))
  
}




### get rsq #####
getR2 <- function(Test_pred, Test_obs) {
  ss_tot <- sum((Test_obs-mean(Test_obs))^2)
  ss_err <- sum((Test_obs-Test_pred)^2)
  return(1-ss_err/ss_tot)
}


extractComparison <- function(allOutput, bd_type=c("infect","PABD"), otu_types = c("PA","count"),replicates=c(1), prop_levels=c(0.2,0.3,0.4,0.5,0.6,0.7,0.8), p_types=c("onlyp","nop","withp")) {
  # allOutput = RF_infect_paramsweep_1  
  # RF_infect_paramsweep_1$PA_nop_infect_0.8_3
  first=TRUE
  # bd_type="infect"
    for (otu_type in otu_types) {
      # otu_type="PA"
      for ( p_type in p_types ) {
        # p_type="nop"
        for ( prop in prop_levels) {
          # prop=0.8
          for (repl in replicates) {
            curr_model <- allOutput[[paste(otu_type,p_type,bd_type,prop,repl, sep="_")]]
            if (first) {
              dat_comparison <- cbind(curr_model$test_train_comparison, p_type=p_type, prop=prop, repl=repl, otu_type=otu_type)
              first <- FALSE
            } else {
              dat_comparison <- rbind(dat_comparison, cbind(curr_model$test_train_comparison, p_type=p_type, prop=prop, repl=repl, otu_type=otu_type))
            }
          }
        }
      }
    }
    return(dat_comparison)
}

extractTrainError <- function(allOutput, bd_type=c("infect","PABD"), otu_types = c("PA","count"), error_type=c("R2","MSE","errorRate"), replicates=c(1:3), prop_levels=c(0.2,0.3,0.4,0.5,0.6,0.7,0.8), p_types=c("onlyp","nop","withp")) {
  # prop_levels <- c(0.2,0.3,0.4,0.5,0.6,0.7,0.8)
  # p_types <- c("onlyp","nop","withp")
  # otu_types=c("PA","count")
  # bd_type="PABD"
  # error_type="comparison"
  # replicates=c(1:3)
  # allOutput= RF_PABD_parasweep_1
  first=TRUE
  for (otu_type in otu_types) {
    for ( p_type in p_types ) {
      for ( prop in prop_levels) {
        for (repl in replicates) {
          # What row to insert
          curr_model <- allOutput[[paste(otu_type,p_type,bd_type,prop,repl, sep="_")]]
          
          if ( bd_type == "infect") {
            if ( error_type == "R2") {
              toAdd <- data.frame(error=as.numeric(curr_model$RF$rsq[length(curr_model$RF$rsq)]), otu_type=otu_type, p_type=p_type, prop=as.numeric(prop), repl=as.numeric(repl), error_type=error_type)
            } else if ( error_type == "MSE") {
              toAdd <- data.frame(error=as.numeric(curr_model$RF$mse[length(curr_model$RF$mse)]), otu_type=otu_type, p_type=p_type, prop=as.numeric(prop), repl=as.numeric(repl), error_type=error_type)
            } else {
              print("Invalid error type")
              toAdd<- NA
              toAdd<- NA
            }
          } else {
            toAdd <- data.frame(error=as.numeric(curr_model$RF$err.rate[,"OOB"][nrow(curr_model$RF$err.rate)]), otu_type=otu_type, p_type=p_type, prop=as.numeric(prop), repl=as.numeric(repl), error_type=error_type)
          }
          
          if (first) {
            dat_train <- toAdd
            first <-FALSE
          } else {
            dat_train <- rbind(dat_train, toAdd)
          }
          
        }
        
      }
    }
  }
  
  return(dat_train)
  
    

}

extractTrainErrorLOO <- function(allOutput, bd_type, otu_type, p_type) {
  first <- TRUE
  for ( curr_model in allOutput ) {
    indivID <- curr_model$test_train_comparison$indivID
    
    if ( bd_type == "infect") {
      error_type = "MSE"
      toAdd <- data.frame(error=as.numeric(curr_model$RF$mse[length(curr_model$RF$mse)]), otu_type=otu_type, p_type=p_type, error_type=error_type, indiv_LOO = indivID)
    } else {
      error_type = "errorRate"
      toAdd <- data.frame(train_error=as.numeric(curr_model$RF$err.rate[,"OOB"][nrow(curr_model$RF$err.rate)]), otu_type=otu_type, p_type=p_type, error_type=error_type, indiv_LOO = indivID)
    }
    
    if (first) {
      dat_train <- toAdd
      first <-FALSE
    } else {
      dat_train <- rbind(dat_train, toAdd)
    }
    
  }
  return(dat_train)
  
}


### Extract
extractItemLOO<- function(allOutput=list(), bd_type, otu_type, p_type, extract="test_train_comparison")  {
  # allOutput <- RF_infect_onlyp_wsp_LOO
  first <- TRUE
  for ( item in allOutput ) {
    indivID <- item$test_train_comparison$indivID
    if (first) {
      temp <- item[[extract]]
      temp$indivID <- indivID
      
      all <- temp
      first <- FALSE
    } else {
      temp <- item[[extract]]
      temp$indivID <- indivID
      all <- rbind(all,temp)
    }
  }
  all <- cbind(all, bd_type=bd_type, otu_type=otu_type, p_type=p_type)
  return(all)
}

extractItem  <- function(allOutput=list(), bd_type=c("infect","PABD"), otu_type = c("rank","PA","count"), p_types=c("onlyp","nop","withp"), reps=c(1:10)
         , extract="test_train_comparison", repstart=1)  {
  i <- repstart
  first <- TRUE
  for ( R in allOutput) {
    for ( r in reps ) {
      # r=1
      # otu_type="PA"
      # p_types="onlyp"
      # bd_type="PABD"
      new <- cbind(R[[paste0(otu_type,"_",p_types,"_", bd_type,"_0.8","_",r)]][[paste0(extract)]],rep=i)
      i <- i+1
      if (first) {
        all <- new
      } else {
        all <- rbind(all,new)
      }
      first <- FALSE
    }
  }
  return(all)
}


# 
# extractImportance <- function(allOutput, bd_type=c("infect","PABD"), otu_type = c("rank","PA","count"), replicates=c(1), top=5) {
#   prop_levels <- c(0.2,0.3,0.4,0.5,0.6,0.7,0.8)
#   p_types <- c("onlyp","nop","withp")
#   
#   if ( bd_type == "infect") {
#     metric <- "X.IncMSE"
#   } else {
#     metric <- "MeanDecreaseAccuracy"
#   }
#   all_tax <- c()
#   for ( p_type in p_types ) {
#     for ( prop in prop_levels) {
#       for (repl in replicates) {
#         curr_model <- allOutput[[paste(otu_type,p_type,bd_type,prop,repl, sep="_")]]
#         all_tax <- c(all_tax,curr_model$importance[1:top,c("taxonomy")])
#       }
#     }
#   }
#   all_tax <- unique(all_tax)
#   
#   dat_importance <- data.frame(taxonomy=rep(all_tax, length(p_types)*length(prop_levels)*length(replicates))
#                                , importance=NA
#                                , model=rep(p_types, each=length(all_tax)*length(prop_levels)*length(replicates))
#                                , prop=rep(rep(prop_levels,each=length(all_tax)*length(replicates)),length(p_types))
#                                , replicates = 1)
#   for ( p_type in p_types ) {
#     for ( prop in prop_levels) {
#       for (repl in replicates) {
#         current_rows <- which(dat_importance$model==p_type& dat_importance$prop==prop)
#         curr_model <- allOutput[[paste(otu_type,p_type,bd_type,prop,repl, sep="_")]]
#         dat_importance[current_rows,"importance"] <- curr_model$importance[match(all_tax,curr_model$importance$taxonomy),metric]
#       }
#     }
#   }
#   return(dat_importance)
# }
