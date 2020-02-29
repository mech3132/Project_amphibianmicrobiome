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


# Random forest loop to find best mtry
randomForest_findmtry <- function(full_set, type, response, frac_train=0.7, nTrainSets=1, nmtry=10) {
  # First, estimate number of mtrys to start with
  if ( type=="class") {
    ave <- sqrt(ncol(full_set)-1)
    lwr <- mean(c(0,ave))
    upr <- mean(c(ave, ncol(full_set)-1 ))
  } else if ( type == "regression") {
    ave <- (ncol(full_set)-1)/2
    lwr <- mean(c(0,ave))
    upr <- mean(c(ave, ncol(full_set)-1 ))
  }
  # Rename response
  full_set <- full_set %>%
    rename(response=response)
  for ( s in 1:nTrainSets) {
    rTrain <- sample(nrow(full_set), size=0.7*nrow(full_set), replace=FALSE)
    trainSet <- full_set[rTrain,]
    testSet <- full_set[-rTrain,]
    a=data.frame()
    for (i in round(seq(lwr,upr, length.out=nmtry))) {
      print(paste0("Currently doing i=",i))
      tempModel <- randomForest(response ~ ., data = trainSet, mtry=i)
      predValid <- predict(tempModel, testSet, type = type)
      if ( type=="class") {
        error <- mean(predValid != testSet$response)
      } else if (type == "regression") {
        error <-  sum((predValid-testSet$response)^2)
      }
      a <- rbind(a, data.frame(error=error, mtry=i, set=s))
    }
  }
  return(a)
}



###### Match inhibitory to non-inhibitory #####



