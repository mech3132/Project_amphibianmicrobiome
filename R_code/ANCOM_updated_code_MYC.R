#### ANCOM ######
# Original code sourced from: https://github.com/mka2136
# Adjustments are coded with "###MYC" in comments.
# I have NOT checked to see if adjusted and repeated designs will work-- I have ONLY modified it so that it will accept one main variable.


### Required libraries ####
library(exactRankTests)
library(nlme)
library(ggplot2)
library(cowplot) # For extracting legend ###MYC

ancom.W.myc = function(otu_data,var_data,
                       adjusted,repeated,
                       main.var,adj.formula,
                       repeat.var,long,rand.formula,
                       multcorr,sig){
  
  # number of otus
  n_otu=dim(otu_data)[2]-1
  # extract otu ids
  otu_ids=colnames(otu_data)[-1]
  
  if(repeated==F){
    data_comp=data.frame(merge(otu_data,var_data,by="Sample.ID",all.y=T),row.names=NULL)
    # colnames(data_comp)[2:(length(otu_ids)+1)] = otu_ids    # EDIT!!!!!! Need to reformat names ###MYC
    # The above edit is necessary when colname characters are not accepted by data.frame; R automatically turns certain characters into
    # "."; so we need to make sure colnames are identical between otu_ids and data_comp. Currently I have commented it out because one should really
    # practice good naming technique... but alas. ###MYC
  }else if(repeated==T){
    data_comp=data.frame(merge(otu_data,var_data,by="Sample.ID"),row.names=NULL)
    # data_comp=data.frame(merge(otu_data,var_data[,c("Sample.ID",main.var,repeat.var)],by="Sample.ID"),row.names=NULL)
  }
  
  base.formula = paste0("lr ~ ",main.var)
  if(repeated==T){
    repeat.formula = paste0(base.formula," | ", repeat.var)
  }
  if(adjusted==T){
    adjusted.formula = paste0(base.formula," + ", adj.formula)
  }
  
  if( adjusted == F & repeated == F ){
    fformula  <- formula(base.formula)
  } else if( adjusted == F & repeated == T & long == T ){
    fformula  <- formula(base.formula)   
  }else if( adjusted == F & repeated == T & long == F ){
    fformula  <- formula(repeat.formula)   
  }else if( adjusted == T & repeated == F  ){
    fformula  <- formula(adjusted.formula)   
  }else if( adjusted == T & repeated == T  ){
    fformula  <- formula(adjusted.formula)   
  }else{
    stop("Problem with data. Dataset should contain OTU abundances, groups, 
         and optionally an ID for repeated measures.")
  }
  
  
  
  if( repeated==FALSE & adjusted == FALSE){
    if( length(unique(data_comp[,which(colnames(data_comp)==main.var)]))==2 ){
      tfun <- exactRankTests::wilcox.exact
    } else{
      tfun <- stats::kruskal.test
    }
  }else if( repeated==FALSE & adjusted == TRUE){
    tfun <- stats::aov
  }else if( repeated== TRUE & adjusted == FALSE & long == FALSE){
    tfun <- stats::friedman.test
  }else if( repeated== TRUE & adjusted == FALSE & long == TRUE){
    tfun <- nlme::lme
  }else if( repeated== TRUE & adjusted == TRUE){
    tfun <- nlme::lme
  }
  
  logratio.mat <- matrix(NA, nrow=n_otu, ncol=n_otu)
  # For progress bar
  print("Beginning logratio calculations...") ### MYC-- added a progress bar so we can see progress in large datasets
  nreps <- (n_otu)*(n_otu-1)/2 ### MYC
  # create progress bar
  pb <- txtProgressBar(min = 0, max = nreps, style = 3)### MYC
  current_n <- 1 ### MYC
  for(ii in 1:(n_otu-1)){
    for(jj in (ii+1):n_otu){
      data.pair <- data_comp[,which(colnames(data_comp)%in%otu_ids[c(ii,jj)])]
      lr <- log((1+as.numeric(data.pair[,1]))/(1+as.numeric(data.pair[,2])))
      lr_dat <- data.frame( lr=lr, data_comp,row.names=NULL )
      
      if(adjusted==FALSE&repeated==FALSE){  ## Wilcox, Kruskal Wallis
        logratio.mat[ii,jj] <- tfun( formula=fformula, data = lr_dat)$p.value
      }else if(adjusted==FALSE&repeated==TRUE&long==FALSE){ ## Friedman's 
        logratio.mat[ii,jj] <- tfun( formula=fformula, data = lr_dat)$p.value
      }else if(adjusted==TRUE&repeated==FALSE){ ## ANOVA
        model=tfun(formula=fformula, data = lr_dat,na.action=na.omit)   
        picker=which(gsub(" ","",row.names(summary(model)[[1]]))==main.var)  
        logratio.mat[ii,jj] <- summary(model)[[1]][["Pr(>F)"]][picker]
      }else if(repeated==TRUE&long==TRUE){ ## GEE
        model=tfun(fixed=fformula,data = lr_dat,
                   random = formula(rand.formula),
                   correlation=corAR1(),
                   na.action=na.omit)   
        picker=which(gsub(" ","",row.names(anova(model)))==main.var)
        logratio.mat[ii,jj] <- anova(model)[["p-value"]][picker]
      }
      # MYC note::
      # I also added a progress bar because not knowing how much longer it will take is INFURIATING!
      current_n <- current_n +1
      setTxtProgressBar(pb, current_n) ### MYC
    }
  } 
  close(pb) ### MYC
  ind <- lower.tri(logratio.mat)
  logratio.mat[ind] <- t(logratio.mat)[ind]
  
  
  logratio.mat[which(is.finite(logratio.mat)==FALSE)] <- 1
  
  mc.pval <- t(apply(logratio.mat,1,function(x){
    s <- p.adjust(x, method = "BH")
    return(s)
  }))
  
  a <- logratio.mat[upper.tri(logratio.mat,diag=FALSE)==TRUE]
  
  b <- matrix(0,ncol=n_otu,nrow=n_otu)
  b[upper.tri(b)==T] <- p.adjust(a, method = "BH")
  diag(b)  <- NA
  ind.1    <- lower.tri(b)
  b[ind.1] <- t(b)[ind.1]
  
  
  #########################################
  ### Code to extract surrogate p-value
  surr.pval <- apply(mc.pval,1,function(x){
    s0=quantile(x[which(as.numeric(as.character(x))<sig)],0.95)
    # s0=max(x[which(as.numeric(as.character(x))<alpha)])
    return(s0)
  })
  #########################################
  ### Conservative
  if(multcorr==1){
    W <- apply(b,1,function(x){
      subp <- length(which(x<sig))
    })
    ### Moderate
  } else if(multcorr==2){
    W <- apply(mc.pval,1,function(x){
      subp <- length(which(x<sig))
    })
    ### No correction
  } else if(multcorr==3){
    W <- apply(logratio.mat,1,function(x){
      subp <- length(which(x<sig))
    })
  }
  
  # ###MYC note::
  # For the section below, I calculate the centered log ratio of OTU table and calculate F-values. 
  # This is for comparison of log ratios between the treatment groups. This is necesary
  # to plot volcano plots in the output. 
  #### BEGIN MYC ADDITION HERE ####
  # Calculate CLR (Centered Log Ratio) of OTU table
  OGmat <- t(otu_data[,-1])
  colnames(OGmat) <- otu_data[,1]
  
  class(OGmat) <- "numeric"
  OGmat_1 <- OGmat+1
  dim.mat <- dim(OGmat_1)
  
  clr.myc <- function(mat, dim) {
    mat=OGmat_1
    dim=dim.mat
    tempMat <- matrix(ncol=ncol(mat), nrow=ncol(mat), data = -1)
    diag(tempMat) <- ncol(mat)-1
    tempMat2 <- log(mat)/ncol(mat)
    finalMat <- tempMat2 %*% tempMat
    dim(finalMat)
    return(finalMat)
  }
  # clr.myc <- function(mat) {
  #   mat=OGmat_1
  #   dim=dim.mat
  #   g <- function(x) {
  #     prod(x)^(1/length(x))
  #   }
  #   
  #   finalMat <- apply(mat, MARGIN=1, FUN=function(x) {
  #     log(x)/g(x)
  #   })
  #   
  #   return(finalMat)
  # }
  
  otu_clr <- clr.myc(OGmat_1)
  colnames(otu_clr) <- colnames(OGmat_1)
  otu_clr_t <- t(otu_clr)
  
  ### NOTE: I have confirmed that these values are the SAME as the clr function in the compositions package
  
  #####################
  
  # get mean abundance in each group
  otu_data_relAbund <- t(OGmat)/rowSums(t(OGmat))
  # order mapping file by samples
  data_comp_sorted <- data_comp[match(rownames(otu_data_relAbund),data_comp[["Sample.ID"]]),]
  
  aveRelAbund <- matrix(ncol=length(unique(data_comp_sorted[[main.var]])), nrow=length(otu_ids), dimnames = list(otu_ids, paste0(levels(factor(data_comp[[main.var]])), "_relAbund"))) 
  f.stat <- matrix(ncol=1, nrow=length(otu_ids), dimnames = list(otu_ids, "W_f")) 
  for ( i in otu_ids) {
    aveRelAbund[i,] <- aggregate(otu_data_relAbund[,i], by=list(data_comp_sorted[[main.var]]), FUN=median)[,"x"]
    # clr (center log ratio transformation) anova
    f.stat[i,] <- anova(lm(otu_clr_t[,i] ~ as.vector(data_comp_sorted[[main.var]])))[["F value"]][1]
  }
  #### END MYC ADDITION HERE ####
  
  return(list(W
              , f.stat
              , aveRelAbund))
}



ANCOM.main.myc = function(OTUdat,Vardat,
                          adjusted,repeated,
                          main.var,adj.formula,
                          repeat.var,longitudinal,
                          random.formula,
                          multcorr,sig,
                          prev.cut){
  print("Removing high zero samples...") ### MYC
  p.zeroes=apply(OTUdat[,-1],2,function(x){
    s=length(which(x==0))/length(x)
  })
  
  zeroes.dist=data.frame(colnames(OTUdat)[-1],p.zeroes,row.names=NULL)
  colnames(zeroes.dist)=c("Taxon","Proportion_zero")
  
  zero.plot = ggplot(zeroes.dist, aes(x=Proportion_zero)) + 
    geom_histogram(binwidth=0.1,colour="black",fill="white") + 
    xlab("Proportion of zeroes") + ylab("Number of taxa") +
    theme_bw() +
    geom_vline(aes(xintercept=prev.cut), col="red", lty=2)
  
  #print(zero.plot)
  
  OTUdat.thinned=OTUdat
  OTUdat.thinned=OTUdat.thinned[,c(1,1+which(p.zeroes<prev.cut))]
  
  otu.names=colnames(OTUdat.thinned)[-1]
  
  W.detected   <- ancom.W.myc(OTUdat.thinned,Vardat,
                              adjusted,repeated,
                              main.var,adj.formula,
                              repeat.var,longitudinal,random.formula,
                              multcorr,sig)
  
  W_stat       <- W.detected[[1]] ### MYC
  W_f <- W.detected[[2]] ### MYC
  W_ave <- W.detected[[3]] ### MYC
  
  ### Bubble plot
  print("Calculating critical statistic values...") ### MYC
  W_frame = data.frame(otu.names,W_stat,W_f,W_ave,row.names=NULL)
  W_frame = W_frame[order(-W_frame$W_stat, -W_frame$W_f),]
  
  ### MYC: plot a CDF for w and statistic to estimate w_crit and stat_crit
  ### MYC note::
  # The section below takes the description of the original ANCOM in qiime and translates it into R
  # (without looking at python source code; this is based on description in supplementary data of paper)
  # What it does is calculates the CDF (Cumulative Distribution Function) of the value "W" and the F-score
  # of the ANOVA ran above. Then, it iterates through different combinations of tau ("drops" in the CDF) to find
  # the point where there is a plateau (restricted to the inner 80th percentile of the distribution). The text
  # describes a theoretical bimodel CDF due to the two groups of "OTU changes between treatments" and "OTU doesn't change
  # between treatments"-- thus, the mid-point plateau is the theoretical barrier between those two groups.
  # This part of the script also outputs the CDF plot for manual inspection.
  # The reality is, I have no idea if this is how this calculation is done in the python version of ANCOM. So please be aware that
  # the results here WILL DIFFER from the original ANCOM code.
  
  ## Note that I make a small addendum to this: instead of setting the bins at a strict 0.05 bin, I actually
  # shift it over a little each iteration in case the bins just HAPPEN not to line up where you want them to.
  
  ### new w_crit 
  ## w_crit
  cdf.w <- ecdf(W_frame$W_stat) # cdf
  bins <- seq(min(W_frame$W_stat),max(W_frame$W_stat), length.out=100)
  iter_bins <- seq(0,2)
  w_crit_mat <- matrix(ncol=5, nrow=length(iter_bins)*20*(100/3), dimnames = list(NULL, c("start","tau","d","dfw","cond_satisfied")))
  r <- 1
  for ( start in iter_bins ) {
    temp_bins <- bins[seq(100-start, 1, by = -3)]
    fw.cdf <- cdf.w(temp_bins) # get F(w)
    dfw.cdf <- -1*diff(fw.cdf) # get dF(w)
    
    tau_seq <- seq(0,max(dfw.cdf), length.out=20) # set drop threshold
    for ( tau in tau_seq ) {
      for ( d in 1:(length(dfw.cdf)-2) ) { # exclude last two
        cond_satisfied <- (dfw.cdf[d] >= tau) & (dfw.cdf[d+1] < tau) & (dfw.cdf[d+2] < tau) #& (bins[d]>quantile(W_frame$W_stat, probs=0.05)) & (bins[d]<quantile(W_frame$W_stat, probs=0.95))
        w_crit_mat[r, ] <- c(start, tau, d, dfw.cdf[d], cond_satisfied)
        r <- r+1
      }
    }
    
  }
  top_w_crit_mat <- w_crit_mat %>%
    as.data.frame() %>%
    filter(cond_satisfied==1) %>%
    arrange((tau)) 
  final_w<- top_w_crit_mat[1,c("start","d", "tau")]
  # Get W that we are looking for
  w_crit <- bins[seq(100-as.numeric(final_w[1]),1,by=-3)][as.numeric(final_w[2]+1)]
  tau_w <- as.numeric(final_w[3])
  
  w.pred <- seq(min(W_frame$W_stat),max(W_frame$W_stat), length.out=100)
  cdf.w.pred <- cdf.w(w.pred)
  plot.cdf.w <- data.frame(w.pred, cdf.w.pred) %>%
    ggplot() +
    geom_line(aes(x=w.pred, y=cdf.w.pred)) +
    xlab("W") +
    ylab("Cumulative Probability (F(W))")
  if (!is.na(w_crit)) {
    plot.cdf.w <- plot.cdf.w + geom_vline(aes(xintercept=w_crit), col="red", lty=2)
  }
  
  
  ####
  ### MYC notes:
  # In this part, I added wcrit (the percentile that is the "critical value" of W) and a bunch of cutoffs for effSize
  # so that you can easily bin "high W, high F-score" in the output.
  
  # Make cutoffs
  W_frame$detected_wcrit=rep(FALSE, dim(W_frame)[1])
  W_frame$detected_0.9=rep(FALSE,dim(W_frame)[1])
  W_frame$detected_0.8=rep(FALSE,dim(W_frame)[1])
  W_frame$detected_0.7=rep(FALSE,dim(W_frame)[1])
  W_frame$detected_0.6=rep(FALSE,dim(W_frame)[1])
  
  W_frame$effSize_statcrit=rep(FALSE,dim(W_frame)[1])
  W_frame$effSize_0.9=rep(FALSE,dim(W_frame)[1])
  W_frame$effSize_0.8=rep(FALSE,dim(W_frame)[1])
  W_frame$effSize_0.7=rep(FALSE,dim(W_frame)[1])
  W_frame$effSize_0.6=rep(FALSE,dim(W_frame)[1])
  
  if ( !is.na(w_crit)) {
    W_frame$detected_wcrit[which(W_frame$W_stat>w_crit)]=TRUE
  }
  W_frame$detected_0.9[which(W_frame$W_stat>0.9*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  W_frame$detected_0.8[which(W_frame$W_stat>0.8*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  W_frame$detected_0.7[which(W_frame$W_stat>0.7*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  W_frame$detected_0.6[which(W_frame$W_stat>0.6*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  
  W_frame$effSize_0.9[which(W_frame$W_f>quantile(W_frame$W_f,0.9))]=TRUE
  W_frame$effSize_0.8[which(W_frame$W_f>quantile(W_frame$W_f,0.8))]=TRUE
  W_frame$effSize_0.7[which(W_frame$W_f>quantile(W_frame$W_f,0.7))]=TRUE
  W_frame$effSize_0.6[which(W_frame$W_f>quantile(W_frame$W_f,0.6))]=TRUE
  
  ## Make volcano plot using both manual selection and automatic 0.9 ### MYC
  print("Generating plots...") ### MYC
  ### MYC note:
  # I plot 2 volcano plots; one using the script-generated wcrit and wstat (if any)
  # and other with a strict 0.9 cutoff. 
  ## Script-selected 
  volcano.plot <- ggplot(W_frame) +
    geom_point(aes(x=W_f, y=W_stat)) +
    xlab("CLR F statistic (Effect size)") +
    ylab("W statistic")
  
  volcano.plot.top <- volcano.plot
  
  if ( !is.na(w_crit) ) {
    onlyTop <- W_frame %>%
      filter(detected_wcrit==TRUE) 
    volcano.plot.top <- volcano.plot +
      geom_point(data=onlyTop, mapping=aes(x=W_f, y=W_stat, col=otu.names), show.legend = FALSE)
    if ( nrow(onlyTop)>0 ) {
      v.leg <- get_legend(volcano.plot +
                            geom_point(data=onlyTop, mapping=aes(x=W_f, y=W_stat, col=otu.names)))
    } else {
      v.leg <- NA
    }
    volcano.plot.top <- volcano.plot.top +
      geom_hline(mapping=aes(yintercept=w_crit), col="red", lty=2)
  }
  
  ### Automatic 0.9
  volcano.plot.auto <- volcano.plot
  volcano.plot.auto <- volcano.plot.auto +
    geom_vline(mapping=aes(xintercept=quantile(W_frame$W_f, probs = 0.9)), col="red", lty=2) +
    geom_hline(mapping=aes(yintercept=0.9*(dim(OTUdat.thinned[,-1])[2]-1)), col="red", lty=2)
  onlyautotop <- W_frame %>%
    filter(detected_0.9==TRUE, effSize_0.9==TRUE)
  if ( nrow(onlyautotop) >0) {
    volcano.plot.auto.top <- volcano.plot.auto +
      geom_point(data=onlyautotop, mapping=aes(x=W_f, y=W_stat, col=otu.names), show.legend = FALSE)
    v.leg.auto <- get_legend(volcano.plot.auto +
                               geom_point(data=onlyautotop, mapping=aes(x=W_f, y=W_stat, col=otu.names)))
  } else {
    volcano.plot.auto.top <- volcano.plot.auto
    v.leg.auto <- NA
  }
  
  final_results=list(W_frame,zero.plot,volcano.plot.top, volcano.plot.auto.top, plot.cdf.w,data.frame(w_crit, tau_w))
  names(final_results)=c("W.taxa","PLot.zeroes","Plot.volcano", "Plot.volcano.0.9", "Plot.critical.w", "Critical.values")
  print("Done")
  return(final_results)
}