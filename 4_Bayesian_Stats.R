#!bin/bash

##### Bayesian statistics ######
library(lme4) # for lmer
library(rstanarm) # For bayesian estimates of alpha and beta diversity
library(car) # for Anova
library(MASS) # for fitdistr
library(tidyverse) # for data manipulation
library(bayestestR)
library(rstan)
library(gridExtra)
library(grid)
library(ggpubr)
##### Load data ######
load("./3_5sp_mapping_otu_downstream/mf_alt_filt_final.RData")
load("./3_5sp_mapping_otu_downstream/otu_filt.RData")
load("./3_5sp_mapping_otu_downstream/otu_filt_inhibOnly.RData")
# load("./3_5sp_mapping_otu_downstream/braycurtis_filt.RData")
# load("./3_5sp_mapping_otu_downstream/unweighted_unifrac_filt.RData")
# load("./3_5sp_mapping_otu_downstream/weighted_unifrac_filt.RData")
# load("./3_5sp_mapping_otu_downstream/alpha_metrics.RData")
# load("./3_5sp_mapping_otu_downstream/beta_metrics.RData")

dir.create("4_Bayesian_Stats")
dir.create("4_Bayesian_Stats/fittedDistributions")
dir.create("4_Bayesian_Stats/BayesianSamplePosteriors")
dir.create("4_Bayesian_Stats/DataPlots")
dir.create("4_Bayesian_Stats/summarystats")
###### Set up mf ######
mf_con <- mf_alt_filt_final %>%
  filter(Bd_exposure == "Control")

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

allVar <- c("observed_otus","shannon"
            ,"inhibRich","percInhib"
            ,"dist_braycurtis","dist_unweighted_unifrac","dist_weighted_unifrac"
            ,"disper_braycurtis","disper_unweighted_unifrac","disper_weighted_unifrac")
allVarNames <- c("ASV richness", "Shannon diversity"
                 , "Richness of \nputative inhibitory bacteria", "Proportion \nputative inhibitory bacteria"
                 , "Instability through time \n(Bray-curtis)", "Instability through time \n(unweighted Unifrac)", "Instability through time \n(weighted Unifrac)"
                 , "Dispersion from centroid \n(Bray-curtis)", "Dispersion from centroid \n(unweighted Unifrac)", "Dispersion from centroid \n(weighted Unifrac)")
names(allVarNames) <- allVar

allVarFilt <- c("observed_otus","shannon","inhibRich","percInhib","dist_weighted_unifrac","disper_weighted_unifrac")

##### Fiting distributions to "uninfected" data ######
# Fit function for each set of values
for (var_temp in allVar) {
  val_temp <- pull(mf_con[,var_temp] %>% drop_na())
  x.fit.temp <- seq(min(val_temp)-sd(val_temp)
                    , max(val_temp)+sd(val_temp)
                    , length.out = 100)
  ## Fitting all types of functions
  # NORMAL
  param.norm <- fitdistr(val_temp, densfun="normal")
  y.pred.norm <- dnorm(x.fit.temp, mean = param.norm$estimate[1], sd = param.norm$estimate[2])
  # LOGNORMAL
  param.lnorm <- fitdistr(val_temp, densfun="lognormal")
  y.pred.lnorm <- dlnorm(x.fit.temp, meanlog=param.lnorm$estimate[1], sdlog = param.lnorm$estimate[2])
  # GAMMA AND BETA
  if (min(val_temp)>0) {
    param.gamma <- fitdistr(val_temp, densfun="gamma")
    y.pred.gamma <- dgamma(x.fit.temp, shape=param.gamma$estimate[1], rate = param.gamma$estimate[2])
    if (max(val_temp)<1) {
      param.beta <- fitdistr(val_temp, densfun="beta", start=list(shape1=6,shape2=6))
      y.pred.beta <- dbeta(x.fit.temp, shape1=param.beta$estimate[1], shape2 = param.beta$estimate[2])
    } else {
      # beta if can't fit
      y.pred.beta <- NA
    }
  } else {
    # gamma if can't fit
    y.pred.gamma <- NA
  }
  # POISSON
  if ( !any(val_temp%%1!=0) ) {
    param.pois <- fitdistr(val_temp, densfun="poisson")
    y.pred.pois <- dpois(round(x.fit.temp), lambda=param.pois$estimate[1])
  } else {
    y.pred.pois <- NA
  }
  # Combine all fit models
  allFits <- data.frame(x = x.fit.temp, Fit_Normal=y.pred.norm, Fit_LogNormal=y.pred.lnorm, Fit_Gamma = y.pred.gamma, Fit_Beta =y.pred.beta, Fit_Poisson=y.pred.pois) %>%
    pivot_longer(starts_with("Fit_"), names_to = "Model", values_to = "y.fit")
  
  ggsave(filename=paste0("4_Bayesian_Stats/fittedDistributions/fitDistr_",var_temp,".png"), height=3, width=4
         ,ggplot(data=data.frame(val_temp) %>% rename_at(vars(val_temp),~paste0(var_temp)), aes(x=get(var_temp))) +
           geom_histogram(aes(y=..density..), bins=30) +
           geom_line(data=allFits, aes(x=x, y=y.fit, col=Model))+
           xlab(paste0(allVarNames[var_temp])) + ylab("Density"))
  
}

###### CHOSEN MODELS
chosen_models <- c(observed_otus = "lognormal", chao1 = "lognormal", shannon = "normal"
                   , disper_braycurtis = "lognormal", disper_unweighted_unifrac = "lognormal", disper_weighted_unifrac = "lognormal"
                   , dist_braycurtis = "beta", dist_unweighted_unifrac = "beta", dist_weighted_unifrac = "beta"
                   , percInhib = "beta", inhibRich = "poisson")


##### @~~~~~~~~~~~~~~~Bayesian models- Control Only~~~~~~~~~~~~~~~@ #####

# mf_nobd_centred <- mf_nonbd %>% mutate(time_centred = time-mean(time, na.rm=TRUE))
set.seed(437745)
seeds <- runif(length(chosen_models), min=0, max = 10000)
if ( !file.exists("4_Bayesian_Stats/all_stan_output_control.RData") ) {
  all_stan_output_control <- list()
  for ( var_temp in allVar) {
    # var_temp <- "disper_weighted_unifrac"
    ind <- match(var_temp, allVar)
    # fit model
    model_temp <- chosen_models[var_temp]
    dat_temp <- mf_con %>% rename_at(vars(paste0(var_temp)), ~c("response"))
    if ( model_temp == "normal" ) {
      stan_temp <- stan_lmer(response ~ species*time + (1|indivID), data=dat_temp
                             , prior_intercept = normal(location = 0,scale = 5, autoscale = TRUE)
                             , prior = normal(location=0, scale=5, autoscale=TRUE)
                             , seed = seeds[ind]
                             , adapt_delta = 0.999
                             , iter = 5000
      )
    } else if (model_temp == "lognormal" ) {
      stan_temp <- stan_lmer(log(response) ~ species*time + (1|indivID), data=dat_temp
                             , prior_intercept = normal(location = 0,scale = 5, autoscale = TRUE)
                             , prior = normal(location=0, scale=5, autoscale=TRUE)
                             , seed = seeds[ind]
                             , adapt_delta = 0.999
                             , iter = 5000
      )
      # } else if (model_temp == "gamma" ) {
      #   stan_temp <- stan_glmer(response ~ species*time + (1|indivID), data=dat_temp
      #                           , prior = normal(0, 5, autoscale = TRUE)
      #                           , family = Gamma(link="identity")
      #                           , seed = seeds[ind]
      #                           , adapt_delta = 0.999
      #                           , iter = 5000
      #   )
    } else if (model_temp == "beta") {
      stan_temp <- stan_glmer(response ~ species*time + (1|indivID) , data=dat_temp
                              , family =mgcv::betar
                              , prior_intercept = normal(location = 0.5,scale = 2.5, autoscale = TRUE)
                              , prior = normal(location=0.5, scale=2.5, autoscale=TRUE)
                              , seed= seeds[ind]
                              , adapt_delta = 0.999
                              , iter = 5000
      )
    } else if (model_temp == "poisson") {
      # SampleID random effect here for overdispersion
      stan_temp <- stan_glmer(response ~ species*time + (1|indivID) + (1|SampleID), data=dat_temp
                              , prior = normal(0, 5, autoscale = TRUE)
                              , family= poisson(link="identity")
                              , seed = seeds[ind]
                              , adapt_delta = 0.999
                              , iter = 5000
      )
    }
    all_stan_output_control[[var_temp]] <- stan_temp
  }
  save(all_stan_output_control,file="4_Bayesian_Stats/all_stan_output_control.RData")
} else {
  load("4_Bayesian_Stats/all_stan_output_control.RData")
}


#### Get samps ####
main_samps_control <- data.frame()
random_samps_control <- data.frame()
if (!file.exists("4_Bayesian_Stats/main_samps_control.RData")) {
  for ( var_temp in names(all_stan_output_control)) {
    # var_temp = "inhibRich"
    stan_temp <- all_stan_output_control[[var_temp]]
    samps_temp <- rstan::extract(stan_temp$stanfit)
    # Rename
    allNames_temp <- names(stan_temp$coefficients)
    main_eff_names <- allNames_temp[-grep("Intercept|b\\[", allNames_temp)]
    main_eff_names <- main_eff_names[grep("species|time",main_eff_names)]
    # Get intercept species
    allSp <- unique(mf_con$species)
    intercept_name <- paste0("species",allSp[which(is.na(match(paste0("species",allSp), main_eff_names)))])
    # Get random effects
    b_eff_names <- allNames_temp[grep("b\\[", allNames_temp)]
    b_eff_sampleid_names <- gsub("^.* |\\]", "", b_eff_names[grep("SampleID", b_eff_names)])
    b_eff_individ_names <- gsub("^.* |\\]", "", b_eff_names[grep("indivID", b_eff_names)])
    
    # Main sample distributions
    long_samps_temp <- data.frame(intercept = samps_temp$alpha) %>% rename_at(vars("intercept"), ~intercept_name) %>%
      cbind(samps_temp$beta) %>% rename_at(vars(paste0(seq(1:ncol(samps_temp$beta)))), ~main_eff_names) %>%
      rownames_to_column(var = "iter") %>% rowwise() %>%
      mutate(Anbo = speciesAnbo, Rhma = speciesAnbo + speciesRhma, Osse = speciesAnbo + speciesOsse, Raca = speciesAnbo + speciesRaca, Rapi = speciesAnbo + speciesRapi) %>%
      mutate(timeAnbo = time, timeRhma = time + `speciesRhma:time`, timeOsse = time + `speciesOsse:time`, timeRaca = time + `speciesRaca:time`, timeRapi = time + `speciesRapi:time`) %>%
      mutate(timeMean = mean(c(timeAnbo, timeRhma, timeOsse, timeRaca, timeRapi)), Mean = mean(c(Anbo, Rhma, Osse, Raca, Rapi))) %>%
      ungroup() %>% select(iter, Anbo, Rhma, Osse, Raca, Rapi, timeAnbo, timeRhma, timeOsse, timeRaca, timeRapi, timeMean, Mean) %>%
      pivot_longer(-iter, names_to = "Predictor", values_to = "samp") %>% rowwise() %>%
      mutate(Response = var_temp, Type = ifelse(length(grep("time",Predictor))>0, "Time", "SpeciesAverage")) %>%
      mutate(Species = gsub("time","",Predictor)) %>%
      ungroup() %>% mutate(Species = factor(Species, levels = c("Anbo","Rhma","Osse","Raca","Rapi","Mean")))
    # Random effects distributions 
    if ( length(b_eff_sampleid_names)==0) {
      long_ranef_samps_temp <- data.frame(samps_temp$b) %>% 
        rename_at(vars(paste0("X",seq(1,ncol(samps_temp$b)))), ~c(b_eff_individ_names,"variation")) %>%
        rownames_to_column(var = "iter") %>%
        pivot_longer(-iter, names_to = "RandomEffect", values_to = "samp") %>%
        separate(RandomEffect, into = c("Type", "Predictor"), sep = ":") %>%
        mutate(Response = var_temp, Species = gsub("_.*$","",Predictor))
    } else {
      long_ranef_samps_temp <- data.frame(samps_temp$b) %>% 
        rename_at(vars(paste0("X",seq(1,ncol(samps_temp$b)))), ~c(b_eff_sampleid_names,"variation1",b_eff_individ_names,"variation2")) %>%
        rownames_to_column(var = "iter") %>%
        pivot_longer(-iter, names_to = "RandomEffect", values_to = "samp") %>%
        separate(RandomEffect, into = c("Type", "Predictor"), sep = ":") %>%
        mutate(Response = var_temp, Species = gsub("_.*$","",Predictor))
    }
    
    ggsave(filename = paste0("4_Bayesian_Stats/BayesianSamplePosteriors/",var_temp,"_MainEffects_posteriors.pdf"), height=8, width=6
           ,ggplot(long_samps_temp) +
             geom_histogram(aes(x=samp)) + geom_vline(aes(xintercept=0), col="red")+
             facet_grid(Species ~ Type, scales="free"))
    
    main_samps_control <- rbind(main_samps_control, long_samps_temp)
    random_samps_control <- rbind(random_samps_control, long_ranef_samps_temp)
  } 
  save(main_samps_control, file="4_Bayesian_Stats/main_samps_control.RData")
  save(random_samps_control, file="4_Bayesian_Stats/random_samps_control.RData")
} else {
  load("4_Bayesian_Stats/main_samps_control.RData")
  load("4_Bayesian_Stats/random_samps_control.RData")
}
#### Summarize main samps ####
summary_main_effects_control <- main_samps_control %>%
  group_by(Response, Type, Predictor, Species) %>%
  summarize(median = median(samp), lower95 = ci(samp, method="HDI")$CI_low, upper95 = ci(samp, method="HDI")$CI_high, postprob = max(c(sum(samp>0)/n(), 1-sum(samp>0)/n())), comp_postprob = 1-postprob) %>%
  ungroup() %>% arrange(rev(Type), Response, Predictor, Species)

summary_randeff_control <- random_samps_control %>%
  filter(Type == "indivID") %>%
  group_by(Response, Predictor, Species) %>%
  summarize(median = median(samp), lower95 = ci(samp, method="HDI")$CI_low, upper95 = ci(samp, method="HDI")$CI_high) %>%
  ungroup()

#### Generate posterior####
# mockData_treat <- mf_treat %>% select(SampleID, species, time, indivID) %>% distinct()
mockData_new <- mf_alt_filt_final %>% select(species, time) %>% distinct() %>%
  unite(species, time, col="SampleID", remove=FALSE) %>% mutate(indivID = paste0(species, "_new"))
if (!file.exists("4_Bayesian_Stats/allPosteriorPredictions.RData")) {
  print("Generating posteriors from control model")
  allPosteriorPredictions <- data.frame()
  total <- length(names(all_stan_output_control))
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  i <- 1
  for ( var_temp in names(all_stan_output_control)) {
    # var_temp = "disper_weighted_unifrac"
    # var_temp = "observed_otus"
    temp_model <- all_stan_output_control[[var_temp]]
    model_temp <- chosen_models[var_temp]
    ### Get model
    if (model_temp == "lognormal") {
      trans_func <- log
      inv_trans_func <- exp
      alt_samps <- function(x) {return(x)} 
      inv_alt_samps <- function(x) {return(x)}
    } else {
      trans_func <- function(x) {return(x)} ### Transform raw values (e.g. log)
      inv_trans_func <- function(x) {return(x)}
      alt_samps <- function(x) {return(x)} ### Link function
      inv_alt_samps <- function(x) {return(x)} ### Inverse link function
    } 
    post_draws_new <- as.data.frame(posterior_predict(temp_model, newdata=mockData_new)) %>%t()%>% as.data.frame() %>% 
      cbind(mockData_new) %>% 
      pivot_longer(-c(SampleID, species, time, indivID), names_to = "iter", values_to="samp") %>%
      mutate(samp_alt = inv_trans_func(samp)) %>%
      mutate(Response=var_temp)
    
    mf_temp <- mf_alt_filt_final %>% select(SampleID, species, indivID, time, Bd_exposure, prepost,one_of(var_temp)) %>%
      left_join(mf_collapsed_bd) %>%
      mutate(maxBd_outcome = ifelse(is.na(maxBd_outcome), 0, maxBd_outcome)
             , PABD = ifelse(is.na(PABD), 2, PABD)) %>%
      mutate(PABD = ifelse(PABD==1, "Yes", ifelse(PABD==2, "Control","No"))) %>%
      arrange(maxBd_outcome) %>% mutate(indivID = factor(indivID, levels=unique(indivID)))
    gg_throughtime <- post_draws_new %>% mutate(Bd_exposure = "Control") %>% full_join(post_draws_new %>% mutate(Bd_exposure = "Bd-exposed")) %>% 
      ggplot() + geom_violin(aes(x=factor(time), y=samp_alt)) +
      geom_jitter(data=mf_temp %>% filter(Bd_exposure == "Control"), aes(x=factor(time), y=get(var_temp), pch= factor(PABD)), col = "grey", height=0, width=0.2) +
      geom_jitter(data=mf_temp %>% filter(Bd_exposure=="Bd-exposed"), aes(x=factor(time), y=get(var_temp), col=maxBd_outcome, pch = factor(PABD)), height=0, width=0.2) +
      facet_grid(species~., scales="free") +
      scale_color_gradient(low="black", high="red") + 
      scale_shape_manual(values=c(8 ,21, 19)) + 
      labs(col = "Maximum Bd load", pch = "Did individual\nbecome infected?", title = paste0("Bayes posterior distribution fits of ", var_temp)) +
      ylab("Posterior distribution") + xlab("Time (Sampling point)") 
    ggsave(filename = paste0("4_Bayesian_Stats/DataPlots/throughTime_",var_temp,".png"), height=8, width=8
           ,gg_throughtime
    )
    gg_throughtime_linepre <- post_draws_new %>% mutate(Bd_exposure = "Control") %>% full_join(post_draws_new %>% mutate(Bd_exposure = "Bd-exposed")) %>% 
      ggplot() + geom_violin(aes(x=factor(time), y=samp_alt)) +
      geom_line(data=mf_temp %>% filter(Bd_exposure=="Control"), aes(x=factor(time), y=get(var_temp), group=indivID), col="grey", lty=1) +
      geom_point(data=mf_temp %>% filter(Bd_exposure=="Bd-exposed", prepost=="Pre"), aes(x=factor(time), y=get(var_temp), group=indivID, col=maxBd_outcome, pch = factor(PABD)), position=position_dodge(width=0.8)) +
      geom_line(data=mf_temp %>% filter(Bd_exposure=="Bd-exposed", prepost=="Pre"), aes(x=factor(time), y=get(var_temp), group=indivID, col=maxBd_outcome, lty=factor(PABD)), position=position_dodge(width=0.8)) +
      facet_grid(species~., scales="free") +
      scale_color_gradient(low="black", high="red") + 
      scale_shape_manual(values=c(21, 19)) + 
      scale_linetype_manual(values=c(3,1)) +
      labs(col = "Maximum Bd load", pch = "Did individual\nbecome infected?",lty = "Did individual\nbecome infected?", title = paste0("Bayes posterior distribution fits of ", var_temp)) +
      ylab("Posterior distribution") + xlab("Time (Sampling point)") 
    
    ggsave(filename = paste0("4_Bayesian_Stats/DataPlots/throughTime_",var_temp,"_withlinepre.png"), height=8, width=8
           ,gg_throughtime_linepre
    )
    
    gg_throughtime_lineall <- post_draws_new %>% mutate(Bd_exposure = "Control") %>% full_join(post_draws_new %>% mutate(Bd_exposure = "Bd-exposed")) %>% 
      ggplot() + geom_violin(aes(x=factor(time), y=samp_alt)) +
      geom_line(data=mf_alt_filt_final %>% filter(Bd_exposure=="Control"), aes(x=factor(time), y=get(var_temp), group=indivID), col="grey", lty=1) +
      geom_line(data=mf_alt_filt_final %>% filter(Bd_exposure=="Bd-exposed"), aes(x=factor(time), y=get(var_temp), group=indivID), position=position_dodge(width=0.8)) +
      geom_point(data=mf_alt_filt_final %>% filter(Bd_exposure=="Bd-exposed"), aes(x=factor(time), y=get(var_temp), group=indivID, col=Bd_load), cex=2,position=position_dodge(width=0.8)) +
      facet_grid(species~., scales="free") +
      scale_color_gradient(low="grey", high="red") + 
      # scale_shape_manual(values=c(21, 19)) + 
      scale_linetype_manual(values=c(3,1)) +
      labs(col = "Maximum Bd load", pch = "Did individual\nbecome infected?",lty = "Did individual\nbecome infected?", title = paste0("Bayes posterior distribution fits of ", var_temp)) +
      ylab("Posterior distribution") + xlab("Time (Sampling point)") 
    ggsave(filename = paste0("4_Bayesian_Stats/DataPlots/throughTime_",var_temp,"_withlineall.png"), height=8, width=8
           ,gg_throughtime_lineall
    )
    
    
    gg_throughtime_lineall_spCollapse <- mf_alt_filt_final %>% filter(Bd_exposure=="Bd-exposed") %>%
      mutate(PABD = ifelse(PABD==1, "Yes", "No") ) %>% ggplot() + 
      # geom_violin(aes(x=factor(time), y=samp_alt, group=species)) +
      # geom_line(data=mf_alt_filt_final %>% filter(Bd_exposure=="Control"), aes(x=factor(time), y=get(var_temp), group=indivID), col="grey", lty=1) +
      geom_line(aes(x=(time), y=get(var_temp), group=indivID, col=species), position=position_dodge(width=0.8)) +
      geom_point(aes(x=(time), y=get(var_temp), group=indivID, fill=Bd_load, cex=factor(PABD)), position=position_dodge(width=0.8), pch=21) +
      # facet_grid(species~., scales="free") +
      # scale_color_manual(values=c("red","green","blue","orange","purple"))+
      scale_fill_gradient(low="white", high="darkred") + 
      scale_size_manual(values=c(2,5)) +
      # scale_shape_manual(values=c(21, 19)) + 
      # scale_linetype_manual(values=c(3,1)) +
      labs(col="Amphibian species", fill = "Maximum Bd load", cex = "Infection detected") +
      ylab(paste0(allVarNames[var_temp])) + xlab("Time (Sampling point)") 

    ggsave(filename = paste0("4_Bayesian_Stats/DataPlots/throughTime_",var_temp,"_withlineall_spCollapse.png"), height=5, width=8
           ,gg_throughtime_lineall_spCollapse
    )
    
    gg_throughtime_lineall_sp <- mf_alt_filt_final %>% filter(Bd_exposure=="Bd-exposed") %>%
      mutate(PABD = ifelse(PABD==1, "Yes", "No") ) %>% ggplot() + 
      # geom_violin(aes(x=factor(time), y=samp_alt, group=species)) +
      # geom_line(data=mf_alt_filt_final %>% filter(Bd_exposure=="Control"), aes(x=factor(time), y=get(var_temp), group=indivID), col="grey", lty=1) +
      geom_line(aes(x=(time), y=get(var_temp), group=indivID, col=species), position=position_dodge(width=0.8)) +
      geom_point(aes(x=(time), y=get(var_temp), group=indivID, fill=Bd_load, cex=factor(PABD)), position=position_dodge(width=0.8), pch=21) +
      facet_grid(species~., scales="free") +
      # scale_color_manual(values=c("red","green","blue","orange","purple"))+
      scale_fill_gradient(low="white", high="darkred") + 
      scale_size_manual(values=c(2,5)) +
      # scale_shape_manual(values=c(21, 19)) + 
      # scale_linetype_manual(values=c(3,1)) +
      labs(col="Amphibian species", fill = "Maximum Bd load", cex = "Infection detected") +
      ylab(paste0(allVarNames[var_temp])) + xlab("Time (Sampling point)") 
    gg_throughtime_lineall_sp
    
    ggsave(filename = paste0("4_Bayesian_Stats/DataPlots/throughTime_",var_temp,"_withlineall_sp.png"), height=8, width=8
           ,gg_throughtime_lineall_sp
    )
    
    gg_raw_PABD <- mf_temp %>% filter(prepost=="Pre", Bd_exposure == "Bd-exposed") %>% 
      mutate(PABD = factor(PABD, levels=c("No","Yes"))) %>%
      ggplot() + 
      geom_boxplot(aes(x=PABD, y=get(var_temp), group=indivID, col=species), position = position_dodge(width=0.4) ) +
      geom_point(aes(x=PABD, y=get(var_temp), group=indivID, col=species, pch=factor(time)), position=position_dodge(width=0.4)) +
      coord_flip() + xlab("Did individuals become infected?") + ylab(paste0(allVarNames[var_temp])) +
      labs(col="Species", pch="Time point")
    ggsave(filename = paste0("4_Bayesian_Stats/DataPlots/rawData_",var_temp,"_PABD.png"), height=4, width=6
           ,gg_raw_PABD
    )
    gg_raw_bdload <- mf_temp %>% filter(prepost=="Pre", Bd_exposure == "Bd-exposed") %>% 
      filter(maxBd_outcome>0) %>%
      ggplot() + 
      geom_boxplot(aes(x=maxBd_outcome, y=get(var_temp), group=indivID, col=species), position = position_dodge(width=0.4) ) +
      geom_point(aes(x=maxBd_outcome, y=get(var_temp), group=indivID, col=species, pch=factor(time)), position=position_dodge(width=0.4)) +
      coord_flip() + xlab("Bd infection load (log)") + ylab(paste0(allVarNames[var_temp])) +
      labs(col="Species", pch="Time point")
    ggsave(filename = paste0("4_Bayesian_Stats/DataPlots/rawData_",var_temp,"_Bdload.png"), height=4, width=6
           ,gg_raw_bdload
    )
    # 
    # ggsave(filename = paste0("4_Bayesian_Stats/DataPlots/throughTime_",var_temp,"_withlineall_sp.png"), height=8, width=8
    #        ,gg_correlations
    # )
    # rank_temp <- mf_treat %>% select(SampleID, species, indivID, one_of(var_temp))%>% left_join(post_draws_treat) %>%
    #   filter(!is.na(get(var_temp))) %>%
    #   # mutate(var_trans = inv_alt_samps(trans_func(get(var_temp)))) %>%
    #   mutate(cump = get(var_temp)>samp_alt) %>% group_by(species, indivID, time, Response) %>% 
    #   summarise(prank = sum(cump)/n()) %>%
    #   ungroup()  %>% #group_by(Response,species, indivID) %>%
    #   # summarise(medrank=median(prank, na.rm=TRUE)) %>%
    #   left_join(mf_collapsed_bd)
    
    #### Save aggregate data
    # allPosteriorPredictions <- rbind(allPosteriorPredictions,rbind(post_draws_con, post_draws_new))
    allPosteriorPredictions <- rbind(allPosteriorPredictions,post_draws_new)
    # allRanks <- rbind(allRanks, rank_temp)
    #------- PROGRESS BAR ------- #
    Sys.sleep(0.1)
    i <- i+1
    # update progress bar
    setTxtProgressBar(pb, i)
    #-----------------------------#
  }
  close(pb)
  save(allPosteriorPredictions, file="4_Bayesian_Stats/allPosteriorPredictions.RData")
  # save(allRanks, file="4_Bayesian_Stats/allRanks.RData")
} else {
  load("4_Bayesian_Stats/allPosteriorPredictions.RData")
  # load("4_Bayesian_Stats/allRanks.RData")
}

#### Get ranks of treatment from control model ####
if (!file.exists("4_Bayesian_Stats/allRanks.RData")) {
  print("Getting ranks for analysis and plotting")
  allRanks <- data.frame()
  for ( var_temp in allVar) {
    # var_temp = "percInhib"
   tempSamps <- allPosteriorPredictions %>% filter(Response==var_temp)
   mf_temp <- mf_alt_filt_final %>% select(SampleID, species, indivID, time, prepost, Bd_exposure, one_of(var_temp)) %>%
     rename_at(vars(all_of(var_temp)), ~"Observed")
   allPValues <- data.frame()

   print(var_temp)
   total <- nrow(mf_temp)
   pb <- txtProgressBar(min = 0, max = total, style = 3)
   i <- 1
   for ( r in 1:nrow(mf_temp)) {
     # r <- 1
     t_temp <- pull(mf_temp[r,"time"])
     sp_temp <- pull(mf_temp[r,"species"])
     # indiv_temp <- pull(mf_temp[r,"indivID"])
     val_temp <- pull(mf_temp[r,"Observed"])
     ## Get distr
     samp_temp_distr <- tempSamps %>% filter(time==t_temp, species==sp_temp) %>% pull(samp_alt)
     p_temp <- sum(val_temp>samp_temp_distr)/length(samp_temp_distr)
     # Save
     allPValues <- rbind(allPValues, data.frame(mf_temp[r,], Response = var_temp, prank = p_temp))
     
     #------- PROGRESS BAR ------- #
     Sys.sleep(0.1)
     i <- i+1
     # update progress bar
     setTxtProgressBar(pb, i)
     #-----------------------------#
     }
   close(pb)
   
    #### Save aggregate data
    allRanks <- rbind(allRanks, allPValues)

  }
  save(allRanks, file="4_Bayesian_Stats/allRanks.RData")
} else {
  load("4_Bayesian_Stats/allRanks.RData")
}


#### Summarize ranks, plot ####
rankSummary_nonInfected <- allRanks %>% filter(prepost=="Pre"| Bd_exposure=="Control") %>%
  group_by(species, indivID, Response) %>%
  summarize(average_rank = mean(prank, na.rm=TRUE), sd = sd(prank, na.rm=TRUE)) %>% ungroup()
allRanks_nonInfected <- allRanks %>% filter(prepost=="Pre"| Bd_exposure=="Control")
mockData_single <- mf_alt_filt_final %>% select(species) %>% distinct() %>%
  rownames_to_column(var="SampleID") %>% mutate(indivID = paste0(species, "_new"), time=0) 
allRealVal <- data.frame()
allRealVal_raw <- data.frame()
allSampPred <- data.frame()
total <- length(allVar)
pb <- txtProgressBar(min = 0, max = total, style = 3)
i <- 1
for ( var_temp in allVar) {
  # var_temp="inhibRich"
  model_temp <- all_stan_output_control[[var_temp]]
  
  samp_pred <- posterior_predict(model_temp, newdata = mockData_single) %>%
    as.data.frame() %>% 
    rownames_to_column(var="iter") %>%
    pivot_longer(-iter,names_to = "SampleID", values_to = "samp") %>%
    left_join(mockData_single)
  temp_ranks <- rankSummary_nonInfected %>% filter(Response==var_temp) 
  temp_ranks_raw <- allRanks_nonInfected %>% filter(Response==var_temp)
  #### For summaries
  realVal <- data.frame()
  for ( r in 1:nrow(temp_ranks)) {
    sp_temp <- pull(temp_ranks[r,"species"])
    indiv_temp <- pull(temp_ranks[r,"indivID"])
    val <- pull(temp_ranks[r,"average_rank"])
    samp_temp <- samp_pred%>% filter(species ==sp_temp) %>% pull(samp)
    est_val <- sort(samp_temp)[length(samp_temp)*val]
    
    realVal <- rbind(realVal, data.frame(species=sp_temp, indivID=indiv_temp, realVal = est_val))
  }
  realVal <- realVal %>% left_join(mf_collapsed_bd) %>%
    mutate(PABD = ifelse(is.na(PABD),"Control", ifelse(PABD==1, "Yes","No"))) %>%
    mutate(Response=var_temp)%>%
    arrange(-maxBd_outcome) %>% mutate(indivID = factor(indivID, levels=unique(indivID)))
  #### For raw values
  realVal_all <- data.frame()
  for ( r in 1:nrow(temp_ranks_raw)) {
    sp_temp <- as.character(temp_ranks_raw[r,"species"])
    indiv_temp <- as.character(temp_ranks_raw[r,"indivID"])
    t_temp <- as.numeric(temp_ranks_raw[r,"time"])
    val <- as.numeric(temp_ranks_raw[r,"prank"])
    samp_temp <- samp_pred%>% filter(species ==sp_temp) %>% pull(samp)
    est_val <- sort(samp_temp)[length(samp_temp)*val+1]
    
    realVal_all <- rbind(realVal_all, data.frame(species=sp_temp, indivID=indiv_temp, time = t_temp, realVal = est_val))
  }
  realVal_all <- realVal_all %>% left_join(mf_collapsed_bd) %>%
    mutate(PABD = ifelse(is.na(PABD),"Control", ifelse(PABD==1, "Yes","No"))) %>%
    mutate(Response=var_temp) %>%
    arrange(-maxBd_outcome) %>% mutate(indivID = factor(indivID, levels=unique(indivID)))
  
  #------- PROGRESS BAR ------- #
  Sys.sleep(0.1)
  i <- i+1
  # update progress bar
  setTxtProgressBar(pb, i)
  #-----------------------------#
  allRealVal <- rbind(allRealVal, realVal)
  allRealVal_raw <- rbind(allRealVal_raw, realVal_all)
  allSampPred <- rbind(allSampPred, data.frame(samp_pred, Response=var_temp))
}
close(pb)

##### Plotting #####

print("Generating posteriors from control model")
total <- length(allVar)
pb <- txtProgressBar(min = 0, max = total, style = 3)
i <- 1
for ( var_temp in allVar) {
  # var_temp = "disper_weighted_unifrac"
  print(var_temp)
  #### Get data
  post_draws_new <- allPosteriorPredictions %>% filter(Response==var_temp)
  allPValues <- allRanks %>% filter(Response==var_temp)
  mf_temp <- mf_alt_filt_final %>% select(SampleID, species, indivID, time, Bd_exposure, prepost,one_of(var_temp)) %>%
    left_join(mf_collapsed_bd) %>%
    mutate(maxBd_outcome = ifelse(is.na(maxBd_outcome), 0, maxBd_outcome)
           , PABD = ifelse(is.na(PABD), 2, PABD)) %>%
    mutate(PABD = ifelse(PABD==1, "Yes", ifelse(PABD==2, "Control","No"))) %>%
    arrange(maxBd_outcome) %>% mutate(indivID = factor(indivID, levels=unique(indivID)))
  samp_pred <- allSampPred %>% filter(Response==var_temp)
  realVal_all <- allRealVal_raw %>% filter(Response==var_temp)
  realVal <- allRealVal %>% filter(Response==var_temp)
  #### through time; violin JUST controls
  
  gg_throughtime_con <- post_draws_new %>% mutate(Bd_exposure = "Control") %>% full_join(post_draws_new %>% mutate(Bd_exposure = "Bd-exposed")) %>% 
    ggplot() + geom_violin(aes(x=factor(time), y=samp_alt)) +
    geom_point(data=mf_con %>% filter(Bd_exposure == "Control"), aes(x=factor(time), y=get(var_temp), col=species)) +
    geom_line(data=mf_con %>% filter(Bd_exposure == "Control"), aes(x=factor(time), y=get(var_temp), group=indivID, col=species)) +
    facet_grid(species~., scales="free") +
    labs(title = paste0("Bayes posterior distribution fits of ", var_temp), col="Species") +
    ylab("Predictive Posterior Distribution") + xlab("Time point") 
  gg_throughtime_con
  ggsave(filename = paste0("4_Bayesian_Stats/DataPlots/throughTime_",var_temp,"_controlsONly.png"), height=8, width=8
         ,gg_throughtime_con
  )
  #### through time; violin
  gg_throughtime <- post_draws_new %>% mutate(Bd_exposure = "Control") %>% full_join(post_draws_new %>% mutate(Bd_exposure = "Bd-exposed")) %>% 
    ggplot() + geom_violin(aes(x=factor(time), y=samp_alt)) +
    geom_jitter(data=mf_temp %>% filter(Bd_exposure == "Control"), aes(x=factor(time), y=get(var_temp), pch= factor(PABD)), col = "grey", height=0, width=0.2) +
    geom_jitter(data=mf_temp %>% filter(Bd_exposure=="Bd-exposed"), aes(x=factor(time), y=get(var_temp), col=maxBd_outcome, pch = factor(PABD)), height=0, width=0.2) +
    facet_grid(species~., scales="free") +
    scale_color_gradient(low="black", high="red") + 
    scale_shape_manual(values=c(8 ,21, 19)) + 
    labs(col = "Maximum Bd load", pch = "Did individual\nbecome infected?", title = paste0("Bayes posterior distribution fits of ", var_temp)) +
    ylab("Posterior distribution") + xlab("Time point)") 
  ggsave(filename = paste0("4_Bayesian_Stats/DataPlots/throughTime_",var_temp,".png"), height=8, width=8
         ,gg_throughtime
  )
  #### through time; violin with line, pre treatment only
  gg_throughtime_linepre <- post_draws_new %>% mutate(Bd_exposure = "Control") %>% full_join(post_draws_new %>% mutate(Bd_exposure = "Bd-exposed")) %>% 
    ggplot() + geom_violin(aes(x=factor(time), y=samp_alt)) +
    geom_line(data=mf_temp %>% filter(Bd_exposure=="Control"), aes(x=factor(time), y=get(var_temp), group=indivID), col="grey", lty=1) +
    geom_point(data=mf_temp %>% filter(Bd_exposure=="Bd-exposed", prepost=="Pre"), aes(x=factor(time), y=get(var_temp), group=indivID, col=maxBd_outcome, pch = factor(PABD)), position=position_dodge(width=0.8)) +
    geom_line(data=mf_temp %>% filter(Bd_exposure=="Bd-exposed", prepost=="Pre"), aes(x=factor(time), y=get(var_temp), group=indivID, col=maxBd_outcome, lty=factor(PABD)), position=position_dodge(width=0.8)) +
    facet_grid(species~., scales="free") +
    scale_color_gradient(low="black", high="red") + 
    scale_shape_manual(values=c(21, 19)) + 
    scale_linetype_manual(values=c(3,1)) +
    labs(col = "Maximum Bd load", pch = "Did individual\nbecome infected?",lty = "Did individual\nbecome infected?", title = paste0("Bayes posterior distribution fits of ", var_temp)) +
    ylab("Posterior distribution") + xlab("Time point") 
  
  ggsave(filename = paste0("4_Bayesian_Stats/DataPlots/throughTime_",var_temp,"_withlinepre.png"), height=8, width=8
         ,gg_throughtime_linepre
  )
  #### through time; violin with line, all samples
  gg_throughtime_lineall <- post_draws_new %>% mutate(Bd_exposure = "Control") %>% full_join(post_draws_new %>% mutate(Bd_exposure = "Bd-exposed")) %>% 
    ggplot() + geom_violin(aes(x=factor(time), y=samp_alt)) +
    geom_line(data=mf_alt_filt_final %>% filter(Bd_exposure=="Control"), aes(x=factor(time), y=get(var_temp), group=indivID), col="grey", lty=1) +
    geom_line(data=mf_alt_filt_final %>% filter(Bd_exposure=="Bd-exposed"), aes(x=factor(time), y=get(var_temp), group=indivID), position=position_dodge(width=0.8)) +
    geom_point(data=mf_alt_filt_final %>% filter(Bd_exposure=="Bd-exposed"), aes(x=factor(time), y=get(var_temp), group=indivID, col=Bd_load), cex=2,position=position_dodge(width=0.8)) +
    facet_grid(species~., scales="free") +
    scale_color_gradient(low="grey", high="red") + 
    # scale_shape_manual(values=c(21, 19)) + 
    scale_linetype_manual(values=c(3,1)) +
    labs(col = "Maximum Bd load", pch = "Did individual\nbecome infected?",lty = "Did individual\nbecome infected?", title = paste0("Bayes posterior distribution fits of ", var_temp)) +
    ylab("Posterior distribution") + xlab("Time point") 
  ggsave(filename = paste0("4_Bayesian_Stats/DataPlots/throughTime_",var_temp,"_withlineall.png"), height=8, width=8
         ,gg_throughtime_lineall
  )
  
  #### through time; line with all species
  gg_throughtime_lineall_spCollapse <- mf_alt_filt_final %>% filter(Bd_exposure=="Bd-exposed") %>%
    mutate(PABD = ifelse(PABD==1, "Yes", "No") ) %>% ggplot() + 
    # geom_violin(aes(x=factor(time), y=samp_alt, group=species)) +
    # geom_line(data=mf_alt_filt_final %>% filter(Bd_exposure=="Control"), aes(x=factor(time), y=get(var_temp), group=indivID), col="grey", lty=1) +
    geom_line(aes(x=(time), y=get(var_temp), group=indivID, col=species), position=position_dodge(width=0.8)) +
    geom_point(aes(x=(time), y=get(var_temp), group=indivID, fill=Bd_load, cex=factor(PABD)), position=position_dodge(width=0.8), pch=21) +
    # facet_grid(species~., scales="free") +
    # scale_color_manual(values=c("red","green","blue","orange","purple"))+
    scale_fill_gradient(low="white", high="darkred") + 
    scale_size_manual(values=c(2,5)) +
    # scale_shape_manual(values=c(21, 19)) + 
    # scale_linetype_manual(values=c(3,1)) +
    labs(col="Amphibian species", fill = "Maximum Bd load", cex = "Infection detected") +
    ylab(paste0(allVarNames[var_temp])) + xlab("Time point") 
  
  ggsave(filename = paste0("4_Bayesian_Stats/DataPlots/throughTime_",var_temp,"_withlineall_spCollapse.png"), height=5, width=8
         ,gg_throughtime_lineall_spCollapse
  )
  #### through time; line with soecues facet
  gg_throughtime_lineall_sp <- mf_alt_filt_final %>% filter(Bd_exposure=="Bd-exposed") %>%
    mutate(PABD = ifelse(PABD==1, "Yes", "No") ) %>% ggplot() + 
    # geom_violin(aes(x=factor(time), y=samp_alt, group=species)) +
    # geom_line(data=mf_alt_filt_final %>% filter(Bd_exposure=="Control"), aes(x=factor(time), y=get(var_temp), group=indivID), col="grey", lty=1) +
    geom_line(aes(x=(time), y=get(var_temp), group=indivID, col=species), position=position_dodge(width=0.8)) +
    geom_point(aes(x=(time), y=get(var_temp), group=indivID, fill=Bd_load, cex=factor(PABD)), position=position_dodge(width=0.8), pch=21) +
    facet_grid(species~., scales="free") +
    # scale_color_manual(values=c("red","green","blue","orange","purple"))+
    scale_fill_gradient(low="white", high="darkred") + 
    scale_size_manual(values=c(2,5)) +
    # scale_shape_manual(values=c(21, 19)) + 
    # scale_linetype_manual(values=c(3,1)) +
    labs(col="Amphibian species", fill = "Maximum Bd load", cex = "Infection detected") +
    ylab(paste0(allVarNames[var_temp])) + xlab("Time point") 
  gg_throughtime_lineall_sp
  
  ggsave(filename = paste0("4_Bayesian_Stats/DataPlots/throughTime_",var_temp,"_withlineall_sp.png"), height=8, width=8
         ,gg_throughtime_lineall_sp
  )
  
  #### raw values of treatment pre; PABD
  gg_raw_PABD <- mf_temp %>% filter(prepost=="Pre", Bd_exposure == "Bd-exposed") %>% 
    mutate(PABD = factor(PABD, levels=c("No","Yes"))) %>%
    ggplot() + 
    geom_boxplot(aes(x=PABD, y=get(var_temp), group=indivID, col=species), position = position_dodge(width=0.4) ) +
    geom_point(aes(x=PABD, y=get(var_temp), group=indivID, col=species, pch=factor(time)), position=position_dodge(width=0.4)) +
    coord_flip() + xlab("Did individuals become infected?") + ylab(paste0(allVarNames[var_temp])) +
    labs(col="Species", pch="Time point")
  ggsave(filename = paste0("4_Bayesian_Stats/DataPlots/rawData_",var_temp,"_PABD.png"), height=4, width=6
         ,gg_raw_PABD
  )
  
  #### raw values of treatment pre; Bd load
  gg_raw_bdload <- mf_temp %>% filter(prepost=="Pre", Bd_exposure == "Bd-exposed") %>% 
    filter(maxBd_outcome>0) %>%
    ggplot() + 
    geom_boxplot(aes(x=maxBd_outcome, y=get(var_temp), group=indivID, col=species), position = position_dodge(width=0.4) ) +
    geom_point(aes(x=maxBd_outcome, y=get(var_temp), group=indivID, col=species, pch=factor(time)), position=position_dodge(width=0.4)) +
    coord_flip() + xlab("Bd infection load (log)") + ylab(paste0(allVarNames[var_temp])) +
    labs(col="Species", pch="Time point")
  ggsave(filename = paste0("4_Bayesian_Stats/DataPlots/rawData_",var_temp,"_Bdload.png"), height=4, width=6
         ,gg_raw_bdload
  )
  ## Test out p over time
  gg_linethroughtime <- allPValues %>%
    left_join(mf_collapsed_bd) %>%
    filter(prepost=="Pre" | Bd_exposure == "Control") %>%
    ggplot() + geom_point(aes(x=time, y=prank, group=indivID, col=maxBd_outcome)) +
    geom_line(aes(x=time, y=prank, group=indivID, lty=Bd_exposure, col=maxBd_outcome)) +
    facet_grid(species~.) +
    scale_color_gradient(low="black", high="red") +
    labs(lty = "Treatment group", col = "Maximum future\nBd load", title=paste0(var_temp)) +
    ylab("Rank in cumulative distribution curve") + xlab("Time")
  ggsave(filename = paste0("4_Bayesian_Stats/DataPlots/Ranks_",var_temp,".png"), width=5, height=4
         ,gg_linethroughtime)
  # Get means
  allPMeans <- allPValues %>% filter(prepost=="Pre", Bd_exposure=="Bd-exposed")  %>%
    group_by(species, indivID) %>%
    summarize(prank = median(prank)) %>% 
    right_join(mf_collapsed_bd) %>%
    mutate(Type="Average")
  gg_rankCorr_pabd <- allPValues %>% filter(prepost=="Pre", Bd_exposure=="Bd-exposed") %>%
    right_join(mf_collapsed_bd) %>% 
    mutate(Type="Individual values") %>%
    full_join(allPMeans) %>%
    mutate(PABD = factor(ifelse(PABD==1, "Yes", "No"), levels=c("No","Yes"))) %>%
    ggplot() + 
    geom_line(aes(x=PABD, y=prank, col=species, group=indivID), position=position_dodge(width=0.4)) +
    geom_point(aes(x=PABD, y=prank, col=species, group=indivID, pch=Type, cex=Type),position=position_dodge(width=0.4)) +
    coord_flip() +
    scale_shape_manual(values=c(19,21))+
    scale_size_manual(values=c(4,1))+
    labs(col="Species", pch ="Datapoint", cex="Datapoint") + 
    xlab("Did individual become infected?") + ylab(paste0(allVarNames[var_temp], "\n(Percentile within species)"))
  ggsave(filename = paste0("4_Bayesian_Stats/DataPlots/RankCorr_",var_temp,"_pabd.png"), width=6.5, height=4
         ,gg_rankCorr_pabd)
  
  gg_rankCorr_bdload <- allPValues %>% filter(prepost=="Pre", Bd_exposure=="Bd-exposed") %>%
    right_join(mf_collapsed_bd) %>% 
    mutate(Type="Individual values") %>%
    full_join(allPMeans) %>%
    filter(maxBd_outcome>0) %>%
    # mutate(PABD = factor(ifelse(PABD==1, "Yes", "No"), levels=c("No","Yes"))) %>%
    ggplot() + 
    geom_line(aes(x=maxBd_outcome, y=prank, col=species, group=indivID), position=position_dodge(width=0.4)) +
    geom_point(aes(x=maxBd_outcome, y=prank, col=species, group=indivID, pch=Type, cex=Type),position=position_dodge(width=0.4)) +
    coord_flip() +
    scale_shape_manual(values=c(19,21))+
    scale_size_manual(values=c(4,1))+
    labs(col="Species", pch ="Datapoint", cex="Datapoint") + 
    xlab("Maximum future Bd load (log)") + ylab(paste0(allVarNames[var_temp], "\n(Percentile within species)"))
  ggsave(filename = paste0("4_Bayesian_Stats/DataPlots/RankCorr_",var_temp,"_bdload.png"), width=6.5, height=4
         ,gg_rankCorr_bdload)
  
  
  #### violin individual 
  gg_violin_indiv <- ggplot(samp_pred) +
    geom_violin(aes(x=species, y=samp)) +
    geom_point(data=realVal, aes(x=species, y=realVal, group=indivID, col=maxBd_outcome, pch=PABD), cex=2, position=position_dodge(width=0.8)) +
    scale_color_gradient(low="black", high="red") +
    scale_shape_manual(values=c(8,21,19)) +
    labs(pch = "Did individual\nbecome infected", col="Maximum future\nBd load") +
    xlab("Species") + ylab(allVarNames[var_temp])
  ggsave(filename = paste0("4_Bayesian_Stats/DataPlots/collapsedTime_",var_temp,"_violinAverages.png"), height=4, width=5
         ,gg_violin_indiv)
  
  #### Violin all data
  gg_violin_indiv_all <- ggplot(samp_pred) +
    geom_violin(aes(x=species, y=samp)) +
    geom_point(data=realVal_all, aes(x=species, y=realVal, col=maxBd_outcome, pch=PABD, group=indivID), cex=2, position=position_dodge(width = 0.8)) +
    geom_line(data=realVal_all, aes(x=species, y=realVal, col=maxBd_outcome, group=indivID), position=position_dodge(width = 0.8)) +
    scale_color_gradient(low="black", high="red") +
    scale_shape_manual(values=c(8,21,19)) +
    labs(pch = "Did individual\nbecome infected", col="Maximum future\nBd load") +
    xlab("Species") + ylab(allVarNames[var_temp])
  ggsave(filename = paste0("4_Bayesian_Stats/DataPlots/collapsedTime_",var_temp,"_violinAllData.png"), height=4, width=8
         ,gg_violin_indiv_all)
  ### violin_violin
  gg_violin_indiv_violinboxplot <- ggplot(samp_pred) +
    geom_violin(aes(x=species, y=samp)) +
    geom_boxplot(data=realVal_all, aes(x=species, y=realVal, fill=maxBd_outcome, group=indivID), position=position_dodge(width = 0.8)) +
    scale_fill_gradient(low="black", high="red", na.value = "grey") +
    # scale_color_gradient(low="black", high="red", na.value = "grey", guide=FALSE) +
    # scale_color_manual(values=c("grey","black","red")) +
    labs(fill="Maximum future\nBd load", col=NA) +
    xlab("Species") + ylab(allVarNames[var_temp])
  ggsave(filename = paste0("4_Bayesian_Stats/DataPlots/collapsedTime_",var_temp,"_violinboxplot.png"), height=4, width=8
         ,gg_violin_indiv_violinboxplot)
  
  #------- PROGRESS BAR ------- #
  Sys.sleep(0.1)
  i <- i+1
  # update progress bar
  setTxtProgressBar(pb, i)
  #-----------------------------#
}
close(pb)


##### @~~~~~~~~~~~~~~~Bayesian models- post~~~~~~~~~~~~~~~@ ##### 
set.seed(65334587)
seeds <- runif(length(chosen_models), min=0, max = 10000)
if ( !file.exists("4_Bayesian_Stats/all_stan_output_Bdload.RData") ) {
  all_stan_output_Bdload <- list()
  for ( var_temp in allVar) {
    # var_temp <- "inhibRich"
    ind <- match(var_temp, allVar)
    # fit model
    model_temp <- chosen_models[var_temp]
    dat_temp <- mf_alt_filt_final %>% rename_at(vars(paste0(var_temp)), ~c("response")) %>%
      mutate(Bd_exposure = factor(Bd_exposure, levels=c("Control","Bd-exposed")), prepost = factor(prepost, levels=c("Pre","Post")))
    if ( model_temp == "normal" ) {
      stan_temp <- stan_lmer(response ~ species*time + prepost + prepost:Bd_exposure + PABD + PABD:Bd_load + (1|indivID), data=dat_temp
                             , prior_intercept = normal(location = 0,scale = 5, autoscale = TRUE)
                             , prior = normal(location=0, scale=5, autoscale=TRUE)
                             , seed = seeds[ind]
                             , adapt_delta = 0.999
                             , iter = 5000
      )

    } else if (model_temp == "lognormal" ) {
      stan_temp <- stan_lmer(log(response) ~ species*time + prepost + prepost:Bd_exposure +  PABD + PABD:Bd_load + (1|indivID), data=dat_temp
                             , prior_intercept = normal(location = 0,scale = 5, autoscale = TRUE)
                             , prior = normal(location=0, scale=5, autoscale=TRUE)
                             , seed = seeds[ind]
                             , adapt_delta = 0.999
                             , iter = 5000
      )
    } else if (model_temp == "gamma" ) {
      stan_temp <- stan_glmer(response ~ species*time + prepost + prepost:Bd_exposure +  PABD + PABD:Bd_load + (1|indivID), data=dat_temp
                              , prior = normal(0, 5, autoscale = TRUE)
                              , family = Gamma(link="identity")
                              , seed = seeds[ind]
                              , adapt_delta = 0.999
                              , iter = 5000
      )
    } else if (model_temp == "beta") {
      stan_temp <- stan_glmer(response ~ species*time + prepost + prepost:Bd_exposure +  PABD + PABD:Bd_load + (1|indivID), data=dat_temp
                              , family =mgcv::betar
                              , prior_intercept = normal(location = 0.5,scale = 2.5, autoscale = TRUE)
                              , prior = normal(location=0.5, scale=2.5, autoscale=TRUE)
                              , seed= seeds[ind]
                              , adapt_delta = 0.999
                              , iter = 5000
      )
    } else if (model_temp == "poisson") {
      # SampleID random effect here for overdispersion
      stan_temp <- stan_glmer(response ~ species*time + prepost + prepost:Bd_exposure +  PABD + PABD:Bd_load + (1|indivID) + (1|SampleID), data=dat_temp
                              , prior = normal(0, 5, autoscale = TRUE)
                              , family= poisson(link="identity")
                              , seed = seeds[ind]
                              , adapt_delta = 0.999
                              , iter = 5000
      )
    }
    all_stan_output_Bdload[[var_temp]] <- stan_temp
    }
  # save(all_stan_output_Bdload,file="4_Bayesian_Stats/all_stan_output_Bdload.RData")
  save(all_stan_output_Bdload,file="4_Bayesian_Stats/all_stan_output_Bdload2.RData")
  
} else {
  # load("4_Bayesian_Stats/all_stan_output_Bdload.RData")
  load("4_Bayesian_Stats/all_stan_output_Bdload2.RData")
}

#### Get samps ####

if (!file.exists("4_Bayesian_Stats/main_samps_Bdload.RData")) {
  main_samps_Bdload <- data.frame()
  random_samps_Bdload <- data.frame()
  for ( var_temp in names(all_stan_output_Bdload)) {
    # var_temp = "inhibRich"
    stan_temp <- all_stan_output_Bdload[[var_temp]]
    samps_temp <- rstan::extract(stan_temp$stanfit)
    # Rename
    allNames_temp <- names(stan_temp$coefficients)
    main_eff_names <- allNames_temp[-grep("Intercept|b\\[", allNames_temp)]
    main_eff_names <- main_eff_names[grep("species|time|PABD|prepost|Bd_exposure",main_eff_names)]
    # Get intercept species
    allSp <- unique(mf_nonbd$species)
    intercept_name <- paste0("species",allSp[which(is.na(match(paste0("species",allSp), main_eff_names)))])
    # Get random effects
    b_eff_names <- allNames_temp[grep("b\\[", allNames_temp)]
    b_eff_sampleid_names <- gsub("^.* |\\]", "", b_eff_names[grep("SampleID", b_eff_names)])
    b_eff_individ_names <- gsub("^.* |\\]", "", b_eff_names[grep("indivID", b_eff_names)])
    
    # Main sample distributions
    long_samps_temp <- data.frame(intercept = samps_temp$alpha) %>% rename_at(vars("intercept"), ~intercept_name) %>%
      cbind(samps_temp$beta) %>% rename_at(vars(paste0(seq(1:ncol(samps_temp$beta)))), ~main_eff_names) %>% 
      rownames_to_column(var = "iter") %>% rowwise() %>% 
      mutate(Anbo = speciesAnbo, Rhma = speciesAnbo + speciesRhma, Osse = speciesAnbo + speciesOsse, Raca = speciesAnbo + speciesRaca, Rapi = speciesAnbo + speciesRapi) %>%
      mutate(timeAnbo = time, timeRhma = time + `speciesRhma:time`, timeOsse = time + `speciesOsse:time`, timeRaca = time + `speciesRaca:time`, timeRapi = time + `speciesRapi:time`) %>%
      mutate(effLiquid = prepostPost, effExposure = `prepostPost:Bd_exposureBd-exposed`, effPABD = PABD, effLoad = `PABD:Bd_load`, diffConTreat = `prepostPre:Bd_exposureBd-exposed`) %>%
      mutate(timeMean = mean(c(timeAnbo, timeRhma, timeOsse, timeRaca, timeRapi)), speciesMean = mean(c(Anbo, Rhma, Osse, Raca, Rapi))) %>%
      ungroup() %>% select(iter, Anbo, Rhma, Osse, Raca, Rapi, timeAnbo, timeRhma, timeOsse, timeRaca, timeRapi, timeMean, speciesMean, effLiquid, effExposure, effPABD, effLoad, diffConTreat) %>%
      pivot_longer(-iter, names_to = "Predictor", values_to = "samp") %>% rowwise() %>%
      mutate(Response = var_temp, Type = ifelse(length(grep("time",Predictor))>0, "Time", ifelse(Predictor %in% c("Anbo","Rhma","Osse","Raca","Rapi"), "SpeciesAverage", "ExperimentalEffects"))) 
      # mutate(Species = gsub("time","",Predictor)) %>%
      # ungroup() %>% mutate(Species = factor(Species, levels = c("Anbo","Rhma","Osse","Raca","Rapi","Mean")))
    # Random effects distributions 
    if ( length(b_eff_sampleid_names)==0) {
      long_ranef_samps_temp <- data.frame(samps_temp$b) %>% 
        rename_at(vars(paste0("X",seq(1,ncol(samps_temp$b)))), ~c(b_eff_individ_names,"variation")) %>%
        rownames_to_column(var = "iter") %>%
        pivot_longer(-iter, names_to = "RandomEffect", values_to = "samp") %>%
        separate(RandomEffect, into = c("Type", "Predictor"), sep = ":") %>%
        mutate(Response = var_temp, Species = gsub("_.*$","",Predictor))
    } else {
      long_ranef_samps_temp <- data.frame(samps_temp$b) %>% 
        rename_at(vars(paste0("X",seq(1,ncol(samps_temp$b)))), ~c(b_eff_sampleid_names,"variation1",b_eff_individ_names,"variation2")) %>%
        rownames_to_column(var = "iter") %>%
        pivot_longer(-iter, names_to = "RandomEffect", values_to = "samp") %>%
        separate(RandomEffect, into = c("Type", "Predictor"), sep = ":") %>%
        mutate(Response = var_temp, Species = gsub("_.*$","",Predictor))
    }
    
    ggsave(filename = paste0("4_Bayesian_Stats/BayesianSamplePosteriors/Bd_load_",var_temp,"_MainEffects_posteriors.pdf"), height=10, width=10
           ,ggplot(long_samps_temp) +
             geom_histogram(aes(x=samp)) + geom_vline(aes(xintercept=0), col="red")+
             facet_wrap(Type ~ Predictor, scales="free", nrow=4))
    
    main_samps_Bdload <- rbind(main_samps_Bdload, long_samps_temp)
    random_samps_Bdload <- rbind(random_samps_Bdload, long_ranef_samps_temp)
  } 
  save(main_samps_Bdload, file="4_Bayesian_Stats/main_samps_Bdload.RData")
  save(random_samps_Bdload, file="4_Bayesian_Stats/random_samps_Bdload.RData")
} else {
  load("4_Bayesian_Stats/main_samps_Bdload.RData")
  load("4_Bayesian_Stats/random_samps_Bdload.RData")
}


#### Summarize main samps ####
summary_main_effects_Bdload <- main_samps_Bdload %>%
  group_by(Response, Type, Predictor) %>%
  summarize(median = median(samp), lower95 = ci(samp, method="HDI")$CI_low, upper95 = ci(samp, method="HDI")$CI_high, postprob = max(c(sum(samp>0)/n(), 1-sum(samp>0)/n()))) %>%
  mutate(sig = postprob>0.95) %>%
  mutate(sigCI = lower95*upper95>0) %>%
  ungroup() %>% arrange(rev(Type), Response, Predictor) 

summary_randeff_Bdload <- random_samps_Bdload %>%
  filter(Type == "indivID") %>%
  group_by(Response, Predictor, Species) %>%
  summarize(median = median(samp), lower95 = ci(samp, method="HDI")$CI_low, upper95 = ci(samp, method="HDI")$CI_high)

#### @~~~~~~~~~~~~~~~Plotting effects~~~~~~~~~~~~~~~@ ####

speciesMeans <- summary_main_effects_control %>% filter(Response %in% allVarFilt, Species =="Mean", Type == "SpeciesAverage") %>%
  select(Response, median, lower95, upper95) %>% rename(speciesMean=median, specieslower=lower95, speciesupper = upper95)
allMeans <- summary_main_effects_control %>% filter(Response %in% allVarFilt, Species =="Mean", Type%in% c("SpeciesAverage","Time")) %>%
  select(Response, Type, median, lower95, upper95) %>% rename(meanall=median, allLower=lower95, allUpper = upper95)
allMeansTimeZero <- allMeans %>%
  mutate(meanAll=ifelse(Type=="Time", 0, NA))

colorForSig <- c("blue", "purple", "black")
names(colorForSig) <- c(">97.5%", ">95%", "<=94%")

fullLabelLevels <- paste0(rep(c("Effect of species on\n","Effect of time on\n"),each=length(allVarFilt)), as.vector(rep(allVarNames[allVarFilt], 2)))

gg_speciesandtime <- summary_main_effects_control %>% filter(Response %in% allVarFilt, Species != "Mean", Type%in%c("SpeciesAverage","Time")) %>%
  left_join(allMeansTimeZero) %>%
  rowwise() %>%
  mutate(postprob = ifelse(Type=="SpeciesAverage", NA, postprob)) %>% ungroup() %>%
  # mutate(sig = postprob>0.95) %>%
  mutate(PosteriorProbability = ifelse(postprob>0.975,">97.5%", ifelse(postprob>0.95, ">95%", "<=94%"))) %>% 
  mutate(PosteriorProbability = factor(ifelse(is.na(PosteriorProbability), "Not applicable", PosteriorProbability), levels=c(">97.5%", ">95%", "<=94%","Not applicable"))) %>%
  ungroup() %>% mutate(newLabels = paste0(ifelse(Type=="SpeciesAverage", "Effect of species on\n", "Effect of time on\n"), allVarNames[Response])) %>%
  mutate(newLabels = factor(newLabels, levels=fullLabelLevels)) %>%
  # mutate(Response = factor(Response, levels=allVarFilt)) %>%
  ggplot() + 
  geom_point(aes(x=Species, y=median, col=PosteriorProbability, cex=PosteriorProbability), position = position_dodge(width=0.8)) +
  geom_linerange(aes(x=Species, ymin=lower95, ymax=upper95, col=PosteriorProbability), position = position_dodge(width=0.8)) +
  # geom_hline(aes(yintercept=speciesMean), col="red", alpha=0.5) +
  # scale_color_manual(values=c("grey","blue")) +
  scale_color_manual(values=c(colorForSig, `Not applicable` = "darkgrey"))+
  scale_size_manual(values=c(3,2,1,1))+
  geom_hline(aes(yintercept = meanAll), col="red", alpha=0.5) +
  # facet_wrap(Response~Type, scales = "free", nrow=2,  dir = "v", labeller = as_labeller(c(allVarNames, c(SpeciesAverage="Species coefficient", Time = "Effect of time by species"))), 
  #            , strip.position = "left") +  
  facet_wrap(newLabels~., scales = "free", nrow=2,  dir = "h",strip.position = "top") +
  ylab("Coefficient estimates for each amphibian species") + xlab("Amphibian Species") +#+labs(col="95% credible interval\noverlaps dataset average") 
  labs(col="One-sided\nprobability of direction", cex="One-sided\nprobability of direction") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5), strip.text.y = element_text(size = 10), legend.position = "bottom")
gg_speciesandtime
ggsave("4_Bayesian_Stats/BayesianSamplePosteriors/effectSpeciesAndTime_control.png"
       , height=6, width=13
       , gg_speciesandtime)


gg_effectspecies <- summary_main_effects_control %>% filter(Response %in% allVarFilt, Species != "Mean", Type=="SpeciesAverage") %>%
  left_join(speciesMeans) %>% 
  rowwise() %>%
  mutate(sig = !between(speciesMean, lower95, upper95)) %>%
  mutate(sig2 = !between(median, specieslower, speciesupper)) %>% ungroup() %>%
  mutate(Response = factor(Response, levels=allVarFilt)) %>%
  ggplot() + 
  geom_point(aes(x=Species, y=median), position = position_dodge(width=0.8)) +
  geom_linerange(aes(x=Species, ymin=lower95, ymax=upper95), position = position_dodge(width=0.8)) +
  geom_hline(aes(yintercept=speciesMean), col="red", alpha=0.5) +
  # scale_color_manual(values=c("grey","blue")) +
  facet_wrap(.~Response, scales = "free", labeller = as_labeller(allVarNames), nrow=1) +
  # coord_flip() +
  ylab("Species coefficient estimates for: ") + xlab("Species predictor") +#+labs(col="95% credible interval\noverlaps dataset average") 
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
gg_effectspecies
ggsave("4_Bayesian_Stats/BayesianSamplePosteriors/effectSpecies_control.png"
       , height=4, width=13
       , gg_effectspecies)

gg_effecttime <- summary_main_effects_control %>% filter(Response %in% allVarFilt, Species != "Mean", Type=="Time",Species !="Mean") %>%
  mutate(Response = factor(Response, levels=allVarFilt)) %>%
  mutate(sig = postprob>0.95) %>%
  mutate(PosteriorProbability = factor(ifelse(postprob>0.975,">97.5%", ifelse(postprob>0.95, ">95%", "<=94%")), levels=c(">97.5%", ">95%", "<=94%"))) %>%
  ggplot() + 
  geom_point(aes(x=Species, y=median, col=PosteriorProbability, cex=PosteriorProbability), position = position_dodge(width=0.8)) +
  geom_linerange(aes(x=Species, ymin=lower95, ymax=upper95, col=PosteriorProbability), position = position_dodge(width=0.8)) +
  # scale_shape_manual(values=c(21,19)) +
  # scale_linetype_manual(values=c(3,1))+
  scale_color_manual(values=c("blue","purple","black"))+
  scale_size_manual(values=c(3,2,1))+
  geom_hline(aes(yintercept=0), col="red", alpha=0.5) +
  scale_shape_manual(values=c(21,19)) +
  facet_wrap(.~Response, scales = "free", labeller = as_labeller(allVarNames[allVarFilt]), nrow=1, strip.position = "left") +
  # coord_flip() +
  ylab("Effect of time on:") +labs(col="One-sided\nprobability of direction", cex="One-sided\nprobability of direction") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
  
gg_effecttime
ggsave("4_Bayesian_Stats/BayesianSamplePosteriors/effectTime_control.png"
       , height=4, width=13
       , gg_effecttime)

predictorNamed <- c(effPABD = "Effect of\npresence/absence infection", effLoad= "Effect of\n Bd load (after accounting\nfor presence/absence)", effLiquid="Effect of\ncontrol (exposure to liquid)", effExposure="Effect of\nBd exposure (after accounting\nfor control)", diffConTreat = "Difference between\ncontrol and treatment groups\n(should be zero)")

# gg_experimentaleffects <- summary_main_effects_Bdload %>% filter(Response %in% allVarFilt, Type == "ExperimentalEffects", Predictor !="speciesMean") %>%
#   mutate(Response = factor(Response, levels=allVarFilt)) %>%
#   rowwise() %>%
#   mutate(Predictors = predictorNamed[Predictor]) %>% ungroup() %>%
#   mutate(Predictors = factor(Predictors, levels=sort(predictorNamed)[c(2,5,3,4,1)])) %>%
#   mutate(PosteriorProbability = factor(ifelse(postprob>0.975,">97.5%", ifelse(postprob>0.95, ">95%", "<=94%")), levels=c(">97.5%", ">95%", "<=94%"))) %>%
#   ggplot() + 
#   geom_point(aes(x=Predictors, y=median, col=PosteriorProbability, cex=PosteriorProbability), position = position_dodge(width=0.8)) +
#   geom_linerange(aes(x=Predictors, ymin=lower95, ymax=upper95, col=PosteriorProbability), position = position_dodge(width=0.8)) +
#   geom_hline(aes(yintercept=0), col="red", alpha=0.5) +
#   scale_color_manual(values=c("blue","purple","black")) +
#   scale_size_manual(values=c(3,2,1)) +
#   theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
#   facet_grid(.~Response, scales = "free", labeller = as_labeller(allVarNames[allVarFilt])) +
#   coord_flip() +ylab("Effect (Coefficient posterior)") +labs(col="One-sided\nposterior probability",cex="One-sided\nposterior probability") +
#   theme_bw()
gg_experimentaleffects <- summary_main_effects_Bdload %>% filter(Response %in% allVarFilt, Type == "ExperimentalEffects", Predictor !="speciesMean") %>%
  mutate(Response = factor(Response, levels=allVarFilt)) %>%
  rowwise() %>%
  mutate(Predictors = predictorNamed[Predictor]) %>% ungroup() %>%
  mutate(Predictors = factor(Predictors, levels=rev(sort(predictorNamed)[c(2,5,3,4,1)]))) %>%
  mutate(PosteriorProbability = factor(ifelse(postprob>0.975,">97.5%", ifelse(postprob>0.95, ">95%", "<=94%")), levels=c(">97.5%", ">95%", "<=94%"))) %>%
  ggplot() + 
  geom_point(aes(x=Predictors, y=median, col=PosteriorProbability, cex=PosteriorProbability), position = position_dodge(width=0.8)) +
  geom_linerange(aes(x=Predictors, ymin=lower95, ymax=upper95, col=PosteriorProbability), position = position_dodge(width=0.8)) +
  geom_hline(aes(yintercept=0), col="red", alpha=0.5) +
  scale_color_manual(values=c("blue","purple","black")) +
  scale_size_manual(values=c(3,2,1)) +
  facet_wrap(.~Response, scales = "free", labeller = as_labeller(allVarNames[allVarFilt]), nrow=2) +
  # coord_flip() +
  ylab("Effect of predictor on each microbial community trait") +labs(col="One-sided\nprobability of direction",cex="One-sided\nprobability of direction") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) 
gg_experimentaleffects
ggsave("4_Bayesian_Stats/BayesianSamplePosteriors/effectExperiment.png"
       , height=8, width=12
       , gg_experimentaleffects)

#### @~~~~~~~~~~~~~~~Printing/saving~~~~~~~~~~~~~~~@ ####

save(summary_main_effects_control, file="4_Bayesian_Stats/summarystats/summary_main_effects_control.RData")
summary_indivID_vs_bd_control <- summary_randeff_control %>% rename(indivID = Predictor) %>% full_join(mf_collapsed_bd) 
save(summary_indivID_vs_bd_control, file="4_Bayesian_Stats/summarystats/summary_indivID_vs_bd_control.RData")
# 
# summary_main_effects_control %>% filter(Type !="SpeciesAverage", Predictor=="timeMean") %>%
#   filter(Response %in% c("shannon","observed_otus","inhibRich","percInhib","dist_braycurtis","disper_braycurtis"))

save(summary_main_effects_Bdload, file = "4_Bayesian_Stats/summarystats/summary_main_effects_Bdload.RData")
# summary_indivID_vs_bd_baseline <- mf_collapsed_bd %>% full_join(summary_randeff_basline) 
save(summary_randeff_Bdload, file="4_Bayesian_Stats/summarystats/summary_randeff_Bdload.RData")

# save(summary_main_effects_baseline, file="4_Bayesian_Stats/summary_main_effects_baseline.RData")
# summary_indivID_vs_bd_baseline <- summary_randeff_baseline %>% rename(indivID = Predictor) %>% full_join(mf_collapsed_bd) 
# save(summary_indivID_vs_bd_baseline, file="4_Bayesian_Stats/summary_indivID_vs_bd_baseline.RData")

# Text versions
toWrite_summary_main_control <- summary_main_effects_control %>% 
  mutate(postprobcomp = 1-postprob, significant = postprobcomp<0.025)
toWrite_summary_rand_control <- summary_randeff_control %>% 
  mutate(significant = sign(lower95)==sign(upper95))  
write.table(toWrite_summary_main_control, file = "4_Bayesian_Stats/summarystats/summary_main_effects_control.txt", quote=FALSE, sep = "\t", row.names = FALSE)
write.table(toWrite_summary_rand_control, file = "4_Bayesian_Stats/summarystats/summary_randeff_control.txt", quote=FALSE, sep = "\t", row.names = FALSE)
# 
# toWrite_summary_main_baseline <- summary_main_effects_baseline %>% 
#   mutate(postprobcomp = 1-postprob, significant = postprobcomp<0.025)
# toWrite_summary_rand_baseline <- summary_randeff_baseline %>% 
#   mutate(significant = sign(lower95)==sign(upper95))  
# write.table(toWrite_summary_main_baseline, file = "4_Bayesian_Stats/summary_main_effects_baseline.txt", quote=FALSE, sep = "\t", row.names = FALSE)
# write.table(toWrite_summary_rand_baseline, file = "4_Bayesian_Stats/summary_randeff_baseline.txt", quote=FALSE, sep = "\t", row.names = FALSE)

toWrite_summary_main_Bdload <- summary_main_effects_Bdload %>% 
  mutate(postprobcomp = 1-postprob, significant = postprobcomp<0.025)
toWrite_summary_rand_Bdload <- summary_randeff_Bdload %>% 
  mutate(significant = sign(lower95)==sign(upper95))  
write.table(toWrite_summary_main_Bdload, file = "4_Bayesian_Stats/summarystats/summary_main_effects_Bdload.txt", quote=FALSE, sep = "\t", row.names = FALSE)
write.table(toWrite_summary_rand_Bdload, file = "4_Bayesian_Stats/summarystats/summary_randeff_Bdload.txt", quote=FALSE, sep = "\t", row.names = FALSE)

##### @~~~~~~~~~~~~~~~Follow-up correlations with load~~~~~~~~~~~~~~~@ #####

#### Custom Bayes models for predictor uncertainty #####

#### Find sd of ranks
# View(mf_collapsed_bd)
allRanks_summary <- allRanks %>%
  filter(prepost=="Pre", Bd_exposure == "Bd-exposed") %>%
  group_by(species, indivID, Response) %>%
  summarise(medrank = median(prank, na.rm = TRUE), sdrank = sd(prank, na.rm=TRUE), serank = sd(prank, na.rm=TRUE)/sqrt(n()), nindiv=n()) %>% ungroup() %>%
  left_join(mf_collapsed_bd) #%>%
  # mutate(maxBd_outcome = ifelse(maxBd_outcome<log(2),0,maxBd_outcome)) %>%
  # mutate(PABD = ifelse(maxBd_outcome>0,1,0)) 
### Fix Osse? Missing individual-level variation because timelines are crunched together.
# meanOsseSD <- allRanks %>%
#   filter(prepost=="Pre", Bd_exposure == "Bd-exposed", Response=="dist_weighted_unifrac") %>%
#   filter(species=="Osse") %>% group_by(indivID) %>%
#   summarise(sdrank = sd(prank, na.rm=TRUE)) %>% pull(sdrank) %>% mean(na.rm=TRUE)
# allRanks_summary <- allRanks_summary %>% mutate(sdrank = ifelse(is.na(sdrank) & species=="Osse", meanOsseSD, sdrank)) %>%
#   mutate(serank = ifelse(is.na(serank) & species=="Osse", sdrank/sqrt(nindiv), serank))

meanOsseSE <- allRanks %>%
  filter(prepost=="Pre", Bd_exposure == "Bd-exposed", Response=="dist_weighted_unifrac") %>%
  filter(species=="Osse") %>% group_by(indivID) %>%
  summarise(serank = sd(prank, na.rm=TRUE)/n()) %>% pull(serank) %>% mean(na.rm=TRUE)
meanOsseSD <- allRanks %>%
  filter(prepost=="Pre", Bd_exposure == "Bd-exposed", Response=="dist_weighted_unifrac") %>%
  filter(species=="Osse") %>% group_by(indivID) %>%
  summarise(sdrank = sd(prank, na.rm=TRUE)) %>% pull(sdrank) %>% mean(na.rm=TRUE)

allRanks_summary <- allRanks_summary %>% mutate(serank = ifelse(is.na(serank) & species=="Osse", meanOsseSE, serank)) %>% 
  mutate(sdrank = ifelse(is.na(sdrank) & species=="Osse", meanOsseSD, sdrank)) 

# View(allRanks_summary)
indivList_treat <- unique((mf_treat$indivID))
if ( !file.exists("4_Bayesian_Stats/allModels_predictBd.RData") ) {
  # set.seed(93845)
  # set.seed(105293)
  # set.seed(6734)
  
  total <- length(allVarFilt)*2
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  i <- 1
  allModels_predictBd <- list()
  aggregatedData_baseline <- data.frame()
  for (bd_type in c("maxBd_outcome","PABD")) {
    # bd_type <- "PABD"

    allModels_predictBd[[paste0(bd_type)]] <- list()
    
    for (var_temp in allVarFilt) {
      
      if ( bd_type == "PABD") {
        data_stan_temp <- allRanks_summary %>% filter(Response==var_temp) %>%
          rename_at(vars(paste0("PABD")), ~c("y")) %>%
          filter(!is.na(y)) %>% 
          mutate(x_meas = medrank) %>%
          mutate(x_meas = (medrank-mean(medrank))/sd(medrank)) %>%
          mutate(sdMean = serank/sd(medrank))
        modelUse <- "4_Bayesian_CustomStan/stan_model_binomial_ME.stan"
        ## Get data formatted
        data_stan_format <- list(N=nrow(data_stan_temp)
                                 , y=data_stan_temp$y
                                 , x_meas = data_stan_temp$medrank
                                 , tau = data_stan_temp$serank
                                 , x_Anbo = as.numeric(data_stan_temp$species=="Anbo")
                                 , x_Rhma = as.numeric(data_stan_temp$species=="Rhma")
                                 , x_Osse = as.numeric(data_stan_temp$species=="Osse")
                                 , x_Raca = as.numeric(data_stan_temp$species=="Raca")
                                 , x_Rapi = as.numeric(data_stan_temp$species=="Rapi")
        )
        fit_temp <- stan(file = modelUse
                         , data=data_stan_format
                         , warmup=5000, iter=10000, chains=4
                         , control = list(adapt_delta = 0.9999, max_treedepth=10)
                         , verbose = FALSE)
        
      } else {
        data_stan_temp <- allRanks_summary %>% filter(Response==var_temp) %>%
          rename_at(vars(paste0("maxBd_outcome")), ~c("y_raw")) %>%
          filter(!is.na(y_raw)) %>%  filter(y_raw>0) %>%
          filter(species!="Rapi") %>%
          mutate( y = (y_raw-mean(y_raw))) %>%
          # mutate( y = y_raw) %>%
          mutate(x_meas = (medrank-mean(medrank))/sd(medrank)) %>%
          mutate(sdMean = serank/sd(medrank))
        # View(data_stan_temp)
        # ggplot(data_stan_temp) + geom_point(aes(x=medrank, y=y, col=species))
        modelUse <- "4_Bayesian_CustomStan/stan_model_normal_ME_norapi.stan"
        # modelUse <- "4_Bayesian_CustomStan/stan_model_normal_ME_norapi_ysd.stan"
        ## Get data formatted
        data_stan_format <- list(N= nrow(data_stan_temp)
                                 , y= data_stan_temp$y
                                 , x_meas = data_stan_temp$medrank
                                 , tau = data_stan_temp$serank
                                 # , ysd = data_stan_temp$sdBdMax
                                 , x_Anbo = as.numeric(data_stan_temp$species=="Anbo")
                                 , x_Rhma = as.numeric(data_stan_temp$species=="Rhma")
                                 , x_Osse = as.numeric(data_stan_temp$species=="Osse")
                                 , x_Raca = as.numeric(data_stan_temp$species=="Raca")
        )
        fit_temp <- stan(file = modelUse
                         , data=data_stan_format
                         , warmup=1500, iter=5000, chains=4
                         , control = list(adapt_delta = 0.9999999, max_treedepth=20)
                         , verbose = FALSE)
        # plot(fit_temp)
        # fit_temp <- stan(file = modelUse
        #                  , data=data_stan_format
        #                  , warmup=5000, iter=10000, chains=4
        #                  , control = list(adapt_delta = 0.99999, max_treedepth=20)
        #                  , verbose = FALSE)
        # plot(fit_temp)
        
      }
      # 
      # data_stan_temp %>% 
      #   ggplot() + geom_point(aes(x=x_meas, y=y, col=species)) +
      #   geom_segment(aes(x=x_meas-sdMean, xend = x_meas+sdMean, y=y, yend=y))
      # 
      
      ### Save models and data
      allModels_predictBd[[paste0(bd_type)]][[paste0(var_temp)]] <- fit_temp
      
      ## Aggregate data
      aggregatedData_baseline <- data_stan_temp %>% select(species, indivID, species, x_meas, sdMean, medrank, sdrank, serank) %>%
        mutate(micro = var_temp) %>%
        rbind(aggregatedData_baseline)
      
      #------- PROGRESS BAR ------- #
      Sys.sleep(0.1)
      i <- i+1
      # update progress bar
      setTxtProgressBar(pb, i)
      #-----------------------------#
    }
    
    
  }
  close(pb)
  
  aggregatedData_baseline <- aggregatedData_baseline %>% distinct()
  # save(allModels_predictBd, file = "4_Bayesian_Stats/allModels_predictBd.RData")
  # save(aggregatedData_baseline, file="4_Bayesian_Stats/aggregatedData_baseline.RData")
} else {
  load("4_Bayesian_Stats/allModels_predictBd.RData")
  load("4_Bayesian_Stats/aggregatedData_baseline.RData")
}

#### Aggregate #####
dir.create("4_Bayesian_Stats/UsedRanksPlots")

for (var_temp in allVarFilt) {
  # var_temp <- "inhibRich"
  
  agg_temp <- aggregatedData_baseline %>% filter(micro==var_temp) %>% rowwise() %>%
    mutate(indivSamp = medrank, indivSamp_lowerse = medrank-serank, indivSamp_uperse = medrank + serank) %>%
    select(species, indivID, indivSamp, indivSamp_lowerse, indivSamp_uperse) %>%
    left_join(mf_collapsed_bd) %>% 
    arrange(-maxBd_outcome) %>% distinct() %>%
    mutate(maxBd_outcome_allzero = maxBd_outcome) %>%
    mutate(maxBd_outcome = ifelse(maxBd_outcome==0,NA,maxBd_outcome)) %>%
    mutate(indivID = factor(indivID, levels=(unique(indivID))))
  
  ggused_ranks <- ggplot() +
    geom_linerange(data=agg_temp, aes(x=species, ymin=indivSamp_lowerse, ymax=indivSamp_uperse),  position=position_dodge2(width=0.6)) +
    geom_point(data=agg_temp, aes(x=species, y=indivSamp, fill=maxBd_outcome, pch=factor(PABD)), cex = 3, position=position_dodge2(width=0.6)) +
    geom_vline(data=data.frame(divide=c(1.5,2.5,3.5,4.5)), aes(xintercept = divide), lwd=3, col="white") +
    # scale_color_manual(values=c(NA, "black")) +
    scale_fill_gradient(low="grey", high="red", na.value = "white") + 
    scale_shape_manual(values=c(19,21)) +
    labs(fill="Future maximum\n Bd load (log)", pch="Did individual\nbecome infected?", col="Did individual\nbecome infected?") +
    xlab("Maximum Bd load") + ylab(paste0("Percentile rank of\n",allVarNames[var_temp])) +
    theme(panel.grid.major.x  = element_blank())
  # ggused_ranks
  ggsave(filename = paste0("4_Bayesian_Stats/UsedRanksPlots/usedRanks_",var_temp,".png"), height=4, width=6
         ,ggused_ranks)
  
  ggused_corr <- ggplot() +
    geom_linerange(data=agg_temp %>% filter(!is.na(maxBd_outcome)), aes(x=maxBd_outcome, ymin=indivSamp_lowerse, ymax=indivSamp_uperse)) +
    geom_point(data=agg_temp %>% filter(!is.na(maxBd_outcome)), aes(x=maxBd_outcome, y=indivSamp, col=species),cex = 3) +
    coord_flip()+
    # geom_vline(data=data.frame(divide=c(1.5,2.5,3.5,4.5)), aes(xintercept = divide), lwd=3, col="white") +
    # scale_color_manual(values=c("white", "black")) +
    # scale_fill_gradient(low="grey", high="red", na.value = "white") + 
    # scale_shape_manual(values=c(21,19)) +
    labs(col="Species") +
    xlab("Maximum Bd load") + ylab(paste0("Percentile rank of\n",allVarNames[var_temp])) +
    theme(panel.grid.major.x  = element_blank())
  # ggused_corr
  ggsave(filename = paste0("4_Bayesian_Stats/UsedRanksPlots/usedRanksCorr_",var_temp,".png"), height=4, width=6
         ,ggused_corr)
  
  ggused_bin <- ggplot() +
    geom_linerange(data=agg_temp, aes(x=PABD, ymin=indivSamp_lowerse, ymax=indivSamp_uperse, group=factor(indivID)),  position=position_dodge2(width=0.3)) +
    geom_point(data=agg_temp, aes(x=PABD, y=indivSamp, fill=species, col=factor(PABD), group=factor(indivID)),pch=21, cex = 3, position=position_dodge2(width=0.3)) +
    coord_flip()+
    # geom_vline(data=data.frame(divide=c(1.5,2.5,3.5,4.5)), aes(xintercept = divide), lwd=3, col="white") +
    scale_color_manual(values=c("white", "black")) +
    # scale_fill_gradient(low="grey", high="red", na.value = "white") + 
    # scale_shape_manual(values=c(21,19)) +
    labs(fill="Species", pch="Did individual\nbecome infected?", col="Did individual\nbecome infected?") +
    xlab("Presence/absence infection") + ylab(paste0("Percentile rank of\n",allVarNames[var_temp])) +
    theme(panel.grid.major.x  = element_blank())
  ggused_bin
  ggsave(filename = paste0("4_Bayesian_Stats/UsedRanksPlots/usedRanksBin_",var_temp,".png"), height=4, width=6
         ,ggused_bin)
  
}
#### Extract effect of micro #####
allSamps_bdresponse <- data.frame()
for (bd_type in c("PABD", "maxBd_outcome")) {
  # bd_type <- "PABD"
  for (var_temp in allVarFilt) {
    # var_temp <- "inhibRich"
    model_temp <- allModels_predictBd[[paste0(bd_type)]][[paste0(var_temp)]]
    samps_temp <- rstan::extract(model_temp)
    if (bd_type == "PABD") {
      allsamps_temp <- data.frame(Anbo = samps_temp$b_Anbo, Rhma = samps_temp$b_Rhma, Osse = samps_temp$b_Osse, Raca = samps_temp$b_Raca, Rapi = samps_temp$b_Rapi, micro = samps_temp$b_micro) %>%
        rownames_to_column(var = "iter") %>%
        pivot_longer(-iter, names_to = "Predictor", values_to = "samp") %>%
        mutate(MicrobiomePredictor = var_temp, Response = bd_type)
    } else {
      allsamps_temp <- data.frame(Anbo = samps_temp$b_Anbo, Rhma = samps_temp$b_Rhma, Osse = samps_temp$b_Osse, Raca = samps_temp$b_Raca, micro = samps_temp$b_micro) %>%
        rownames_to_column(var = "iter") %>%
        pivot_longer(-iter, names_to = "Predictor", values_to = "samp") %>%
        mutate(MicrobiomePredictor = var_temp, Response = bd_type)
    }
    
    allSamps_bdresponse <- rbind(allSamps_bdresponse,allsamps_temp)
  }
}

allSamps_bdresponse_adj <- allSamps_bdresponse %>% 
  mutate(Predictor = factor(Predictor, levels=c("Anbo","Rhma","Osse","Raca","Rapi","micro"))) %>%
  rowwise() %>% mutate(MP = allVarNames[MicrobiomePredictor]) %>%
  mutate(MP = factor(MP, levels=rev(allVarNames))) %>%
  group_by(Predictor, MP, Response) %>%
  summarise(medsamp = median(samp), lower95 = ci(samp, method="HDI")$CI_low, upper95 = ci(samp, method="HDI")$CI_high, postprob = max(c(sum(samp>0), sum(samp<0)))/n()) %>% 
  ungroup() %>%
  mutate(sig = postprob>0.95)
#### @~~~~~~~~~~~~~~~Plotting effects~~~~~~~~~~~~~~~@ ####

## Look at species-level
colors_of_traits <- c("red","darkred","purple","pink","green","blue")
names(colors_of_traits) <- allVarNames[allVarFilt]
# gg_species <- allSamps_bdresponse_adj %>%
#   filter(Predictor !="micro") %>%
#   ggplot() + 
#   geom_linerange(aes(x=Predictor, ymin=lower95, ymax=upper95, col=MP), position = position_dodge(width=0.8))+
#   geom_point(aes(x=Predictor, y=medsamp, col=MP), position = position_dodge(width=0.8)) +
#   geom_vline(data=data.frame(divide=c(1.5,2.5,3.5,4.5)),aes(xintercept=divide), col="darkgrey", lwd=0.5)+
#   facet_grid(. ~ Response, labeller = as_labeller(c(maxBd_outcome = "Maximum future Bd load\n (log)", PABD = "Probability of infection\n(Presence/absence)"))) +
#   coord_flip() + labs(col="Microbiome community trait") +
#   scale_color_manual(values=c(colors_of_traits)) + xlab("Species") + ylab("Effect on infection load or probability") +
#   theme_bw()
gg_species <- allSamps_bdresponse_adj %>%
  filter(Predictor !="micro") %>%
  ggplot() + 
  geom_linerange(aes(x=Predictor, ymin=lower95, ymax=upper95, col=MP), position = position_dodge(width=0.8))+
  geom_point(aes(x=Predictor, y=medsamp, col=MP), position = position_dodge(width=0.8)) +
  geom_vline(data=data.frame(divide=c(1.5,2.5,3.5,4.5)),aes(xintercept=divide), col="darkgrey", lwd=0.5)+
  facet_wrap(. ~ Response, ncol =1, labeller = as_labeller(c(maxBd_outcome = "Maximum future Bd load\n (log)", PABD = "Probability of infection\n(Presence/absence)"))
             , strip.position="left") +
  # coord_flip() + 
  labs(col="Microbial community trait") +
  scale_color_manual(values=c(colors_of_traits)) + xlab("Species") + ylab("Effect of species on") +
  theme_bw()+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
gg_species
ggsave(filename = "4_Bayesian_Stats/BayesianSamplePosteriors/bayes_custom_species_effects.png", height=10, width=7
, gg_species)

colorForSig <- c("blue", "purple", "black")
names(colorForSig) <- c(">97.5%", ">95%", "<=94%")
# gg_microbiome <- allSamps_bdresponse_adj %>%
#   filter(Predictor == "micro")  %>%
#   mutate(Response = factor(Response, levels=c("PABD","maxBd_outcome"))) %>%
#   mutate(sigCI = lower95*upper95>0) %>%
#   mutate(PosteriorProbability = factor(ifelse(postprob>0.975,">97.5%", ifelse(postprob>0.95, ">95%", "<=94%")), levels=c(">97.5%", ">95%", "<=94%"))) %>%
#   ggplot() + 
#   geom_point(aes(x=MP, y=medsamp, col=PosteriorProbability, cex=PosteriorProbability)) +
#   geom_segment(aes(x=MP, xend = MP, y=lower95, yend=upper95, col=PosteriorProbability))+
#   geom_hline(aes(yintercept=0), col="red",alpha=0.5)+
#   facet_grid(. ~ Response, labeller = as_labeller(c(maxBd_outcome = "Maximum future Bd load\n (log)", PABD = "Probability of infection\n(Presence/absence)"))) +
#   coord_flip() +labs(col="One-sided\nposterior probability", cex="One-sided\nposterior probability") +
#   scale_color_manual(values=colorForSig) +
#   scale_size_manual(values=c(3,2,1)) +
#   xlab("Microbial community trait") + ylab("Effect on infection load or probability\n(Coefficient estimate)")+
#   theme_bw()
gg_microbiome <- allSamps_bdresponse_adj %>%
  filter(Predictor == "micro")  %>%
  mutate(MP = factor(MP, levels=c(allVarNames[allVarFilt]))) %>%
  mutate(Response = factor(Response, levels=c("PABD","maxBd_outcome"))) %>%
  mutate(sigCI = lower95*upper95>0) %>%
  mutate(PosteriorProbability = factor(ifelse(postprob>0.975,">97.5%", ifelse(postprob>0.95, ">95%", "<=94%")), levels=c(">97.5%", ">95%", "<=94%"))) %>%
  ggplot() + 
  geom_point(aes(x=MP, y=medsamp, col=PosteriorProbability, cex=PosteriorProbability)) +
  geom_segment(aes(x=MP, xend = MP, y=lower95, yend=upper95, col=PosteriorProbability))+
  geom_hline(aes(yintercept=0), col="red",alpha=0.5)+
  facet_wrap(. ~ Response, ncol=1, strip.position = "left", labeller = as_labeller(c(maxBd_outcome = "Maximum future Bd load\n (log)", PABD = "Probability of infection\n(Presence/absence)"))) +
  # coord_flip() +
  labs(col="One-sided\nprobability of direction", cex="One-sided\nprobability of direction") +
  scale_color_manual(values=colorForSig) +
  scale_size_manual(values=c(3,2,1)) +
  xlab("Microbial community trait predictor") + ylab("Effect of microbial community trait on:")+
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5), legend.position = "right")

gg_microbiome
ggsave(filename = "4_Bayesian_Stats/BayesianSamplePosteriors/bayes_custom_microbiome_effects.png", height=7, width=6
,gg_microbiome)

##### Multipanel figure ####
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
spCol <- gg_color_hue(5)
names(spCol) <- c("Anbo","Rhma","Osse","Raca","Rapi")
agg_inhibRich <- aggregatedData_baseline %>% filter(micro=="inhibRich") %>% rowwise() %>%
  mutate(indivSamp = medrank, indivSamp_lowerse = medrank-serank, indivSamp_uperse = medrank + serank) %>%
  select(species, indivID, indivSamp, indivSamp_lowerse, indivSamp_uperse, micro) %>%
  left_join(mf_collapsed_bd) %>% 
  arrange(-maxBd_outcome) %>% distinct() %>%
  mutate(maxBd_outcome_allzero = maxBd_outcome) %>%
  mutate(maxBd_outcome = ifelse(maxBd_outcome==0,NA,maxBd_outcome)) %>%
  mutate(indivID = factor(indivID, levels=(unique(indivID)))) 

agg_percInhib<- aggregatedData_baseline %>% filter(micro %in% c("percInhib")) %>% rowwise() %>%
  mutate(indivSamp = medrank, indivSamp_lowerse = medrank-serank, indivSamp_uperse = medrank + serank) %>%
  select(species, indivID, indivSamp, indivSamp_lowerse, indivSamp_uperse, micro) %>%
  left_join(mf_collapsed_bd) %>% 
  arrange(-maxBd_outcome) %>% distinct() %>%
  mutate(maxBd_outcome_allzero = maxBd_outcome) %>%
  mutate(maxBd_outcome = ifelse(maxBd_outcome==0,NA,maxBd_outcome)) %>%
  mutate(indivID = factor(indivID, levels=(unique(indivID))))

agg_dist<- aggregatedData_baseline %>% filter(micro %in% c("dist_weighted_unifrac")) %>% rowwise() %>%
  mutate(indivSamp = medrank, indivSamp_lowerse = medrank-serank, indivSamp_uperse = medrank + serank) %>%
  select(species, indivID, indivSamp, indivSamp_lowerse, indivSamp_uperse, micro) %>%
  left_join(mf_collapsed_bd) %>% 
  arrange(-maxBd_outcome) %>% distinct() %>%
  mutate(maxBd_outcome_allzero = maxBd_outcome) %>%
  mutate(maxBd_outcome = ifelse(maxBd_outcome==0,NA,maxBd_outcome)) %>%
  mutate(indivID = factor(indivID, levels=(unique(indivID))))

forLegend <- ggplot() +
  geom_linerange(data=agg_inhibRich %>% mutate(PABD = ifelse(PABD==1,"Yes","No")), aes(x=factor(PABD), ymin=indivSamp_lowerse, ymax=indivSamp_uperse, group=factor(indivID)),  position=position_dodge2(width=0.3)) +
  geom_point(data=agg_inhibRich %>%  mutate(PABD = ifelse(PABD==1,"Yes","No")), aes(x=factor(PABD), y=indivSamp,col=species, group=factor(indivID)),cex = 3, position=position_dodge2(width=0.3)) +
  coord_flip()+
  labs(col="Amphibian Species") +
  xlab("Became infected?") + ylab("Richness of putative inhibitory bacteria\n(Percentile rank within speciesxtime)") +
  theme(panel.grid.major.x  = element_blank()) +
  scale_color_manual(values=spCol)
g_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  legend
} 
gg_specieslegend <- g_legend(forLegend)
forLegend2a <- ggplot() +
  geom_linerange(data=agg_inhibRich %>% mutate(PABD = ifelse(PABD==1,"Yes","No")), aes(x=factor(PABD), ymin=indivSamp_lowerse, ymax=indivSamp_uperse, group=factor(indivID)),  position=position_dodge2(width=0.3)) +
  geom_point(data=agg_inhibRich %>%  mutate(PABD = ifelse(PABD==1,"Yes","No")), aes(x=factor(PABD), y=indivSamp,col=species, group=factor(indivID)),cex = 3, position=position_dodge2(width=0.3)) +
  coord_flip()+
  labs(col="Amphibian Species") +
  xlab("Became infected?") + ylab("Richness of putative inhibitory bacteria\n(Percentile rank within speciesxtime)") +
  theme(panel.grid.major.x  = element_blank(), legend.position = "bottom")+
  scale_color_manual(values=spCol)

gg_specieslegend2a <- g_legend(forLegend2a)


gg_inhibRich <- ggplot() +
  geom_point(data=agg_inhibRich %>%  mutate(PABD = ifelse(PABD==1,"Yes","No")), aes(x=factor(PABD), y=indivSamp,col=species, group=factor(indivID)),cex = 3, position=position_dodge2(width=0.3), show.legend = FALSE) +
  geom_linerange(data=agg_inhibRich %>% mutate(PABD = ifelse(PABD==1,"Yes","No")), aes(x=factor(PABD), ymin=indivSamp_lowerse, ymax=indivSamp_uperse, group=factor(indivID)),  position=position_dodge2(width=0.3)) +
  coord_flip()+
  labs(col="Amphibian Species") +
  xlab("Became infected?") + ylab("Richness of putative inhibitory bacteria\n(Percentile rank within speciesxtime)") +
  theme(panel.grid.major.x  = element_blank())+
  scale_color_manual(values=spCol)
gg_inhibRich
# ggsave(filename = paste0("4_Bayesian_Stats/UsedRanksPlots/custom_inhibRich.png"), height=4, width=6
#        ,gg_inhibRich)


gg_percInhib <- ggplot() +
  geom_point(data=agg_percInhib %>% filter(!is.na(maxBd_outcome)), aes(x=maxBd_outcome, y=indivSamp, col=species),cex = 3, show.legend = FALSE) +
  geom_linerange(data=agg_percInhib %>% filter(!is.na(maxBd_outcome)), aes(x=maxBd_outcome, ymin=indivSamp_lowerse, ymax=indivSamp_uperse)) +
  coord_flip()+
  labs(col="Amphibian Species") +
  xlab(expression(paste("Maximum ",italic("Bd")," load (log+1)"))) + ylab("Proportion putative inhibitory bacteria\n(Percentile rank within speciesxtime)") +
  theme(panel.grid.major.x  = element_blank(),) +
  scale_color_manual(values=spCol)
gg_percInhib

# ggsave(filename = paste0("4_Bayesian_Stats/UsedRanksPlots/custom_percInhib.png"), height=4, width=6
#        ,gg_percInhib)

gg_dist <- ggplot() +
  geom_point(data=agg_dist %>% filter(!is.na(maxBd_outcome)), aes(x=maxBd_outcome, y=indivSamp, col=species),cex = 3, show.legend = FALSE) +
  geom_linerange(data=agg_dist %>% filter(!is.na(maxBd_outcome)), aes(x=maxBd_outcome, ymin=indivSamp_lowerse, ymax=indivSamp_uperse)) +
  coord_flip()+
  labs(col="Amphibian Species") +
  xlab(expression(paste("Maximum ",italic("Bd")," load (log+1)"))) + ylab("Instability between time points\n(Percentile rank within speciesxtime)") +
  theme(panel.grid.major.x  = element_blank()) +
  scale_color_manual(values=spCol)
gg_dist
# ggsave(filename = paste0("4_Bayesian_Stats/UsedRanksPlots/custom_dist.png"), height=4, width=6
#        ,gg_dist)

gg_microbiome_forpanel <- allSamps_bdresponse_adj %>%
  filter(Predictor == "micro")  %>%
  mutate(MP = factor(MP, levels=c(allVarNames[allVarFilt]))) %>%
  mutate(Response = factor(Response, levels=c("PABD","maxBd_outcome"))) %>%
  mutate(sigCI = lower95*upper95>0) %>%
  mutate(PosteriorProbability = factor(ifelse(postprob>0.975,">97.5%", ifelse(postprob>0.95, ">95%", "<=94%")), levels=c(">97.5%", ">95%", "<=94%"))) %>%
  ggplot() + 
  geom_point(aes(x=MP, y=medsamp, col=PosteriorProbability, cex=PosteriorProbability), show.legend = FALSE) +
  geom_segment(aes(x=MP, xend = MP, y=lower95, yend=upper95, col=PosteriorProbability), show.legend = FALSE)+
  geom_hline(aes(yintercept=0), col="red",alpha=0.5)+
  facet_wrap(. ~ Response, ncol=1, strip.position = "left", labeller = as_labeller(c(maxBd_outcome = "Maximum future Bd load\n (log)", PABD = "Probability of infection\n(Presence/absence)"))) +
  # coord_flip() +
  labs(col="One-sided\nprobability of direction", cex="One-sided\nprobability of direction") +
  scale_color_manual(values=colorForSig) +
  scale_size_manual(values=c(3,2,1)) +
  xlab("Microbial community trait predictor") + ylab("Effect of microbial community trait on:")+
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5), legend.position = "right")

forLegend2b <- allSamps_bdresponse_adj %>%
  filter(Predictor == "micro")  %>%
  mutate(MP = factor(MP, levels=c(allVarNames[allVarFilt]))) %>%
  mutate(Response = factor(Response, levels=c("PABD","maxBd_outcome"))) %>%
  mutate(sigCI = lower95*upper95>0) %>%
  mutate(PosteriorProbability = factor(ifelse(postprob>0.975,">97.5%", ifelse(postprob>0.95, ">95%", "<=94%")), levels=c(">97.5%", ">95%", "<=94%"))) %>%
  ggplot() + 
  geom_point(aes(x=MP, y=medsamp, col=PosteriorProbability, cex=PosteriorProbability)) +
  geom_segment(aes(x=MP, xend = MP, y=lower95, yend=upper95, col=PosteriorProbability))+
  geom_hline(aes(yintercept=0), col="red",alpha=0.5)+
  facet_wrap(. ~ Response, ncol=1, strip.position = "left", labeller = as_labeller(c(maxBd_outcome = "Maximum future Bd load\n (log)", PABD = "Probability of infection\n(Presence/absence)"))) +
  # coord_flip() +
  labs(col="One-sided\nprobability of direction", cex="One-sided\nprobability of direction") +
  scale_color_manual(values=colorForSig) +
  scale_size_manual(values=c(3,2,1)) +
  xlab("Microbial community trait predictor") + ylab("Effect of microbial community trait on:")+
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5), legend.position = "bottom")
gg_specieslegend2b <- g_legend(forLegend2b)

### ONE GRID

# lay <- rbind(c(NA,5,1,1,4,4,4)
#              ,c(NA,5,1,1,4,4,4)
#              ,c(2,2,3,3,4,4,4)
#              ,c(2,2,3,3,4,4,4)
#              ,c(NA,NA,NA,NA,4,4,4)
# )
lay <- rbind(c(1,1,5,5,4,4,4)
             ,c(1,1,5,5,4,4,4)
             ,c(2,2,3,3,4,4,4)
             ,c(2,2,3,3,4,4,4)
             ,c(NA,NA,NA,NA,4,4,4)
)
lay2<- rbind(c(5,5,5,5,6,6,6)
             ,c(1,1,1,NA,4,4,4)
             ,c(1,1,1,NA,4,4,4)
             ,c(2,2,2,NA,4,4,4)
             ,c(2,2,2,NA,4,4,4)
             ,c(3,3,3,NA,4,4,4)
             ,c(3,3,3,NA,4,4,4)
             # ,c(NA,NA,NA,4,4,4)
)
leg_temp <- as_ggplot(gg_specieslegend)
leg_temp2a <- as_ggplot(gg_specieslegend2a)
leg_temp2b <- as_ggplot(gg_specieslegend2b)

# grid.arrange(test[[1]], test[[2]])
gg_all <- as_ggplot(grid.arrange(grobs=list(gg_inhibRich, gg_percInhib, gg_dist, gg_microbiome, leg_temp), layout_matrix=lay)) 
ggsave(filename = "4_Bayesian_Stats/BayesianSamplePosteriors/all_microbiome_effects.png", height=6, width=13
       ,gg_all)


gg_all2 <- as_ggplot(grid.arrange(grobs=list(gg_inhibRich, gg_percInhib, gg_dist, gg_microbiome_forpanel, leg_temp2a, leg_temp2b), layout_matrix=lay2)) 
ggsave(filename = "4_Bayesian_Stats/BayesianSamplePosteriors/all_microbiome_effects_v2.png", height=10, width=10
       ,gg_all2)
# dev.off()
#### @~~~~~~~~~~~~~~~Printing/saving~~~~~~~~~~~~~~~@ ####

allSamps_bdresponse_toWrite <- allSamps_bdresponse %>% 
  mutate(Predictor = factor(Predictor, levels=c("Anbo","Rhma","Osse","Raca","Rapi","micro"))) %>%
  # rowwise() %>% mutate(MP = allVarNames[MicrobiomePredictor]) %>%
  # mutate(MP = factor(MP, levels=rev(allVarNames))) %>%
  group_by(Predictor, MicrobiomePredictor, Response) %>%
  summarise(medsamp = median(samp), lower95 = ci(samp, method="HDI")$CI_low, upper95 = ci(samp, method="HDI")$CI_high, postprob = max(c(sum(samp>0), sum(samp<0)))/n()) %>% 
  ungroup() %>%
  mutate(sig = postprob>0.95) %>% arrange(Response)

toWriteCustomStats <- allSamps_bdresponse_toWrite %>%
  filter(Predictor == "micro") 
write.table(toWriteCustomStats, file="4_Bayesian_Stats/summarystats/microbiome_effects_stats.txt", sep="\t", row.names = FALSE, quote = FALSE)

#### @ ~~~~~~~~~~~~~~~~ Predictor correlations ~~~~~~~~~~~~~~~ @ ####

allRanks_raw <- allRanks %>% select(-prank) %>% pivot_wider(values_from=Observed, names_from=Response)%>%
  mutate(prepost = factor(prepost, levels=c("Pre","Post"))) %>%
  mutate(Bd_exposure = factor(Bd_exposure, levels=c("Control","Bd-exposed")))%>%
  mutate(inhibRich_stand = (inhibRich-mean(inhibRich, na.rm=TRUE))/sd(inhibRich, na.rm=TRUE)
         , observed_otus_stand = (log(observed_otus)-mean(log(observed_otus), na.rm=TRUE))/sd(log(observed_otus), na.rm=TRUE)
         , percInhib_stand = (percInhib-mean(percInhib, na.rm=TRUE))/sd(percInhib, na.rm=TRUE), disper_weighted_unifrac_stand =  (disper_weighted_unifrac-mean(disper_weighted_unifrac, na.rm=TRUE))/sd(disper_weighted_unifrac, na.rm=TRUE))

stan_richcorr_basic <- stan_glm(inhibRich_stand ~ observed_otus_stand, data=allRanks_raw, iter=5000)
# stan_inhibdisper_basic <- stan_glm(percInhib_stand ~ disper_weighted_unifrac_stand, data=allRanks_raw, iter=5000)
# stan_inhibcorr_basic <- stan_glm(percInhib_stand ~ inhibRich_stand, data=allRanks_raw, iter=5000)

sampsinhibrich <- rstan::extract(stan_richcorr_basic$stanfit)$beta
# sampsinhibdisper <- rstan::extract(stan_inhibdisper_basic$stanfit)$beta
# sampsinhibcorr <- rstan::extract(stan_inhibcorr_basic$stanfit)$beta

# Calculate R^2
getRsqBayes <- function(model) {
  predictedy <- as.data.frame(posterior_predict(model))
  residuals <- apply(predictedy,1,function(x) x-model$y )
  var_u <- apply(predictedy, 1, var)
  var_res <- apply(residuals, 2, var)
  rsq <- var_u/(var_u+var_res)
  return(rsq)
}

allCorr_summary <- as.data.frame(rbind(data.frame(metric = "inhibrich_vs_rich", samps= sampsinhibrich, R2 = getRsqBayes(stan_richcorr_basic))
                                       # ,data.frame(metric = "percinhib_vs_disper", samps= sampsinhibdisper, R2 = getRsqBayes(stan_inhibdisper_basic))
                                       # ,data.frame(metric = "percinhib_vs_inhibrich", samps= sampsinhibcorr, R2 = getRsqBayes(stan_inhibcorr_basic))
)) %>%
  mutate(samps=as.numeric(samps)) %>%
  group_by(metric) %>%
  summarise(slope = median(samps)
            , lower95 = ci(samps, method="HDI")$CI_low
            , upper95 = ci(samps, method="HDI")$CI_high
            , postprob = max(c(sum(samps<0), sum(samps>0)))/n()
            , R2_median = median(R2)
            , lower95_R2 = ci(R2, method="HDI")$CI_low
            , upper95_R2 = ci(R2, method="HDI")$CI_high)
allCorr_summary

gg_corrich <- allRanks_raw %>%
  # filter(prepost=="Pre", Bd_exposure =="Bd-exposed") %>%
  ggplot() + geom_point(aes(x=log(observed_otus), y=inhibRich, col=species), cex=2) +
  geom_smooth(aes(x=log(observed_otus), y=inhibRich), method="lm", col="black") +
  xlab("Overall bacterial richness (log)") + ylab("Richness of putative\nBd-inhibitory bacteria") + labs(col="Amphibian\nspecies")
ggsave("4_Bayesian_Stats/DataPlots/CORR_inhibRich_observedotus.png", height=4, width=6
       ,gg_corrich)


allRanks_prank %>%
  ggplot() + geom_point(aes(x=inhibRich, y=percInhib, col=species), cex=2) +
  geom_smooth(aes(x=inhibRich, y=percInhib), method="lm", col="black") +
  ylab("Proportion of putative\nBd-inhibitory bacteria") + xlab("Richness of putative\nBd-inhibitory bacteria") + labs(col="Amphibian\nspecies")

allRanks_prank %>%
  ggplot() + geom_point(aes(x=disper_weighted_unifrac, y=percInhib, col=species), cex=2) +
  geom_smooth(aes(x=disper_weighted_unifrac, y=percInhib), method="lm", col="black") +
  ylab("Proportion of putative\nBd-inhibitory bacteria") + xlab("Dispersion between individuals\n(by time point)") + labs(col="Amphibian\nspecies")


allRanks_prank <- allRanks %>% select(-Observed) %>% pivot_wider(values_from=prank, names_from=Response)%>%
  mutate(prepost = factor(prepost, levels=c("Pre","Post"))) %>%
  mutate(Bd_exposure = factor(Bd_exposure, levels=c("Control","Bd-exposed")))
allRanks_raw%>%
  mutate(prepost = factor(ifelse(prepost=="Pre","Pre-exposure", "Post-exposure"), levels=c("Pre-exposure","Post-exposure"))) %>%
  mutate(Bd_exposure = ifelse(Bd_exposure=="Control","Control (H-broth)", "Treatment (Bd solution)")) %>%
  ggplot(aes(x=observed_otus, y=inhibRich, col=species)) + geom_point(cex=2) + geom_smooth(method="lm", se=FALSE) +
  facet_grid(prepost~Bd_exposure) +ylab("Richness of putative\nBd-inhibitory bacteria") + xlab("Overall bacterial richness") +
  labs(col="Amphibian\nspecies")
if (!file.exists("4_Bayesian_Stats/stan_richcorr3.RData")) {
  stan_richcorr3 <- stan_glmer(inhibRich ~ observed_otus*species*prepost*Bd_exposure + (1|indivID)
                               , data=allRanks_raw
                               , iter=5000, adapt_delta = 0.99)
  save(stan_richcorr3, file = "4_Bayesian_Stats/stan_richcorr3.RData")
} else {
  load("4_Bayesian_Stats/stan_richcorr3.RData")
}

if (!file.exists("4_Bayesian_Stats/stan_richcorr.RData")) {
  stan_richcorr <- stan_glmer(inhibRich ~ observed_otus*species*prepost*Bd_exposure + (1|indivID)
                               , data=allRanks_prank
                               , iter=5000, adapt_delta = 0.99)
  save(stan_richcorr, file = "4_Bayesian_Stats/stan_richcorr.RData")
} else {
  load("4_Bayesian_Stats/stan_richcorr.RData")
}

allpred <- names(fixef(stan_richcorr))
allpred[grep("observed_otus", allpred)]
plot(stan_richcorr, pars=allpred[grep("observed_otus", allpred)])

samps <- rstan::extract(stan_richcorr$stanfit)
long_samps_renamed <- cbind(samps$alpha,samps$beta) %>% as.data.frame() %>% rename_at(vars(all_of(paste0("V",seq(1:(length(allpred)))))), ~all_of(allpred)) %>%
  rownames_to_column(var="iter") %>%
  pivot_longer(-iter, names_to="Coefficient", values_to="samp") %>%
  rowwise() %>%
  mutate(prepost = ifelse(length(grep("Post",Coefficient))>0, "Post", "Pre" )) %>%
  mutate(Bd_exposure = ifelse(length(grep("Bd-exposed",Coefficient))>0, "Treatment", "Control")) %>%
  mutate(Predictor = ifelse(length(grep("observed_otus",Coefficient))>0, "Slope", "Intercept")) %>%
  mutate(issp = ifelse(length(grep("species", Coefficient))>0, TRUE,FALSE)) %>%
  mutate(species = ifelse(issp, gsub(":.*$", "", gsub("^.*species","",Coefficient)), "Anbo")) %>% select(-issp) %>% ungroup()

long_samps_effectsadj <- long_samps_renamed %>% 
  filter(Predictor == "Slope") %>% 
  select(-Coefficient, -Predictor) %>% 
  # mutate(issp = ifelse(species=="Baseline","Baseline","Species")) %>%
  # pivot_wider(names_from=issp, values_from=samp)
  pivot_wider(names_from=Bd_exposure, values_from=samp) %>%
  pivot_wider(names_from=prepost, values_from=c("Control", "Treatment")) %>%
  pivot_wider(names_from=species, values_from = c("Control_Pre", "Control_Post", "Treatment_Pre", "Treatment_Post")) %>%
  mutate(ANBO_Control_Pre = Control_Pre_Anbo
         # , ANBO_Control_Pre = BASELINE_Control_Pre + Control_Pre_Anbo
         , RHMA_Control_Pre = ANBO_Control_Pre + Control_Pre_Rhma
         , OSSE_Control_Pre = ANBO_Control_Pre + Control_Pre_Osse
         , RACA_Control_Pre = ANBO_Control_Pre + Control_Pre_Raca
         , RAPI_Control_Pre = ANBO_Control_Pre + Control_Pre_Rapi
         , ANBO_Treatment_Pre = ANBO_Control_Pre + Treatment_Pre_Anbo
         # , ANBO_Treatment_Pre = Treatment_Pre_Baseline + ANBO_Control_Pre + Treatment_Pre_Anbo
         , RHMA_Treatment_Pre = ANBO_Treatment_Pre + RHMA_Control_Pre + Treatment_Pre_Rhma
         , OSSE_Treatment_Pre = ANBO_Treatment_Pre + OSSE_Control_Pre + Treatment_Pre_Osse
         , RACA_Treatment_Pre = ANBO_Treatment_Pre + RACA_Control_Pre + Treatment_Pre_Raca
         , RAPI_Treatment_Pre = ANBO_Treatment_Pre + RAPI_Control_Pre + Treatment_Pre_Rapi
         
         , ANBO_Control_Post = Control_Pre_Anbo + Control_Post_Anbo
         # , ANBO_Control_Post = BASELINE_Control_Post + Control_Pre_Anbo + Control_Post_Anbo
         , RHMA_Control_Post = ANBO_Control_Post + Control_Pre_Rhma + Control_Post_Rhma
         , OSSE_Control_Post = ANBO_Control_Post+ Control_Pre_Osse + Control_Post_Osse
         , RACA_Control_Post = ANBO_Control_Post + Control_Pre_Raca + Control_Post_Raca
         , RAPI_Control_Post = ANBO_Control_Post + Control_Pre_Rapi + Control_Post_Rapi
         , ANBO_Treatment_Post = ANBO_Control_Post + Treatment_Pre_Anbo + Treatment_Post_Anbo
         # , ANBO_Treatment_Post = ANBO_Control_Post + Treatment_Pre_Baseline + Treatment_Post_Baseline + Treatment_Post_Anbo
         , RHMA_Treatment_Post = RHMA_Control_Post + Treatment_Pre_Anbo + Treatment_Post_Anbo + Treatment_Post_Rhma
         , OSSE_Treatment_Post = OSSE_Control_Post + Treatment_Pre_Anbo + Treatment_Post_Anbo + Treatment_Post_Osse
         , RACA_Treatment_Post = RACA_Control_Post + Treatment_Pre_Anbo + Treatment_Post_Anbo + Treatment_Post_Raca
         , RAPI_Treatment_Post = RAPI_Control_Post + Treatment_Pre_Anbo + Treatment_Post_Anbo + Treatment_Post_Rapi
  ) %>% 
  select(starts_with("ANBO"), starts_with("RHMA"), starts_with("OSSE"), starts_with("RACA"), starts_with("RAPI")) %>%
  rownames_to_column(var="iter") %>%
  pivot_longer(-iter, names_to="Effect", values_to="samps") %>%
  separate(Effect, into = c("Species", "Treatment", "Prepost"), sep="_") 

# SUMMARISE EFFECTS
summary_effects_observedotus_on_inhibRich <- long_samps_effectsadj %>%
  group_by(Species, Treatment, Prepost) %>% 
  summarize(median = median(samps), lower95 = ci(samps,method="HDI")$CI_low, upper95 = ci(samps, method="HDI")$CI_high, postprob = max(c(sum(samps>0), sum(samps<0)))/n()) %>%
  mutate(sig = lower95*upper95>0)

summary_effects_observedotus_on_inhibRich %>%
  mutate(PosteriorProbability = ifelse(postprob>0.975,">97.5%", ifelse(postprob>0.95, ">95%", "<=94%"))) %>% 
  mutate(PosteriorProbability = factor(ifelse(is.na(PosteriorProbability), "Not applicable", PosteriorProbability), levels=c(">97.5%", ">95%", "<=94%","Not applicable"))) %>%
  ggplot() + geom_linerange(aes(x=Species, ymin=lower95, ymax=upper95, col=PosteriorProbability)) +
  geom_point(aes(x=Species, y=median, col=PosteriorProbability)) +
  scale_color_manual(values=c("blue","purple","black")) +
  labs(col="One-sided\nprobability of direction") +
  facet_grid(Prepost ~ Treatment) + ylab("Correlation coefficient between\noverall bacterial richness and putative Bd-inhibitory richness")



#### Correlation overall and richness- prank and control/treat/pre merged####

allRanks_prank_groups <- allRanks_prank %>% mutate(Group = ifelse(prepost=="Pre", "Pre-exposure", ifelse(Bd_exposure=="Control","H-broth","Bd-solution"))) %>%
  mutate(Group = factor(Group, levels=c("Pre-exposure", "H-broth", "Bd-solution")))

if (!file.exists("4_Bayesian_Stats/stan_richcorrstan_richcorr_merged.RData")) {
  stan_richcorr_merged <- stan_glmer(inhibRich ~ Group*species*observed_otus + (1|indivID)
                              , data=allRanks_prank_groups
                              , iter=5000, adapt_delta = 0.99)
  save(stan_richcorr_merged, file = "4_Bayesian_Stats/stan_richcorr_merged.RData")
} else {
  load("4_Bayesian_Stats/stan_richcorrstan_richcorr_merged.RData")
}

allpred <- names(fixef(stan_richcorr_merged))
allpred[grep("observed_otus", allpred)]
plot(stan_richcorr_merged, pars=allpred[grep("observed_otus", allpred)])

samps <- rstan::extract(stan_richcorr_merged$stanfit)
long_samps_renamed <- cbind(samps$alpha,samps$beta) %>% as.data.frame() %>% rename_at(vars(all_of(paste0("V",seq(1:(length(allpred)))))), ~all_of(allpred)) %>%
  rownames_to_column(var="iter") %>%
  pivot_longer(-iter, names_to="Coefficient", values_to="samp") %>%
  rowwise() %>%
  mutate(Predictor = ifelse(length(grep("observed_otus",Coefficient))>0, "Slope", "Intercept")) %>%
  mutate(prepost = ifelse(length(grep("GroupH-broth",Coefficient))>0, "H-broth", ifelse(length(grep("GroupBd-solution",Coefficient))>0, "Bd-solution", "Pre") )) %>%
  mutate(issp = ifelse(length(grep("species", Coefficient))>0, TRUE,FALSE)) %>%
  mutate(species = ifelse(issp, gsub(":.*$", "", gsub("^.*species","",Coefficient)), "Anbo")) %>% select(-issp) %>% ungroup()

long_samps_effectsadj <- long_samps_renamed %>% 
  filter(Predictor == "Slope") %>% 
  select(-Coefficient, -Predictor) %>% 
  # mutate(issp = ifelse(species=="Baseline","Baseline","Species")) %>%
  # pivot_wider(names_from=issp, values_from=samp)
  pivot_wider(names_from=species, values_from=samp) %>%
  pivot_wider(names_from=prepost, values_from=c("Anbo","Rhma","Osse","Raca","Rapi")) %>%
  # pivot_wider(names_from=species, values_from = c("Control_Pre", "Control_Post", "Treatment_Pre", "Treatment_Post")) %>%
  mutate(
    Pre_Anbo = Anbo_Pre
    , Pre_Rhma = Anbo_Pre + Rhma_Pre
    , Pre_Osse = Anbo_Pre + Osse_Pre
    , Pre_Raca = Anbo_Pre + Raca_Pre
    , Pre_Rapi = Anbo_Pre + Rapi_Pre
    , Hbroth_Anbo = Pre_Anbo + `Anbo_H-broth`
    , Hbroth_Rhma = Pre_Rhma + `Rhma_H-broth`
    , Hbroth_Osse = Pre_Osse + `Osse_H-broth`
    , Hbroth_Raca = Pre_Raca + `Raca_H-broth`
    , Hbroth_Rapi = Pre_Rapi + `Rapi_H-broth`
    , Bdsolution_Anbo = Pre_Anbo + `Anbo_Bd-solution`
    , Bdsolution_Rhma = Pre_Rhma + `Rhma_Bd-solution`
    , Bdsolution_Osse = Pre_Osse + `Osse_Bd-solution`
    , Bdsolution_Raca = Pre_Raca + `Raca_Bd-solution`
    , Bdsolution_Rapi = Pre_Rapi + `Rapi_Bd-solution`
  ) %>% 
  select(starts_with("Pre"), starts_with("Hbroth"), starts_with("Bdsolution")) %>%
  rownames_to_column(var="iter") %>%
  pivot_longer(-iter, names_to="Effect", values_to="samps") %>%
  separate(Effect, into = c("TreatmentGroup", "Species"), sep="_") 

# SUMMARISE EFFECTS
summary_effects_observedotus_on_inhibRich <- long_samps_effectsadj %>%
  group_by(Species, TreatmentGroup) %>% 
  summarize(median = median(samps), lower95 = ci(samps,method="HDI")$CI_low, upper95 = ci(samps, method="HDI")$CI_high, postprob = max(c(sum(samps>0), sum(samps<0)))/n()) %>%
  mutate(sig = lower95*upper95>0)

summary_effects_observedotus_on_inhibRich %>%
  mutate(PosteriorProbability = ifelse(postprob>0.975,">97.5%", ifelse(postprob>0.95, ">95%", "<=94%"))) %>% 
  mutate(PosteriorProbability = factor(ifelse(is.na(PosteriorProbability), "Not applicable", PosteriorProbability), levels=c(">97.5%", ">95%", "<=94%","Not applicable"))) %>%
  ggplot() + geom_linerange(aes(x=Species, ymin=lower95, ymax=upper95, col=PosteriorProbability)) +
  geom_point(aes(x=Species, y=median, col=PosteriorProbability)) +
  scale_color_manual(values=c("blue","purple","black")) +
  labs(col="One-sided\nprobability of direction") +
  facet_grid(. ~ TreatmentGroup) + ylab("Correlation coefficient between\noverall bacterial richness and putative Bd-inhibitory richness")



#### Correlation overall and richness- raw and control/treat/pre merged STANDARDIZED ####

allRanks_raw_groups <- allRanks_raw %>% mutate(Group = ifelse(prepost=="Pre", "Pre-exposure", ifelse(Bd_exposure=="Control","H-broth","Bd-solution"))) %>%
  mutate(Group = factor(Group, levels=c("Pre-exposure", "H-broth", "Bd-solution"))) %>%
  mutate(inhibRich_stand = (inhibRich-mean(inhibRich, na.rm=TRUE))/sd(inhibRich, na.rm=TRUE), observed_otus_stand = (observed_otus-mean(observed_otus, na.rm=TRUE))/sd(observed_otus, na.rm=TRUE) )

if (!file.exists("4_Bayesian_Stats/stan_richcorrstan_richcorr_merged_raw.RData")) {
  stan_richcorr_merged_raw <- stan_glmer(inhibRich_stand ~ Group*species*observed_otus_stand + (1|indivID)
                                     , data=allRanks_raw_groups
                                     , iter=5000, adapt_delta = 0.99)
  save(stan_richcorr_merged_raw, file = "4_Bayesian_Stats/stan_richcorr_merged_raw.RData")
  save(stan_richcorr_basic, file = "4_Bayesian_Stats/stan_richcorr_basic.RData")
} else {
  load("4_Bayesian_Stats/stan_richcorrstan_richcorr_merged_raw.RData")
  load("4_Bayesian_Stats/stan_richcorr_basic.RData")
}

allpred <- names(fixef(stan_richcorr_merged_raw))
allpred[grep("observed_otus", allpred)]
plot(stan_richcorr_merged_raw, pars=allpred[grep("observed_otus", allpred)])


samps <- rstan::extract(stan_richcorr_merged_raw$stanfit)
long_samps_renamed <- cbind(samps$alpha,samps$beta) %>% as.data.frame() %>% rename_at(vars(all_of(paste0("V",seq(1:(length(allpred)))))), ~all_of(allpred)) %>%
  rownames_to_column(var="iter") %>%
  pivot_longer(-iter, names_to="Coefficient", values_to="samp") %>%
  rowwise() %>%
  mutate(Predictor = ifelse(length(grep("observed_otus",Coefficient))>0, "Slope", "Intercept")) %>%
  mutate(prepost = ifelse(length(grep("GroupH-broth",Coefficient))>0, "H-broth", ifelse(length(grep("GroupBd-solution",Coefficient))>0, "Bd-solution", "Pre") )) %>%
  mutate(issp = ifelse(length(grep("species", Coefficient))>0, TRUE,FALSE)) %>%
  mutate(species = ifelse(issp, gsub(":.*$", "", gsub("^.*species","",Coefficient)), "Anbo")) %>% select(-issp) %>% ungroup()

long_samps_effectsadj <- long_samps_renamed %>% 
  filter(Predictor == "Slope") %>% 
  select(-Coefficient, -Predictor) %>% 
  # mutate(issp = ifelse(species=="Baseline","Baseline","Species")) %>%
  # pivot_wider(names_from=issp, values_from=samp)
  pivot_wider(names_from=species, values_from=samp) %>%
  pivot_wider(names_from=prepost, values_from=c("Anbo","Rhma","Osse","Raca","Rapi")) %>%
  # pivot_wider(names_from=species, values_from = c("Control_Pre", "Control_Post", "Treatment_Pre", "Treatment_Post")) %>%
  mutate(
    Pre_Anbo = Anbo_Pre
    , Pre_Rhma = Anbo_Pre + Rhma_Pre
    , Pre_Osse = Anbo_Pre + Osse_Pre
    , Pre_Raca = Anbo_Pre + Raca_Pre
    , Pre_Rapi = Anbo_Pre + Rapi_Pre
    , Hbroth_Anbo = Pre_Anbo + `Anbo_H-broth`
    , Hbroth_Rhma = Pre_Rhma + `Rhma_H-broth`
    , Hbroth_Osse = Pre_Osse + `Osse_H-broth`
    , Hbroth_Raca = Pre_Raca + `Raca_H-broth`
    , Hbroth_Rapi = Pre_Rapi + `Rapi_H-broth`
    , Bdsolution_Anbo = Pre_Anbo + `Anbo_Bd-solution`
    , Bdsolution_Rhma = Pre_Rhma + `Rhma_Bd-solution`
    , Bdsolution_Osse = Pre_Osse + `Osse_Bd-solution`
    , Bdsolution_Raca = Pre_Raca + `Raca_Bd-solution`
    , Bdsolution_Rapi = Pre_Rapi + `Rapi_Bd-solution`
  ) %>% 
  select(starts_with("Pre"), starts_with("Hbroth"), starts_with("Bdsolution")) %>%
  rownames_to_column(var="iter") %>%
  pivot_longer(-iter, names_to="Effect", values_to="samps") %>%
  separate(Effect, into = c("TreatmentGroup", "Species"), sep="_") 

# SUMMARISE EFFECTS
summary_effects_observedotus_on_inhibRich <- long_samps_effectsadj %>%
  group_by(Species, TreatmentGroup) %>% 
  summarize(median = median(samps), lower95 = ci(samps,method="HDI")$CI_low, upper95 = ci(samps, method="HDI")$CI_high, postprob = max(c(sum(samps>0), sum(samps<0)))/n()) %>%
  mutate(sig = lower95*upper95>0)

summary_effects_observedotus_on_inhibRich %>%
  mutate(PosteriorProbability = ifelse(postprob>0.975,">97.5%", ifelse(postprob>0.95, ">95%", "<=94%"))) %>% 
  mutate(PosteriorProbability = factor(ifelse(is.na(PosteriorProbability), "Not applicable", PosteriorProbability), levels=c(">97.5%", ">95%", "<=94%","Not applicable"))) %>%
  ggplot() + geom_linerange(aes(x=Species, ymin=lower95, ymax=upper95, col=PosteriorProbability)) +
  geom_point(aes(x=Species, y=median, col=PosteriorProbability)) +
  scale_color_manual(values=c("blue","purple","black")) +
  labs(col="One-sided\nprobability of direction") +
  facet_grid(. ~ TreatmentGroup) + ylab("Correlation coefficient between\noverall bacterial richness and putative Bd-inhibitory richness")







allRanks_summarisedall_prank %>%
  ggplot() + 
  geom_segment(aes(x=m_observed_otus-s_observed_otus, xend=m_observed_otus+s_observed_otus, y=m_inhibRich, yend=m_inhibRich)) +  
  geom_segment(aes(y=m_inhibRich-s_inhibRich, yend=m_inhibRich+s_inhibRich, x=m_observed_otus, xend=m_observed_otus)) +
  geom_point(aes(x=m_observed_otus, y=m_inhibRich, col=species), cex=3) +
  facet_grid(prepost ~ Bd_exposure)



allRanks_summarisedall%>%
  ggplot() + 
  geom_segment(aes(x=m_observed_otus-s_observed_otus, xend=m_observed_otus+s_observed_otus, y=m_inhibRich, yend=m_inhibRich)) +  
  geom_segment(aes(y=m_inhibRich-s_inhibRich, yend=m_inhibRich+s_inhibRich, x=m_observed_otus, xend=m_observed_otus)) +
  geom_point(aes(x=m_observed_otus, y=m_inhibRich, col=species), cex=3) +
  facet_grid(prepost ~ Bd_exposure)
  
allRanks_summarisedall %>%
  filter(Bd_exposure=="Control" ) %>%
  ggplot() + 
  geom_segment(aes(x=m_observed_otus-s_observed_otus, xend=m_observed_otus+s_observed_otus, y=m_inhibRich, yend=m_inhibRich)) +  
  geom_segment(aes(y=m_inhibRich-s_inhibRich, yend=m_inhibRich+s_inhibRich, x=m_observed_otus, xend=m_observed_otus)) +
  geom_point(aes(x=m_observed_otus, y=m_inhibRich, col=species), cex=3) +
  geom_smooth(aes(x=m_observed_otus, y=m_inhibRich, col=species), method="lm")

allRanks_summarisedall%>%
  filter(Bd_exposure=="Control" ) %>%
  ggplot() + 
  geom_segment(aes(x=m_observed_otus-s_observed_otus, xend=m_observed_otus+s_observed_otus, y=m_inhibRich, yend=m_inhibRich)) +  
  geom_segment(aes(y=m_inhibRich-s_inhibRich, yend=m_inhibRich+s_inhibRich, x=m_observed_otus, xend=m_observed_otus)) +
  geom_point(aes(x=m_observed_otus, y=m_inhibRich, col=species), cex=3) +
  geom_smooth(aes(x=m_observed_otus, y=m_inhibRich, col=species), method="lm")


allRanks_summarisedall %>%
  filter(prepost=="Pre" ) %>%
  ggplot() + 
  geom_segment(aes(x=m_observed_otus-s_observed_otus, xend=m_observed_otus+s_observed_otus, y=m_inhibRich, yend=m_inhibRich)) +  
  geom_segment(aes(y=m_inhibRich-s_inhibRich, yend=m_inhibRich+s_inhibRich, x=m_observed_otus, xend=m_observed_otus)) +
  geom_point(aes(x=m_observed_otus, y=m_inhibRich, col=species), cex=3) +
  geom_smooth(aes(x=m_observed_otus, y=m_inhibRich, col=species), method="lm")

allRanks %>% group_by(species, indivID, prepost, Bd_exposure, Response) %>%
  filter(prepost=="Pre" ) %>%
  summarize(m=mean(prank, na.rm=TRUE), s=sd(prank,na.rm=TRUE)/sqrt(n())) %>%
  pivot_wider(names_from=Response, values_from=c(m, s)) %>%
  ggplot() + 
  geom_segment(aes(x=m_observed_otus-s_observed_otus, xend=m_observed_otus+s_observed_otus, y=m_inhibRich, yend=m_inhibRich)) +  
  geom_segment(aes(y=m_inhibRich-s_inhibRich, yend=m_inhibRich+s_inhibRich, x=m_observed_otus, xend=m_observed_otus)) +
  geom_point(aes(x=m_observed_otus, y=m_inhibRich, col=species), cex=3) +
  geom_smooth(aes(x=m_observed_otus, y=m_inhibRich, col=species), method="lm", se = FALSE)


allRanks %>% group_by(species, indivID, prepost, Bd_exposure, Response) %>%
  filter(Bd_exposure=="Bd_exposed" | prepost=="Pre") %>%
  summarize(m=mean(Observed, na.rm=TRUE), s=sd(Observed,na.rm=TRUE)/sqrt(n())) %>%
  pivot_wider(names_from=Response, values_from=c(m, s)) %>%
  ggplot() + 
  geom_segment(aes(x=m_observed_otus-s_observed_otus, xend=m_observed_otus+s_observed_otus, y=m_inhibRich, yend=m_inhibRich)) +  
  geom_segment(aes(y=m_inhibRich-s_inhibRich, yend=m_inhibRich+s_inhibRich, x=m_observed_otus, xend=m_observed_otus)) +
  geom_point(aes(x=m_observed_otus, y=m_inhibRich, col=species), cex=3) +
  geom_smooth(aes(x=m_observed_otus, y=m_inhibRich, col=species), method="lm")

allRanks %>% group_by(species, indivID, prepost, Bd_exposure, Response) %>%
  filter(Bd_exposure=="Bd-exposed" | prepost=="Pre") %>%
  summarize(m=mean(prank, na.rm=TRUE), s=sd(prank,na.rm=TRUE)/sqrt(n())) %>%
  pivot_wider(names_from=Response, values_from=c(m, s)) %>%
  ggplot() + 
  geom_segment(aes(x=m_observed_otus-s_observed_otus, xend=m_observed_otus+s_observed_otus, y=m_inhibRich, yend=m_inhibRich)) +  
  geom_segment(aes(y=m_inhibRich-s_inhibRich, yend=m_inhibRich+s_inhibRich, x=m_observed_otus, xend=m_observed_otus)) +
  geom_point(aes(x=m_observed_otus, y=m_inhibRich, col=species), cex=3) +
  geom_smooth(aes(x=m_observed_otus, y=m_inhibRich, col=species), method="lm")
