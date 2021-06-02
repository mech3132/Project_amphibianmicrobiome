#!bin/bash

##### Bayesian statistics ######
# library(lme4) # for lmer
# library(rstanarm) # For bayesian estimates of alpha and beta diversity
# library(car) # for Anova
# library(MASS) # for fitdistr
library(bayestestR) # for CI
# library(brms) # for modelling
library(rstan)
library(tidyverse) # for data manipulation
library(rstanarm)
# library(rstanarm)
##### Load data ######
load("./3_5sp_mapping_otu_downstream/mf_alt_filt_final.RData")
# load("./3_5sp_mapping_otu_downstream/otu_filt.RData")
# load("./3_5sp_mapping_otu_downstream/otu_filt_inhibOnly.RData")
# load("./3_5sp_mapping_otu_downstream/braycurtis_filt.RData")
# load("./3_5sp_mapping_otu_downstream/unweighted_unifrac_filt.RData")
# load("./3_5sp_mapping_otu_downstream/weighted_unifrac_filt.RData")
# load("./3_5sp_mapping_otu_downstream/alpha_metrics.RData")
# load("./3_5sp_mapping_otu_downstream/beta_metrics.RData")
load("./4_Bayesian_Stats/all_stan_output_baseline.RData")
load("./4_Bayesian_Stats/allRanks.RData")

# setwd("..")
dir.create("4_Bayesian_CustomStan")
setwd("4_Bayesian_CustomStan")
###### Set up mf ######
mf_con <- mf_alt_filt_final %>%
  filter(Bd_exposure == "Control")

mf_treat <- mf_alt_filt_final %>%
  filter(Bd_exposure == "Bd-exposed")

mf_pre <- mf_alt_filt_final %>%
  filter(prepost=="Pre")

mf_nonbd <- mf_treat %>% filter(prepost == "Pre") %>% full_join(mf_con)

mf_collapsed_bd <- mf_alt_filt_final %>%
  filter(Bd_exposure == "Bd-exposed") %>% 
  group_by(species, indivID, prepost) %>%
  summarize(maxBd = max(Bd_load, na.rm = TRUE), medBd = median(Bd_load, na.rm = TRUE), fracBd = sum(PABD, na.rm=FALSE)/n()) %>%
  group_by(species, indivID) %>% mutate(maxBd_outcome = max(maxBd, na.rm=TRUE), fracBd_outcome = max(fracBd, na.rm=TRUE)) %>%
  filter(prepost=="Pre") %>%
  mutate(PABD = ifelse(fracBd_outcome>0,1,0)) %>%
  select(species, indivID, maxBd_outcome, fracBd_outcome, PABD) %>% #rename(Species = species) %>%
  ungroup()

allVar <- c("observed_otus","shannon"
            ,"inhibRich","percInhib"
            ,"dist_braycurtis","dist_unweighted_unifrac","dist_weighted_unifrac"
            ,"disper_braycurtis","disper_unweighted_unifrac","disper_weighted_unifrac")
allVarNames <- c("ASV richness", "Shannon diversity"
                 , "Richness of \nputative inhibitory bacteria", "Proportion \nputative inhibitory bacteria"
                 , "Turnover through time \n(Bray-curtis)", "Turnover through time \n(unweighted Unifrac)", "Turnover through time \n(weighted Unifrac)"
                 , "Dispersion from centroid \n(Bray-curtis)", "Dispersion from centroid \n(unweighted Unifrac)", "Dispersion from centroid \n(weighted Unifrac)")
names(allVarNames) <- allVar

#### Find sd of ranks ####

allRanks_summary <- allRanks %>%
  filter(prepost=="Pre", Bd_exposure == "Bd-exposed") %>%
  group_by(species, indivID, Response) %>%
  summarise(medrank = median(prank, na.rm = TRUE), sdrank = sd(prank, na.rm=TRUE), serank = sd(prank, na.rm=TRUE)/sqrt(n()), nindiv=n()) %>% ungroup() %>%
  left_join(mf_collapsed_bd) 

### Fix Osse? Missing individual-level variation because timelines are crunched together.
maxOsseSD <- allRanks %>%
  filter(prepost=="Pre", Bd_exposure == "Bd-exposed", Response==var_temp) %>%
  filter(species=="Osse") %>% group_by(indivID) %>%
  summarise(sdrank = sd(prank, na.rm=TRUE)) %>% pull(sdrank) %>% max(na.rm=TRUE)

allRanks_summary <- allRanks_summary %>% mutate(sdrank = ifelse(is.na(sdrank) & species=="Osse", maxOsseSD, sdrank)) %>%
  mutate(serank = ifelse(is.na(serank) & species=="Osse", sdrank/sqrt(nindiv), serank))
###### CHOSEN MODELS
chosen_models <- c(observed_otus = "lognormal", chao1 = "lognormal", shannon = "normal"
                   , disper_braycurtis = "lognormal", disper_unweighted_unifrac = "lognormal", disper_weighted_unifrac = "lognormal"
                   , dist_braycurtis = "beta", dist_unweighted_unifrac = "beta", dist_weighted_unifrac = "beta"
                   , percInhib = "beta", inhibRich = "poisson")

#### Custom Stan Models with uncertainty ####

write("//Stan bernoulli model with uncertainty
data {
  int<lower=0> N;       // number of cases
  int<lower=0, upper=1> y[N];          // outcome (variate)
  vector[N] x_meas;     // measurement of x
  vector[N] tau;       // sd of measurment of x
  vector[N] x_Anbo;       // species
  vector[N] x_Rhma;       // species
  vector[N] x_Osse;       // species
  vector[N] x_Raca;       // species
  vector[N] x_Rapi;       // species
}
parameters {
  real b_Anbo;            // slope for species
  real b_Rhma;            // slope for species
  real b_Osse;            // slope for species
  real b_Raca;            // slope for species
  real b_Rapi;            // slope for species
  real b_micro;          // slope for microbiome predictor
  vector[N] x;          // unknown true value
}
model {
  b_Anbo ~ normal(0, 10);
  b_Rhma ~ normal(0, 10);
  b_Osse ~ normal(0, 10);
  b_Raca ~ normal(0, 10);
  b_Rapi ~ normal(0, 10);
  b_micro ~ normal(0,10);
  x ~ normal(0, 10);  // prior
  x_meas ~ normal(x, tau);    // measurement model
  y ~ bernoulli_logit(b_Anbo * x_Anbo + b_Rhma * x_Rhma + b_Osse * x_Osse + b_Raca * x_Raca+ b_Rapi * x_Rapi + b_micro * x);
  
}",
"stan_model_binomial_ME.stan")
# stanc("stan_model_binomial_ME.stan")
modelstan_binME <- "stan_model_binomial_ME.stan"


write("//Stan normal model with uncertainty
data {
  int<lower=0> N;       // number of cases
  vector[N] y;          // outcome (variate)
  vector[N] x_meas;     // measurement of x
  vector[N] tau;       // sd of measurment of x
  vector[N] x_Anbo;       // species
  vector[N] x_Rhma;       // species
  vector[N] x_Osse;       // species
  vector[N] x_Raca;       // species
  vector[N] x_Rapi;       // species
}
parameters {
  real b_Anbo;            // slope for species
  real b_Rhma;            // slope for species
  real b_Osse;            // slope for species
  real b_Raca;            // slope for species
  real b_Rapi;            // slope for species
  real b_micro;          // slope for microbiome predictor
  vector[N] x;          // unknown true value
  real<lower=0> sigma;  // sample-level variation
}
model {
  b_Anbo ~ normal(0, 10);
  b_Rhma ~ normal(0, 10);
  b_Osse ~ normal(0, 10);
  b_Raca ~ normal(0, 10);
  b_Rapi ~ normal(0, 10);
  b_micro ~ normal(0,10);
  sigma ~ cauchy(0, 5);
  x ~ normal(0, 10);  // prior
  x_meas ~ normal(x, tau);    // measurement model
  y ~ normal(b_Anbo * x_Anbo + b_Rhma * x_Rhma + b_Osse * x_Osse + b_Raca * x_Raca + b_Rapi * x_Rapi + b_micro * x, sigma);
  
}",
"stan_model_normal_ME.stan")
# stanc("stan_model_normal_ME.stan")
modelstan_normME <- "stan_model_normal_ME.stan"



write("//Stan normal model with uncertainty NO RAPI
data {
  int<lower=0> N;       // number of cases
  real y[N];          // outcome (variate)
  vector[N] x_meas;     // measurement of x
  vector[N] tau;       // sd of measurment of x
  vector[N] x_Anbo;       // species
  vector[N] x_Rhma;       // species
  vector[N] x_Osse;       // species
  vector[N] x_Raca;       // species
  //vector[N] x_Rapi;       // species
}
parameters {
  real b_Anbo;            // slope for species
  real b_Rhma;            // slope for species
  real b_Osse;            // slope for species
  real b_Raca;            // slope for species
  //real b_Rapi;            // slope for species
  real b_micro;          // slope for microbiome predictor
  vector[N] x;          // unknown true value
  real<lower=0> sigma;  // sample-level variation
}
model {
  b_Anbo ~ normal(0, 10);
  b_Rhma ~ normal(0, 10);
  b_Osse ~ normal(0, 10);
  b_Raca ~ normal(0, 10);
  //b_Rapi ~ normal(0, 10);
  b_micro ~ normal(0,10);
  sigma ~ cauchy(0, 5);
  x ~ normal(0, 10);  // prior
  x_meas ~ normal(x, tau);    // measurement model
  y ~ normal(b_Anbo * x_Anbo + b_Rhma * x_Rhma + b_Osse * x_Osse + b_Raca * x_Raca + b_micro * x, sigma);
  
}",
"stan_model_normal_ME_norapi.stan")
# stanc("stan_model_normal_ME_norapi.stan")
modelstan_normME_norapi <- "stan_model_normal_ME_norapi.stan"




write("//Stan normal model with uncertainty NO RAPI
data {
  int<lower=0> N;       // number of cases
  real<lower=0> y[N];          // outcome (variate)
  //vector[N] y;          // outcome (variate)
  vector[N] ysd;        //outcome uncertainty
  vector[N] x_meas;     // measurement of x
  vector[N] tau;       // sd of measurment of x
  vector[N] x_Anbo;       // species
  vector[N] x_Rhma;       // species
  vector[N] x_Osse;       // species
  vector[N] x_Raca;       // species
  //vector[N] x_Rapi;       // species
}
parameters {
  real b_Anbo;            // slope for species
  real b_Rhma;            // slope for species
  real b_Osse;            // slope for species
  real b_Raca;            // slope for species
  //real b_Rapi;            // slope for species
  real b_micro;          // slope for microbiome predictor
  vector[N] x;          // unknown true value
  //real<lower=0> sigma;  // sample-level variation
  //vector[N] y2;
}
model {
  b_Anbo ~ normal(0, 10);
  b_Rhma ~ normal(0, 10);
  b_Osse ~ normal(0, 10);
  b_Raca ~ normal(0, 10);
  //b_Rapi ~ normal(0, 10);
  b_micro ~ normal(0,10);
  //sigma ~ cauchy(0, 5);
  x ~ normal(0, 10);  // prior
  x_meas ~ normal(x, tau);    // measurement model
  y ~ normal(b_Anbo * x_Anbo + b_Rhma * x_Rhma + b_Osse * x_Osse + b_Raca * x_Raca + b_micro * x, ysd);
  //y ~ normal(y2, ysd);
  
}",
"stan_model_normal_ME_norapi_ysd.stan")
stanc("stan_model_normal_ME_norapi_ysd.stan")
modelstan_normME_norapi_ysd <- "stan_model_normal_ME_norapi_ysd.stan"


write("//Stan normal model with uncertainty NO RAPI with species interaction
data {
  int<lower=0> N;       // number of cases
  int<lower=1> K;       // number of species
  real y[N];          // outcome (variate)
  vector[N] x_meas;     // measurement of x
  vector[N] tau;       // sd of measurment of x
  matrix[N,K] x_species; //all species data matrix
  //vector[N] x_Anbo;       // species
  //vector[N] x_Rhma;       // species
  //vector[N] x_Osse;       // species
  //vector[N] x_Raca;       // species

}
parameters {
  vector[K] b_species;
  //real b_Anbo;            // slope for species
  //real b_Rhma;            // slope for species
  //real b_Osse;            // slope for species
  //real b_Raca;            // slope for species
  vector[K] r_micro;
  //real b_Anbo_inter;            // slope for interaction
  //real b_Rhma_inter;            // slope for interaction
  //real b_Osse_inter;            // slope for interaction
  //real b_Raca_inter;            // slope for interaction
  real b_micro;          // overall slope for microbiome predictor
  vector[N] x;          // unknown true value
  real<lower=0> sigma;  // sample-level variation
  real<lower=0> ksigma; // random effect variation
}
model {
  b_micro ~ normal(0,10);
  sigma ~ cauchy(0, 5);
  ksigma ~ cauchy(0,5);
  for ( k in 1:K ) 
      r_micro[k] ~ normal(b_micro, ksigma);
  
  x ~ normal(0, 10);  // prior
  x_meas ~ normal(x, tau);    // measurement model
  y ~ normal(x_species*b_species + x_species*r_micro + b_micro * x, sigma);
  
}",
"stan_model_normal_ME_norapi_interact.stan")
stanc("stan_model_normal_ME_norapi_interact.stan")
modelstan_normME_norapi_interact <- "stan_model_normal_ME_norapi_interact.stan"


##### Run iterative models #####
allVarFilt <- c("observed_otus","shannon","inhibRich","percInhib","dist_weighted_unifrac","disper_weighted_unifrac")

indivList_treat <- unique((mf_treat$indivID))
if ( !file.exists("allModels_predictBd.RData") ) {
  # set.seed(93845)
  set.seed(105293)
  
  total <- length(allVarFilt)*2
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  i <- 1
  allModels_predictBd <- list()
  aggregatedData_baseline <- data.frame()
  for (bd_type in c("PABD", "maxBd_outcome")) {
    # bd_type <- "maxBd_outcome"
    allModels_predictBd[[paste0(bd_type)]] <- list()
    
    for (var_temp in allVarFilt) {
      # var_temp <- allVar[8]
      #tempBaseline <- all_stan_output_baseline[[paste0(var_temp)]]
      if ( bd_type == "PABD") {
        # data_stan_temp <- tempBaseline$stan_summary %>% as.data.frame() %>% rownames_to_column(var="Coef") %>%
        #   rowwise() %>% filter(length(grep("b[(Intercept) indivID:", Coef, fixed=TRUE))>0) %>%
        #   mutate(Coef = gsub("b[(Intercept) indivID:", "", Coef, fixed = TRUE)) %>%
        #   mutate(Coef = gsub("]","", Coef, fixed = TRUE)) %>% 
        #   filter(Coef %in% indivList_treat) %>%
        #   select(Coef, mean, se_mean, sd) %>%
        #   rename(indivID=Coef, x_meas_raw = mean, seMean_raw = se_mean, sdMean_raw = sd) %>%
        #   separate(indivID, into=c("x_sp", "indiv"), remove=FALSE, sep="_") %>%
        #   full_join(mf_collapsed_bd) %>%
        #   rename_at(vars(paste0(bd_type)), ~c("y")) %>%
        #   filter(!is.na(y)) %>%
        #   #mutate( y = (y_raw-mean(y_raw))/sd(y_raw)) %>%
        #   mutate(x_meas = (x_meas_raw-mean(x_meas_raw))/sd(x_meas_raw)) %>%
        #   mutate(sdMean = (sdMean_raw/sd(x_meas_raw)))
        data_stan_temp <- allRanks_summary %>% filter(Response==var_temp) %>%
          rename_at(vars(paste0(bd_type)), ~c("y")) %>%
          filter(!is.na(y)) %>% 
          mutate(x_meas = medrank) %>%
          mutate(x_meas = (medrank-mean(medrank))/sd(medrank)) %>%
          mutate(sdMean = serank/sd(medrank))
          # mutate(sdMean = serank)
        # View(data_stan_temp)
        modelUse <- "stan_model_binomial_ME.stan"
        ## Get data formatted
        data_stan_format <- list(N=nrow(data_stan_temp)
                                 , y=data_stan_temp$y
                                 # , x_meas = data_stan_temp$x_meas
                                 # , tau = data_stan_temp$sdMean                                 , x_meas = data_stan_temp$x_meas
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
        # ranks and uncertainty
        # data_stan_temp <- tempBaseline$stan_summary %>% as.data.frame() %>% rownames_to_column(var="Coef") %>%
        #   rowwise() %>% filter(length(grep("b[(Intercept) indivID:", Coef, fixed=TRUE))>0) %>%
        #   mutate(Coef = gsub("b[(Intercept) indivID:", "", Coef, fixed = TRUE)) %>%
        #   mutate(Coef = gsub("]","", Coef, fixed = TRUE)) %>% 
        #   filter(Coef %in% indivList_treat) %>%
        #   select(Coef, mean, se_mean, sd) %>%
        #   rename(indivID=Coef, x_meas_raw = mean, seMean_raw = se_mean, sdMean_raw = sd) %>%
        #   separate(indivID, into=c("x_sp", "indiv"), remove=FALSE, sep="_") %>%
        #   full_join(mf_collapsed_bd) %>%
        #   rename_at(vars(paste0(bd_type)), ~c("y_raw")) %>%
        #   filter(!is.na(y_raw)) %>% filter(y_raw>0) %>%
        #   mutate( y = (y_raw-mean(y_raw))/sd(y_raw)) %>%
        #   mutate(x_meas = (x_meas_raw-mean(x_meas_raw))/sd(x_meas_raw)) %>%
        #   mutate(sdMean = (sdMean_raw/sd(x_meas_raw)))
        # 
        # var_temp = "dist_braycurtis"
        data_stan_temp <- allRanks_summary %>% filter(Response==var_temp) %>%
          rename_at(vars(paste0(bd_type)), ~c("y_raw")) %>%
          filter(!is.na(y_raw)) %>% filter(y_raw>0) %>%
          filter(species!="Rapi") %>%
          # mutate( y = (y_raw-mean(y_raw))/sd(y_raw)) %>%
          mutate( y = y_raw) %>%
          mutate(x_meas = (medrank-mean(medrank))/sd(medrank)) %>%
          # mutate(x_meas = medrank) %>%
          mutate(sdMean = serank/sd(medrank))
          # mutate(sdMean = serank)
        # modelUse <- "stan_model_normal_ME.stan"
        modelUse <- "stan_model_normal_ME_norapi.stan"
        # View(data_stan_temp)
        ## Get data formatted
        data_stan_format <- list(N=nrow(data_stan_temp)
                                 , y=data_stan_temp$y
                                 # , x_meas = data_stan_temp$x_meas
                                 # , tau = data_stan_temp$sdMean
                                 , x_meas = data_stan_temp$medrank
                                 , tau = data_stan_temp$serank
                                 , x_Anbo = as.numeric(data_stan_temp$species=="Anbo")
                                 , x_Rhma = as.numeric(data_stan_temp$species=="Rhma")
                                 , x_Osse = as.numeric(data_stan_temp$species=="Osse")
                                 , x_Raca = as.numeric(data_stan_temp$species=="Raca")
                                 # , x_Rapi = as.numeric(data_stan_temp$species=="Rapi")
        )
        fit_temp <- stan(file = modelUse
                         , data=data_stan_format
                         , warmup=1500, iter=3000, chains=4
                         , control = list(adapt_delta = 0.99999, max_treedepth=20)
                         , verbose = FALSE)
        
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
  save(allModels_predictBd, file = "allModels_predictBd.RData")
  save(aggregatedData_baseline, file="aggregatedData_baseline.RData")
} else {
  load("allModels_predictBd.RData")
  load("aggregatedData_baseline.RData")
  
  # load("allModels_predictBd_PABD.RData")
  # load("aggregatedData_baseline_PABD.RData")
}


#### Aggregate #####
dir.create("UsedRanksPlots")
for (var_temp in allVarFilt) {
  # var_temp <- "inhibRich"

  agg_temp <- aggregatedData_baseline %>% filter(micro==var_temp) %>% rowwise() %>%
    mutate(indivSamp = medrank, indivSamp_lowerse = medrank-serank, indivSamp_uperse = medrank + serank) %>%
    select(species, indivID, indivSamp, indivSamp_lowerse, indivSamp_uperse) %>%
    left_join(mf_collapsed_bd) %>% 
    arrange(maxBd_outcome) %>% distinct() %>%
    mutate(maxBd_outcome_allzero = maxBd_outcome) %>%
    mutate(maxBd_outcome = ifelse(maxBd_outcome==0,NA,maxBd_outcome))

    ggused_ranks <- ggplot() +
      geom_linerange(data=agg_temp, aes(x=species, ymin=indivSamp_lowerse, ymax=indivSamp_uperse),  position=position_dodge2(width=0.6)) +
      geom_point(data=agg_temp, aes(x=species, y=indivSamp, fill=maxBd_outcome, pch=factor(PABD)), cex = 3, position=position_dodge2(width=0.6)) +
      geom_vline(data=data.frame(divide=c(1.5,2.5,3.5,4.5)), aes(xintercept = divide), lwd=3, col="white") +
      # scale_color_manual(values=c(NA, "black")) +
      scale_fill_gradient(low="grey", high="red", na.value = "white") + 
      scale_shape_manual(values=c(19,21)) +
      labs(fill="Future maximum\n Bd load (log)", pch="Did individual\nbecome infected?", col="Did individual\nbecome infected?") +
      xlab("Species") + ylab(paste0("Percentile rank of\n",allVarNames[var_temp])) +
      theme(panel.grid.major.x  = element_blank())
    # ggused_ranks
    ggsave(filename = paste0("./UsedRanksPlots/usedRanks_",var_temp,".png"), height=4, width=6
           ,ggused_ranks)
  
    ggused_corr <- ggplot() +
      geom_linerange(data=agg_temp, aes(x=maxBd_outcome_allzero, ymin=indivSamp_lowerse, ymax=indivSamp_uperse),  position=position_dodge2(width=0.6)) +
      geom_point(data=agg_temp, aes(x=maxBd_outcome_allzero, y=indivSamp, fill=species, col=factor(PABD)),pch=21, cex = 3, position=position_dodge2(width=0.6)) +
      coord_flip()+
      # geom_vline(data=data.frame(divide=c(1.5,2.5,3.5,4.5)), aes(xintercept = divide), lwd=3, col="white") +
      scale_color_manual(values=c("white", "black")) +
      # scale_fill_gradient(low="grey", high="red", na.value = "white") + 
      # scale_shape_manual(values=c(21,19)) +
      labs(fill="Species", pch="Did individual\nbecome infected?", col="Did individual\nbecome infected?") +
      xlab("Species") + ylab(paste0("Percentile rank of\n",allVarNames[var_temp])) +
      theme(panel.grid.major.x  = element_blank())
    # ggused_corr
    ggsave(filename = paste0("./UsedRanksPlots/usedRanksCorr_",var_temp,".png"), height=4, width=6
           ,ggused_corr)
  
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

## Look at species-level
colors_of_traits <- c("red","darkred","purple","pink","green","blue")
names(colors_of_traits) <- allVarNames[allVarFilt]
gg_species <- allSamps_bdresponse_adj %>%
  filter(Predictor !="micro") %>%
  ggplot() + 
  geom_linerange(aes(x=Predictor, ymin=lower95, ymax=upper95, col=MP), position = position_dodge(width=0.8))+
  geom_point(aes(x=Predictor, y=medsamp, col=MP), position = position_dodge(width=0.8)) +
  facet_grid(. ~ Response, labeller = as_labeller(c(maxBd_outcome = "Maximum future Bd load\n (log)", PABD = "Probability of infection\n(Presence/absence)"))) +
  coord_flip() + labs(col="Microbiome community trait") +
  scale_color_manual(values=c(colors_of_traits)) + xlab("Species") + ylab("Effect on infection load or probability")
ggsave(filename = "bayes_custom_species_effects.png", height=6, width=8
       , gg_species)

gg_microbiome <- allSamps_bdresponse_adj %>%
  filter(Predictor == "micro")  %>%
  mutate(Response = factor(Response, levels=c("PABD","maxBd_outcome"))) %>%
  ggplot() + 
  geom_point(aes(x=MP, y=medsamp, col=sig)) +
  geom_segment(aes(x=MP, xend = MP, y=lower95, yend=upper95, col=sig))+
  geom_hline(aes(yintercept=0), col="red",alpha=0.5)+
  facet_grid(. ~ Response, labeller = as_labeller(c(maxBd_outcome = "Maximum future Bd load\n (log)", PABD = "Probability of infection\n(Presence/absence)"))) +
  coord_flip() +labs(col="One-sided posterior\nprobability >95%") +
  scale_color_manual(values=c("grey","blue")) +
  xlab("Microbial community trait") + ylab("Effect on infection load or probability\n(Coefficient estimate)")
ggsave(filename = "bayes_custom_microbiome_effects.png", height=4, width=7
       ,gg_microbiome)

toWriteCustomStats <- allSamps_bdresponse_adj %>%
  filter(Predictor == "micro") 
write.table(toWriteCustomStats, file="microbiome_effects_stats.txt", sep="\t", row.names = FALSE, quote = FALSE)

setwd("..")
