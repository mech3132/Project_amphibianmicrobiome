# Codebits

anova_log_lmer_3way <- function(mf, dep, indep1, indep2, indep3) {
  # change names in mapping file
  temp_mf <- mf %>%
    rename(indep1=(indep1), indep2=(indep2), indep3=(indep3), dep=(dep)) %>%
    dplyr::select(dep, indep1, indep2, indep3, indivID)
  
  inter_3way <- Anova(lmer(log(dep) ~ indep1*indep2*indep3 + (1|indivID), data=temp_mf)) # First do a sequential test
  inter_3way_sig <- inter_3way$`Pr(>Chisq)`[7] < 0.05
  if ( inter_3way_sig ) { # if indep1:time:Bd_exposure IS significant, we should do an ANOVA III to test all 2-way effects
    # Test each 2-way interaction
    inter_2way <- Anova(lmer(log(dep) ~ indep1*indep2*indep3 + (1|indivID), data=temp_mf), type=3) 
    indep1_indep2_sig <- (inter_2way$`Pr(>Chisq)`[5] < 0.05)
    indep1_indep3_sig <- (inter_2way$`Pr(>Chisq)`[6] < 0.05)
    indep2_indep3_sig <- (inter_2way$`Pr(>Chisq)`[7] < 0.05)
    
    if ( indep1_indep2_sig ) {
      if ( indep1_indep3_sig ) {
        if ( indep2_indep3_sig ) {
          # if all 3 are sig
          main_1way <- Anova(lmer(log(dep) ~ indep1*indep2*indep3 + (1|indivID), data=temp_mf), type=3)
        } else {
          # if 1,2 and 1,3 are sig, NOT 2,3
          main_1way <- Anova(lmer(log(dep) ~ indep1*indep2*indep3-indep2:indep3 + (1|indivID), data=temp_mf), type=3)
        }
      } else {
        if ( indep2_indep3_sig ) {
          #if 1,2 and 2,3 are sig ( NOT 1,3)
          main_1way <- Anova(lmer(log(dep) ~ indep1*indep2*Bd_exposure-indep1:indep3 + (1|indivID), data=temp_mf), type=3)
        } else {
          # if 1,2 is sig (NOT 1,3 or 2,3)
          main_1way <- Anova(lmer(log(dep) ~ indep1*indep2*Bd_exposure-indep1:indep3-indep2:indep3 + (1|indivID), data=temp_mf), type=3)
        }
      }
    } else {
      if ( indep1_indep3_sig ) {
        if ( indep2_indep3_sig ) {
          # 1,3 and 2,3 is significant (NOT 1,2)
          main_1way <- Anova(lmer(log(dep) ~ indep1*indep2*indep3-indep1:indep2 + (1|indivID), data=temp_mf), type=3)
        } else {
          # 1,3 is sig (NOT 1,2 and 2,3)
          main_1way <- Anova(lmer(log(dep) ~ indep1*indep2*indep3-indep1:indep2-indep2:indep3 + (1|indivID), data=temp_mf), type=3)
        }
      } else {
        if ( indep2_indep3_sig ) {
          # 2,3 is sig (NOT 1,2 and 1,3)
          main_1way <- Anova(lmer(log(dep) ~ indep1*indep2*indep3-indep1:indep2-species:indep3 + (1|indivID), data=temp_mf), type=3)
        } else {
          # Nothing is significant
          main_1way <- Anova(lmer(log(dep) ~ indep1*indep2*indep3-indep1:indep2-indep1:indep3-indep2:indep3 + (1|indivID), data=temp_mf), type=3)
        }
      }
    }
    
    
  } else { # if there is NO 3-way significant interaction
    
    # Test each 2-way interaction
    inter_2way <- Anova(lmer(log(dep) ~ indep1*indep2*indep3-indep1:indep2:indep3 + (1|indivID), data=temp_mf), type=3)
    indep1_indep2_sig <- (inter_2way$`Pr(>Chisq)`[5] < 0.05)
    indep1_indep3_sig <- (inter_2way$`Pr(>Chisq)`[6] < 0.05)
    indep2_indep3_sig <- (inter_2way$`Pr(>Chisq)`[7] < 0.05)
    
    if ( indep1_indep2_sig ) {
      if ( indep1_indep3_sig ) {
        if ( indep2_indep3_sig ) {
          # if all 3 are sig
          main_1way <- Anova(lmer(log(dep) ~ indep1*indep2*indep3-indep1:indep2:indep3 +(1|indivID), data=temp_mf), type=3)
        } else {
          # if 1,2 and 1,3 are sig, NOT 2,3
          main_1way <- Anova(lmer(log(dep) ~ indep1*indep2*indep3-indep1:indep2:indep3-indep2:indep3+(1|indivID), data=temp_mf), type=3)
        }
      } else {
        if ( indep2_indep3_sig ) {
          #if 1,2 and 2,3 are sig ( NOT 1,3)
          main_1way <- Anova(lmer(log(dep) ~ indep1*indep2*indep3-indep1:indep2:indep3-indep1:indep3+(1|indivID), data=temp_mf), type=3)
        } else {
          # if 1,2 is sig (NOT 1,3 or 2,3)
          main_1way <- Anova(lmer(log(dep) ~ indep1*indep2*indep3-indep1:indep2:indep3-indep1:indep3-indep2:indep3+(1|indivID), data=temp_mf), type=3)
        }
      }
    } else {
      if ( indep1_indep3_sig ) {
        if ( indep2_indep3_sig ) {
          # 1,3 and 2,3 is significant (NOT 1,2)
          main_1way <- Anova(lmer(log(dep) ~ indep1*indep2*indep3-indep1:indep2:indep3-indep1:indep2+(1|indivID), data=temp_mf), type=3)
        } else {
          # 1,3 is sig (NOT 1,2 and 2,3)
          main_1way <- Anova(lmer(log(dep) ~ indep1*indep2*indep3-indep1:indep2:indep3-indep1:indep2-indep2:indep3+(1|indivID), data=temp_mf), type=3)
        }
      } else {
        if ( indep2_indep3_sig ) {
          # 2,3 is sig (NOT 1,2 and 1,3)
          main_1way <- Anova(lmer(log(dep) ~ indep1*indep2*indep3-indep1:indep2:indep3-indep1:indep2-indep1:indep3+(1|indivID), data=temp_mf), type=3)
        } else {
          # Nothing is significant
          main_1way <- Anova(lmer(log(dep) ~ indep1+indep2+indep3+(1|indivID), data=temp_mf), type=3)
        }
      }
    }
  }
  
  ## Now, spit out results
  print("3-way interaction")
  print(paste0(indep1,":", indep2,":", indep3, " p = ",inter_3way$`Pr(>Chisq)`[7], ", Chisq(",inter_3way$Df[7],",",inter_3way$Df[length(inter_3way$Df)],") = ",inter_3way$`Chisq`[7]))
  
  print("2-way interaction")
  print(paste0(indep1,":", indep2, " p = ",inter_2way$`Pr(>Chisq)`[5], ", Chisq(",inter_2way$Df[5],",",inter_2way$Df[length(inter_2way$Df)],") = ",inter_2way$`Chisq`[5]))
  print(paste0(indep1,":", indep3, " p = ",inter_2way$`Pr(>Chisq)`[6], ", Chisq(",inter_2way$Df[6],",",inter_2way$Df[length(inter_2way$Df)],") = ",inter_2way$`F value`[6]))
  print(paste0(indep2,":", indep3, " p = ",inter_2way$`Pr(>Chisq)`[7], ", Chisq(",inter_2way$Df[7],",",inter_2way$Df[length(inter_2way$Df)],") = ",inter_2way$`F value`[7]))
  
  print("1-way main interactions")
  print(paste0(indep1, " p = ",main_1way$`Pr(>Chisq)`[2], ", Chisq(",main_1way$Df[2],",",inter_3way$Df[length(inter_3way$Df)],") = ",main_1way$`Chisq`[2]))
  print(paste0(indep2, " p = ",main_1way$`Pr(>Chisq)`[3], ", Chisq(",main_1way$Df[3],",",inter_3way$Df[length(inter_3way$Df)],") = ",main_1way$`Chisq`[3]))
  print(paste0(indep3, " p = ",main_1way$`Pr(>Chisq)`[4], ", Chisq(",main_1way$Df[4],",",inter_3way$Df[length(inter_3way$Df)],") = ",main_1way$`Chisq`[4]))
  
  return(list(inter_3way, inter_2way, main_1way))
  
}

anova_log_lm_3way <- function(mf, dep, indep1, indep2, indep3) {
  # change names in mapping file
  temp_mf <- mf %>%
    rename(indep1=(indep1), indep2=(indep2), indep3=(indep3), dep=(dep)) %>%
    dplyr::select(dep, indep1, indep2, indep3)
  
  inter_3way <- anova(lm(log(dep) ~ indep1*indep2*indep3, data=temp_mf)) # First do a sequential test
  inter_3way_sig <- inter_3way$`Pr(>F)`[7] < 0.05
  if ( inter_3way_sig ) { # if indep1:time:Bd_exposure IS significant, we should do an ANOVA III to test all 2-way effects
    # Test each 2-way interaction
    inter_2way <- Anova(lm(log(dep) ~ indep1*indep2*indep3, data=temp_mf), type=3) 
    indep1_indep2_sig <- (inter_2way$`Pr(>F)`[5] < 0.05)
    indep1_indep3_sig <- (inter_2way$`Pr(>F)`[6] < 0.05)
    indep2_indep3_sig <- (inter_2way$`Pr(>F)`[7] < 0.05)
    
    if ( indep1_indep2_sig ) {
      if ( indep1_indep3_sig ) {
        if ( indep2_indep3_sig ) {
          # if all 3 are sig
          main_1way <- Anova(lm(log(dep) ~ indep1*indep2*indep3, data=temp_mf), type=3)
        } else {
          # if 1,2 and 1,3 are sig, NOT 2,3
          main_1way <- Anova(lm(log(dep) ~ indep1*indep2*indep3-indep2:indep3, data=temp_mf), type=3)
        }
      } else {
        if ( indep2_indep3_sig ) {
          #if 1,2 and 2,3 are sig ( NOT 1,3)
          main_1way <- Anova(lm(log(dep) ~ indep1*indep2*Bd_exposure-indep1:indep3, data=temp_mf), type=3)
        } else {
          # if 1,2 is sig (NOT 1,3 or 2,3)
          main_1way <- Anova(lm(log(dep) ~ indep1*indep2*Bd_exposure-indep1:indep3-indep2:indep3, data=temp_mf), type=3)
        }
      }
    } else {
      if ( indep1_indep3_sig ) {
        if ( indep2_indep3_sig ) {
          # 1,3 and 2,3 is significant (NOT 1,2)
          main_1way <- Anova(lm(log(dep) ~ indep1*indep2*indep3-indep1:indep2, data=temp_mf), type=3)
        } else {
          # 1,3 is sig (NOT 1,2 and 2,3)
          main_1way <- Anova(lm(log(dep) ~ indep1*indep2*indep3-indep1:indep2-indep2:indep3, data=temp_mf), type=3)
        }
      } else {
        if ( indep2_indep3_sig ) {
          # 2,3 is sig (NOT 1,2 and 1,3)
          main_1way <- Anova(lm(log(dep) ~ indep1*indep2*indep3-indep1:indep2-species:indep3, data=temp_mf), type=3)
        } else {
          # Nothing is significant
          main_1way <- Anova(lm(log(dep) ~ indep1*indep2*indep3-indep1:indep2-indep1:indep3-indep2:indep3, data=temp_mf), type=3)
        }
      }
    }
    
    
  } else { # if there is NO 3-way significant interaction
    
    # Test each 2-way interaction
    inter_2way <- Anova(lm(log(dep) ~ indep1*indep2*indep3-indep1:indep2:indep3, data=temp_mf), type=3)
    indep1_indep2_sig <- (inter_2way$`Pr(>F)`[5] < 0.05)
    indep1_indep3_sig <- (inter_2way$`Pr(>F)`[6] < 0.05)
    indep2_indep3_sig <- (inter_2way$`Pr(>F)`[7] < 0.05)
    
    if ( indep1_indep2_sig ) {
      if ( indep1_indep3_sig ) {
        if ( indep2_indep3_sig ) {
          # if all 3 are sig
          main_1way <- Anova(lm(log(dep) ~ indep1*indep2*indep3-indep1:indep2:indep3, data=temp_mf), type=3)
        } else {
          # if 1,2 and 1,3 are sig, NOT 2,3
          main_1way <- Anova(lm(log(dep) ~ indep1*indep2*indep3-indep1:indep2:indep3-indep2:indep3, data=temp_mf), type=3)
        }
      } else {
        if ( indep2_indep3_sig ) {
          #if 1,2 and 2,3 are sig ( NOT 1,3)
          main_1way <- Anova(lm(log(dep) ~ indep1*indep2*indep3-indep1:indep2:indep3-indep1:indep3, data=temp_mf), type=3)
        } else {
          # if 1,2 is sig (NOT 1,3 or 2,3)
          main_1way <- Anova(lm(log(dep) ~ indep1*indep2*indep3-indep1:indep2:indep3-indep1:indep3-indep2:indep3, data=temp_mf), type=3)
        }
      }
    } else {
      if ( indep1_indep3_sig ) {
        if ( indep2_indep3_sig ) {
          # 1,3 and 2,3 is significant (NOT 1,2)
          main_1way <- Anova(lm(log(dep) ~ indep1*indep2*indep3-indep1:indep2:indep3-indep1:indep2, data=temp_mf), type=3)
        } else {
          # 1,3 is sig (NOT 1,2 and 2,3)
          main_1way <- Anova(lm(log(dep) ~ indep1*indep2*indep3-indep1:indep2:indep3-indep1:indep2-indep2:indep3, data=temp_mf), type=3)
        }
      } else {
        if ( indep2_indep3_sig ) {
          # 2,3 is sig (NOT 1,2 and 1,3)
          main_1way <- Anova(lm(log(dep) ~ indep1*indep2*indep3-indep1:indep2:indep3-indep1:indep2-indep1:indep3, data=temp_mf), type=3)
        } else {
          # Nothing is significant
          main_1way <- Anova(lm(log(dep) ~ indep1+indep2+indep3, data=temp_mf), type=3)
        }
      }
    }
  }
  
  ## Now, spit out results
  print("3-way interaction")
  print(paste0(indep1,":", indep2,":", indep3, " p = ",inter_3way$`Pr(>F)`[7], ", F(",inter_3way$Df[7],",",inter_3way$Df[length(inter_3way$Df)],") = ",inter_3way$`F value`[7]))
  
  print("2-way interaction")
  print(paste0(indep1,":", indep2, " p = ",inter_2way$`Pr(>F)`[5], ", F(",inter_2way$Df[5],",",inter_2way$Df[length(inter_2way$Df)],") = ",inter_2way$`F value`[5]))
  print(paste0(indep1,":", indep3, " p = ",inter_2way$`Pr(>F)`[6], ", F(",inter_2way$Df[6],",",inter_2way$Df[length(inter_2way$Df)],") = ",inter_2way$`F value`[6]))
  print(paste0(indep2,":", indep3, " p = ",inter_2way$`Pr(>F)`[7], ", F(",inter_2way$Df[7],",",inter_2way$Df[length(inter_2way$Df)],") = ",inter_2way$`F value`[7]))
  
  print("1-way main interactions")
  print(paste0(indep1, " p = ",main_1way$`Pr(>F)`[2], ", F(",main_1way$Df[2],",",inter_3way$Df[length(inter_3way$Df)],") = ",main_1way$`F value`[2]))
  print(paste0(indep2, " p = ",main_1way$`Pr(>F)`[3], ", F(",main_1way$Df[3],",",inter_3way$Df[length(inter_3way$Df)],") = ",main_1way$`F value`[3]))
  print(paste0(indep3, " p = ",main_1way$`Pr(>F)`[4], ", F(",main_1way$Df[4],",",inter_3way$Df[length(inter_3way$Df)],") = ",main_1way$`F value`[4]))
  
  return(list(inter_3way, inter_2way, main_1way))
  
}


# Or do manual comparisons of control and Bd_exposure
# Make a function to calculate this so we don't have to type it out every time
anova_log_lm_2way <- function(mf, dep, indep1, indep2) {
  # change names in mapping file
  temp_mf <- mf %>%
    rename(indep1=(indep1), indep2=(indep2), dep=(dep)) %>%
    dplyr::select(dep, indep1, indep2)
  
  inter_2way <- anova(lm(log(dep) ~ indep1*indep2, data=temp_mf)) 
  inter_2way_sig <- inter_2way$`Pr(>F)`[3] < 0.05
  if ( inter_2way_sig ) {
    main_1way <- Anova(lm(log(dep) ~ indep1*indep2, data=temp_mf), type=3)
  } else {
    main_1way <- Anova(lm(log(dep) ~ indep1+indep2, data=temp_mf), type=3)
  }
  
  print("2-way interaction")
  print(paste0(indep1,":", indep2, " p = ",inter_2way$`Pr(>F)`[3], ", F(",inter_2way$Df[3],",",inter_2way$Df[length(inter_2way$Df)],") = ",inter_2way$`F value`[3]))
  
  print("1-way main interactions")
  print(paste0(indep1, " p = ",main_1way$`Pr(>F)`[2], ", F(",main_1way$Df[2],",",main_1way$Df[length(main_1way$Df)],") = ",main_1way$`F value`[2]))
  print(paste0(indep2, " p = ",main_1way$`Pr(>F)`[3], ", F(",main_1way$Df[3],",",main_1way$Df[length(main_1way$Df)],") = ",main_1way$`F value`[3]))
  
  return(list(inter_2way, main_1way))
  
}



# Set up a function to calculate ANOVAs
# Make a function to calculate this so we don't have to type it out every time
# mf=mf_alt_filt_final
# dep = paste0("dist_",b)
# indep1="species"
# indep2="time"
# indep3="Bd_exposure"
anova_betareg_3way <- function(mf, dep, indep1, indep2, indep3) {
  # change names in mapping file
  temp_mf <- mf %>%
    rename(indep1=(indep1), indep2=(indep2), indep3=(indep3), dep=(dep)) %>%
    dplyr::select(dep, indep1, indep2, indep3, indivID)

  inter_3way <- Anova(betareg((dep) ~ indep1*indep2*indep3, data=temp_mf), type=3) # First do a sequential test
  # inter_3way <- Anova(betareg((dep) ~ indep1*indep2*indep3 + (1|indivID), data=temp_mf), type=3) # First do a sequential test
  inter_3way_sig <- inter_3way$`Pr(>Chisq)`[8] < 0.05
  if ( inter_3way_sig ) { # if indep1:time:Bd_exposure IS significant, we should do an ANOVA III to test all 2-way effects
    # Test each 2-way interaction
    inter_2way <- Anova(betareg((dep) ~ indep1*indep2*indep3, data=temp_mf), type=3) 
    indep1_indep2_sig <- (inter_2way$`Pr(>Chisq)`[5] < 0.05)
    indep1_indep3_sig <- (inter_2way$`Pr(>Chisq)`[6] < 0.05)
    indep2_indep3_sig <- (inter_2way$`Pr(>Chisq)`[7] < 0.05)
    
    if ( indep1_indep2_sig ) {
      if ( indep1_indep3_sig ) {
        if ( indep2_indep3_sig ) {
          # if all 3 are sig
          main_1way <- Anova(betareg((dep) ~ indep1*indep2*indep3, data=temp_mf), type=3)
        } else {
          # if 1,2 and 1,3 are sig, NOT 2,3
          main_1way <- Anova(betareg((dep) ~ indep1*indep2*indep3-indep2:indep3, data=temp_mf), type=3)
        }
      } else {
        if ( indep2_indep3_sig ) {
          #if 1,2 and 2,3 are sig ( NOT 1,3)
          main_1way <- Anova(betareg((dep) ~ indep1*indep2*Bd_exposure-indep1:indep3, data=temp_mf), type=3)
        } else {
          # if 1,2 is sig (NOT 1,3 or 2,3)
          main_1way <- Anova(betareg((dep) ~ indep1*indep2*Bd_exposure-indep1:indep3-indep2:indep3, data=temp_mf), type=3)
        }
      }
    } else {
      if ( indep1_indep3_sig ) {
        if ( indep2_indep3_sig ) {
          # 1,3 and 2,3 is significant (NOT 1,2)
          main_1way <- Anova(betareg((dep) ~ indep1*indep2*indep3-indep1:indep2, data=temp_mf), type=3)
        } else {
          # 1,3 is sig (NOT 1,2 and 2,3)
          main_1way <- Anova(betareg((dep) ~ indep1*indep2*indep3-indep1:indep2-indep2:indep3, data=temp_mf), type=3)
        }
      } else {
        if ( indep2_indep3_sig ) {
          # 2,3 is sig (NOT 1,2 and 1,3)
          main_1way <- Anova(betareg((dep) ~ indep1*indep2*indep3-indep1:indep2-species:indep3, data=temp_mf), type=3)
        } else {
          # Nothing is significant
          main_1way <- Anova(betareg((dep) ~ indep1*indep2*indep3-indep1:indep2-indep1:indep3-indep2:indep3, data=temp_mf), type=3)
        }
      }
    }
    
    
  } else { # if there is NO 3-way significant interaction
    
    # Test each 2-way interaction
    inter_2way <- Anova(betareg((dep) ~ indep1*indep2*indep3-indep1:indep2:indep3, data=temp_mf), type=3)
    indep1_indep2_sig <- (inter_2way$`Pr(>Chisq)`[5] < 0.05)
    indep1_indep3_sig <- (inter_2way$`Pr(>Chisq)`[6] < 0.05)
    indep2_indep3_sig <- (inter_2way$`Pr(>Chisq)`[7] < 0.05)
    
    if ( indep1_indep2_sig ) {
      if ( indep1_indep3_sig ) {
        if ( indep2_indep3_sig ) {
          # if all 3 are sig
          main_1way <- Anova(betareg((dep) ~ indep1*indep2*indep3-indep1:indep2:indep3, data=temp_mf), type=3)
        } else {
          # if 1,2 and 1,3 are sig, NOT 2,3
          main_1way <- Anova(betareg((dep) ~ indep1*indep2*indep3-indep1:indep2:indep3-indep2:indep3, data=temp_mf), type=3)
        }
      } else {
        if ( indep2_indep3_sig ) {
          #if 1,2 and 2,3 are sig ( NOT 1,3)
          main_1way <- Anova(betareg((dep) ~ indep1*indep2*indep3-indep1:indep2:indep3-indep1:indep3, data=temp_mf), type=3)
        } else {
          # if 1,2 is sig (NOT 1,3 or 2,3)
          main_1way <- Anova(betareg((dep) ~ indep1*indep2*indep3-indep1:indep2:indep3-indep1:indep3-indep2:indep3, data=temp_mf), type=3)
        }
      }
    } else {
      if ( indep1_indep3_sig ) {
        if ( indep2_indep3_sig ) {
          # 1,3 and 2,3 is significant (NOT 1,2)
          main_1way <- Anova(betareg((dep) ~ indep1*indep2*indep3-indep1:indep2:indep3-indep1:indep2, data=temp_mf), type=3)
        } else {
          # 1,3 is sig (NOT 1,2 and 2,3)
          main_1way <- Anova(betareg((dep) ~ indep1*indep2*indep3-indep1:indep2:indep3-indep1:indep2-indep2:indep3, data=temp_mf), type=3)
        }
      } else {
        if ( indep2_indep3_sig ) {
          # 2,3 is sig (NOT 1,2 and 1,3)
          main_1way <- Anova(betareg((dep) ~ indep1*indep2*indep3-indep1:indep2:indep3-indep1:indep2-indep1:indep3, data=temp_mf), type=3)
        } else {
          # Nothing is significant
          main_1way <- Anova(betareg((dep) ~ indep1+indep2+indep3, data=temp_mf), type=3)
        }
      }
    }
  }
  ## Now, spit out results
  print("3-way interaction")
  print(paste0(indep1,":", indep2,":", indep3, " p = ",inter_3way$`Pr(>Chisq)`[8], ", Chisq(",inter_3way$Df[8],") = ",inter_3way$`Chisq`[8]))
  
  print("2-way interaction")
  print(paste0(indep1,":", indep2, " p = ",inter_2way$`Pr(>Chisq)`[5], ", Chisq(",inter_2way$Df[5],") = ",inter_2way$`Chisq`[5]))
  print(paste0(indep1,":", indep3, " p = ",inter_2way$`Pr(>Chisq)`[6], ", Chisq(",inter_2way$Df[6],") = ",inter_2way$`Chisq`[6]))
  print(paste0(indep2,":", indep3, " p = ",inter_2way$`Pr(>Chisq)`[7], ", Chisq(",inter_2way$Df[7],") = ",inter_2way$`Chisq`[7]))
  
  print("1-way main interactions")
  print(paste0(indep1, " p = ",main_1way$`Pr(>Chisq)`[2], ", Chisq(",main_1way$Df[2],") = ",main_1way$`Chisq`[2]))
  print(paste0(indep2, " p = ",main_1way$`Pr(>Chisq)`[3], ", Chisq(",main_1way$Df[3],") = ",main_1way$`Chisq`[3]))
  print(paste0(indep3, " p = ",main_1way$`Pr(>Chisq)`[4], ", Chisq(",main_1way$Df[4],") = ",main_1way$`Chisq`[4]))
  
  return(list(inter_3way, inter_2way, main_1way))
  
}



anova_log_lmer_2way <- function(mf, dep, indep1, indep2) {
  # change names in mapping file
  temp_mf <- mf %>%
    rename(indep1=(indep1), indep2=(indep2), dep=(dep)) %>%
    dplyr::select(dep, indep1, indep2, indivID)
  
  inter_2way <- Anova(lmer(log(dep) ~ indep1*indep2 + (1|indivID), data=temp_mf)) 
  inter_2way_sig <- inter_2way$`Pr(>Chisq)`[3] < 0.05
  if ( inter_2way_sig ) {
    main_1way <- Anova(lmer(log(dep) ~ indep1*indep2+ (1|indivID), data=temp_mf), type=3)
  } else {
    main_1way <- Anova(lmer(log(dep) ~ indep1+indep2+ (1|indivID), data=temp_mf), type=3)
  }
  
  print("2-way interaction")
  print(paste0(indep1,":", indep2, " p = ",inter_2way$`Pr(>Chisq)`[3], ", Chisq(",inter_2way$Df[3],",",inter_2way$Df[length(inter_2way$Df)],") = ",inter_2way$Chisq[3]))
  
  print("1-way main interactions")
  print(paste0(indep1, " p = ",main_1way$`Pr(>Chisq)`[2], ", Chisq(",main_1way$Df[2],",",main_1way$Df[length(main_1way$Df)],") = ",main_1way$Chisq[2]))
  print(paste0(indep2, " p = ",main_1way$`Pr(>Chisq)`[3], ", Chisq(",main_1way$Df[3],",",main_1way$Df[length(main_1way$Df)],") = ",main_1way$Chisq[3]))
  
  return(list(inter_2way, main_1way))
  
}

# Or do manual comparisons of control and Bd_exposure
# Make a function to calculate this so we don't have to type it out every time
anova_betareg_2way <- function(mf, dep, indep1, indep2) {
  # change names in mapping file
  temp_mf <- mf %>%
    rename(indep1=(indep1), indep2=(indep2), dep=(dep)) %>%
    dplyr::select(dep, indep1, indep2)
  
  inter_2way <- Anova(betareg((dep) ~ indep1*indep2, data=temp_mf), type=3) 
  inter_2way_sig <- inter_2way$`Pr(>Chisq)`[4] < 0.05
  if ( inter_2way_sig ) {
    main_1way <- Anova(betareg((dep) ~ indep1*indep2, data=temp_mf), type=3)
  } else {
    main_1way <- Anova(betareg((dep) ~ indep1*indep2, data=temp_mf), type=2)
  }
  
  print("2-way interaction")
  print(paste0(indep1,":", indep2, " p = ",inter_2way$`Pr(>Chisq)`[4], ", Chisq(",inter_2way$Df[4],") = ",inter_2way$`Chisq`[4]))
  
  print("1-way main interactions")
  print(paste0(indep1, " p = ",main_1way$`Pr(>Chisq)`[1], ", Chisq(",main_1way$Df[1],") = ",main_1way$`Chisq`[1]))
  print(paste0(indep2, " p = ",main_1way$`Pr(>Chisq)`[2], ", Chisq(",main_1way$Df[2],") = ",main_1way$`Chisq`[2]))
  
  return(list(inter_2way, main_1way))
  
}

# Or do manual comparisons of control and Bd_exposure
# Make a function to calculate this so we don't have to type it out every time
anova_betareg_2way_percInhib <- function(mf, dep, indep1, indep2) {
  # change names in mapping file
  temp_mf <- mf %>%
    rename(indep1=(indep1), indep2=(indep2), dep=(dep)) %>%
    dplyr::select(dep, indep1, indep2)
  
  inter_2way <- Anova(betareg((dep) ~ indep1*indep2, data=temp_mf), type=3) 
  inter_2way_sig <- inter_2way$`Pr(>Chisq)`[4] < 0.05
  if ( inter_2way_sig ) {
    main_1way <- Anova(betareg((dep) ~ indep1*indep2, data=temp_mf), type=3)
  } else {
    main_1way <- Anova(betareg((dep) ~ indep1*indep2, data=temp_mf), type=2)
  }
  
  print("2-way interaction")
  print(paste0(indep1,":", indep2, " p = ",inter_2way$`Pr(>Chisq)`[4], ", Chisq(",inter_2way$Df[4],") = ",inter_2way$`Chisq`[4]))
  
  print("1-way main interactions")
  print(paste0(indep1, " p = ",main_1way$`Pr(>Chisq)`[2], ", Chisq(",main_1way$Df[2],") = ",main_1way$`Chisq`[2]))
  print(paste0(indep2, " p = ",main_1way$`Pr(>Chisq)`[3], ", Chisq(",main_1way$Df[3],") = ",main_1way$`Chisq`[3]))
  
  return(list(inter_2way, main_1way))
  
}



# Set up a function to calculate ANOVAs
# Make a function to calculate this so we don't have to type it out every time
anova_poisson_3way <- function(mf, dep, indep1, indep2, indep3) {
  # change names in mapping file
  temp_mf <- mf %>%
    rename(indep1=(indep1), indep2=(indep2), indep3=(indep3), dep=(dep)) %>%
    dplyr::select(dep, indep1, indep2, indep3)
  
  inter_3way <- Anova(glm((dep) ~ indep1*indep2*indep3, data=temp_mf, family=poisson(link="identity")), type=3) # First do a sequential test
  print(inter_3way)
  inter_3way_sig <- inter_3way$`Pr(>Chisq)`[8] < 0.05
  if ( inter_3way_sig ) { # if indep1:time:Bd_exposure IS significant, we should do an ANOVA III to test all 2-way effects
    # Test each 2-way interaction
    inter_2way <- Anova(glm((dep) ~ indep1*indep2*indep3, data=temp_mf, family=poisson(link="identity")), type=3) 
    indep1_indep2_sig <- (inter_2way$`Pr(>Chisq)`[5] < 0.05)
    indep1_indep3_sig <- (inter_2way$`Pr(>Chisq)`[6] < 0.05)
    indep2_indep3_sig <- (inter_2way$`Pr(>Chisq)`[7] < 0.05)
    
    if ( indep1_indep2_sig ) {
      if ( indep1_indep3_sig ) {
        if ( indep2_indep3_sig ) {
          # if all 3 are sig
          main_1way <- Anova(glm((dep) ~ indep1*indep2*indep3, data=temp_mf, family=poisson(link="identity")), type=3)
        } else {
          # if 1,2 and 1,3 are sig, NOT 2,3
          main_1way <- Anova(glm((dep) ~ indep1*indep2*indep3-indep2:indep3, data=temp_mf, family=poisson(link="identity")), type=3)
        }
      } else {
        if ( indep2_indep3_sig ) {
          #if 1,2 and 2,3 are sig ( NOT 1,3)
          main_1way <- Anova(glm((dep) ~ indep1*indep2*Bd_exposure-indep1:indep3, data=temp_mf, family=poisson(link="identity")), type=3)
        } else {
          # if 1,2 is sig (NOT 1,3 or 2,3)
          main_1way <- Anova(glm((dep) ~ indep1*indep2*Bd_exposure-indep1:indep3-indep2:indep3, data=temp_mf, family=poisson(link="identity")), type=3)
        }
      }
    } else {
      if ( indep1_indep3_sig ) {
        if ( indep2_indep3_sig ) {
          # 1,3 and 2,3 is significant (NOT 1,2)
          main_1way <- Anova(glm((dep) ~ indep1*indep2*indep3-indep1:indep2, data=temp_mf, family=poisson(link="identity")), type=3)
        } else {
          # 1,3 is sig (NOT 1,2 and 2,3)
          main_1way <- Anova(glm((dep) ~ indep1*indep2*indep3-indep1:indep2-indep2:indep3, data=temp_mf, family=poisson(link="identity")), type=3)
        }
      } else {
        if ( indep2_indep3_sig ) {
          # 2,3 is sig (NOT 1,2 and 1,3)
          main_1way <- Anova(glm((dep) ~ indep1*indep2*indep3-indep1:indep2-species:indep3, data=temp_mf, family=poisson(link="identity")), type=3)
        } else {
          # Nothing is significant
          main_1way <- Anova(glm((dep) ~ indep1*indep2*indep3-indep1:indep2-indep1:indep3-indep2:indep3, data=temp_mf, family=poisson(link="identity")), type=3)
        }
      }
    }
    
    
  } else { # if there is NO 3-way significant interaction
    
    # Test each 2-way interaction
    inter_2way <- Anova(glm((dep) ~ indep1*indep2*indep3-indep1:indep2:indep3, data=temp_mf, family=poisson(link="identity")), type=3)
    indep1_indep2_sig <- (inter_2way$`Pr(>Chisq)`[5] < 0.05)
    indep1_indep3_sig <- (inter_2way$`Pr(>Chisq)`[6] < 0.05)
    indep2_indep3_sig <- (inter_2way$`Pr(>Chisq)`[7] < 0.05)
    
    if ( indep1_indep2_sig ) {
      if ( indep1_indep3_sig ) {
        if ( indep2_indep3_sig ) {
          # if all 3 are sig
          main_1way <- Anova(glm((dep) ~ indep1*indep2*indep3-indep1:indep2:indep3, data=temp_mf, family=poisson(link="identity")), type=3)
        } else {
          # if 1,2 and 1,3 are sig, NOT 2,3
          main_1way <- Anova(glm((dep) ~ indep1*indep2*indep3-indep1:indep2:indep3-indep2:indep3, data=temp_mf, family=poisson(link="identity")), type=3)
        }
      } else {
        if ( indep2_indep3_sig ) {
          #if 1,2 and 2,3 are sig ( NOT 1,3)
          main_1way <- Anova(glm((dep) ~ indep1*indep2*indep3-indep1:indep2:indep3-indep1:indep3, data=temp_mf, family=poisson(link="identity")), type=3)
        } else {
          # if 1,2 is sig (NOT 1,3 or 2,3)
          main_1way <- Anova(glm((dep) ~ indep1*indep2*indep3-indep1:indep2:indep3-indep1:indep3-indep2:indep3, data=temp_mf, family=poisson(link="identity")), type=3)
        }
      }
    } else {
      if ( indep1_indep3_sig ) {
        if ( indep2_indep3_sig ) {
          # 1,3 and 2,3 is significant (NOT 1,2)
          main_1way <- Anova(glm((dep) ~ indep1*indep2*indep3-indep1:indep2:indep3-indep1:indep2, data=temp_mf, family=poisson(link="identity")), type=3)
        } else {
          # 1,3 is sig (NOT 1,2 and 2,3)
          main_1way <- Anova(glm((dep) ~ indep1*indep2*indep3-indep1:indep2:indep3-indep1:indep2-indep2:indep3, data=temp_mf, family=poisson(link="identity")), type=3)
        }
      } else {
        if ( indep2_indep3_sig ) {
          # 2,3 is sig (NOT 1,2 and 1,3)
          main_1way <- Anova(glm((dep) ~ indep1*indep2*indep3-indep1:indep2:indep3-indep1:indep2-indep1:indep3, data=temp_mf, family=poisson(link="identity")), type=3)
        } else {
          # Nothing is significant
          main_1way <- Anova(glm((dep) ~ indep1+indep2+indep3, data=temp_mf, family=poisson(link="identity")), type=3)
        }
      }
    }
  }
  ## Now, spit out results
  print("3-way interaction")
  print(paste0(indep1,":", indep2,":", indep3, " p = ",inter_3way$`Pr(>Chisq)`[8], ", Chisq(",inter_3way$Df[8],") = ",inter_3way$`Chisq`[8]))
  
  print("2-way interaction")
  print(paste0(indep1,":", indep2, " p = ",inter_2way$`Pr(>Chisq)`[5], ", Chisq(",inter_2way$Df[5],") = ",inter_2way$`Chisq`[5]))
  print(paste0(indep1,":", indep3, " p = ",inter_2way$`Pr(>Chisq)`[6], ", Chisq(",inter_2way$Df[6],") = ",inter_2way$`Chisq`[6]))
  print(paste0(indep2,":", indep3, " p = ",inter_2way$`Pr(>Chisq)`[7], ", Chisq(",inter_2way$Df[7],") = ",inter_2way$`Chisq`[7]))
  
  print("1-way main interactions")
  print(paste0(indep1, " p = ",main_1way$`Pr(>Chisq)`[2], ", Chisq(",main_1way$Df[2],") = ",main_1way$`Chisq`[2]))
  print(paste0(indep2, " p = ",main_1way$`Pr(>Chisq)`[3], ", Chisq(",main_1way$Df[3],") = ",main_1way$`Chisq`[3]))
  print(paste0(indep3, " p = ",main_1way$`Pr(>Chisq)`[4], ", Chisq(",main_1way$Df[4],") = ",main_1way$`Chisq`[4]))
  
  return(list(inter_3way, inter_2way, main_1way))
  
}

# Or do manual comparisons of control and Bd_exposure
# Make a function to calculate this so we don't have to type it out every time
anova_pois_2way <- function(mf, dep, indep1, indep2) {
  # change names in mapping file
  temp_mf <- mf %>%
    rename(indep1=(indep1), indep2=(indep2), dep=(dep)) %>%
    dplyr::select(dep, indep1, indep2)
  
  inter_2way <- Anova(glm((dep) ~ indep1*indep2, data=temp_mf, family=poisson(link="identity")), type=3) 
  inter_2way_sig <- inter_2way$`Pr(>Chisq)`[3] < 0.05
  if ( inter_2way_sig ) {
    main_1way <- Anova(glm((dep) ~ indep1*indep2, data=temp_mf, family=poisson(link="identity")), type=3)
  } else {
    main_1way <- Anova(glm((dep) ~ indep1*indep2, data=temp_mf, family=poisson(link="identity")), type=2)
  }
  print("2-way interaction")
  print(paste0(indep1,":", indep2, " p = ",inter_2way$`Pr(>Chisq)`[3], ", LR Chisq(",inter_2way$Df[3],") = ",inter_2way$`LR Chisq`[3]))
  
  print("1-way main interactions")
  print(paste0(indep1, " p = ",main_1way$`Pr(>Chisq)`[1], ", LR Chisq(",main_1way$Df[1],") = ",main_1way$`LR Chisq`[1]))
  print(paste0(indep2, " p = ",main_1way$`Pr(>Chisq)`[2], ", LR Chisq(",main_1way$Df[2],") = ",main_1way$`LR Chisq`[2]))
  
  return(list(inter_2way, main_1way))
  
}





process_glmer <- function(g_lmer, dep, name_dep, transform_raw_func, intercept_present, fit_distr, time_factor=FALSE, overdispersion=FALSE, mf_con, mf_treat) {
  
  # ### FOR TESTING
  # g_lmer = glmer_inhibRich
  # dep= "inhibRich"
  # name_dep="inhibRich"
  # transform_raw_func="None"
  # intercept_present=T
  # fit_distr="Poisson"
  # overdispersion=F
  # time_factor=T
  
  
  # Need overdispersion term for poisson
  
  if (transform_raw_func == "log") {
    trans_func <- log
    alt_samps <- function(x) {return(x)} 
  } else if (transform_raw_func == "None") {
    trans_func <- function(x) {return(x)} 
    alt_samps <- function(x) {return(x)} 
  } else if (transform_raw_func == "beta") {
    trans_func <- function(x) {return(x)}
    a <- function(mu,phi){
      mu*phi
    }
    b <- function(mu,phi) {
      phi-mu*phi
    }
    mu <- function(a,phi) {
      a/phi
    }
    alt_samps <- function(x) {
      exp(x)/(exp(x)+1)
    }
    logit <- function(p) {
      log(p/(1-p))
    }
    mu_ab <- function(a,b) {
      a/(a+b)
    }
  } 
  
  
  # Alter mapping files
  mf_con_temp <- mf_con %>%
    rename(dep=paste0(dep)) %>%
    mutate(dep=trans_func(dep)) 
  mf_treat_temp <- mf_treat %>%
    rename(dep=paste0(dep)) %>%
    mutate(dep=trans_func(dep)) 
  # Look at distributions according to models
  samps <- rstan::extract(g_lmer$stanfit)
  pre_test_set <- mf_treat_temp %>%
    filter(time<6)
  # If time factor present
  if ( time_factor & intercept_present ) {
    t_var <- alt_samps(mean(samps$beta[,5]))
  } else if (time_factor) {
    t_var <- alt_samps(mean(samps$beta[,6]))
  } else {
    t_var <- 0
  }
  
  ## PLOT OF EXPECTED FOR EACH SPECIES
  if ( intercept_present ) {
    samps$beta <- samps$beta %>%
      as.data.frame() %>%
      mutate(Anbo=samps$alpha) %>%
      mutate(Rhma=alt_samps(Anbo+V1), Osse=alt_samps(Anbo+V2), Raca=alt_samps(Anbo+V3), Rapi=alt_samps(Anbo+V4)) %>%
      dplyr::select(Anbo, Rhma, Osse, Raca, Rapi) %>%
      rename(V1=Anbo, V2=Rhma, V3=Osse, V4=Raca, V5=Rapi) %>%
      as.matrix()
    gg_distr <- samps$beta %>%
      as.data.frame() %>%
      rename(Anbo=V1, Rhma=V2, Osse=V3, Raca=V4, Rapi=V5) %>%
      dplyr::select(Anbo,Rhma,Osse,Raca,Rapi) %>%
      gather(key=species, value=dependent_var) %>%
      mutate(species=factor(species, levels=c("Anbo","Rhma","Osse","Raca","Rapi"))) %>%
      ggplot(mapping=aes(x=species, y=dependent_var))+
      geom_violin() +
      geom_point(data=mf_con_temp, aes(y=dep-time*(t_var), x=species), position = position_jitter(width = 0.1, height=0), col="blue") +
      geom_point(data=pre_test_set, aes(y=dep-time*(t_var), x=species), position=position_jitter(width = 0.1, height=0), col="red") +
      ylab(label=paste0(name_dep))
    
  } else {
    gg_distr <- samps$beta %>%
      as.data.frame() %>%
      rename(Anbo=V1, Rhma=V2, Osse=V3, Raca=V4, Rapi=V5) %>%
      mutate(Anbo=alt_samps(Anbo), Rhma=alt_samps(Rhma), Osse=alt_samps(Osse), Raca=alt_samps(Raca), Rapi=alt_samps(Rapi)) %>%
      dplyr::select(Anbo,Rhma,Osse,Raca,Rapi) %>%
      gather(key=species, value=depenent_var) %>%
      mutate(species=factor(species, levels=c("Anbo","Rhma","Osse","Raca","Rapi"))) %>%
      ggplot(mapping=aes(x=species, y=depenent_var))+
      geom_violin() +
      geom_point(data=mf_con_temp, aes(y=dep-time*(t_var), x=species), position = position_jitter(width = 0.1, height=0), col="blue") +
      geom_point(data=pre_test_set, aes(y=dep-time*(t_var), x=species), position=position_jitter(width = 0.1, height=0), col="red") +
      ylab(label=paste0(name_dep))
    # Blue is controls, which the model is base don
    # Red is treatment individuals, which we are testing
  }
  
  ## Get standard deviation between toad individuals and samples
  ## For poisson, need to account for overdispersion
  # sd (variation) between samples
  samp_sigma <- samps$aux
  if ( overdispersion ) { 
    # sd (variation) between individuals of the same species
    indiv_intercept <- ranef(g_lmer)$indivID
    indivID_sigma <- sd(samps$b[,ncol(samps$b)])
    # sd (variation) between samples taken (overdispersion)
    samp_intercept <- ranef(g_lmer)$SampleID
    sampleID_sigma <- sd(samps$b[,nrow(samp_intercept)+1])
    # get toad_samps, which is change per toad
    # Get standard deviation between toad individuals and samples
    samp_indiv <- samps$b[,(nrow(samp_intercept)+2):(nrow(samp_intercept)+1+nrow(indiv_intercept))]
    colnames(samp_indiv) <- rownames(indiv_intercept)
    
  } else {
    # sd (variation) between individuals of the same species
    indivID_sigma <- sd(samps$b[,ncol(samps$b)])
    sampleID_sigma <- 0
  }
  
  
  # Now, we can calculate the probability that the "test" dataset values come from this distribution
  # List of individuals
  treat_indiv <- unique(mf_treat_temp$indivID)
  # List of each species
  species_list <- levels(factor(mf_con_temp$species))
  
  
  species_key <- treat_indiv %>%
    as_tibble() %>%
    rename(indivID=value) %>%
    separate(indivID,into=c("species","indiv"), remove=FALSE)
  
  exp_distr <- as.data.frame(matrix(ncol=length(species_list), nrow=4000, dimnames = list(1:4000, species_list)))
  for ( num_sp in 1:length(species_list)) {
    if ( fit_distr == "normal" ) {
      # Expected value at each species level
      exp_sp <- rnorm(length(samps$beta[,num_sp]), mean=samps$beta[,num_sp], sd=indivID_sigma)
      # Expected value after "sampling"
      exp_distr[,num_sp] <- rnorm(length(samps$beta[,num_sp]), mean=exp_sp, sd=samp_sigma)
    } else if ( fit_distr == "Gamma" ) {
      # Expected value at each species level
      exp_sp <- rnorm(length(samps$beta[,num_sp]), mean=samps$beta[,num_sp], sd=indivID_sigma)
      # Expected value after "sampling"
      exp_distr[,num_sp] <- rgamma(n=length(samps$beta[,num_sp]), shape=exp_sp, scale=1) # pretty sure default is 1
    } else if  ( fit_distr == "Beta" ) {
      
      ## TEST
      mu <- alt_samps(rnorm(4000, mean=samps$beta[,num_sp], sd=indivID_sigma))
      
      exp_distr[,num_sp] <-  rbeta(length(samps$beta[,num_sp])
                                   ,shape1=a(mu,(samps$aux))
                                   ,shape2=b(mu,(samps$aux)))
      
      # # ORIGINAL
      # mu <- alt_samps(rnorm(4000, mean=samps$beta[,num_sp], sd=indivID_sigma))
      # # exp_distr[,nump_sp] <- mu
      # exp_distr[,num_sp] <- rbeta(length(samps$beta[,num_sp])
      #                             ,shape1=a(mu,samps$aux)
      #                             ,shape2=b(mu,samps$aux))
      # 
      
    } else if ( fit_distr == "Poisson" ) {
      mu <- rnorm(length(samps$beta[,num_sp]), mean=samps$beta[,num_sp], sd=sampleID_sigma)
      rpois(length(samps$beta[,num_sp]), lambda=mu)
      exp_distr[,num_sp] <- rpois(length(samps$beta[,num_sp]), lambda=mu)
      
    }
    
  }
  
  # Loop through and calculate probability of having diversity at that level
  pre_exp_indiv <- data.frame(indivID=treat_indiv, exp=rep(NA, length(treat_indiv)), p=rep(NA, length(treat_indiv)), infect=rep(NA, length(treat_indiv)))
  for ( i in treat_indiv ) {
    n_row <- match(i, treat_indiv)
    sp <- unlist(strsplit(i,"_"))
    num_sp <- match(sp[1], levels(factor(mf_con_temp$species)))
    temp_dep <- mf_treat_temp %>%
      filter(indivID==i, time <=5 ) %>%
      dplyr::select(dep) %>%
      pull()
    temp_time <- mf_treat_temp %>%
      filter(indivID==i, time <=5 ) %>%
      dplyr::select(time) %>%
      pull()
    
    x.fit <- temp_dep-temp_time*mean(t_var)
    x.fit <- x.fit[!is.na(x.fit)]
    # I don't use a beta distribution here because often times optimization fails.
    # I think the sample size is too small to adequeately estimate shape1 and shape2?
    if ( fit_distr != "Beta" ) {
      if ( length(x.fit)>1 ) {
        fitted <- fitdistr(x.fit, densfun = paste0(fit_distr))$estimate
        if ( fit_distr!="Gamma") {
          exp <- fitted[1]
        } else {
          exp <-  fitted[1]/fitted[2]
          
        }
      } else if (length(x.fit) == 1 ) {
        exp <- x.fit
      } else {
        exp <- NA
      }
    } else {
      if ( length(x.fit)>1 ) {
        exp <-  fitdistr(x.fit, densfun = "normal")$estimate[1]
      } else if (length(x.fit) == 1) {
        exp <- (x.fit)
      } else {
        exp <- NA
      }
    }
    
    
    # pred_distr <- rnorm(length(samps_lmer_shannon$beta[,num_sp]), mean=rnorm(length(samps_lmer_shannon$beta[,num_sp]), mean=samps_lmer_shannon$beta[,num_sp], sd=toadID_sigma), sd=samp_sigma)
    p <- sum(exp_distr[,sp[1]]<exp)/length(exp_distr[,sp[1]])
    
    ### Did they get infected?
    infect <- max(mf_treat %>%
                    filter(indivID==i) %>%
                    dplyr::select(Bd_load) %>%
                    pull()
    )
    
    pre_exp_indiv[n_row,c("exp","p","infect")] <- c(exp, p, infect)
    
  }
  
  if ( intercept_present ) {
    # Get estimates for control toads
    species_exp <- fixef(g_lmer)[1:5]  %>%
      tibble::enframe() %>%
      dplyr::select(value) %>%
      t() %>%
      as.data.frame() %>%
      mutate("Anbo"=V1, "Rhma"=V1+V2, "Osse"=V1+V3, "Raca"=V1+V4, "Rapi"=V1+V5) %>%
      dplyr::select(Anbo, Rhma, Osse, Raca, Rapi) %>%
      t() %>%
      as.data.frame()  %>%
      mutate(species=c("Anbo","Rhma","Osse","Raca","Rapi")) %>%
      rename(value=V1)
  } else {
    # Get estimates for control toads
    species_exp <- fixef(g_lmer)[1:5]  %>%
      tibble::enframe() %>%
      mutate(species=c("Anbo", "Rhma","Osse","Raca","Rapi"))
  }
  
  # Need to extract ranef manually for beta; not sure why fixef works and ranef doesn't
  if ( fit_distr == "Beta" ) {
    randeffects_dist <- g_lmer$coefficients[grep("Intercept", names(g_lmer$coefficients))]
    names_randef <- names(randeffects_dist)
    names(randeffects_dist) <- gsub("]","",gsub("^.* indivID:", "", x=names_randef), fixed=TRUE)
    temp_randeffects <- data.frame(X=as.numeric(randeffects_dist), row.names = names(randeffects_dist))
    temp_randeffects <- temp_randeffects %>%
      rename(sp_est="X")
  } else {
    temp_randeffects <- ranef(g_lmer)$indivID %>%
      rename(sp_est="(Intercept)")
  }
  
  con_toad_est <- temp_randeffects %>%
    mutate(indivID=rownames(temp_randeffects)) %>%
    separate(indivID, into=c("species","num"), sep="_",remove=FALSE) %>%
    dplyr::select(-num) %>%
    left_join(species_exp, by = "species") %>%
    mutate(est=sp_est+value)
  
  con_exp_indiv <- data.frame(indivID=con_toad_est$indivID, exp=rep(NA, length(con_toad_est$indivID)), p=rep(NA, length(con_toad_est$indivID)))
  for ( i in 1:nrow(con_toad_est) ) {
    s <- con_toad_est[i,"species"]
    
    exp <- alt_samps(con_toad_est[i,"est"])
    p <- sum(exp_distr[,s]<exp)/length(exp_distr[,s])
    
    con_exp_indiv[i,c("exp","p")] <- c(exp, p)
  }
  con_exp_indiv <- con_exp_indiv %>%
    separate(indivID, into=c("species","indiv"), remove=FALSE)
  
  
  # create species column
  pre_exp_indiv <- pre_exp_indiv %>%
    separate(indivID, into=c("species","indiv"), remove = FALSE) %>%
    mutate(species=factor(species, levels=c("Anbo","Rhma","Osse","Raca","Rapi")))
  # Plot results 
  gg_p <- ggplot(pre_exp_indiv, aes(x=p, y=log(infect+1))) +
    geom_point(aes(color=species), cex=4) +
    geom_smooth(aes(color=species),method=lm, se = FALSE) +
    geom_smooth(method=lm, se=FALSE, col="black") +
    xlab(label=paste0("p_",name_dep))
  # if we'd JUST plotted raw values
  gg_raw <- ggplot(pre_exp_indiv, aes(x=exp, y=log(infect+1)))+
    geom_point(aes(color=species), cex=4) +
    geom_smooth(method=lm, se=FALSE, col="black") +
    geom_smooth(aes(color=species),method=lm, se = FALSE) +
    xlab(label=paste0("exp_",name_dep))
  #grid.arrange(gg_p, gg_raw, nrow=1)
  
  gg_exp_distr <- exp_distr %>%
    # add time adjustment to raw data values
    gather(key=species, value=dep_var) %>%
    mutate(species=factor(species,levels=c("Anbo","Rhma","Osse","Raca","Rapi")))%>%
    ggplot(aes(x=species, y=dep_var)) +
    geom_violin() +
    geom_point(data=pre_exp_indiv, aes(x=species, y=exp, col=log(infect+1)), cex=4, position=position_jitter(height=0, width=0.1)) +
    ylab(label=paste0(name_dep))
  
  gg_exp_distr_controls <- exp_distr %>%
    # add time adjustment to raw data values
    gather(key=species, value=dep_var) %>%
    mutate(species=factor(species,levels=c("Anbo","Rhma","Osse","Raca","Rapi")))%>%
    ggplot(aes(x=species, y=dep_var)) +
    geom_violin() +
    geom_point(data=con_exp_indiv, aes(x=species, y=exp), cex=3, position=position_jitter(height=0, width=0.1)) +
    ylab(label=paste0(name_dep))
  
  all_p <- pre_exp_indiv %>%
    dplyr::select(indivID, exp, p,infect) 
  
  # Make list of items
  output <- list()
  output[["Var"]] <- data.frame( dep=dep, name_dep=name_dep
                                 , transform_raw_func=transform_raw_func
                                 , intercept_present=intercept_present
                                 ,fit_distr=fit_distr)
  output[["gg_model_Distribution_of_all_values"]] <- gg_distr
  output[["Control_individuals"]] <- con_exp_indiv
  output[["gg_p"]] <- gg_p
  output[["gg_raw"]] <- gg_raw
  output[["gg_ExpectedDistribution_and_Bd_exposed"]] <- gg_exp_distr
  output[["gg_ExpectedDistribution_controls"]] <- gg_exp_distr_controls
  output[["all_p"]] <- all_p
  
  return(output)
}







process_glmer_withtime <- function(g_lmer, dep, name_dep, transform_raw_func, intercept_present, fit_distr, time_factor=FALSE, time_inter = FALSE, overdispersion=FALSE, mf_con, mf_treat) {
  # ### FOR TESTING
  # g_lmer = glmer_inhibRich
  # dep= "inhibRich"
  # name_dep="inhibRich"
  # transform_raw_func="None"
  # intercept_present=T
  # fit_distr="Poisson"
  # overdispersion=F
  # time_factor=T
  # g_lmer = lmer_log_observed_otus_wtime
  # dep= "observed_otus"
  # name_dep="log_observed_otus"
  # transform_raw_func="log"
  # intercept_present=F
  # fit_distr="normal"
  # mf_con=mf_con
  # mf_treat=mf_treat
  # overdispersion=F
  # time_factor=F
  # time_inter=T
  # 
  # g_lmer = lmer_shannon_wtime
  # dep= "shannon"
  # name_dep="shannon"
  # transform_raw_func="None"
  # intercept_present=F
  # fit_distr="normal"
  # mf_con=mf_con
  # mf_treat=mf_treat
  # time_inter=T
  # time_factor = T
  # g_lmer = glmer_dist_unweighted_unifrac_wtime
  # dep= "dist_unweighted_unifrac"
  # name_dep="dist_unweighted_unifrac"
  # transform_raw_func="beta"
  # intercept_present=F
  # fit_distr="Beta"
  # mf_con=mf_con
  # mf_treat=mf_treat
  # time_factor = T
  # time_inter = F
  # g_lmer =glmer_inhibRich_wtime
  # dep= "inhibRich"
  # name_dep="inhibRich"
  # transform_raw_func="None"
  # intercept_present=T
  # fit_distr="Poisson"
  # mf_con=mf_con
  # mf_treat=mf_treat
  # time_factor=T
  # time_inter=T
  # # Need overdispersion term for poisson
  
  if (transform_raw_func == "log") {
    trans_func <- log
    alt_samps <- function(x) {return(x)} 
    inv_alt_samps <- function(x) {return(x)}
  } else if (transform_raw_func == "None") {
    trans_func <- function(x) {return(x)} 
    alt_samps <- function(x) {return(x)} 
    inv_alt_samps <- function(x) {return(x)}
  } else if (transform_raw_func == "beta") {
    trans_func <- function(x) {return(x)}
    a <- function(mu,phi){
      mu*phi
    }
    b <- function(mu,phi) {
      phi-mu*phi
    }
    mu <- function(a,phi) {
      a/phi
    }
    alt_samps <- function(x) {
      exp(x)/(exp(x)+1)
    }
    logit <- function(p) {
      log(p/(1-p))
    }
    inv_alt_samps <- logit
    mu_ab <- function(a,b) {
      a/(a+b)
    }
  } 
  
  if (fit_distr == "Poisson" ) {
    alt_samps <- function(x) {
      round(x)
    }
  }
  # Alter mapping files
  mf_con_temp <- mf_con %>%
    rename(dep=paste0(dep)) %>%
    mutate(dep=trans_func(dep)) 
  mf_treat_temp <- mf_treat %>%
    rename(dep=paste0(dep)) %>%
    mutate(dep=trans_func(dep)) 
  # Look at distributions according to models
  samps <- rstan::extract(g_lmer$stanfit)
  pre_test_set <- mf_treat_temp %>%
    filter(time<6)
  # If time factor present
  if ( time_factor & intercept_present & time_inter ) {
    # t_var <- alt_samps(mean(samps$beta[,5]) + colMeans(samps$beta[,6:10])) ### NEED TO CHECK THIS
    t_var <- (mean(samps$beta[,5]) + c(0,colMeans(samps$beta[,6:9]))) ### NEED TO CHECK THIS
  } else if ( time_factor & intercept_present & !time_inter ) {
    # t_var <- rep(alt_samps(mean(samps$beta[,5])), 5)
    t_var <- rep((mean(samps$beta[,5])), 5)
  } else if (time_factor  & !intercept_present & time_inter) {
    # t_var <- alt_samps(mean(samps$beta[,6]) + c(0, colMeans(samps$beta[,7:10])))
    t_var <- (mean(samps$beta[,6]) + c(0, colMeans(samps$beta[,7:10])))
  } else if (!time_factor  & intercept_present & time_inter) {
    # t_var <- alt_samps(mean(samps$beta[,5]) + c(0, colMeans(samps$beta[,6:9])))
    t_var <- (mean(samps$beta[,5]) + c(0, colMeans(samps$beta[,6:9])))
    # t_var <- alt_samps(colMeans(samps$beta[,5:9]))
  } else if (!time_factor  & !intercept_present & time_inter) {
    # t_var <- alt_samps(colMeans(samps$beta[,6:10]))
    t_var <- (colMeans(samps$beta[,6:10]))
  } else if (time_factor  & !intercept_present & !time_inter) {
    # t_var <- rep(alt_samps(mean(samps$beta[,6])), 5)
    t_var <- rep((mean(samps$beta[,6])), 5)
  } else {
    t_var <- rep(0,5)
  }
  names(t_var) <- c("Anbo","Rhma","Osse","Raca","Rapi")
  
  ### NEED TO ADD TIME TO CON/TREAT MF
  mf_con_temp <- mf_con_temp %>% mutate(dep_wtime = alt_samps(inv_alt_samps(dep)-time*t_var[species])) # minus time because we want to SUBTRACT effect of time 
  pre_test_set <- pre_test_set %>% mutate(dep_wtime = alt_samps(inv_alt_samps(dep)-time*t_var[species]))

  ## PLOT OF EXPECTED FOR EACH SPECIES
  if ( intercept_present ) {
    samps$beta <- samps$beta %>%
      as.data.frame() %>%
      mutate(Anbo=samps$alpha) %>%
      mutate(Rhma=alt_samps(Anbo+V1), Osse=alt_samps(Anbo+V2), Raca=alt_samps(Anbo+V3), Rapi=alt_samps(Anbo+V4)) %>%
      dplyr::select(Anbo, Rhma, Osse, Raca, Rapi) %>%
      rename(V1=Anbo, V2=Rhma, V3=Osse, V4=Raca, V5=Rapi) %>%
      as.matrix()
    
    gg_distr <- samps$beta %>%
      as.data.frame() %>%
      rename(Anbo=V1, Rhma=V2, Osse=V3, Raca=V4, Rapi=V5) %>%
      dplyr::select(Anbo,Rhma,Osse,Raca,Rapi) %>%
      gather(key=species, value=dependent_var) %>%
      mutate(species=factor(species, levels=c("Anbo","Rhma","Osse","Raca","Rapi"))) %>%
      ggplot(mapping=aes(x=species, y=dependent_var))+
      geom_violin() +
      geom_point(data=mf_con_temp, aes(y=dep_wtime, x=species), position = position_jitter(width = 0.1, height=0), col="blue") +
      geom_point(data=pre_test_set, aes(y=dep_wtime, x=species), position=position_jitter(width = 0.1, height=0), col="red") +
      ylab(label=paste0(name_dep))
    
  } else {
    gg_distr <- samps$beta %>%
      as.data.frame() %>%
      rename(Anbo=V1, Rhma=V2, Osse=V3, Raca=V4, Rapi=V5) %>%
      mutate(Anbo=alt_samps(Anbo), Rhma=alt_samps(Rhma), Osse=alt_samps(Osse), Raca=alt_samps(Raca), Rapi=alt_samps(Rapi)) %>%
      dplyr::select(Anbo,Rhma,Osse,Raca,Rapi) %>%
      gather(key=species, value=depenent_var) %>%
      mutate(species=factor(species, levels=c("Anbo","Rhma","Osse","Raca","Rapi"))) %>%
      ggplot(mapping=aes(x=species, y=depenent_var))+
      geom_violin() +
      geom_point(data=mf_con_temp, aes(y=dep_wtime, x=species), position = position_jitter(width = 0.1, height=0), col="blue") +
      geom_point(data=pre_test_set, aes(y=dep_wtime, x=species), position=position_jitter(width = 0.1, height=0), col="red") +
      ylab(label=paste0(name_dep))
    # Blue is controls, which the model is base don
    # Red is treatment individuals, which we are testing
  }
  
  ## Get standard deviation between toad individuals and samples
  ## For poisson, need to account for overdispersion
  # sd (variation) between samples
  samp_sigma <- samps$aux
  if ( overdispersion ) { 
    # sd (variation) between individuals of the same species
    indiv_intercept <- ranef(g_lmer)$indivID
    indivID_sigma <- sd(samps$b[,ncol(samps$b)])
    # sd (variation) between samples taken (overdispersion)
    samp_intercept <- ranef(g_lmer)$SampleID
    sampleID_sigma <- sd(samps$b[,nrow(samp_intercept)+1])
    # get toad_samps, which is change per toad
    # Get standard deviation between toad individuals and samples
    samp_indiv <- samps$b[,(nrow(samp_intercept)+2):(nrow(samp_intercept)+1+nrow(indiv_intercept))]
    colnames(samp_indiv) <- rownames(indiv_intercept)
    
  } else {
    # sd (variation) between individuals of the same species
    indivID_sigma <- sd(samps$b[,ncol(samps$b)])
    sampleID_sigma <- 0
  }
  
  
  # Now, we can calculate the probability that the "test" dataset values come from this distribution
  # List of individuals
  treat_indiv <- unique(mf_treat_temp$indivID)
  # List of each species
  species_list <- levels(factor(mf_con_temp$species))
  
  
  species_key <- treat_indiv %>%
    as_tibble() %>%
    rename(indivID=value) %>%
    separate(indivID,into=c("species","indiv"), remove=FALSE)
  
  exp_distr <- as.data.frame(matrix(ncol=length(species_list), nrow=4000, dimnames = list(1:4000, species_list)))
  for ( num_sp in 1:length(species_list)) {
    if ( fit_distr == "normal" ) {
      # Expected value at each species level
      exp_sp <- rnorm(length(samps$beta[,num_sp]), mean=samps$beta[,num_sp], sd=indivID_sigma)
      # Expected value after "sampling"
      exp_distr[,num_sp] <- rnorm(length(samps$beta[,num_sp]), mean=exp_sp, sd=samp_sigma)
    } else if ( fit_distr == "Gamma" ) {
      # Expected value at each species level
      exp_sp <- rnorm(length(samps$beta[,num_sp]), mean=samps$beta[,num_sp], sd=indivID_sigma)
      # Expected value after "sampling"
      exp_distr[,num_sp] <- rgamma(n=length(samps$beta[,num_sp]), shape=exp_sp, scale=1) # pretty sure default is 1
    } else if  ( fit_distr == "Beta" ) {
      
      ## TEST
      mu <- alt_samps(rnorm(4000, mean=samps$beta[,num_sp], sd=indivID_sigma))
      
      exp_distr[,num_sp] <-  rbeta(length(samps$beta[,num_sp])
                                   ,shape1=a(mu,(samps$aux))
                                   ,shape2=b(mu,(samps$aux)))
      
      # # ORIGINAL
      # mu <- alt_samps(rnorm(4000, mean=samps$beta[,num_sp], sd=indivID_sigma))
      # # exp_distr[,nump_sp] <- mu
      # exp_distr[,num_sp] <- rbeta(length(samps$beta[,num_sp])
      #                             ,shape1=a(mu,samps$aux)
      #                             ,shape2=b(mu,samps$aux))
      # 
      
    } else if ( fit_distr == "Poisson" ) {
      mu <- rnorm(length(samps$beta[,num_sp]), mean=samps$beta[,num_sp], sd=sampleID_sigma)
      # rpois(length(samps$beta[,num_sp]), lambda=mu)
      exp_distr[,num_sp] <- rpois(length(samps$beta[,num_sp]), lambda=mu)
      
    }
    
  }
  
  # Loop through and calculate probability of having diversity at that level
  pre_exp_indiv <- data.frame(indivID=treat_indiv, exp=rep(NA, length(treat_indiv)), p=rep(NA, length(treat_indiv)), infect=rep(NA, length(treat_indiv)))
  for ( i in treat_indiv ) {
    n_row <- match(i, treat_indiv)
    sp <- unlist(strsplit(i,"_"))
    num_sp <- match(sp[1], levels(factor(mf_con_temp$species)))
    temp_dep <- mf_treat_temp %>%
      filter(indivID==i, time <=5 ) %>%
      dplyr::select(dep) %>%
      pull()
    temp_time <- mf_treat_temp %>%
      filter(indivID==i, time <=5 ) %>%
      dplyr::select(time) %>%
      pull()
    
    x.fit <- alt_samps(inv_alt_samps(temp_dep)-temp_time*t_var[sp[1]])
    x.fit <- x.fit[!is.na(x.fit)]
    # I don't use a beta distribution here because often times optimization fails.
    # I think the sample size is too small to adequeately estimate shape1 and shape2?
    if ( fit_distr != "Beta" ) {
      if ( length(x.fit)>1 ) {
        fitted <- fitdistr(x.fit, densfun = paste0(fit_distr))$estimate
        if ( fit_distr!="Gamma") {
          exp <- fitted[1]
        } else {
          exp <-  fitted[1]/fitted[2]
          
        }
      } else if (length(x.fit) == 1 ) {
        exp <- x.fit
      } else {
        exp <- NA
      }
    } else {
      if ( length(x.fit)>1 ) {
        exp <-  fitdistr(x.fit, densfun = "normal")$estimate[1]
      } else if (length(x.fit) == 1) {
        exp <- (x.fit)
      } else {
        exp <- NA
      }
    }
    
    
    # pred_distr <- rnorm(length(samps_lmer_shannon$beta[,num_sp]), mean=rnorm(length(samps_lmer_shannon$beta[,num_sp]), mean=samps_lmer_shannon$beta[,num_sp], sd=toadID_sigma), sd=samp_sigma)
    p <- sum(exp_distr[,sp[1]]<exp)/length(exp_distr[,sp[1]])
    
    ### Did they get infected?
    infect <- max(mf_treat %>%
                    filter(indivID==i) %>%
                    dplyr::select(Bd_load) %>%
                    pull()
    )
    
    pre_exp_indiv[n_row,c("exp","p","infect")] <- c(exp, p, infect)
    
  }
  
  if ( intercept_present ) {
    # Get estimates for control toads
    species_exp <- fixef(g_lmer)[1:5]  %>%
      tibble::enframe() %>%
      dplyr::select(value) %>%
      t() %>%
      as.data.frame() %>%
      mutate("Anbo"=V1, "Rhma"=V1+V2, "Osse"=V1+V3, "Raca"=V1+V4, "Rapi"=V1+V5) %>%
      dplyr::select(Anbo, Rhma, Osse, Raca, Rapi) %>%
      t() %>%
      as.data.frame()  %>%
      mutate(species=c("Anbo","Rhma","Osse","Raca","Rapi")) %>%
      rename(value=V1)
  } else {
    # Get estimates for control toads
    species_exp <- fixef(g_lmer)[1:5]  %>%
      tibble::enframe() %>%
      mutate(species=c("Anbo", "Rhma","Osse","Raca","Rapi"))
  }
  
  # Need to extract ranef manually for beta; not sure why fixef works and ranef doesn't
  if ( fit_distr == "Beta" ) {
    randeffects_dist <- g_lmer$coefficients[grep("Intercept", names(g_lmer$coefficients))]
    names_randef <- names(randeffects_dist)
    names(randeffects_dist) <- gsub("]","",gsub("^.* indivID:", "", x=names_randef), fixed=TRUE)
    temp_randeffects <- data.frame(X=as.numeric(randeffects_dist), row.names = names(randeffects_dist))
    temp_randeffects <- temp_randeffects %>%
      rename(sp_est="X")
  } else {
    temp_randeffects <- ranef(g_lmer)$indivID %>%
      rename(sp_est="(Intercept)")
  }
  
  con_toad_est <- temp_randeffects %>%
    mutate(indivID=rownames(temp_randeffects)) %>%
    separate(indivID, into=c("species","num"), sep="_",remove=FALSE) %>%
    dplyr::select(-num) %>%
    left_join(species_exp, by = "species") %>%
    mutate(est=sp_est+value)
  
  con_exp_indiv <- data.frame(indivID=con_toad_est$indivID, exp=rep(NA, length(con_toad_est$indivID)), p=rep(NA, length(con_toad_est$indivID)))
  for ( i in 1:nrow(con_toad_est) ) {
    s <- con_toad_est[i,"species"]
    
    exp <- alt_samps(con_toad_est[i,"est"])
    p <- sum(exp_distr[,s]<exp)/length(exp_distr[,s])
    
    con_exp_indiv[i,c("exp","p")] <- c(exp, p)
  }
  con_exp_indiv <- con_exp_indiv %>%
    separate(indivID, into=c("species","indiv"), remove=FALSE)
  
  
  # create species column
  pre_exp_indiv <- pre_exp_indiv %>%
    separate(indivID, into=c("species","indiv"), remove = FALSE) %>%
    mutate(species=factor(species, levels=c("Anbo","Rhma","Osse","Raca","Rapi")))
  # Plot results 
  gg_p <- ggplot(pre_exp_indiv, aes(x=p, y=log(infect+1))) +
    geom_point(aes(color=species), cex=4) +
    geom_smooth(aes(color=species),method=lm, se = FALSE) +
    geom_smooth(method=lm, se=FALSE, col="black") +
    xlab(label=paste0("p_",name_dep))
  # if we'd JUST plotted raw values
  gg_raw <- ggplot(pre_exp_indiv, aes(x=exp, y=log(infect+1)))+
    geom_point(aes(color=species), cex=4) +
    geom_smooth(method=lm, se=FALSE, col="black") +
    geom_smooth(aes(color=species),method=lm, se = FALSE) +
    xlab(label=paste0("exp_",name_dep))
  #grid.arrange(gg_p, gg_raw, nrow=1)
  
  gg_exp_distr <- exp_distr %>%
    # add time adjustment to raw data values
    gather(key=species, value=dep_var) %>%
    mutate(species=factor(species,levels=c("Anbo","Rhma","Osse","Raca","Rapi")))%>%
    ggplot(aes(x=species, y=dep_var)) +
    geom_violin() +
    geom_point(data=pre_exp_indiv, aes(x=species, y=exp, col=log(infect+1)), cex=4, position=position_jitter(height=0, width=0.1)) +
    ylab(label=paste0(name_dep))
  
  gg_exp_distr_controls <- exp_distr %>%
    # add time adjustment to raw data values
    gather(key=species, value=dep_var) %>%
    mutate(species=factor(species,levels=c("Anbo","Rhma","Osse","Raca","Rapi")))%>%
    ggplot(aes(x=species, y=dep_var)) +
    geom_violin() +
    geom_point(data=con_exp_indiv, aes(x=species, y=exp), cex=3, position=position_jitter(height=0, width=0.1)) +
    ylab(label=paste0(name_dep))
  
  all_p <- pre_exp_indiv %>%
    dplyr::select(indivID, exp, p,infect) 
  
  # Make list of items
  output <- list()
  output[["Var"]] <- data.frame( dep=dep, name_dep=name_dep
                                 , transform_raw_func=transform_raw_func
                                 , intercept_present=intercept_present
                                 ,fit_distr=fit_distr)
  output[["gg_model_Distribution_of_all_values"]] <- gg_distr
  output[["Control_individuals"]] <- con_exp_indiv
  output[["gg_p"]] <- gg_p
  output[["gg_raw"]] <- gg_raw
  output[["gg_ExpectedDistribution_and_Bd_exposed"]] <- gg_exp_distr
  output[["gg_ExpectedDistribution_controls"]] <- gg_exp_distr_controls
  output[["all_p"]] <- all_p
  
  return(output)
}



process_glmer_all <- function(g_lmer, dep, name_dep, transform_raw_func, intercept_present, fit_distr, time_factor=FALSE, overdispersion=FALSE, mf_con, mf_treat) {
  
  # Need overdispersion term for poisson
  #  g_lmer = lmer_log_observed_otus_all
  #  dep= "observed_otus"
  #   name_dep="log_observed_otus"
  #   transform_raw_func="log"
  #   intercept_present=F
  # fit_distr="normal"
  #   mf_con=mf_con
  #   mf_treat=mf_treat
  #   time_factor=FALSE
  #   overdispersion=FALSE
  
  if (transform_raw_func == "log") {
    trans_func <- log
    alt_samps <- function(x) {return(x)} 
  } else if (transform_raw_func == "None") {
    trans_func <- function(x) {return(x)} 
    alt_samps <- function(x) {return(x)} 
  } else if (transform_raw_func == "beta") {
    trans_func <- function(x) {return(x)}
    a <- function(mu,phi){
      mu*phi
    }
    b <- function(mu,phi) {
      phi-mu*phi
    }
    mu <- function(a,phi) {
      a/phi
    }
    alt_samps <- function(x) {
      exp(x)/(exp(x)+1)
    }
    logit <- function(p) {
      log(p/(1-p))
    }
    mu_ab <- function(a,b) {
      a/(a+b)
    }
  } 
  
  # Alter mapping files
  mf_con_temp <- mf_con %>%
    rename(dep=paste0(dep)) %>%
    mutate(dep=trans_func(dep)) 
  mf_treat_temp <- mf_treat %>%
    rename(dep=paste0(dep)) %>%
    mutate(dep=trans_func(dep)) 
  # combine pre and post mapping files
  mf_uninfected <- mf_treat_temp %>%
    filter(time<6) %>%
    rbind(mf_con_temp)
  # Look at distributions according to models
  samps <- rstan::extract(g_lmer$stanfit)
  post_test_set <- mf_treat_temp %>%
    filter(time>=6) 
  
  # If time factor present
  if ( time_factor & intercept_present ) {
    t_var <- alt_samps(mean(samps$beta[,5]))
  } else if (time_factor) {
    t_var <- alt_samps(mean(samps$beta[,6]))
  } else {
    t_var <- 0
  }
  
  ## PLOT OF EXPECTED FOR EACH SPECIES
  if ( intercept_present ) {
    samps$beta <- samps$beta %>%
      as.data.frame() %>%
      mutate(Anbo=samps$alpha) %>%
      mutate(Rhma=alt_samps(Anbo+V1), Osse=alt_samps(Anbo+V2), Raca=alt_samps(Anbo+V3), Rapi=alt_samps(Anbo+V4)) %>%
      dplyr::select(Anbo, Rhma, Osse, Raca, Rapi) %>%
      rename(V1=Anbo, V2=Rhma, V3=Osse, V4=Raca, V5=Rapi) %>%
      as.matrix()
    gg_distr <- samps$beta %>%
      as.data.frame() %>%
      rename(Anbo=V1, Rhma=V2, Osse=V3, Raca=V4, Rapi=V5) %>%
      dplyr::select(Anbo,Rhma,Osse,Raca,Rapi) %>%
      gather(key=species, value=dependent_var) %>%
      mutate(species=factor(species, levels=c("Anbo","Rhma","Osse","Raca","Rapi"))) %>%
      ggplot(mapping=aes(x=species, y=dependent_var))+
      geom_violin() +
      geom_point(data=mf_uninfected, aes(y=dep-time*(t_var), x=species), position = position_jitter(width = 0.1, height=0), col="blue") +
      geom_point(data=post_test_set, aes(y=dep-time*(t_var), x=species), position=position_jitter(width = 0.1, height=0), col="red") +
      ylab(label=paste0(name_dep))
    
  } else {
    gg_distr <- samps$beta %>%
      as.data.frame() %>%
      rename(Anbo=V1, Rhma=V2, Osse=V3, Raca=V4, Rapi=V5) %>%
      mutate(Anbo=alt_samps(Anbo), Rhma=alt_samps(Rhma), Osse=alt_samps(Osse), Raca=alt_samps(Raca), Rapi=alt_samps(Rapi)) %>%
      dplyr::select(Anbo,Rhma,Osse,Raca,Rapi) %>%
      gather(key=species, value=depenent_var) %>%
      mutate(species=factor(species, levels=c("Anbo","Rhma","Osse","Raca","Rapi"))) %>%
      ggplot(mapping=aes(x=species, y=depenent_var))+
      geom_violin() +
      geom_point(data=mf_uninfected, aes(y=dep-time*(t_var), x=species), position = position_jitter(width = 0.1, height=0), col="blue") +
      geom_point(data=post_test_set, aes(y=dep-time*(t_var), x=species), position=position_jitter(width = 0.1, height=0), col="red") +
      ylab(label=paste0(name_dep))
    # Blue is controls, which the model is base don
    # Red is treatment individuals, which we are testing
  }
  
  ## Get standard deviation between toad individuals and samples
  ## For poisson, need to account for overdispersion
  # sd (variation) between samples
  samp_sigma <- samps$aux
  if ( overdispersion ) { 
    # sd (variation) between individuals of the same species
    randeffects_dist <- g_lmer$coefficients[grep("Intercept", names(g_lmer$coefficients))]
    names_randef <- names(randeffects_dist)
    names(randeffects_dist) <- gsub("]","",gsub("^.* indivID:", "", x=names_randef), fixed=TRUE)
    
    indiv_intercept <- data.frame(Intercept=randeffects_dist)
    indivID_sigma <- sd(samps$b[,ncol(samps$b)])
    
    # sd (variation) between samples taken (overdispersion)
    samprandeffects_dist <- g_lmer$coefficients[grep("SampleID", names(g_lmer$coefficients))]
    names_samprandef <- names(samprandeffects_dist)
    names(samprandeffects_dist) <- gsub("]","",gsub("^.* indivID:", "", x=names_samprandef), fixed=TRUE)
    
    samp_intercept <- data.frame(Intercept_SampleID=samprandeffects_dist)
    sampleID_sigma <- sd(samps$b[,nrow(samp_intercept)+1])
    # get toad_samps, which is change per toad
    # Get standard deviation between toad individuals and samples
    samp_indiv <- samps$b[,(nrow(samp_intercept)+2):(nrow(samp_intercept)+1+nrow(indiv_intercept))]
    colnames(samp_indiv) <- rownames(indiv_intercept)
    
  } else {
    # Get standard deviation between toad individuals and samples
    randeffects_dist <- g_lmer$coefficients[grep("Intercept", names(g_lmer$coefficients))]
    names_randef <- names(randeffects_dist)
    names(randeffects_dist) <- gsub("]","",gsub("^.* indivID:", "", x=names_randef), fixed=TRUE)
    
    indiv_intercept <- data.frame(Intercept=randeffects_dist)
    samp_indiv <- samps$b[,1:nrow(indiv_intercept)]
    colnames(samp_indiv) <- rownames(indiv_intercept)
    
    # sd (variation) between individuals of the same species
    indivID_sigma <- sd(samps$b[,ncol(samps$b)])
    sampleID_sigma <- 0
  }
  
  
  # Now, we can calculate the probability that the "test" dataset values come from this distribution
  # List of individuals-- including everything
  all_indiv <- unique(rbind(mf_con_temp, mf_treat_temp)$indivID)
  # List of each species
  species_list <- levels(factor(mf_con_temp$species))
  
  
  species_key <- all_indiv %>%
    as_tibble() %>%
    rename(indivID=value) %>%
    separate(indivID,into=c("species","indiv"), remove=FALSE)
  species_order <- levels(as.factor(mf_uninfected$species))
  
  exp_distr <- as.data.frame(matrix(ncol=length(all_indiv), nrow=4000, dimnames = list(1:4000, all_indiv)))
  for ( num_indiv in 1:length(all_indiv)) {
    indiv <- all_indiv[num_indiv]
    sp <- pull(species_key[num_indiv,"species"])
    num_sp <- match(sp, species_order)
    
    if ( fit_distr == "normal" ) {
      # Expected value at each species level
      exp_sp <- rnorm(length(samps$beta[,num_sp]), mean=samps$beta[,num_sp], sd=indivID_sigma)
      # Expected value after "sampling"
      exp_distr[,num_indiv] <- rnorm(length(samps$beta[,num_sp]), mean=exp_sp+samp_indiv[,indiv], sd=samp_sigma)
    } else if ( fit_distr == "Gamma" ) {
      # Expected value at each species level
      exp_sp <- rnorm(length(samps$beta[,num_sp]), mean=samps$beta[,num_sp], sd=indivID_sigma)
      # Expected value after "sampling"
      exp_distr[,num_indiv] <- rgamma(n=length(samps$beta[,num_sp]), shape=exp_sp+samp_indiv[,indiv], scale=1) # pretty sure default is 1
    } else if  ( fit_distr == "Beta" ) {
      
      mu <- alt_samps(rnorm(4000, mean=samps$beta[,num_sp], sd=indivID_sigma) +samp_indiv[,indiv])
      
      exp_distr[,num_indiv] <-  rbeta(length(samps$beta[,num_sp])
                                      ,shape1=a(mu,(samps$aux))
                                      ,shape2=b(mu,(samps$aux)))
      
      
    } else if ( fit_distr == "Poisson" ) {
      mu <- rnorm(length(samps$beta[,num_sp]), mean=samps$beta[,num_sp], sd=sampleID_sigma) + +samp_indiv[,indiv]
      exp_distr[,num_indiv] <- rpois(length(samps$beta[,num_sp]), lambda=mu)
      
    }
    
  }
  
  # Loop through and calculate probability of having diversity at that level
  pos_exp_indiv <- mf_treat_temp %>%
    as.data.frame() %>%
    filter(time>5, Bd_exposure=="Bd-exposed") %>%
    dplyr::select(indivID, species, Bd_load) %>%
    group_by(indivID, species) %>%
    summarize(max_Bd_load=max(Bd_load)) %>%
    mutate(p=NA, exp=NA)
  ##################### PASTED ########
  for ( i in pos_exp_indiv$indivID ) {
    #n_row <- match(i, treat_indiv)
    sp <- unlist(strsplit(i,"_"))
    num_sp <- match(sp[1], levels(factor(mf_treat_temp$species)))
    temp_dep <- mf_treat_temp %>%
      filter(indivID==i, time >5 ) %>%
      dplyr::select(dep) %>%
      pull()
    temp_time <- mf_treat_temp %>%
      filter(indivID==i, time >5 ) %>%
      dplyr::select(time) %>%
      pull()
    
    x.fit <- temp_dep-temp_time*mean(t_var)
    x.fit <- x.fit[!is.na(x.fit)]
    # I don't use a beta distribution here because often times optimization fails.
    # I think the sample size is too small to adequeately estimate shape1 and shape2?
    if ( fit_distr != "Beta" ) {
      if ( length(x.fit)>1 ) {
        fitted <- fitdistr(x.fit, densfun = paste0(fit_distr))$estimate
        if ( fit_distr!="Gamma") {
          exp <- fitted[1]
        } else {
          exp <-  fitted[1]/fitted[2]
          
        }
      } else if (length(x.fit) == 1 ) {
        exp <- x.fit
      } else {
        exp <- NA
      }
    } else {
      if ( length(x.fit)>1 ) {
        exp <-  fitdistr(x.fit, densfun = "normal")$estimate[1]
      } else if (length(x.fit) == 1) {
        exp <- (x.fit)
      } else {
        exp <- NA
      }
    }
    
    
    # pred_distr <- rnorm(length(samps_lmer_shannon$beta[,num_sp]), mean=rnorm(length(samps_lmer_shannon$beta[,num_sp]), mean=samps_lmer_shannon$beta[,num_sp], sd=toadID_sigma), sd=samp_sigma)
    p <- sum(exp_distr[,i]<exp)/length(exp_distr[,i])
    
    ### Did they get infected?
    # infect <- max(mf_treat %>%
    #                 filter(indivID==i) %>%
    #                 dplyr::select(Bd_load) %>%
    #                 pull()
    # )
    
    pos_exp_indiv[match(i, pos_exp_indiv$indivID),c("exp","p")] <- c(exp, p)
    
  }
  
  # Get all control individuals and individuals before exposure
  con_exp_indiv <- data.frame(indivID=mf_uninfected$indivID, exp=rep(NA, length(mf_uninfected$indivID)), p=rep(NA, length(mf_uninfected$indivID)), time=mf_uninfected$time)
  for ( i in 1:nrow(mf_uninfected) ) {
    indiv <- pull(mf_uninfected[i,"indivID"])
    exp <- pull(mf_uninfected[i,"dep"])-pull(mf_uninfected[i,"time"])*t_var
    
    p <- sum(exp_distr[,indiv]<exp)/length(exp_distr[,indiv])
    
    con_exp_indiv[i,c("exp","p")] <- c(exp, p)
  }
  con_exp_indiv <- con_exp_indiv %>%
    separate(indivID, into=c("species","indiv"), remove=FALSE)
  
  # Plot results 
  gg_p <- ggplot(pos_exp_indiv, aes(x=p, y=max_Bd_load)) +
    geom_point(aes(color=species), cex=4) +
    geom_smooth(aes(color=species),method=lm, se = FALSE) +
    geom_smooth(method=lm, se=FALSE, col="black") +
    xlab(label=paste0("p_",name_dep))
  # if we'd JUST plotted raw values
  gg_raw <- ggplot(pos_exp_indiv, aes(x=exp, y=max_Bd_load))+
    geom_point(aes(color=species), cex=4) +
    geom_smooth(method=lm, se=FALSE, col="black") +
    geom_smooth(aes(color=species),method=lm, se = FALSE) +
    xlab(label=paste0("exp_",name_dep))
  #grid.arrange(gg_p, gg_raw, nrow=1)
  
  gg_exp_distr <- exp_distr %>%
    # add time adjustment to raw data values
    gather(key=indivID, value=dep_var) %>%
    separate(indivID, into=c("species","indiv"), remove=FALSE) %>%
    mutate(species=factor(species,levels=c("Anbo","Rhma","Osse","Raca","Rapi")))%>%
    ggplot(aes(x=species, y=dep_var)) +
    geom_violin() +
    geom_point(data=pos_exp_indiv, aes(x=species, y=exp, col=max_Bd_load),  cex=2, position=position_jitter(height=0, width=0.1)) +
    ylab(label=paste0(name_dep))
  
  gg_exp_distr_controls <- exp_distr %>%
    # add time adjustment to raw data values
    gather(key=indivID, value=dep_var) %>%
    separate(indivID, into=c("species","indiv"), remove=FALSE) %>%
    mutate(species=factor(species,levels=c("Anbo","Rhma","Osse","Raca","Rapi")))%>%
    ggplot(aes(x=species, y=dep_var)) +
    geom_violin() +
    geom_point(data=con_exp_indiv, aes(x=species, y=exp),  cex=2, position=position_jitter(height=0, width=0.1)) +
    ylab(label=paste0(name_dep))
  
  all_p <- pos_exp_indiv %>%
    rename(infect=max_Bd_load) %>%
    dplyr::select(indivID, exp, p, infect) 
  
  # Make list of items
  output <- list()
  output[["Var"]] <- data.frame( dep=dep, name_dep=name_dep
                                 , transform_raw_func=transform_raw_func
                                 , intercept_present=intercept_present
                                 ,fit_distr=fit_distr)
  output[["gg_model_Distribution_of_all_values"]] <- gg_distr
  output[["gg_p"]] <- gg_p
  output[["gg_raw"]] <- gg_raw
  output[["gg_ExpectedDistribution_and_Bd_exposed"]] <- gg_exp_distr
  output[["gg_ExpectedDistribution_controls"]] <- gg_exp_distr_controls
  output[["all_p"]] <- all_p
  
  return(output)
}


process_glmer_all_wtime <- function(g_lmer, dep, name_dep, transform_raw_func, intercept_present, fit_distr, time_factor=FALSE, overdispersion=FALSE, mf_con, mf_treat, time_inter = FALSE) {
  # g_lmer = lmer_log_observed_otus_all_wtime
  # dep= "observed_otus"
  # name_dep="log_observed_otus"
  # transform_raw_func="log"
  # intercept_present=F
  # fit_distr="normal"
  # mf_con=mf_con
  # mf_treat=mf_treat
  # time_factor = T
  # time_inter = T
  # 
  # g_lmer = glmer_faith_pd_all_wtime
  # dep= "faith_pd"
  # name_dep="faith_pd"
  # transform_raw_func="None"
  # intercept_present=T
  # fit_distr="Gamma"
  # mf_con = mf_con
  # mf_treat=mf_treat
  # time_factor=T
  # time_inter = T

  if (transform_raw_func == "log") {
    trans_func <- log
    alt_samps <- function(x) {return(x)} 
    inv_alt_samps <- function(x) {return(x)}
  } else if (transform_raw_func == "None") {
    trans_func <- function(x) {return(x)} 
    alt_samps <- function(x) {return(x)} 
    inv_alt_samps <- function(x) {return(x)}
  } else if (transform_raw_func == "beta") {
    trans_func <- function(x) {return(x)}
    a <- function(mu,phi){
      mu*phi
    }
    b <- function(mu,phi) {
      phi-mu*phi
    }
    mu <- function(a,phi) {
      a/phi
    }
    alt_samps <- function(x) {
      exp(x)/(exp(x)+1)
    }
    logit <- function(p) {
      log(p/(1-p))
    }
    inv_alt_samps <- logit
    mu_ab <- function(a,b) {
      a/(a+b)
    }
  } 
  
  if (fit_distr == "Poisson" ) {
    alt_samps <- function(x) {
      round(x)
    }
  }
  
  if (fit_distr == "Gamma") {
    alt_samps <- function(x) {
      for ( i in 1:length(x)) {
        if (x[i]<0) {
          x[i] <- 0.01
        }
      }
      return(x)
    }
    
  }
  
  # Alter mapping files
  mf_con_temp <- mf_con %>%
    rename(dep=paste0(dep)) %>%
    mutate(dep=trans_func(dep)) 
  mf_treat_temp <- mf_treat %>%
    rename(dep=paste0(dep)) %>%
    mutate(dep=trans_func(dep)) 
  # combine pre and post mapping files
  mf_uninfected <- mf_treat_temp %>%
    filter(time<6) %>%
    rbind(mf_con_temp)
  # Look at distributions according to models
  samps <- rstan::extract(g_lmer$stanfit)
  post_test_set <- mf_treat_temp %>%
    filter(time>=6) 
  
  # If time factor present
  if ( time_factor & intercept_present & time_inter ) {
    # t_var <- alt_samps(mean(samps$beta[,5]) + colMeans(samps$beta[,6:10])) ### NEED TO CHECK THIS
    t_var <- (mean(samps$beta[,5]) + c(0,colMeans(samps$beta[,6:9]))) ### NEED TO CHECK THIS
  } else if ( time_factor & intercept_present & !time_inter ) {
    # t_var <- rep(alt_samps(mean(samps$beta[,5])), 5)
    t_var <- rep((mean(samps$beta[,5])), 5)
  } else if (time_factor  & !intercept_present & time_inter) {
    # t_var <- alt_samps(mean(samps$beta[,6]) + c(0, colMeans(samps$beta[,7:10])))
    t_var <- (mean(samps$beta[,6]) + c(0, colMeans(samps$beta[,7:10])))
  } else if (!time_factor  & intercept_present & time_inter) {
    # t_var <- alt_samps(mean(samps$beta[,5]) + c(0, colMeans(samps$beta[,6:9])))
    t_var <- (mean(samps$beta[,5]) + c(0, colMeans(samps$beta[,6:9])))
    # t_var <- alt_samps(colMeans(samps$beta[,5:9]))
  } else if (!time_factor  & !intercept_present & time_inter) {
    # t_var <- alt_samps(colMeans(samps$beta[,6:10]))
    t_var <- (colMeans(samps$beta[,6:10]))
  } else if (time_factor  & !intercept_present & !time_inter) {
    # t_var <- rep(alt_samps(mean(samps$beta[,6])), 5)
    t_var <- rep((mean(samps$beta[,6])), 5)
  } else {
    t_var <- rep(0,5)
  }
  names(t_var) <- c("Anbo","Rhma","Osse","Raca","Rapi")
  
  ### NEED TO ADD TIME TO CON/TREAT MF
  mf_uninfected <- mf_uninfected %>% mutate(dep_wtime = alt_samps(inv_alt_samps(dep)-time*t_var[species])) # minus time because we want to SUBTRACT effect of time 
  post_test_set <- post_test_set %>% mutate(dep_wtime = alt_samps(inv_alt_samps(dep)-time*t_var[species]))
  
  
  ## PLOT OF EXPECTED FOR EACH SPECIES
  if ( intercept_present ) {
    samps$beta <- samps$beta %>%
      as.data.frame() %>%
      mutate(Anbo=samps$alpha) %>%
      mutate(Rhma=alt_samps(Anbo+V1), Osse=alt_samps(Anbo+V2), Raca=alt_samps(Anbo+V3), Rapi=alt_samps(Anbo+V4)) %>%
      dplyr::select(Anbo, Rhma, Osse, Raca, Rapi) %>%
      rename(V1=Anbo, V2=Rhma, V3=Osse, V4=Raca, V5=Rapi) %>%
      as.matrix()
    gg_distr <- samps$beta %>%
      as.data.frame() %>%
      rename(Anbo=V1, Rhma=V2, Osse=V3, Raca=V4, Rapi=V5) %>%
      dplyr::select(Anbo,Rhma,Osse,Raca,Rapi) %>%
      gather(key=species, value=dependent_var) %>%
      mutate(species=factor(species, levels=c("Anbo","Rhma","Osse","Raca","Rapi"))) %>%
      ggplot(mapping=aes(x=species, y=dependent_var))+
      geom_violin() +
      geom_point(data=mf_uninfected, aes(y=dep_wtime, x=species), position = position_jitter(width = 0.1, height=0), col="blue") +
      geom_point(data=post_test_set, aes(y=dep_wtime, x=species), position=position_jitter(width = 0.1, height=0), col="red") +
      ylab(label=paste0(name_dep))
    
  } else {
    gg_distr <- samps$beta %>%
      as.data.frame() %>%
      rename(Anbo=V1, Rhma=V2, Osse=V3, Raca=V4, Rapi=V5) %>%
      mutate(Anbo=alt_samps(Anbo), Rhma=alt_samps(Rhma), Osse=alt_samps(Osse), Raca=alt_samps(Raca), Rapi=alt_samps(Rapi)) %>%
      dplyr::select(Anbo,Rhma,Osse,Raca,Rapi) %>%
      gather(key=species, value=depenent_var) %>%
      mutate(species=factor(species, levels=c("Anbo","Rhma","Osse","Raca","Rapi"))) %>%
      ggplot(mapping=aes(x=species, y=depenent_var))+
      geom_violin() +
      geom_point(data=mf_uninfected, aes(y=dep_wtime, x=species), position = position_jitter(width = 0.1, height=0), col="blue") +
      geom_point(data=post_test_set, aes(y=dep_wtime, x=species), position=position_jitter(width = 0.1, height=0), col="red") +
      ylab(label=paste0(name_dep))
    # Blue is controls, which the model is base don
    # Red is treatment individuals, which we are testing
  }
  
  ## Get standard deviation between toad individuals and samples
  ## For poisson, need to account for overdispersion
  # sd (variation) between samples
  samp_sigma <- samps$aux
  if ( overdispersion ) { 
    # sd (variation) between individuals of the same species
    randeffects_dist <- g_lmer$coefficients[grep("Intercept", names(g_lmer$coefficients))]
    names_randef <- names(randeffects_dist)
    names(randeffects_dist) <- gsub("]","",gsub("^.* indivID:", "", x=names_randef), fixed=TRUE)
    
    indiv_intercept <- data.frame(Intercept=randeffects_dist)
    indivID_sigma <- sd(samps$b[,ncol(samps$b)])
    
    # sd (variation) between samples taken (overdispersion)
    samprandeffects_dist <- g_lmer$coefficients[grep("SampleID", names(g_lmer$coefficients))]
    names_samprandef <- names(samprandeffects_dist)
    names(samprandeffects_dist) <- gsub("]","",gsub("^.* indivID:", "", x=names_samprandef), fixed=TRUE)
    
    samp_intercept <- data.frame(Intercept_SampleID=samprandeffects_dist)
    sampleID_sigma <- sd(samps$b[,nrow(samp_intercept)+1])
    # get toad_samps, which is change per toad
    # Get standard deviation between toad individuals and samples
    samp_indiv <- samps$b[,(nrow(samp_intercept)+2):(nrow(samp_intercept)+1+nrow(indiv_intercept))]
    colnames(samp_indiv) <- rownames(indiv_intercept)
    
  } else {
    # Get standard deviation between toad individuals and samples
    randeffects_dist <- g_lmer$coefficients[grep("Intercept", names(g_lmer$coefficients))]
    names_randef <- names(randeffects_dist)
    names(randeffects_dist) <- gsub("]","",gsub("^.* indivID:", "", x=names_randef), fixed=TRUE)
    
    indiv_intercept <- data.frame(Intercept=randeffects_dist)
    samp_indiv <- samps$b[,1:nrow(indiv_intercept)]
    colnames(samp_indiv) <- rownames(indiv_intercept)
    
    # sd (variation) between individuals of the same species
    indivID_sigma <- sd(samps$b[,ncol(samps$b)])
    sampleID_sigma <- 0
  }
  
  
  # Now, we can calculate the probability that the "test" dataset values come from this distribution
  # List of individuals-- including everything
  all_indiv <- unique(rbind(mf_con_temp, mf_treat_temp)$indivID)
  # List of each species
  species_list <- levels(factor(mf_con_temp$species))
  
  
  species_key <- all_indiv %>%
    as_tibble() %>%
    rename(indivID=value) %>%
    separate(indivID,into=c("species","indiv"), remove=FALSE)
  species_order <- levels(as.factor(mf_uninfected$species))
  
  exp_distr <- as.data.frame(matrix(ncol=length(all_indiv), nrow=4000, dimnames = list(1:4000, all_indiv)))
  for ( num_indiv in 1:length(all_indiv)) {
    indiv <- all_indiv[num_indiv]
    sp <- pull(species_key[num_indiv,"species"])
    num_sp <- match(sp, species_order)
    
    if ( fit_distr == "normal" ) {
      # Expected value at each species level
      exp_sp <- rnorm(length(samps$beta[,num_sp]), mean=samps$beta[,num_sp], sd=indivID_sigma)
      # Expected value after "sampling"
      exp_distr[,num_indiv] <- rnorm(length(samps$beta[,num_sp]), mean=exp_sp+samp_indiv[,indiv], sd=samp_sigma)
    } else if ( fit_distr == "Gamma" ) {
      # Expected value at each species level
      exp_sp <- rnorm(length(samps$beta[,num_sp]), mean=samps$beta[,num_sp], sd=indivID_sigma)
      # Expected value after "sampling"
      exp_distr[,num_indiv] <- rgamma(n=length(samps$beta[,num_sp]), shape=exp_sp+samp_indiv[,indiv], scale=1) # pretty sure default is 1
    } else if  ( fit_distr == "Beta" ) {
      
      mu <- alt_samps(rnorm(4000, mean=samps$beta[,num_sp], sd=indivID_sigma) +samp_indiv[,indiv])
      
      exp_distr[,num_indiv] <-  rbeta(length(samps$beta[,num_sp])
                                      ,shape1=a(mu,(samps$aux))
                                      ,shape2=b(mu,(samps$aux)))
      
      
    } else if ( fit_distr == "Poisson" ) {
      mu <- rnorm(length(samps$beta[,num_sp]), mean=samps$beta[,num_sp], sd=sampleID_sigma) + +samp_indiv[,indiv]
      exp_distr[,num_indiv] <- rpois(length(samps$beta[,num_sp]), lambda=mu)
      
    }
    
  }
  
  # Loop through and calculate probability of having diversity at that level
  pos_exp_indiv <- mf_treat_temp %>%
    as.data.frame() %>%
    filter(time>5, Bd_exposure=="Bd-exposed") %>%
    dplyr::select(indivID, species, Bd_load) %>%
    group_by(indivID, species) %>%
    summarize(max_Bd_load=max(Bd_load)) %>%
    mutate(p=NA, exp=NA) %>%
    mutate(p = as.numeric(p), exp = as.numeric(exp))
  ##################### PASTED ########
  for ( i in pos_exp_indiv$indivID ) {
    # i = "Osse_5"
    #n_row <- match(i, treat_indiv)
    sp <- unlist(strsplit(i,"_"))
    num_sp <- match(sp[1], levels(factor(mf_treat_temp$species)))
    temp_dep <- mf_treat_temp %>%
      filter(indivID==i, time >5 ) %>%
      dplyr::select(dep) %>%
      pull()
    temp_time <- mf_treat_temp %>%
      filter(indivID==i, time >5 ) %>%
      dplyr::select(time) %>%
      pull()
    
    x.fit <- alt_samps(inv_alt_samps(temp_dep)-temp_time*mean(t_var[sp[1]]))
    x.fit <- x.fit[!is.na(x.fit)]
    # I don't use a beta distribution here because often times optimization fails.
    # I think the sample size is too small to adequeately estimate shape1 and shape2?
    if ( fit_distr != "Beta" ) {
      if ( length(x.fit)>1 ) {
        fitted <- fitdistr(x.fit, densfun = paste0(fit_distr))$estimate
        if ( fit_distr!="Gamma") {
          exp <- fitted[1]
        } else {
          exp <-  fitted[1]/fitted[2]
          
        }
      } else if (length(x.fit) == 1 ) {
        exp <- x.fit
      } else {
        exp <- NA
      }
    } else {
      if ( length(x.fit)>1 ) {
        exp <-  fitdistr(x.fit, densfun = "normal")$estimate[1]
      } else if (length(x.fit) == 1) {
        exp <- (x.fit)
      } else {
        exp <- NA
      }
    }
    
    
    # pred_distr <- rnorm(length(samps_lmer_shannon$beta[,num_sp]), mean=rnorm(length(samps_lmer_shannon$beta[,num_sp]), mean=samps_lmer_shannon$beta[,num_sp], sd=toadID_sigma), sd=samp_sigma)
    p <- sum(exp_distr[,i]<exp)/length(exp_distr[,i])
    
    ### Did they get infected?
    # infect <- max(mf_treat %>%
    #                 filter(indivID==i) %>%
    #                 dplyr::select(Bd_load) %>%
    #                 pull()
    # )
    
    pos_exp_indiv[match(i, pos_exp_indiv$indivID),c("exp","p")] <- data.frame(exp =exp, p = p)
    
  }
  
  # Get all control individuals and individuals before exposure
  con_exp_indiv <- data.frame(indivID=mf_uninfected$indivID, exp=rep(NA, length(mf_uninfected$indivID)), p=rep(NA, length(mf_uninfected$indivID)), time=mf_uninfected$time)
  for ( i in 1:nrow(mf_uninfected) ) {
    indiv <- pull(mf_uninfected[i,"indivID"])
    exp <- pull(mf_uninfected[i,"dep"])-pull(mf_uninfected[i,"time"])*t_var
    
    p <- sum(exp_distr[,indiv]<exp)/length(exp_distr[,indiv])
    
    con_exp_indiv[i,c("exp","p")] <- c(exp, p)
  }
  con_exp_indiv <- con_exp_indiv %>%
    separate(indivID, into=c("species","indiv"), remove=FALSE)
  
  # Plot results 
  gg_p <- ggplot(pos_exp_indiv, aes(x=p, y=max_Bd_load)) +
    geom_point(aes(color=species), cex=4) +
    geom_smooth(aes(color=species),method=lm, se = FALSE) +
    geom_smooth(method=lm, se=FALSE, col="black") +
    xlab(label=paste0("p_",name_dep))
  # if we'd JUST plotted raw values
  gg_raw <- ggplot(pos_exp_indiv, aes(x=exp, y=max_Bd_load))+
    geom_point(aes(color=species), cex=4) +
    geom_smooth(method=lm, se=FALSE, col="black") +
    geom_smooth(aes(color=species),method=lm, se = FALSE) +
    xlab(label=paste0("exp_",name_dep))
  #grid.arrange(gg_p, gg_raw, nrow=1)
  
  gg_exp_distr <- exp_distr %>%
    # add time adjustment to raw data values
    gather(key=indivID, value=dep_var) %>%
    separate(indivID, into=c("species","indiv"), remove=FALSE) %>%
    mutate(species=factor(species,levels=c("Anbo","Rhma","Osse","Raca","Rapi")))%>%
    ggplot(aes(x=species, y=dep_var)) +
    geom_violin() +
    geom_point(data=pos_exp_indiv, aes(x=species, y=exp, col=max_Bd_load),  cex=2, position=position_jitter(height=0, width=0.1)) +
    ylab(label=paste0(name_dep))
  
  gg_exp_distr_controls <- exp_distr %>%
    # add time adjustment to raw data values
    gather(key=indivID, value=dep_var) %>%
    separate(indivID, into=c("species","indiv"), remove=FALSE) %>%
    mutate(species=factor(species,levels=c("Anbo","Rhma","Osse","Raca","Rapi")))%>%
    ggplot(aes(x=species, y=dep_var)) +
    geom_violin() +
    geom_point(data=con_exp_indiv, aes(x=species, y=exp),  cex=2, position=position_jitter(height=0, width=0.1)) +
    ylab(label=paste0(name_dep))
  
  all_p <- pos_exp_indiv %>%
    rename(infect=max_Bd_load) %>%
    dplyr::select(indivID, exp, p, infect) 
  
  # Make list of items
  output <- list()
  output[["Var"]] <- data.frame( dep=dep, name_dep=name_dep
                                 , transform_raw_func=transform_raw_func
                                 , intercept_present=intercept_present
                                 ,fit_distr=fit_distr)
  output[["gg_model_Distribution_of_all_values"]] <- gg_distr
  output[["gg_p"]] <- gg_p
  output[["gg_raw"]] <- gg_raw
  output[["gg_ExpectedDistribution_and_Bd_exposed"]] <- gg_exp_distr
  output[["gg_ExpectedDistribution_controls"]] <- gg_exp_distr_controls
  output[["all_p"]] <- all_p
  
  return(output)
}

