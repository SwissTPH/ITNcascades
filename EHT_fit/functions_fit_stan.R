fit_EHT_multinomial_stan=function(data,iter = 6000,
                           warmup= 3000,
                           chains= 4,
                           path, stanpath, stanmodel= "EHT_fitting_model.stan"){
  
  ### building our data for stan
  nb_treat <- dplyr::n_distinct(data$treatment)
  nb_days  <- dplyr::n_distinct(data$day)
  
  
  
  # data for the random effect model
  stan_data<-list(N=dim(data)[1],
                  tr=nb_treat-1,
                  treat=data$treatment,
                  y=matrix(c(data$UA, data$UD,data$FA ,data$FD ), ncol=4)
                  )
  
  
  
  stan_file = file.path(stanpath,stanmodel)
  pars=c("InitialPostprandialkillingEfficacy",
         "KillingDuringHostSeeking",
         "InitialRepellencyRate",
         "InitialPreprandialkillingEfficacy",
         "alpha_0",
         "mu_0")
  gen_inits <- function() {
    list=list(
      InitialPostprandialkillingEfficacy = runif(nb_treat,0,1),
      KillingDuringHostSeeking =runif(nb_treat,0,1),
      InitialRepellencyRate = runif(nb_treat,0,1)
    )
  }
  
  stan_output <-rstan::stan(file=stan_file,
                            data=stan_data,
                            pars=pars,
                            include = TRUE, #only samples for parameters in pars are stored in the fitted results
                            iter=iter,
                            init=gen_inits,
                            warmup = warmup, #burnin
                            chains=chains)
  
  
  saveRDS(stan_output , path)
  
  return(stan_output)
  
}



get_IRSmodel_summary_output=function(results,
                                     decay="weibull",
                                     save=FALSE,
                                     path=NULL, data) {
  
  mypars=c("InitialPostprandialkillingEfficacy",
           "InitialPreprandialkillingEfficacy",
           "InitialRepellencyRate",
           "KillingDuringHostSeeking")
  
  if ("kappa" %in% results@model_pars){
    mypars=c(mypars, "kappa")
  }
  
  if("L" %in% results@model_pars) {
    mypars=c(mypars, "L")
  }
  
  fit_summary <- rstan::summary(results, pars=mypars)
  summary_outputs=fit_summary$summary %>%
    as.data.frame()%>%
    round(digits = 2) %>%
    dplyr::select(c("mean", "2.5%", "97.5%"))%>%
    dplyr::mutate(index=rownames(fit_summary$summary))%>%
    tidyr::separate(index, sep = "\\[", into=c("param", "treatment"))%>%
    dplyr::mutate(treatment = as.numeric(gsub("\\]", "", treatment))-1) %>%
    dplyr::filter(treatment>0)%>%
    dplyr::left_join(unique(data %>% dplyr::select(insecticide_name, treatment)))%>%
    dplyr::select(insecticide_name, param, mean, "2.5%", "97.5%")
  
  
  
  if (save){
    gt(summary_outputs) %>% gtsave(filename = "summary_outputs.rtf", path=path)
  }
  gt::gt(summary_outputs)
  
  return(summary_outputs)
}