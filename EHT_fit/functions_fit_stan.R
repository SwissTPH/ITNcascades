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
         "mu_0",
         "pc_0")
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



get_stan_summary_output=function(results,
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


###############################
# extract_stan_posteriormax
#' @title extract_stan_posteriormax
#'
#' @description  \code{extract_stan_uncertainty}
#' Extract a sample of posterior values from stan output
#'
#' @param res stanfit file containing your stan results
#'
#' @return a list of interventions
#'
#' @export
extract_stan_posteriormax <- function(res, data=NULL){
  
  
  mypars=c("InitialPostprandialkillingEfficacy",
           "InitialPreprandialkillingEfficacy",
           "InitialRepellencyRate",
           "KillingDuringHostSeeking",
           "alpha_0", "mu_0", "lp__")
  
  fit_extract <- rstan::extract(res, pars=mypars, inc_warmup = FALSE)
  
  id_MAP=which.max(fit_extract$lp__)
  full_df_extract=data.frame(sample=id_MAP)
  for(ii in 1:4){
    mydf=as.data.frame(t(fit_extract[[ii]][id_MAP,]))
    names(mydf)=paste0(mypars[[ii]], "_", 1:ncol(mydf))
    full_df_extract=cbind(full_df_extract,mydf )
  }
  
  for(ii in 5:6){
    mydf=as.data.frame(t(fit_extract[[ii]][id_MAP]))
    names(mydf)=paste0(mypars[[ii]], "_", 0)
    full_df_extract=cbind(full_df_extract,mydf )
  }
  
  new_output=full_df_extract %>%
    tidyr::pivot_longer(cols=names(full_df_extract)[-1], names_to = "param", values_to = "value")%>%
    mutate(param=ifelse(str_detect( param,"_0_0"), gsub("_0_0", "_0", param), param))%>%
    tidyr::separate(param, sep="_" ,into =c( "param" ,"treatment"))%>%
    mutate(treatment=as.numeric(treatment)-1)%>%
    filter(treatment!=0)
  
  if(!is.null(data)){
    new_output=new_output  %>%
    dplyr::left_join(unique(data %>% dplyr::select(insecticide_name, treatment)))
  }
  new_output$param[new_output$param=="InitialRepellencyRate"]="InitialRepellentEfficacy"
  new_output$param[new_output$param=="alpha"]="alpha_0"
  new_output$param[new_output$param=="mu"]="mu_0"
  
  return(new_output)
}



