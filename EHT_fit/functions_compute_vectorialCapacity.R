# interventions_vectorial_capacity_wrapper
#' @title interventions_vectorial_capacity_wrapper
#'
#' @description  \code{interventions_vectorial_capacity_wrapper}
#' Wrapper around vecotrla capacity calculation functions, in order to choose whether uncertainty and/or insecticide decay are included
#'
#' @param results stanfit file containing stan results for unwashed nets
#' @param results_washed20 stanfit file containing stan results for washed nets
#' @param decay either "weibull", "linear", "hill" or "smooth-compact", is "weibull" by default
#' @param model_p model built by AnophelesModel package
#' @param insecticides id of the number corresponding to the insecticide you want to get the interventions for
#' @param names list of names corresponding to the insecticide numbers
#' @param L attrition half life
#' @param kappa attrition shape parameter (for Weibull in particular)
#' @param npoints number of interpolation points for AnophelesModel
#' @param washedDecay boolean indicating if insecticide decay should be included
#' @param uncertainty boolean indicating if uncertainty should be included
#'
#' @return a list of interventions
#'
#' @export

interventions_vectorial_capacity_wrapper <- function(results,results_washed20,
                                                     decay="weibull",
                                                     model_p,
                                                     insecticides,
                                                     names, L, kappa,
                                                     inbed_exposure,
                                                     npoints,washedDecay=TRUE, uncertainty=FALSE, cov=NULL, nsamples=100, aged_time=3 ){
  if(!washedDecay){
    results_washed20=NULL
  }

  cov=cov*inbed_exposure
  
  if(uncertainty){
    if(is.null(cov)){
      stop("please provide coverage value when uncertainty is included")
    }
    out=interventions_vectorial_capacity_uncertainty(results=results,results_washed20=results_washed20,
                                                     decay=decay,
                                                     model_p=model_p,
                                                     insecticides=insecticides,
                                                     names=names, L=L, kappa=kappa,
                                                     nsamples=nsamples, seed=1,
                                                     npoints=npoints, cov=cov, aged_time=aged_time)
  } else{
    out=interventions_vectorial_capacity(results=results,results_washed20=results_washed20,
                                         decay=decay,
                                         model_p=model_p,
                                         insecticides=insecticides,
                                         names=names, L=L, kappa=kappa,
                                         npoints=npoints, cov=cov, aged_time=aged_time)
  }

  return(out)

}

# extract_stan_Mean
#' @title extract_stan_Mean
#'
#' @description  \code{extract_stan_Mean}
#' Extract mean value from stan output
#'
#' @param res stanfit file containing stan results
#' @param insecticides id of the number corresponding to the insecticide you want to get the interventions for
#'
#' @return a dataframe with the mean posterior values
#'
#' @export
extract_stan_Mean <- function(res){


  mypars=c("InitialPostprandialkillingEfficacy",
           "InitialPreprandialkillingEfficacy",
           "InitialRepellencyRate",
           "KillingDuringHostSeeking",
           "alpha_0", "mu_0")

  fit_summary <- rstan::summary(res, pars=mypars)
  output <- as.data.frame(fit_summary$summary, optional = TRUE) %>% dplyr::select(mean)
  names(output)=c("mean")

  new_output <- output %>%
    tibble::rownames_to_column( "param") %>%
    mutate(param=ifelse(str_detect( param,"\\["), param, paste0(param, "[0]")))%>%
    tidyr::separate(param, sep="\\[" ,into =c( "param" ,"insecticide")) %>%
    mutate(insecticide=as.numeric(gsub("]", "", insecticide))-1 )
    

  new_output$param[new_output$param=="InitialRepellencyRate"]="InitialRepellentEfficacy"

  return(new_output)
}


# interventions_vectorial_capacity
#' @title interventions_vectorial_capacity
#'
#' @description  \code{interventions_vectorial_capacity}
#' The purpose of this function is get an intervention vector that will then be used to calculate the impact of said intervention.
#' You get this intervention vector from stan results. You have to provide which insecticide you want to get the intervention.
#'
#' @param results stanfit file containing stan results for unwashed nets
#' @param results_washed20 stanfit file containing stan results for washed nets
#' @param decay either "weibull", "linear", "hill" or "smooth-compact", is "weibull" by default
#' @param model_p model built by AnophelesModel package
#' @param insecticides id of the number corresponding to the insecticide you want to get the interventions for
#' @param names list of names corresponding to the insecticide number
#' @param L attrition half life
#' @param kappa attrition shape parameter (for Weibull in particular)
#' @param npoints number of interpolation points for AnophelesModel
#' @param cov effective coverage of the intervention (ITN use)
#'
#' @return a list of interventions
#'
#' @export


interventions_vectorial_capacity <- function(results,results_washed20,
                                             decay="weibull",
                                             model_p,
                                             insecticides,
                                             names, L=NULL, kappa=NULL,
                                             npoints, cov, aged_time){


  #new_output=extract_stan_Mean(res=results)
  new_output=extract_stan_posteriormax(res=results)%>% dplyr::rename(insecticide=treatment, mean=value)
  if(! is.null(results_washed20)){
    #new_output_washed20=extract_stan_Mean(res=results_washed20)
    new_output_washed20=extract_stan_posteriormax(res=results_washed20)%>% dplyr::rename(insecticide=treatment, mean=value)
    
  } else {
    new_output_washed20=NULL
  }


  myInterv <- get_interventionObj_custom(new_output= new_output,new_output_washed20= new_output_washed20,
                                         decay=decay, model_p=model_p, i=insecticides,L=L, kappa=kappa, name=names, npoints=npoints, aged_time=aged_time)

  myImpacts =  calculate_impact(interventions_vec =myInterv,
                                      coverage_vec  =  c(0,cov),
                                      model_p = model_p, Nv0 = 10000, num_ip_points = npoints)
  return(myImpacts$interventions_vec[[1]]$effects$avg_impact[2])
}


######### get_interventionObj_custom #####
#' @title get_interventionObj_custom
#'
#' @description  \code{get_interventionObj_custom}
#'This function is used inside the previous function to obtain the intervention object associated with a desired set of parameters
#'
#' @param new_output an object
#' @param decay either "weibull", "linear", "hill" or "smooth-compact", is "weibull" by default
#' @param model_p for Anopheles Package
#' @param i the number of the corresponding insecticide you want to get the interventions for
#' @param name name corresponding to the insecticide number
#'
#' @return one intervention
#'
get_interventionObj_custom <- function(new_output,new_output_washed20,
                                       decay,  model_p,i, name,npoints, L, kappa, aged_time=3){


  #our_insecticide <- new_output %>% filter(insecticide==i) %>% dplyr::select(-insecticide)
  new_output_filter <- new_output %>% filter(insecticide==i | insecticide==-1) %>% dplyr::select(-insecticide)
  myproba=getproba(stan_output=new_output_filter)
  
  if(! is.null(new_output_washed20)){
    new_output_filter_washed20 <- new_output_washed20 %>% filter(insecticide==i| insecticide==-1) %>% dplyr::select(-insecticide)
    myproba_washed20 <-getproba(stan_output=new_output_filter_washed20)
    
    halflife_insecticide <- aged_time/(1-myproba_washed20/myproba)
    names(halflife_insecticide)=c("Deterrency","PrePrandial","PostPrandial")

  } else {
    halflife_insecticide = NULL
  }

  interv =get_interv(myproba=myproba,
                     L=L, kappa=kappa, npoints = npoints, model_p = model_p, name=name, decay=decay, halflife_insecticide=halflife_insecticide)

  return(interv)
}



######### getproba #####
#' @title getproba
#'
#' @description  \code{getproba}
#' Get deterrence, preprandial and postprandial parameters
#'
#' @param stan_output an object
#'
#' @return updated vector of probabilities
#'
getproba <- function(stan_output){

  mypi=stan_output$mean[stan_output$param=="InitialRepellentEfficacy"]
  myphi=stan_output$mean[stan_output$param=="InitialPreprandialkillingEfficacy"]
  myxi=stan_output$mean[stan_output$param=="InitialPostprandialkillingEfficacy"]
  
  myproba_update=c(mypi,myphi,myxi)
  names(myproba_update)=c("Deterrency","PrePrandial","PostPrandial")
  
  return(myproba_update)
}
#########get_interv#####
#' @title get_interv
#'
#' @description  \code{get_interv}
#' The purpose of this function is to get your intervention for the Anopheles package given the pre prandial
#' post prandial and the two parameters of the weibull decay L and kappa
#'
#' @param myproba list of pre, post prandial and deterency
#' @param interventions_vec an intervention vector that is NULL by default
#' @param which_factor dunno
#' @param L parameter of the decay
#' @param kappa parameter of the decay
#' @param npoints equal to 0 at the beginning, as many as the days in the data
#' @param effect_only dunno
#' @param model_p dunno
#' @param name the nam of the new parameterization
#' @param intervention_type here is by default IRS and can only be IRS
#' @param decay either "weibull", "linear", "hill" or "smooth-compact", is "weibull" by default
#'
#'
#' @return one interventions
#'
#' @export


# function to modify intervention effects, for a Weibull decay
get_interv =function(myproba,L=NULL, kappa=NULL, npoints,
                     model_p=NULL, name=NULL, intervention_type="LLINs", decay="weibull",halflife_insecticide){
  
  if(!is.null(model_p)){
    list_interv =list(intervention_obj_examples$LLINs_example)
    interventions_vec = suppressMessages(def_interventions_effects(list_interv, model_p, num_ip_points =max(npoints,3), verbose=FALSE))
  } else{
    stop("model_p")
  }
  
  interv = interventions_vec[[1]]$effects
  
  # attrition decay
  decay_factor=decay_function(L=L,kappa=kappa,  datej=c(0:(npoints-1)),decay=decay)
  interv$survival=decay_factor
  
  # entomological effects + insecticide decay + activity rhythms
  if("Deterrency" %in%  names(myproba)){
    if(! is.null(halflife_insecticide)){
      #if(halflife_insecticide["Deterrency"]>0){
      insecticide_decay_factor_det=decay_function(L=365*halflife_insecticide["Deterrency"],kappa=NULL,  datej=c(0:(npoints-1)),decay="linear")
      #} else{
      #  insecticide_decay_factor_det=1
      #}
    } else {
      insecticide_decay_factor_det=1
    }
    interv$alphai[,1]=interv$alphai[,2]*(1-myproba["Deterrency"]*insecticide_decay_factor_det)
  } else{
    interv$alphai[,1]=interv$alphai[,2]
  }
  if("PrePrandial" %in%  names(myproba)){
    if(! is.null(halflife_insecticide)){
      #if(halflife_insecticide["PrePrandial"]>0){
      insecticide_decay_factor_pre=decay_function(L=365*halflife_insecticide["PrePrandial"],kappa=NULL,  datej=c(0:(npoints-1)),decay="linear")
      #} else{
      #  insecticide_decay_factor_pre=1
      #}
    } else {
      insecticide_decay_factor_pre=1
    }
    interv$PBi[,1]=interv$PBi[,2]*(1-myproba["PrePrandial"]*insecticide_decay_factor_pre)
  } else{
    interv$PBi[,1]=interv$PBi[,2]
  }
  if("PostPrandial" %in%  names(myproba)){
    if(! is.null(halflife_insecticide)){
      #if(halflife_insecticide["PostPrandial"]>0){
      insecticide_decay_factor_post=decay_function(L=365*halflife_insecticide["PostPrandial"],kappa=NULL,  datej=c(0:(npoints-1)),decay="linear")
      #} else {
      #  insecticide_decay_factor_post=1
      #}
    } else {
      insecticide_decay_factor_post=1
    }
    interv$PCi[,1]=interv$PCi[,2]*(1-myproba["PostPrandial"]*insecticide_decay_factor_post)
  } else{
    interv$PCi[,1]=interv$PCi[,2]
  }
  
  interv$alphai = interv$alphai[1:max(1, npoints), ]
  interv$PBi =  interv$PBi[1:max(1, npoints), ]
  interv$PCi = interv$PCi[1:max(1, npoints), ]
  interv$PDi = interv$PDi[1:max(1, npoints), ]
  interv$PEi = interv$PEi[1:max(1, npoints), ]
  
  
  result=interventions_vec
  result[[1]]$effects=interv
  if(!is.null(name)){
    result[[1]]$parameterisation=name
    result[[1]]$description=name
  }
  return(result)
}


###############################
# extract_stan_uncertainty
#' @title extract_stan_uncertainty
#'
#' @description  \code{extract_stan_uncertainty}
#' Extract a sample of posterior values from stan output
#'
#' @param res stanfit file containing your stan results
#' @param insecticides id of the number corresponding to the insecticide you want to get the interventions for
#' @param nsamples number of stan samples to retain
#' @param seed seed for random selection of samples
#'
#' @return a list of interventions
#'
#' @export
extract_stan_uncertainty <- function(res,nsamples, seed){


  mypars=c("InitialPostprandialkillingEfficacy",
           "InitialPreprandialkillingEfficacy",
           "InitialRepellencyRate",
           "KillingDuringHostSeeking",
           "alpha_0", "mu_0")

  set.seed(seed)

  fit_extract <- rstan::extract(res, pars=mypars, inc_warmup = FALSE)

  mysample=sample(1:nrow(fit_extract[[1]]), nsamples)

  full_df_extract=data.frame(nsample=1:nsamples)
  for(ii in 1:4){
    mydf=data.frame(fit_extract[[ii]][mysample,])
    names(mydf)=paste0(mypars[[ii]], "_", 1:ncol(mydf))
    full_df_extract=cbind(full_df_extract,mydf )
  }

  for(ii in 5:6){
    mydf=data.frame(fit_extract[[ii]][mysample])
    names(mydf)=paste0(mypars[[ii]], "_", 0)
    full_df_extract=cbind(full_df_extract,mydf )
  }

  new_output=full_df_extract %>%
    tidyr::pivot_longer(cols=names(full_df_extract)[-1], names_to = "param", values_to = "mean")%>%
    mutate(param=ifelse(str_detect( param,"_0_0"), gsub("_0_0", "_0", param), param))%>%
    tidyr::separate(param, sep="_" ,into =c( "param" ,"insecticide"))%>%
    mutate(insecticide=as.numeric(insecticide)-1)
  new_output$param[new_output$param=="InitialRepellencyRate"]="InitialRepellentEfficacy"
  new_output$param[new_output$param=="alpha"]="alpha_0"
  new_output$param[new_output$param=="mu"]="mu_0"

  return(new_output)
}



###############################

# interventions_vectorial_capacity with uncertainty
#' @title interventions_vectorial_capacity_uncertainty
#'
#' @description  \code{interventions_vectorial_capacity_uncertainty}
#' The purpose of this function is get an intervention vector that will then be used to calculate the impact of said intervention.
#' You get this intervention vector from stan results. You have to provide which insecticides you want to get the intervention.
#'
#' @param results stanfit file containing your stan results
#' @param decay either "weibull", "linear", "hill" or "smooth-compact", is "weibull" by default
#' @param model_p for Anopheles Package
#' @param insecticides id of the number corresponding to the insecticide you want to get the interventions for
#' @param names list of names corresponding to the insecticide number
#' @param nsamples number of stan samples to retain
#' @param seed seed for random selection of samples
#'
#' @return a list of interventions
#'
#' @export
interventions_vectorial_capacity_uncertainty=function(results,results_washed20,
                                                      decay="weibull",
                                                      model_p,
                                                      insecticides,
                                                      names, L=NULL, kappa=NULL, nsamples=100, seed=1, npoints, cov,aged_time){


  new_output=extract_stan_uncertainty(results,nsamples, seed=seed)

  if(! is.null(results_washed20)){
    new_output_washed20=extract_stan_uncertainty(results_washed20,nsamples, seed=seed)

  } else {
    new_output_washed20=NULL
  }

  VCred = c()

  for(this_sample in 1:nsamples){

      this.new_output=new_output %>% dplyr::filter(nsample ==this_sample) %>% dplyr::select(-nsample)
      if(! is.null(new_output_washed20)){
        this.new_output_washed20 = new_output_washed20 %>% dplyr::filter(nsample ==this_sample) %>% dplyr::select(-nsample)
      } else {
        this.new_output_washed20 = NULL
      }
      myList <- get_interventionObj_custom(
        new_output= this.new_output,
        new_output_washed20=this.new_output_washed20,
        decay=decay,model_p=model_p,i=insecticides,name=paste0(names, "_", this_sample),
        L=L, kappa=kappa,npoints = npoints, aged_time=aged_time)

    this.VCred=   calculate_impact(interventions_vec =myList,
                                   coverage_vec  =  c(0, cov),
                                   model_p =model_p, Nv0 = 10000, num_ip_points = npoints)

    VCred=c(VCred, this.VCred$interventions_vec[[1]]$effects$avg_impact[2])
  }


  return(VCred)
}



#########decay#####
#' @title decay functiom
#'
#' @description  \code{decay function}
#'
#' @param L half life of the decay
#' @param kappa shape parameter of the decay (for weibull, hill and smooth-compact decays only)
#' @param datej vector of time points (in days) for which the decay function needs to be evaluated
#' @param decay can either be "weibull", "hill", "linear" or "smooth-compact"
#'
#' @return summary data
#'

decay_function <- function(L,
                           kappa=NULL,
                           datej,
                           decay){

  if (decay %in% c("weibull" , "Weibull")){
    return( exp( -(datej/(L))^kappa * log(2) ))
  }

  if (decay %in% c("linear","Linear")){
    return( 1-datej/L)
  }

  if (decay %in% c("hill","Hill")){

    return (1/(1+(datej/L)^kappa))
  }

  if (decay %in% c("smooth-compact" , "smoothcompact")){

    return (exp(kappa - kappa/(1-(datej/L)^2)))
  }
  if (decay %in% c("none" , "no decay")){

    return (1)
  }

}

