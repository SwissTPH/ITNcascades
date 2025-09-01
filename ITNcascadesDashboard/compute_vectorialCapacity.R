# interventions_vectorial_capacity
#' @title interventions_vectorial_capacity
#'
#' @description  \code{interventions_vectorial_capacity}
#' The purpose of this function is get an intervention vector that will then be used to calculate the impact of said intervention.
#' You get this intervention vector from stan results. You have to provide which insecticide you want to get the intervention.
#'
#' @param outputs_posteriormax file containing stan results (posterior maximum)
#' @param insecticide_decay boolean indicating if insecticide decay has to be included
#' @param decay either "weibull", "linear", "hill" or "smooth-compact", is "weibull" by default
#' @param model_p model built by AnophelesModel package
#' @param names list of names corresponding to the insecticide number
#' @param L attrition half life
#' @param kappa attrition shape parameter (for Weibull in particular)
#' @param cov effective coverage of the intervention (ITN use)
#'
#' @return a list of interventions
#'
#' @export


interventions_vectorial_capacity <- function(outputs_posteriormax,
                                             insecticide_decay,
                                             decay="weibull",
                                             model_p,
                                              L=NULL, kappa=NULL,
                                             cov){
  
  npoints=3*365
  
  myInterv <- get_interventionObj_custom(outputs_posteriormax= outputs_posteriormax,insecticide_decay= insecticide_decay,
                                         decay=decay, model_p=model_p, L=L, kappa=kappa,  npoints=npoints)
  
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
get_interventionObj_custom <- function(outputs_posteriormax,insecticide_decay,
                                       decay,  model_p,npoints, L, kappa){
  
  myproba=c(outputs_posteriormax$Unwashed[outputs_posteriormax$parameter=="Repellent"],
                   outputs_posteriormax$Unwashed[outputs_posteriormax$parameter=="Preprandialkilling"],
                   outputs_posteriormax$Unwashed[outputs_posteriormax$parameter=="Postprandialkilling"])
  
  names(myproba)=c("Deterrency","PrePrandial","PostPrandial")

  if(insecticide_decay){
    halflife_insecticide=c(outputs_posteriormax$halflife_insecticide[outputs_posteriormax$parameter=="Repellent"],
                     outputs_posteriormax$halflife_insecticide[outputs_posteriormax$parameter=="Preprandialkilling"],
                     outputs_posteriormax$halflife_insecticide[outputs_posteriormax$parameter=="Postprandialkilling"])
    
    names(halflife_insecticide)=c("Deterrency","PrePrandial","PostPrandial")
    
  } else {
    halflife_insecticide = NULL
  }
  
  interv =get_interv(myproba=myproba,
                     L=L, kappa=kappa, npoints = npoints, model_p = model_p, name=NULL, decay=decay, halflife_insecticide=halflife_insecticide)
  
  return(interv)
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
