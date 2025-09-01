return_anophelesParam=function(mosqs, value_effect, param, mosqs_exposures){
  anophelesParams = list()
  for (m in mosqs){
    indoor=str_detect(m, "indoor")
    exposure_m=as.numeric(mosqs_exposures[m])

    if(! is.null(value_effect)){
      value_effect_simple=value_effect$Unwashed
      names(value_effect_simple)=value_effect$parameter
      #value_effect_corrected=round(correct_inbed_exposure_gvi(myproba=value_effect_simple, inbed = exposure_m), digits = 2)
      value_effect_corrected=round(getvalue_gvi(myproba=value_effect_simple), digits = 2)
      myvalue=as.character(as.numeric(value_effect_corrected[param]))
    } else {
      myvalue=0
    }

    this.list=list(
      m_in = list(
        propActive = ifelse(indoor, 1, 0),
        value = ifelse(indoor, myvalue, 0)
      ))
    names(this.list)=m
    anophelesParams=append(anophelesParams, this.list)
  }
  return(anophelesParams)
}

#' #' Calculate phi parameter (preprandial killing effect)
#' #'
#' #' @inheritParams P_UA_new
#' #' @return Updated phi value
#' calculate_phi=function(alpha0, mu0, mypi, mykappa){
#'   mu_i=mu0+alpha0*mykappa
#'   alpha_i=alpha0*(1-mypi)
#' 
#'   p_UD_0=P_UD(alpha0,mu0)
#'   p_UD_i=P_UD(alpha_i,mu_i)
#' 
#'   phi=1-(1-p_UD_i)/(1-p_UD_0)
#'   return(phi)
#' }
#' 
#' #' Probability of remaining Unfed Alive (UA)
#' #'
#' #' @param alpha Transition rate to fed
#' #' @param mu Transition rate to dead
#' #' @param tt Time
#' #' @return Probability of staying UA at time tt
#' P_UA=function(alpha, mu, tt=1){
#'   return(
#'     exp(-(alpha+mu)*tt)
#'   )
#' }
#' 
#' 
#' #' Probability of becoming Unfed Dead (UD)
#' #'
#' #' @inheritParams P_UA
#' #' @return Probability of being UD at time tt
#' P_UD=function(alpha, mu, tt=1){
#' 
#' 
#'   return(
#'     (1-P_UA(alpha,mu, tt))*mu/(alpha+mu)
#'   )
#' }
#' 
#' 
#' ######### correct_inbed_exposure #####
#' #' @title correct_inbed_exposure
#' #'
#' #' @description  \code{correct_inbed_exposure}
#' #' Correct deterrence, preprandial and postprandial parameters for in-bed exposure
#' #'
#' #' @param stan_output an object
#' #' @param inbed exposure coefficient, from AnophelesModel
#' #'
#' #' @return updated vector of probabilities
#' #'
#' correct_inbed_exposure_gvi <- function(myproba, inbed){
#' 
#'   deterrency_update=myproba["Repellent"]*inbed
#'   kappa_update=myproba["KillingDuringHostSeeking"]*inbed
#'   postprandial_update=myproba["Postprandialkilling"]*inbed*(1-myproba["Repellent"])/(1-myproba["Repellent"]*inbed)
#'   preprandial_update=calculate_phi(alpha0=myproba["alpha_0"],
#'                                    mu0=myproba["mu_0"],
#'                                    deterrency_update, kappa_update)
#' 
#'   myproba_update=c(deterrency_update,preprandial_update,postprandial_update)
#'   names(myproba_update)=c("Deterrence","PrePrandial","PostPrandial")
#' 
#'   return(myproba_update)
#' }

######### correct_inbed_exposure #####
#' @title correct_inbed_exposure
#'
#' @description  \code{correct_inbed_exposure}
#' Correct deterrence, preprandial and postprandial parameters for in-bed exposure
#'
#' @param stan_output an object
#' @param inbed exposure coefficient, from AnophelesModel
#'
#' @return updated vector of probabilities
#'
getvalue_gvi <- function(myproba){
  
  deterrency_update=myproba["Repellent"]
  preprandial_update=myproba["Preprandialkilling"]
  postprandial_update=myproba["Postprandialkilling"]
  myproba_update=c(deterrency_update,preprandial_update,postprandial_update)
  names(myproba_update)=c("Deterrence","PrePrandial","PostPrandial")
  
  return(myproba_update)
}

#' Extract GVi values for ITNs from an external file containing fits from EHT data. Overall wrapper
#'
#' @param EHT name of the EHT from which fitted parameters are extracted
#' @param halflife_functionalSurvival Half life of the net representing functional survival. Assumes a Weibull decay
#' @param kappa_functionalSurvival Shape parameter for the decay in functional survival. Assumes a Weibull decay
#' @param exposure_correction 
#' @param myname a name to give to the snippet
#' @param mosqs a mosq object from OpenMalariaUtilities
#' @param parameters_GVI the dataframe containing parameter values
#' @return the input for the function called defineGVI_simple
#' @export
extract_GVI_params=function(EHT, netType, halflife_functionalSurvival, kappa_functionalSurvival, mosqs_exposures, myname, mosqs, parameters_GVI, insecticide_decay=T, decay="weibull"){
  my_GVI_params=update_halflife_insecticideDecay(my_parameters_GVI=parameters_GVI,
                                                 halflife=halflife_functionalSurvival, kappa=kappa_functionalSurvival,
                                                 my_EHT =EHT, my_netType=netType, insecticide_decay=insecticide_decay)
  
  
  if(insecticide_decay){
    deterrency_snippet=create_vectorInterventionParameters(my_GVI_params, param="deterrence",decay= decay, myname=paste0(myname, "_deterrency"), mosqs=mosqs, mosqs_exposures=mosqs_exposures)
    preprandial_snippet=create_vectorInterventionParameters(my_GVI_params, param="preprandial",decay= decay, myname=paste0(myname, "_preprandial"), mosqs=mosqs, mosqs_exposures=mosqs_exposures)
    postprandial_snippet=create_vectorInterventionParameters(my_GVI_params, param="postprandial",decay= decay, myname=paste0(myname, "_postprandial"), mosqs=mosqs, mosqs_exposures=mosqs_exposures)
  
    output=list(deterrency_snippet, preprandial_snippet,postprandial_snippet)
  } else {
    gvi_snippet=create_vectorInterventionParameters(my_GVI_params, param="all",L=halflife_functionalSurvival,kappa=kappa_functionalSurvival,
                                                    decay= decay, myname=myname, mosqs=mosqs, mosqs_exposures=mosqs_exposures)
    output=gvi_snippet
  }

  return(output)
}


#' Create GVI snippet for the point estimate, the lower and the upper bound of the credible intervals
#'
#' @return a list of inputs for the function called defineGVI_simple
#' @export
create_vectorInterventionParameters = function(my_GVI_params, param="all",#deterrency, preprandial, postprandial,
                                               L=NULL, kappa=NULL, decay= "weibull", myname, mosqs, mosqs_exposures){
  
  
  my_GVI_params_det=my_GVI_params
  my_GVI_params_pre=my_GVI_params
  my_GVI_params_post=my_GVI_params
  
  if(param=="deterrence"){
    my_GVI_params_pre=NULL
    my_GVI_params_post=NULL
    L=my_GVI_params$final_hl[my_GVI_params$parameter=="Repellent"]
    kappa=my_GVI_params$final_kappa[my_GVI_params$parameter=="Repellent"]
  }
  if(param=="preprandial"){
    my_GVI_params_det=NULL
    my_GVI_params_post=NULL
    L=my_GVI_params$final_hl[my_GVI_params$parameter=="Preprandialkilling"]
    kappa=my_GVI_params$final_kappa[my_GVI_params$parameter=="Preprandialkilling"]
    
  }
  if(param=="postprandial"){
    my_GVI_params_det=NULL
    my_GVI_params_pre=NULL
    L=my_GVI_params$final_hl[my_GVI_params$parameter=="Postprandialkilling"]
    kappa=my_GVI_params$final_kappa[my_GVI_params$parameter=="Postprandialkilling"]
    
  }
  
  anophParam_deterrency=return_anophelesParam(mosqs, value_effect = my_GVI_params_det, param="Deterrence", mosqs_exposures=mosqs_exposures)
  anophParam_preprandial=return_anophelesParam(mosqs, value_effect = my_GVI_params_pre, param="PrePrandial", mosqs_exposures=mosqs_exposures)
  anophParam_postprandial=return_anophelesParam(mosqs, value_effect = my_GVI_params_post, param="PostPrandial", mosqs_exposures=mosqs_exposures)
  
  myvectorInterventionParameters <- list(
    myname = list(
      deterrency = list(
        decay = list(
          L = as.character(L), k=as.character(kappa),
          "function" = decay
        ),
        anophelesParams = anophParam_deterrency
      ),
      preprandialKillingEffect = list(
        decay = list(
          L = as.character(L), k=as.character(kappa),
          "function" = decay
        ),
        anophelesParams = anophParam_preprandial
      ),
      postprandialKillingEffect = list(
        decay = list(
          L = as.character(L), k=as.character(kappa),
          "function" = decay
        ),
        anophelesParams = anophParam_postprandial
      )
    )
  )
  names(myvectorInterventionParameters)=myname
  return(myvectorInterventionParameters)
}

#' #' Create GVI snippet for the point estimate, the lower and the upper bound of the credible intervals
#' #'
#' #' @return a list of inputs for the function called defineGVI_simple
#' #' @export
#' create_vectorInterventionParameters_old = function(deterrency, preprandial, postprandial,
#'                                                 L, kappa, decay= "weibull", myname, mosqs, mosqs_exposures){
#' 
#'   anophParam_deterrency=return_anophelesParam(mosqs, value_effect = deterrency, mosqs_exposures=mosqs_exposures)
#'   anophParam_preprandial=return_anophelesParam(mosqs, value_effect = preprandial, mosqs_exposures=mosqs_exposures)
#'   anophParam_postprandial=return_anophelesParam(mosqs, value_effect = postprandial, mosqs_exposures=mosqs_exposures)
#' 
#'   myvectorInterventionParameters <- list(
#'     myname = list(
#'       deterrency = list(
#'         decay = list(
#'           L = as.character(L), k=as.character(kappa),
#'           "function" = decay
#'         ),
#'         anophelesParams = anophParam_deterrency
#'       ),
#'       preprandialKillingEffect = list(
#'         decay = list(
#'           L = as.character(L), k=as.character(kappa),
#'           "function" = decay
#'         ),
#'         anophelesParams = anophParam_preprandial
#'       ),
#'       postprandialKillingEffect = list(
#'         decay = list(
#'           L = as.character(L), k=as.character(kappa),
#'           "function" = decay
#'         ),
#'         anophelesParams = anophParam_postprandial
#'       )
#'     )
#'   )
#'   names(myvectorInterventionParameters)=myname
#'   return(myvectorInterventionParameters)
#' }

########
# Add the GVI parameters to the baseList for OpenMalariaUtilities


defineGVI_simple=function (baseList, vectorInterventionParameters, append = TRUE,
                           name = NULL, verbatim = TRUE, hist = FALSE)
{
  assertCol <- checkmate::makeAssertCollection()
  checkmate::assert(checkmate::checkLogical(hist), add = assertCol)
  checkmate::reportAssertions(assertCol)
  if (is.null(baseList$interventions$human)) {
    stop("To append, the baseList needs a child called '$interventions$human'")
  }
  for (intervention in names(vectorInterventionParameters)) {
    for (effect in names(vectorInterventionParameters[[intervention]])) {
      if (!setequal(names(vectorInterventionParameters[[intervention]][[effect]][["anophelesParams"]]),
                    unique(unlist(lapply(baseList$entomology$vector,
                                         function(x) x$mosquito))))) {
        stop("To append, each vector species defined in the entomology section must be the same as in the intervention component.")
      }
    }
  }
  baseList <- openMalariaUtilities:::.defineInterventionsHeader(baseList = baseList)
  for (k in names(vectorInterventionParameters)) {
    componentData <- vectorInterventionParameters[[k]]
    component_id=k
    if (verbatim) {
      message(paste0("Defining intervention with component_id: ",
                     component_id))
    }
    GVIList <- list(decay = if (hist) {
      list(L = 1, `function` = "step")
    } else {
      componentData[[effect]][["decay"]]
    })
    for (vector_species in names(componentData[[effect]]$anophelesParams)) {
      print(paste0("Writing effect values for vector species: ",
                   vector_species))
      values <- c(deterrency = 0, preprandialKillingEffect = 0,
                  postprandialKillingEffect = 0)
      for (effect in names(componentData)) {
        values[effect] <- componentData[[effect]][["anophelesParams"]][[vector_species]][["value"]]
      }
      GVIList <- append(GVIList, list(anophelesParams = list(mosquito = vector_species,
                                                             propActive = componentData[[effect]][["anophelesParams"]][[vector_species]][["propActive"]],
                                                             deterrency = list(value = values[["deterrency"]]),
                                                             preprandialKillingEffect = list(value = values[["preprandialKillingEffect"]]),
                                                             postprandialKillingEffect = list(value = values[["postprandialKillingEffect"]]))))
    }
    baseList <- openMalariaUtilities:::.xmlAddList(data = baseList, sublist = c("interventions",
                                                                                "human"), append = append, entry = "component",
                                                   input = list(id = component_id, name = if (is.null(name)) "your_tag" else name[[k]],
                                                                GVI = GVIList))
    #}
  }
  return(baseList)
}


###########################
# Update decay parameters to include insecticide decay


# function to refit weibull decay
optimise_decay_param=function(param, insecticide, df){
  #if(param[3]<1){
  out=data.frame(
    time=df %>% filter(setting==insecticide) %>% pull(time),
    decay_obs=df %>% filter(setting==insecticide) %>% pull(ITNcov)
  )%>%
    mutate(decay=param[3]*exp( -(time/(param[1]))^param[2] * log(2) ),
           diff=(decay-decay_obs)^2)%>%
    summarise(diff=sum(diff))
  return(out)
}

# function to refit hill decay
optimise_decay_param_hill=function(param, insecticide){
  #if(param[3]<1){
  data.frame(
    time=ITN_cov %>% filter(setting==insecticide) %>% pull(time),
    decay_obs=ITN_cov %>% filter(setting==insecticide) %>% pull(ITNcov)
  )%>%
    mutate(decay=param[3] / (1 + (time/param[1])^param[2]),
           diff=(decay-decay_obs)^2)%>%
    summarise(diff=sum(diff))
  return(out)
}

# optimisation routine to find decay parameters
calculate_param=function(insecticide, opti_fun="weibull", df){
  opti_function=ifelse(opti_fun=="weibull", optimise_decay_param, optimise_decay_param_hill)
  opti=optim(c(0.448, 1.11, 0.7), opti_function, insecticide=insecticide, df=df , method = "L-BFGS-B", lower=c(0, 0, 0), upper = c(Inf, Inf, 1) )$par
  names(opti)=c("L", "k", "a")
  opti$setting=insecticide
  opti$decay=opti_fun
  return(opti)
}


# Update halflife of the net to include insecticide decay from EHT data
update_halflife_insecticideDecay=function(my_parameters_GVI, halflife, kappa, my_EHT, my_netType, insecticide_decay){
  this.parameters_GVI=my_parameters_GVI%>%
    filter(netType %in% c(my_netType, "control"), EHT==my_EHT)%>%
    mutate(L_functionalSurvival=halflife,
           kappa_functionalSurvival=kappa)
  
  if(insecticide_decay){
    
    df=data.frame(time=seq(0, 365*3))%>%
      mutate( weibull=exp( -(time/halflife)^kappa * log(2) ))
    
    for(i in 1:nrow(this.parameters_GVI)){
      this.line=this.parameters_GVI[i,]
      
      if(is.na(this.line$halflife_insecticide) | this.line$halflife_insecticide<0 ){
        this.parameters_GVI$final_hl[i]=halflife
        this.parameters_GVI$final_kappa[i]=kappa
      } else {
        this.df=df %>% mutate(
          linear= 1-time/this.line$halflife_insecticide,
          ITNcov=weibull*linear, setting="test")
        
        product_fit=calculate_param("test", opti_fun="weibull", df=this.df)
        
        this.parameters_GVI$final_hl[i]=round(product_fit$L, digits = 1)
        this.parameters_GVI$final_kappa[i]=round(product_fit$k, digits = 1)
        this.df=NULL
      }
    }
  } else{
    this.parameters_GVI=this.parameters_GVI%>%
      mutate(final_hl=halflife,
             final_kappa=kappa)
  }

  return(this.parameters_GVI)
}

