return_anophelesParam=function(mosqs, value_effect){
  anophelesParams = list()
  for (m in mosqs){
    indoor=str_detect(m, "indoor")
    this.list=list(
      m_in = list(
        propActive = ifelse(indoor, 1, 0),
        value = ifelse(indoor, as.character(value_effect), 0)
      ))
    names(this.list)=m
    anophelesParams=append(anophelesParams, this.list)
  }
  return(anophelesParams)
}



extract_GVI_params=function(EHT, netType, halflife_functionalSurvival, kappa_functionalSurvival, exposure_correction, myname, mosqs, parameters_GVI){
  my_GVI_params=update_halflife_insecticideDecay(my_parameters_GVI=parameters_GVI,
                                                 halflife=halflife_functionalSurvival, kappa=kappa_functionalSurvival,
                                                 my_EHT =EHT, my_netType=netType)
  
  
  create_vectorInterventionParameters_3decays(deterrency = my_GVI_params$mean_Unwashed[my_GVI_params$parameter=="Reduction in host availability"],
                                              preprandial = my_GVI_params$mean_Unwashed[my_GVI_params$parameter=="Pre-prandial killing effect"],
                                              postprandial = my_GVI_params$mean_Unwashed[my_GVI_params$parameter=="Post-prandial killing effect"],
                                              deterrency_inf = my_GVI_params$q025_Unwashed[my_GVI_params$parameter=="Reduction in host availability"],
                                              preprandial_inf = my_GVI_params$q025_Unwashed[my_GVI_params$parameter=="Pre-prandial killing effect"],
                                              postprandial_inf =  my_GVI_params$q025_Unwashed[my_GVI_params$parameter=="Post-prandial killing effect"],
                                              deterrency_sup = my_GVI_params$q975_Unwashed[my_GVI_params$parameter=="Reduction in host availability"],
                                              preprandial_sup = my_GVI_params$q975_Unwashed[my_GVI_params$parameter=="Pre-prandial killing effect"],
                                              postprandial_sup =  my_GVI_params$q975_Unwashed[my_GVI_params$parameter=="Post-prandial killing effect"],
                                              
                                              L_deterrency=my_GVI_params$final_hl[my_GVI_params$parameter=="Reduction in host availability"],
                                              kappa_deterrency=my_GVI_params$final_kappa[my_GVI_params$parameter=="Reduction in host availability"],
                                              L_preprandial=my_GVI_params$final_hl[my_GVI_params$parameter=="Pre-prandial killing effect"],
                                              kappa_preprandial=my_GVI_params$final_kappa[my_GVI_params$parameter=="Pre-prandial killing effect"],
                                              L_postprandial=my_GVI_params$final_hl[my_GVI_params$parameter=="Post-prandial killing effect"],
                                              kappa_postprandial=my_GVI_params$final_kappa[my_GVI_params$parameter=="Post-prandial killing effect"],
                                              decay= "weibull",exposure_correction=exposure_correction,  myname= myname, mosqs = mosqs)
}




create_vectorInterventionParameters_3decays = function(deterrency, preprandial, postprandial,
                                                         deterrency_inf, preprandial_inf, postprandial_inf,
                                                         deterrency_sup, preprandial_sup, postprandial_sup,
                                                         L_deterrency, kappa_deterrency,
                                                         L_preprandial, kappa_preprandial,
                                                         L_postprandial, kappa_postprandial, decay= "weibull",exposure_correction,  myname, mosqs){

  deterrency_snippet=create_vectorInterventionParameters(deterrency=deterrency*exposure_correction, preprandial=0, postprandial=0,
                                                          L=L_deterrency, kappa=kappa_deterrency, decay= decay, myname=paste0(myname, "_deterrency"), mosqs=mosqs)
  preprandial_snippet=create_vectorInterventionParameters(deterrency=0, preprandial=preprandial*exposure_correction, postprandial=0,
                                                           L=L_preprandial, kappa=kappa_preprandial, decay= decay, myname=paste0(myname, "_preprandial"), mosqs=mosqs)
  postprandial_snippet=create_vectorInterventionParameters(deterrency=0, preprandial=0, postprandial=postprandial*exposure_correction,
                                                            L=L_postprandial, kappa=kappa_postprandial, decay= decay, myname=paste0(myname, "_postprandial"), mosqs=mosqs)

  deterrency_snippet_inf=create_vectorInterventionParameters(deterrency=deterrency_inf*exposure_correction, preprandial=0, postprandial=0,
                                                              L=L_deterrency, kappa=kappa_deterrency, decay= decay, myname=paste0(myname, "_inf_deterrency"), mosqs=mosqs)
  preprandial_snippet_inf=create_vectorInterventionParameters(deterrency=0, preprandial=preprandial_inf*exposure_correction, postprandial=0,
                                                               L=L_preprandial, kappa=kappa_preprandial, decay= decay, myname=paste0(myname, "_inf_preprandial"), mosqs=mosqs)
  postprandial_snippet_inf=create_vectorInterventionParameters(deterrency=0, preprandial=0, postprandial=postprandial_inf*exposure_correction,
                                                                L=L_postprandial, kappa=kappa_postprandial, decay= decay, myname=paste0(myname, "_inf_postprandial"), mosqs=mosqs)


  deterrency_snippet_sup=create_vectorInterventionParameters(deterrency=deterrency_sup*exposure_correction, preprandial=0, postprandial=0,
                                                              L=L_deterrency, kappa=kappa_deterrency, decay= decay, myname=paste0(myname, "_sup_deterrency"), mosqs=mosqs)
  preprandial_snippet_sup=create_vectorInterventionParameters(deterrency=0, preprandial=preprandial_sup*exposure_correction, postprandial=0,
                                                               L=L_preprandial, kappa=kappa_preprandial, decay= decay, myname=paste0(myname, "_sup_preprandial"), mosqs=mosqs)
  postprandial_snippet_sup=create_vectorInterventionParameters(deterrency=0, preprandial=0, postprandial=postprandial_sup*exposure_correction,
                                                                L=L_postprandial, kappa=kappa_postprandial, decay= decay, myname=paste0(myname, "_sup_postprandial"), mosqs=mosqs)


  return(list(deterrency_snippet, preprandial_snippet,postprandial_snippet,
              deterrency_snippet_inf, preprandial_snippet_inf,postprandial_snippet_inf,
              deterrency_snippet_sup, preprandial_snippet_sup,postprandial_snippet_sup))
}




create_vectorInterventionParameters = function(deterrency, preprandial, postprandial,
                                                L, kappa, decay= "weibull", myname, mosqs){

  anophParam_deterrency=return_anophelesParam(mosqs, value_effect = deterrency)
  anophParam_preprandial=return_anophelesParam(mosqs, value_effect = preprandial)
  anophParam_postprandial=return_anophelesParam(mosqs, value_effect = postprandial)

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

optimise_decay_param=function(param, insecticide, df){
  #if(param[3]<1){
  out=data.frame(
    time=df %>% filter(setting==insecticide) %>% pull(time),
    decay_obs=df %>% filter(setting==insecticide) %>% pull(ITNcov)
  )%>%
    mutate(decay=param[3]*exp( -(time/(param[1]))^param[2] * log(2) ),
           diff=(decay-decay_obs)^2)%>%
    summarise(diff=sum(diff))
  # } else {
  #   out=Inf
  # }
  return(out)
}
optimise_decay_param_hill=function(param, insecticide){
  #if(param[3]<1){
  data.frame(
    time=ITN_cov %>% filter(setting==insecticide) %>% pull(time),
    decay_obs=ITN_cov %>% filter(setting==insecticide) %>% pull(ITNcov)
  )%>%
    mutate(decay=param[3] / (1 + (time/param[1])^param[2]),
           diff=(decay-decay_obs)^2)%>%
    summarise(diff=sum(diff))
  # } else {
  #   out=Inf
  # }
}

calculate_param=function(insecticide, opti_fun="weibull", df){
  opti_function=ifelse(opti_fun=="weibull", optimise_decay_param, optimise_decay_param_hill)
  opti=optim(c(0.448, 1.11, 0.7), opti_function, insecticide=insecticide, df=df , method = "L-BFGS-B", lower=c(0, 0, 0), upper = c(Inf, Inf, 1) )$par
  names(opti)=c("L", "k", "a")
  opti$setting=insecticide
  opti$decay=opti_fun
  return(opti)
}

update_halflife_insecticideDecay=function(my_parameters_GVI, halflife, kappa, my_EHT, my_netType){
  this.parameters_GVI=my_parameters_GVI%>%
    filter(netType==my_netType, EHT==my_EHT)%>%
    mutate(L_functionalSurvival=halflife,
           kappa_functionalSurvival=kappa)
  
  df=data.frame(time=seq(0, 365*3))%>%
    mutate( weibull=exp( -(time/halflife)^kappa * log(2) ))
  
  for(i in 1:nrow(this.parameters_GVI)){
    this.line=this.parameters_GVI[i,]
    
    if(is.na(this.line$halflife_insecticide)){
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
  return(this.parameters_GVI)
}

