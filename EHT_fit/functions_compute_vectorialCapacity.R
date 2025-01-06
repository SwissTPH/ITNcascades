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
#' @param insecticides list of numbers corresponding to the insecticides you want to get the interventions for
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
                                                     npoints,washedDecay=TRUE, uncertainty=FALSE, cov=NULL, pyrethroid, nsamples=100 ){
  if(!washedDecay){
    results_washed20=NULL
  }
  
  if(uncertainty){
    if(is.null(cov)){
      stop("please provide coverage value when uncertainty is included")
    }
    out=interventions_vectorial_capacity_uncertainty(results=results,results_washed20=results_washed20,
                                                            decay=decay,
                                                            model_p=model_p,
                                                            insecticides=insecticides,
                                                            names=names, L=L, kappa=kappa, nsamples=nsamples, seed=1,  
                                                            npoints=npoints, cov=cov, pyrethroid=pyrethroid)
  } else{
    out=interventions_vectorial_capacity(results=results,results_washed20=results_washed20,
                                                decay=decay,
                                                model_p=model_p,
                                                insecticides=insecticides,
                                                names=names, L=L, kappa=kappa,  
                                                npoints=npoints)
  }

  return(out)
  
}



# interventions_vectorial_capacity
#' @title interventions_vectorial_capacity
#'
#' @description  \code{interventions_vectorial_capacity}
#' The purpose of this function is get an intervention vector that will then be used to calculate the impact of said intervention.
#' You get this intervention vector from stan results. You have to provide which insecticides you want to get the intervention.
#'
#' @param results stanfit file containing stan results for unwashed nets
#' @param results_washed20 stanfit file containing stan results for washed nets
#' @param decay either "weibull", "linear", "hill" or "smooth-compact", is "weibull" by default
#' @param model_p model built by AnophelesModel package
#' @param insecticides list of numbers corresponding to the insecticides you want to get the interventions for
#' @param names list of names corresponding to the insecticide number
#' @param L attrition half life
#' @param kappa attrition shape parameter (for Weibull in particular)
#' @param npoints number of interpolation points for AnophelesModel
#'
#' @return a list of interventions
#'
#' @export


interventions_vectorial_capacity <- function(results,results_washed20,
                                                    decay="weibull",
                                                    model_p,
                                                    insecticides,
                                                    names, L=NULL, kappa=NULL,  
                                                    npoints){
  
  
  mypars=c("InitialPostprandialkillingEfficacy",
           "InitialPreprandialkillingEfficacy",
           "InitialRepellencyRate")
  
  fit_summary <- rstan::summary(results, pars=mypars)
  output <- as.data.frame(fit_summary$summary, optional = TRUE) %>% dplyr::select(mean, sd)
  names(output)=c("mean", "sd")
  
  new_output <- output %>%
    tibble::rownames_to_column( "param") %>%
    tidyr::separate(param, sep="\\[" ,into =c( "param" ,"insecticide")) %>%
    mutate(insecticide=gsub("]", "", insecticide) )
  
  new_output$param[new_output$param=="InitialRepellencyRate"]="InitialRepellentEfficacy"
  
  if(! is.null(results_washed20)){
    fit_summary_washed20 <- rstan::summary(results_washed20, pars=mypars)
    output_washed20 <- as.data.frame(fit_summary_washed20$summary, optional = TRUE) %>% dplyr::select(mean, sd)
    names(output_washed20)=c("mean", "sd")
    
    new_output_washed20 <- output_washed20 %>%
      tibble::rownames_to_column( "param") %>%
      tidyr::separate(param, sep="\\[" ,into =c( "param" ,"insecticide"))  %>%
      mutate(insecticide=gsub("]", "", insecticide) )
    
    new_output_washed20$param[new_output_washed20$param=="InitialRepellencyRate"]="InitialRepellentEfficacy"
    
  } else {
    new_output_washed20=NULL
  }
  
  
  myList <- c()
  for (i in 1:length(insecticides))  {
    myList[i] <- interventions_vectorial_capacity_1insecticide(new_output= new_output,new_output_washed20= new_output_washed20,
                                                                        decay=decay,
                                                                        model_p=model_p,
                                                                        i=insecticides[i],L=L, kappa=kappa,
                                                                        name=names[i], npoints=npoints)
  }
  
  return(myList)
}


######### interventions_vectorial_capacity_1insecticide #####
#' @title interventions_vectorial_capacity_1insecticide
#'
#' @description  \code{interventions_vectorial_capacity_1insecticide}
#'This function is used inside the previous functions "getting_our_values". It gets you the intervention for a specific insecticide i.
#'
#' @param new_output an object
#' @param decay either "weibull", "linear", "hill" or "smooth-compact", is "weibull" by default
#' @param model_p for Anopheles Package
#' @param i the number of the corresponding insecticide you want to get the interventions for
#' @param name name corresponding to the insecticide number
#'
#' @return one intervention
#'
#'
#'
interventions_vectorial_capacity_1insecticide <- function(new_output,new_output_washed20,
                                                                   decay,
                                                                   model_p,
                                                                   i,
                                                                   name,npoints, L, kappa){
  
  
  our_insecticide <- new_output %>% filter(insecticide==i) %>% dplyr::select(-insecticide)
  myproba <-  c((our_insecticide %>%  filter(param=="InitialRepellentEfficacy") %>% dplyr::select(mean))[[1]],
                (our_insecticide %>%  filter(param=="InitialPreprandialkillingEfficacy") %>% dplyr::select(mean))[[1]],
                (our_insecticide %>%  filter(param=="InitialPostprandialkillingEfficacy") %>% dplyr::select(mean))[[1]])
  names(myproba)=c("Deterrency","PrePrandial","PostPrandial")
  
  if(! is.null(new_output_washed20)){
    our_insecticide_washed20 <- new_output_washed20 %>% filter(insecticide==i) %>% dplyr::select(-insecticide)
    myproba_washed20 <-  c((our_insecticide_washed20 %>%  filter(param=="InitialRepellentEfficacy") %>% dplyr::select(mean))[[1]],
                           (our_insecticide_washed20 %>%  filter(param=="InitialPreprandialkillingEfficacy") %>% dplyr::select(mean))[[1]],
                           (our_insecticide_washed20 %>%  filter(param=="InitialPostprandialkillingEfficacy") %>% dplyr::select(mean))[[1]])
    
    halflife_insecticide <- 3/(1-myproba_washed20/myproba)
    names(halflife_insecticide)=c("Deterrency","PrePrandial","PostPrandial")
    
  } else {
    halflife_insecticide = NULL
  }
  
  interv =get_interv(myproba=myproba,
                              L=L, kappa=kappa, npoints = npoints, model_p = model_p, name=name, decay=decay, halflife_insecticide=halflife_insecticide)
  
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
    interv$alphai[,1]=interv$alphai[,2]*(1-myproba["Deterrency"]*get_exposure_multiplier("alphai",model_p, intervention_type)*insecticide_decay_factor_det)
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
    interv$PBi[,1]=interv$PBi[,2]*(1-myproba["PrePrandial"]*get_exposure_multiplier("PBi",model_p, intervention_type)*insecticide_decay_factor_pre)
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
    interv$PCi[,1]=interv$PCi[,2]*(1-myproba["PostPrandial"]*get_exposure_multiplier("PCi",model_p, intervention_type)*insecticide_decay_factor_post)
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
#' @param insecticides list of numbers corresponding to the insecticides you want to get the interventions for
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
                                                             names, L=NULL, kappa=NULL, nsamples=100, seed=1, npoints, cov, pyrethroid){
  
  
  mypars=c("InitialPostprandialkillingEfficacy",
           "InitialPreprandialkillingEfficacy",
            "InitialRepellencyRate")
  set.seed(seed)
  
  fit_extract <- rstan::extract(results, pars=mypars, inc_warmup = FALSE)
  
  mysample=sample(1:nrow(fit_extract[[1]]), nsamples)

  full_df_extract=data.frame(nsample=1:nsamples)
  for(ii in 1:length(mypars)){
    mydf=data.frame(fit_extract[[ii]][mysample,])
    names(mydf)=paste0(mypars[[ii]], "_", 1:ncol(mydf))
    full_df_extract=cbind(full_df_extract,mydf )
  }
  
  new_output=full_df_extract %>%
    tidyr::pivot_longer(cols=names(full_df_extract)[-1], names_to = "param", values_to = "mean")%>%
    tidyr::separate(param, sep="_" ,into =c( "param" ,"insecticide"))
  new_output$param[new_output$param=="InitialRepellencyRate"]="InitialRepellentEfficacy"
  
  if(! is.null(results_washed20)){
    fit_extract_washed20 <- rstan::extract(results_washed20, pars=mypars, inc_warmup = FALSE)
    mysample_washed20=sample(1:nrow(fit_extract_washed20[[1]]), nsamples)
    
    full_df_extract_washed20=data.frame(nsample=1:nsamples)
    for(ii in 1:length(mypars)){
      mydf_washed20=data.frame(fit_extract_washed20[[ii]][mysample_washed20,])
      names(mydf_washed20)=paste0(mypars[[ii]], "_", 1:ncol(mydf))
      full_df_extract_washed20=cbind(full_df_extract_washed20,mydf_washed20 )
    }
    
    
    new_output_washed20=full_df_extract_washed20 %>%
      tidyr::pivot_longer(cols=names(full_df_extract)[-1], names_to = "param", values_to = "mean")%>%
      tidyr::separate(param, sep="_" ,into =c( "param" ,"insecticide"))
    new_output_washed20$param[new_output_washed20$param=="InitialRepellencyRate"]="InitialRepellentEfficacy"
    
  } else {
    new_output_washed20=NULL
  }

  VCred = c()
  
  for(this_sample in 1:nsamples){
    myList <- c()
    for (i in 1:length(insecticides))  {
      this.new_output=new_output %>% dplyr::filter(nsample ==this_sample) %>% dplyr::select(-nsample)
      if(! is.null(new_output_washed20)){
        this.new_output_washed20 = new_output_washed20 %>% dplyr::filter(nsample ==this_sample) %>% dplyr::select(-nsample)
      } else {
        this.new_output_washed20 = NULL
      }
      myList[i] <- interventions_vectorial_capacity_1insecticide(
        new_output= this.new_output,
        new_output_washed20=this.new_output_washed20,
        decay=decay,model_p=model_p,i=insecticides[i],name=paste0(names[i], "_", this_sample),
        L=L, kappa=kappa, npoints = npoints)
    }
    
    this.VCred=   calculate_combined_impact(intervention1=myList[[1]],
                                                intervention2=pyrethroid[[1]],
                                                cov_intervention1 =  cov,
                                                cov_intervention2 = cov,
                                                cov_combination = cov,
                                                N_vec =10000)
    
    VCred=c(VCred, this.VCred)
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

###############################

# Plot vectorial capacity reduction with uncertainty
#' @title plot_VC_uncertaintyIntervention
#'
#' @description  \code{plot_VC_uncertaintyIntervention}
#' Plot vectorial capacity reduction including uncertinaty in GVI parameters
#'
#' @param impacts output of the calculate_impact function in the AnophelesModel package
#' @param list_colors list of colors to be included
#' @param alpha alpha opacity for ribbon in the ggplot2 package
#'
#' @return a list of interventions
#'
#' @export
plot_VC_uncertaintyIntervention=function(impacts, list_colors=c( "red", "orange", "dodgerblue","darkblue"), alpha=0.2){
  
  get_impact_df=function(impacts_obj) {
    plot_df = NULL
    for(interv_imp in impacts_obj$interventions_vec) {
      intervention_name = interv_imp$description
      intervention_coverage = interv_imp$coverages
      intervention_impact = round(interv_imp$effects$avg_impact,
                                  digits = 3)
      plot_tab = cbind.data.frame(intervention_name,
                                  intervention_coverage,
                                  intervention_impact)
      plot_df = rbind.data.frame(plot_df, plot_tab)
    }
    return(plot_df)
  }
  
  impact_df=get_impact_df(impacts) %>%
    tidyr::separate(intervention_name,sep="_" ,into =c( "intervention_name" ,"nsample")) %>%
    group_by(intervention_name, intervention_coverage) %>%
    summarise(mean=mean(intervention_impact),
              min=min(intervention_impact),max=max(intervention_impact))
  
  p = ggplot(impact_df) +
    theme_light() + theme_linedraw() + theme_bw(base_size = 9) +
    theme(legend.position = "bottom") +
    ylim(c(0, 100)) +
    geom_line(aes(x=intervention_coverage, y=mean*100,
                  col=factor(intervention_name)), size=1) +
    geom_ribbon(aes(x=intervention_coverage, ymin=min*100, ymax=max*100,
                    fill=factor(intervention_name)), size=1, alpha=alpha) +
    guides(color = guide_legend(nrow = 2,byrow = TRUE)) +
    labs(fill="",
         x="Coverage", y="Mean reduction in \nvectorial capacity (%)",
         col="")+
    scale_color_manual(values=list_colors)+
    scale_fill_manual(values=list_colors)
  
  return(p)
  
}



