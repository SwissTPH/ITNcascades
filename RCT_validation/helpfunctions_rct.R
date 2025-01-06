compute_CI = function(est,size,crit=1.96){
  sd=sqrt(est*(1-est)/size)
  LCI=est-crit*sd
  UCI=est+crit*sd
  return(list(sd=sd,LCI=LCI,UCI=UCI))
}



calibration=function(my_all_simul, fitdat, select_subset =.5){
  minEIR=min(my_all_simul$EIR)
  maxEIR=max(my_all_simul$EIR)

  simul_fit = fitdat %>%
    dplyr::select(c(setting,age,year,month)) %>%
    #keep only one future net type per setting, in any case calibration is before deployment
    left_join(my_all_simul %>% filter(year %in% fitdat$year) %>% rename("age"="age_group", "PR"="prevalenceRate") %>% mutate(EIR_num=EIR)) %>%
    mutate(sub=setting) %>%
    mutate(futITNcov=0) %>%
    ungroup() %>% unique()

  # for each seed and each setting, get EIR that minimises distance to observed prevalence
  best_EIRs=get_best_EIR(simul_fit,data_to_fit = fitdat)
  # quadratic fit
  quadratic_l_approx=get_quadratic_fit_cf(df=best_EIRs$best_EIR_ML_gaussian
                                          , df_all =  best_EIRs$all_llk, select_subset =select_subset, lambda = 1
                                          ,data_to_fit = fitdat)
  best_EIRs_l_ci_approx=get_profile_llk_ci_approx(quadratic_l_approx,minEIR,maxEIR
                                                  ,lower_threshold = minEIR+2
                                                  ,upper_threshold = maxEIR-10
                                                  ,upper_threshold2 = maxEIR-40,
                                                  inside_window = TRUE)

  calibration_reformat=best_EIRs_l_ci_approx %>%
    ungroup() %>%
    dplyr::select(setting,sub,starts_with("EIR")) %>%
    pivot_longer(cols=starts_with("EIR"),names_to="EIR_type",names_prefix="EIR",values_to="EIR_num") %>%
    dplyr::filter(EIR_type %in% c("_lci","mle_tilde_round","_uci")) %>%
    mutate(EIR_type=ifelse(EIR_type=="mle_tilde_round","middle",ifelse(EIR_type=="_lci","lower","upper"))) %>%
    unique() %>%
    mutate(EIR_num=ifelse(is.na(EIR_num),ifelse(EIR_type=="uci",maxEIR,minEIR),EIR_num))

  return(list("best_EIRs_l_ci_approx"=best_EIRs_l_ci_approx, "quadratic_l_approx"=quadratic_l_approx,
              "calibration_reformat"=calibration_reformat))
}



merge_calibration_simulation=function(calibration_reformat, my_all_simul){

  my_all_simul= my_all_simul%>% mutate(EIR_num=as.numeric(EIR)) %>%rename("PR"="prevalenceRate", "age"="age_group")
  simul_calib=inner_join(calibration_reformat, my_all_simul,by=c("setting","EIR_num")) %>%
    unique()%>%
    dplyr::select(EIR_type,setting, sub, seed,year,month, age,PR,EIR_num )%>%
    tidyr::pivot_wider(names_from=EIR_type, values_from=all_of(c("PR","EIR_num"))
                       , id_cols=c(setting, sub, seed,year,month, age
                       )) %>%
    rowwise()


  for (o in c("PR","EIR_num")){
    simul_calib=simul_calib %>% mutate(!! sym(paste(o,"inf",sep="_")) := min(!! sym(paste(o,"lower",sep="_"))
                                                                             ,!! sym(paste(o,"middle",sep="_"))
                                                                             ,na.rm=T)
                                       ,!! sym(paste(o,"sup",sep="_")) := max(!! sym(paste(o,"upper",sep="_"))
                                                                              ,!! sym(paste(o,"middle",sep="_"))
                                                                              ,na.rm=T))

  }

  simul_calib=simul_calib %>% group_by(setting,year,month,age)

  for (q in c("PR")){
    simul_calib=simul_calib %>%
      mutate(!! sym(paste(q,"mini",sep="_")) := min(!! sym(paste(q,"inf",sep="_")),na.rm=T)
             ,!! sym(paste(q,"maxi",sep="_")) := max(!! sym(paste(q,"sup",sep="_")),na.rm=T))
  }

  #average across seeds
  seed_avg = simul_calib %>% group_by(setting,year,month,age) %>%
    mutate(PR_middle=as.numeric(PR_middle)) %>%
    summarise_if(is.numeric,.funs=mean, na.rm=TRUE) %>%
    mutate(name_month=ifelse(month==1,"Jan",""))

  return(list("simul_calib"=simul_calib, "seed_avg"=seed_avg))
}


merge_calibration_simulation_uncertainty=function(calibration_reformat, my_all_simul){
  simul_calib=inner_join(calibration_reformat, my_all_simul,by=c("setting","EIR", "EIR_type")) %>%
    rename("PR"=prevalenceRate)%>%
    unique()%>%
    tidyr::pivot_wider(names_from=EIR_type, values_from=all_of(c("PR","EIR"))
                       , id_cols=c(setting, seed,year,month, age_group
                       )) %>%
    rowwise()


  for (o in c("PR","EIR")){
    simul_calib=simul_calib %>% mutate(!! sym(paste(o,"inf",sep="_")) := min(!! sym(paste(o,"lower",sep="_"))
                                                                             ,!! sym(paste(o,"middle",sep="_"))
                                                                             ,na.rm=T)
                                       ,!! sym(paste(o,"sup",sep="_")) := max(!! sym(paste(o,"upper",sep="_"))
                                                                              ,!! sym(paste(o,"middle",sep="_"))
                                                                              ,na.rm=T))

  }

  simul_calib=simul_calib %>% group_by(setting,year,month,age_group)

  for (q in c("PR")){
    simul_calib=simul_calib %>%
      mutate(!! sym(paste(q,"mini",sep="_")) := min(!! sym(paste(q,"inf",sep="_")),na.rm=T)
             ,!! sym(paste(q,"maxi",sep="_")) := max(!! sym(paste(q,"sup",sep="_")),na.rm=T))
  }

  #average across seeds
  seed_avg = simul_calib %>% group_by(setting,year,month,age_group) %>%
    mutate(PR_middle=as.numeric(PR_middle)) %>%
    summarise_if(is.numeric,.funs=mean, na.rm=TRUE) %>%
    mutate(name_month=ifelse(month==1,"Jan",""))

  return(list("simul_calib"=simul_calib, "seed_avg"=seed_avg))
}


plot_validation=function(my_seed_avg, fitdat, validat, distribution_time=2019+27/365,
                         myAgePR, mylimits=c(2018.5,2021.5), mybreaks=2018:2021,
                         mylabel="Prevalence in 6 months\nto 14 years old"){
  p=ggplot(my_seed_avg %>% filter(age==myAgePR), aes(x=year+(month-1)/12))+
    geom_ribbon(aes(ymin=PR_mini,ymax=PR_maxi),fill="grey",alpha=0.5)+
    geom_line(aes(y=PR_middle),color="grey",size=1.5)+
    geom_vline(aes(xintercept=distribution_time,linetype="Net distribution"),size=1.5)+
    geom_pointrange(data=fitdat,aes(y=PR_obs,ymin=PR_obs_min,ymax=PR_obs_max,fill="Calibration points"),
                    size=3,stroke=3,shape=21,color=rgb(44,160,137,maxColorValue = 255))+
    geom_pointrange(data=validat,aes(y=PR_obs,ymin=PR_obs_min,ymax=PR_obs_max,fill="Validation points"),
                    size=3,color=rgb(44,160,137,maxColorValue = 255),shape=21)+
    scale_y_continuous(labels=scales::percent)+
    scale_fill_manual(values=c("white",rgb(44,160,137,maxColorValue = 255)))+
    scale_x_continuous(limits=mylimits,breaks=mybreaks,labels=function(x){return(paste0("Jan\n",x))})+
    labs(x="Time",y=mylabel,fill="",linetype="")+
    facet_wrap(.~setting)+
    theme_minimal(base_size = 60)
  return(p)
}


#'Plot the likelihood resulting from ClaraFit, the quadratic approximation and best EIR.
#'
#' @param quadratic_l_approx Result data frame from \emph{get_quadratic_fit_cf} from \emph{openMalariaUtilities}
#' @param list_setting Vector of admin1 names. Better to have only one admin1, otherwise there'll be lot's of subplots. If empty, the default is the name of the first found setting.
#' @param list_sub Vector of admin2 names to include in the plot. If empty, the default are all the districts in the specified list_setting.
#' @param llk_gaussian Boolean variable. Default is TRUE, will show the likelihood curves associated with llk_gaussian. If false, will show the values from llk_laplace.
#' @return Likelihood curve and quadratic approximation, per admin2.
#' @examples
#'plotQuadraticApprox(quadratic_l_approx = quadratic_l_approx
#'                    , list_setting = "maputocidade"
#'                    , list_sub = "kanyaka")
#' @import ggthemes
#' @export
plotQuadraticApprox <- function(quadratic_l_approx
                                , list_setting = NULL
                                , list_sub = NULL
                                , llk_gaussian = TRUE){
  if(is.null(list_setting)) list_setting = quadratic_l_approx$setting[1]
  if(is.null(list_sub)) list_sub = unique(quadratic_l_approx$sub[which(quadratic_l_approx$setting %in% list_setting)])
  if(llk_gaussian){
    g = ggplot2::ggplot(quadratic_l_approx %>% dplyr::filter(sub %in% list_sub & setting %in% list_setting)
    ) +
      ggplot2::geom_line(aes(x=EIR_num, y=llk_gaussian, group=seed), color="grey68")+
      ggplot2::geom_line(aes(x=EIR_num, y=llk_tilde, color= "approx"))+
      ggplot2::geom_line(aes(x=EIR_num, y=smoothed_llk, color= "smoothed"))+
      ggplot2::geom_vline(aes(xintercept = EIRmle_tilde,  color= "approx"))+
      ggplot2::facet_wrap(setting ~sub, scales = "free")+
      ggplot2::theme_minimal()+
      ggplot2::labs(x = "EIR", y = "llk Gaussian", color = element_blank())+
      ggthemes::scale_color_tableau()

  } else {
    g = ggplot2::ggplot(quadratic_l_approx %>% dplyr::filter(sub %in% list_sub & setting %in% list_setting)
    ) +
      ggplot2::geom_line(aes(x=EIR_num, y=llk_laplace, group=seed), color="grey68")+
      ggplot2::geom_line(aes(x=EIR_num, y=llk_tilde, color= "approx"))+
      ggplot2::geom_line(aes(x=EIR_num, y=smoothed_llk, color= "smoothed"))+
      ggplot2::geom_vline(aes(xintercept = EIRmle_tilde,  color= "approx"))+
      ggplot2::facet_wrap(setting ~sub, scales = "free")+
      ggplot2::theme_minimal()+
      ggplot2::labs(x = "EIR", y = "llk Laplace", color = element_blank())+
      ggthemes::scale_color_tableau()
  }
  return(g)
}


#'Plot the likelihood resulting from ClaraFit, the quadratic approximation, best EIR & min and max EIR (C.I).
#'
#' @param best_EIRs_l_ci_apporx Result data frame from \emph{get_profile_llk_ci_approx} from \emph{openMalariaUtilities}
#' @param quadratic_l_approx Result data frame from \emph{get_quadratic_fit_cf} from \emph{openMalariaUtilities}
#' @param list_setting Vector of admin1 names. Better to have only one admin1, otherwise there'll be lot's of subplots. If empty, the default is the name of the first found setting.
#' @param list_sub Vector of admin2 names to include in the plot. If empty, the default are all the districts in the specified list_setting.
#' @param llk_gaussian Boolean variable. Default is TRUE, will show the likelihood curves associated with llk_gaussian. If false, will show the values from llk_laplace.
#' @return Likelihood curve, quadratic approximation and EIR uncertainty range, per admin2.
#' @examples
#'plotQuadraticApproxMaxMinEIR(quadratic_l_approx = quadratic_l_approx
#'                             , best_EIRs_l_ci_approx = best_EIRs_l_ci_approx
#'                             , list_setting = "maputocidade")
#' @import ggplot2
#' @export
plotQuadraticApproxMaxMinEIR <- function(quadratic_l_approx
                                         , best_EIRs_l_ci_approx
                                         , list_setting = NULL
                                         , list_sub = NULL
                                         , llk_gaussian = TRUE){

  if(is.null(list_setting)) list_setting = quadratic_l_approx$setting[1]
  if(is.null(list_sub)) list_sub = unique(quadratic_l_approx$sub[which(quadratic_l_approx$setting %in% list_setting)])
  foo = best_EIRs_l_ci_approx %>% dplyr::filter(sub %in% list_sub & setting %in% list_setting)
  g = plotQuadraticApprox(quadratic_l_approx, list_setting, list_sub, llk_gaussian)
  g = g + ggplot2::geom_vline(data = foo, aes(xintercept = EIR_lci), linetype = "dotted") +
    ggplot2::geom_vline(data = foo, aes(xintercept = EIR_uci), linetype = "dotted")

  return(g)

}

