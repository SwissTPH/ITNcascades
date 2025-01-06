optimise_decay_param=function(param, insecticide){
  #if(param[3]<1){
    out=data.frame(
      time=ITN_cov %>% filter(setting==insecticide) %>% pull(time),
      decay_obs=ITN_cov %>% filter(setting==insecticide) %>% pull(ITNcov)
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

calculate_param=function(insecticide, opti_fun="weibull"){
  opti_function=ifelse(opti_fun=="weibull", optimise_decay_param, optimise_decay_param_hill)
  opti=optim(c(0.448, 1.11, 0.7), opti_function, insecticide=insecticide , method = "L-BFGS-B", lower=c(0, 0, 0), upper = c(Inf, Inf, 1) )$par
  names(opti)=c("L", "k", "a")
  opti$setting=insecticide
  opti$decay=opti_fun
  return(opti)
}

