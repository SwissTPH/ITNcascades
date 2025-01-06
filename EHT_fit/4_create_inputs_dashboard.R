library(dplyr)
library(rstan)
library(readr)
mainDir="."
dirInputs=file.path(mainDir, "stan_outputs")
scriptDir=file.path(mainDir,"ITNcascadesdashboard")
dirDashboard=file.path(scriptDir, "Inputs/")

mypars=c("InitialPostprandialkillingEfficacy",
         "InitialPreprandialkillingEfficacy",
         "InitialRepellencyRate")

save_param=function(results_EHT, name_EHT, write=T){
  fit_summary <- rstan::summary(results_EHT, pars=mypars)
  output <- as.data.frame(fit_summary$summary, optional = TRUE) %>% dplyr::select(mean, sd)
  names(output)=c("mean", "sd")
  
  new_output <- output %>%
    tibble::rownames_to_column( "param") %>%
    tidyr::separate(param, sep="\\[" ,into =c( "param" ,"insecticide")) %>%
    tidyr::separate(insecticide , sep="\\]" , into= c("insecticide"))
  
  new_output$param[new_output$param=="InitialRepellencyRate"]="InitialRepellentEfficacy"
  
  #return(new_output)
  
  if(write){
    write_rds(new_output, file = file.path(dirDashboard,paste0("/output_",name_EHT,".rds")))
  } else {
      return(new_output)
    }
  
}


results_kibondo_w<- readRDS(file.path(dirInputs,"/stan_kibondo_washed_controlUnw_72.rds"))
results_kibondo_unw<- readRDS(file.path(dirInputs,"/stan_kibondo_unwashed_72.rds"))

results_bit055_w<- readRDS(file.path(dirInputs,"/stan_BIT055_washed_controlUnw_24.rds"))
results_bit055_unw<- readRDS(file.path(dirInputs,"/stan_BIT055_unwashed_24.rds"))

results_bit059_w<- readRDS(file.path(dirInputs,"/stan_BIT059_washed_controlUnw_24.rds"))
results_bit059_unw<- readRDS(file.path(dirInputs,"/stan_BIT059_unwashed_24.rds"))

results_bit103_unw<- readRDS(file.path(dirInputs,"/stan_BIT103_unwashed_72_ifakara.rds"))
results_bit103_w<- readRDS(file.path(dirInputs,"/stan_BIT103_washed_controlUnw_72_ifakara.rds"))

results_bit080_unw<- readRDS(file.path(dirInputs,"/stan_BIT080_unwashed_72.rds"))
results_bit080_w<- readRDS(file.path(dirInputs,"/stan_BIT080_washed_controlUnw_72.rds"))


save_param(results_kibondo_unw, "kibondo_unwashed")
save_param(results_kibondo_w, "kibondo_washed")
save_param(results_bit055_unw, "bit055_unwashed")
save_param(results_bit055_w, "bit055_washed")
save_param(results_bit103_unw, "bit103_unwashed")
save_param(results_bit103_w, "bit103_washed")
save_param(results_bit059_unw, "bit059_unwashed")
save_param(results_bit059_w, "bit059_washed")
save_param(results_bit080_unw, "bit080_unwashed")
save_param(results_bit080_w, "bit080_washed")
