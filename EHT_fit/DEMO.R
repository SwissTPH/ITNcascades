rm(list=ls())
library(dplyr)
library(stringr)
library(tidyr)
library(rstan)

mainDir="."
scriptDir=file.path(".","EHT_fit")

stanDir=file.path(mainDir, "results/stan_outputs")
outputsDir=file.path(mainDir, "results/csv_outputs")
source(file.path(scriptDir, "functions_fit_stan.R"))

niter=1000
nwarmup=500
nchains=3
rerun=FALSE
######################################################
# simulated data
######################################################

mu_0=0.14
alpha_0=0.15
pc_0=0.92

pi=0.69
kappa=3.8
xi=0.52

return_proba=function(alpha, mu, pc){
  proba=c(exp(-alpha-mu), 
                  (1-exp(-alpha-mu))*mu/(alpha+mu), 
                  (1-exp(-alpha-mu))*alpha/(alpha+mu)*pc,
                  (1-exp(-alpha-mu))*alpha/(alpha+mu)*(1-pc))
  names(proba)=c("UA", "UD", "FA", "FD")
  return(proba)
}

proba_control=return_proba(alpha_0, mu_0, pc_0)
proba_i=return_proba(alpha_0*(1-pi), mu_0+alpha_0*kappa, pc_0*(1-xi))
n_samples=1000



sample_control=data.frame(t(rmultinom(n_samples,4, proba_control)))
sample_i=data.frame(t(rmultinom(n_samples,4, proba_i)))

simulated_data=rbind(sample_control %>% mutate(treatment=0, insecticide_name="control"),
      sample_i %>% mutate(treatment=1, insecticide_name="newITN"))%>%
  mutate(fed=FA+FD, total=UA+UD+FA+FD)


fit_EHT_multinomial_stan(data=simulated_data,
                           iter=niter,
                           warmup=nwarmup,
                           chains=nchains,
                           path=file.path(stanDir,"/stan_simulated_data.rds"),
                           stanpath = scriptDir, stanmodel = "EHT_fitting_model.stan")
  

results_simulated<- readRDS(file.path(stanDir,"stan_simulated_data.rds"))

estimates_simulated=get_stan_summary_output(results_simulated, decay="none",
                                                     save=FALSE,
                                                     path=NULL, data=simulated_data)

map_simulated=extract_stan_posteriormax(results_simulated,  data=simulated_data)


