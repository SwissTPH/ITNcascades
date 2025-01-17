library(devtools)
#devtools::install_github("SwissTPH/r-openMalariaUtilities",
#                         ref = "master")
#devtools::install("~/Git_repo/omuaddons") #run this if you have pulled in your local directory

library(openMalariaUtilities)
library(tidyverse)
library(ggplot2)
library(OMAddons)
library(omuslurm)
library(omucompat)

rm(list=ls())

## Demographic data for the following countries included:
## BEN, CMR, GHA, HTI, MOZ, TZA, UGA
COUNTRY <- "TZA"
## Start date of monitoring. OM's recommendation is to start 20 to 50 years
## before the first measures are deployed.
SIMSTART <- "1918-01-01"
## openMalariaUtilities supports Open Malaria version 43 and up
versionnum <- 44L



source("helpfunctions_rct.R")
source("helpfunction_omu_gvi.R")
main_folder=(".")
parameters_GVI=read.csv(file.path(main_folder, "fitted_parameters_all_new_VCred.csv"))

expName <- "RCT_protopopoff_validation"
setwd(main_folder)
experimentDir <- paste0(main_folder,expName)
experiment_folder=file.path(main_folder, expName)

## Define mosquitoes
mosqs <- c("arabiensis_indoor","funestus_indoor", "gambiae_indoor")#personal communication
contrib <- c("@ain@","@fin@", "@gin@" )
cbind( mosqs, contrib)

## Basic skeleton
baseList <- list(
  ## Mandatory
  expName = expName,
  ## Optional
  ## subDir = "sub",
  ## Mandatory
  OMVersion = versionnum,
  ## Optional
  ## analysisNo = 1L,
  ## Optional
  ## xmlBasename = "gha_base",
  ## Mandatory
  demography = list(),
  monitoring = list(),
  interventions = list(),
  healthSystem = list(),
  entomology = list(),
  ## These are optional for OM
  ## parasiteGenetics = list(),
  ## pharmacology = list(),
  ## diagnostics = list(),
  model = list()
)

## Create demography
baseList <- write_demography_compat(
  baseList = baseList, pop = "@pop@", country = COUNTRY
)

## Create monitoring
baseList <- write_monitoring_compat(
  baseList = baseList, y1 = 2011, y2 = 2019, detect = 100, SIMSTART = SIMSTART,
  upperbounds = c(.5, 14,100)
)

## Entomology section: MANDATORY
baseList <- make_ento_compat(
  baseList = baseList, mosqs, contrib,
  EIR = "@EIR@",
  seasonality = paste0("@m", 1:12, "@"),
  propInfected = .078,
  propInfectious = .021
)

## Define historical ITN
baseList <- define_ITN(
  baseList = baseList, component = "dummyITN", mosquitos = mosqs, historical = TRUE,
  resist = TRUE, halflife = 0, strong=FALSE
)

exposure=.69
halflife_PBO=2.4
kappa_PBO=1.7

histITN= create_vectorInterventionParameters(deterrency = 0,
                                                preprandial =0.05,
                                                postprandial = 0,
                                                kappa=3.1,
                                                L=1, myname = "histITN", decay = "step", mosqs =mosqs)

## Define standard ITN(setting 13.1 has longer half-life)
futITN= create_vectorInterventionParameters(deterrency = 0,
                                               preprandial =  0.05,
                                               postprandial = 0,
                                               kappa=3.1,
                                               L=2.7, myname = "futITN", mosqs = mosqs)



GVI_PBO_BIT055= extract_GVI_params(EHT="BIT055",
                                   netType="Olyset Plus",
                                   halflife_functionalSurvival=halflife_PBO, kappa_functionalSurvival=kappa_PBO,
                                   exposure_correction=exposure, myname="futPBO_BIT055", mosqs=mosqs, parameters_GVI=parameters_GVI)

GVI_PBO_BIT059= extract_GVI_params(EHT="BIT059",
                                   netType="Olyset Plus",
                                   halflife_functionalSurvival=halflife_PBO, kappa_functionalSurvival=kappa_PBO,
                                   exposure_correction=exposure, myname="futPBO_BIT059", mosqs=mosqs, parameters_GVI=parameters_GVI)

GVI_PBO_BIT103= extract_GVI_params(EHT="BIT103",
                                   netType="Olyset Plus",
                                   halflife_functionalSurvival=halflife_PBO, kappa_functionalSurvival=kappa_PBO,
                                   exposure_correction=exposure, myname="futPBO_BIT103", mosqs=mosqs,  parameters_GVI=parameters_GVI)


list_GVI_snippets=list(histITN, futITN)

list_GVI_snippets=append(list_GVI_snippets, GVI_PBO_BIT103)
list_GVI_snippets=append(list_GVI_snippets, GVI_PBO_BIT055)
list_GVI_snippets=append(list_GVI_snippets, GVI_PBO_BIT059)

for(GVI_param in list_GVI_snippets){
  baseList= defineGVI_simple(baseList = baseList,
                             vectorInterventionParameters=GVI_param,
                             append = TRUE,
                             verbatim = TRUE,
                             hist = FALSE)
}

## Deployment section

## Historical ITN
baseList <- deploy_it_compat(
  baseList = baseList, component = "histITN", coverage = "@histITNcov@",
  ## Allowing for different hist coverage levels
  byyear = FALSE,
  ## Annual deployments
  y1 = 2012, y2 = 2014, every = 1, interval = "year",
  m1 = 2, m2 = 2, d1 = 25, d2 = 25, #to end just when the successive distribution starts
  SIMSTART = SIMSTART
)

## Future ITN
baseList <- deploy_it_compat(
  baseList = baseList, component = "futITN",
  coverage = "@futITNcov@",
  ## One future level of coverage
  byyear = FALSE,
  y1 = 2015, y2 = 2015, every = 3, interval = "year",
  m1 = 2, m2 = 2, d1 = 25, d2 = 25, #to have June-July survey really 4 months after distribution
  SIMSTART = SIMSTART
)


## Future IG2
list_GVI_snippets_id=c("futPBO_BIT055","futPBO_BIT055_sup","futPBO_BIT055_inf",
                       "futPBO_BIT059","futPBO_BIT059_sup","futPBO_BIT059_inf",
                       "futPBO_BIT103","futPBO_BIT103_sup","futPBO_BIT103_inf")


for(GVI_snippet_id in list_GVI_snippets_id){
  baseList <- deploy_it_compat(
    baseList = baseList, component = paste0(GVI_snippet_id, "_deterrency") ,
    coverage = paste0("@", GVI_snippet_id,"_cov@"),
    ## allowing for different hist coverage levels
    byyear = FALSE,
    ##  annual deployments
    y1 = 2015, y2 = 2015, every = 3, interval = "year",
    m1 = 2, m2 = 2, d1 = 25, d2 = 25, #to have June-July survey really 4 months after distribution
    SIMSTART = SIMSTART
  )

  # Append additional component
  # Get index
  idx <- extractList(baseList, value = paste0(GVI_snippet_id, "_deterrency"), onlyIndex = TRUE)
  ## The "component" is the second entry, the first one is the "name" with the
  ## same value
  ## Because the order of the entries matters, we need to splice the new
  ## component id into the list at the correct position
  ## The "-2" means we go two levels up which results in full entry of interest.
  ## in this list we need to add the new component after the old one (which is
  ## the second entry) and then append the rest again.

  tmp <- c(baseList[[head(idx[[2]], -1)]][1:2],
           list(component = list(id =  paste0(GVI_snippet_id, "_preprandial"))),
           list(component = list(id = paste0(GVI_snippet_id, "_postprandial"))),
           list(component = list(id = "futITN")),
           baseList[[head(idx[[2]], -1)]][3:length(baseList[[head(idx[[2]], -1)]])])
  baseList[[head(idx[[2]], -1)]] <- tmp
}



## Importation: MANDATORY
baseList <- define_importedInfections_compat(baseList = baseList, 10, time = 0)


## Write health system: MANDATORY
baseList <- write_healthsys_compat(
  baseList = baseList, access = "@EffCovconv@"
)


## Specify seed: MANDATORY
baseList <- write_end_compat(
  baseList = baseList, seed = "@seed@", modelname = "base"
)


setwd(main_folder)
## Create base XML file

setupDirs(experimentName = expName)
#createBaseXml(baseList, replace = TRUE)


## Copy necessary Open Malaria files. Currently needs to be called after
## createBaseXml, otherwise the cache is not set up.
setupOM(version = versionnum, dir=experiment_folder)


#########################################
### seasonality
season_dir="/scicore/home/pothin/champa0000/IG2modelling/"
seasonality_chirps=read.csv(file.path(season_dir, "seasonality_Muleba.csv"))

seasonality_full=rbind(seasonality_chirps %>% mutate(season="chirps"))
names(seasonality_full)=c("m4", "m5", "m6","m7","m8","m9" ,"m10", "m11", "m12", "m1", "m2", "m3", "season")


##-- creating the experiment
full             = list()
full $ seed      = 1:10
full $ setting   = c("pyrethroid",
                     "PBO_BIT055","PBO_BIT055_sup","PBO_BIT055_inf",
                     "PBO_BIT059","PBO_BIT059_sup","PBO_BIT059_inf",
                     "PBO_BIT103","PBO_BIT103_sup","PBO_BIT103_inf" )
#c("pyrethroid","PBO_BIT055_briet", "PBO_BIT055_denz")
full $ pop       = 10000
full$ EffCovconv = convert_cm(0.3782992) #DHS 2016, effective coverage, code by Katya an Rémi
full $ EIR =seq(10,200,2)
# species composition: personal communication
full$ain = .04 
full$fin = .04 
full$gin = .92 
full$season = c("chirps")

#### 'scens' will do all possible combinations of those factors
scens = expand.grid( full )
scens = scens %>% left_join(seasonality_full) #%>%



futPBO_cov=0.81
NetCov <- data.frame(setting=c("pyrethroid","PBO"),
                     histITNcov=c(902/2996,810/3078),
                     futITNcov=c(0.77 , 0)
)

NetCov_full=NetCov
for(mysetting in full$setting[-1]){
  print(mysetting)
  myNetCov=NetCov %>%
    filter(setting %in% ("PBO"))%>%
    mutate(futcov=futPBO_cov, setting=mysetting)
  names(myNetCov)[4]=paste0("fut", mysetting, "_cov")
  NetCov_full=merge(NetCov_full, myNetCov, all=TRUE)
}
NetCov_full[is.na(NetCov_full)]=0

scens <- left_join(scens,NetCov_full)

scens <- data.frame(scens)



scens <- finalize_scenarios(scens)
storeScenarios(scens)

dim(scens)

length(unique(scens$setting))

# small check on the scenarios, to check if the scenario has the right placeholder names
# setup_scenarios(scenarios=scens[])

## Generate scenarios

## SLURM
slurmPrepareScenarios(expName = expName, scenarios = scens, nCPU = 1, qos = "30min", time = "00:30:00", 
                      rModule = "R/4.4.1-gfbf-2023b")
#slurmCreateScenarios()
check_scenarios_created(experiment_folder)

validateXML(xmlfile = file.path(experiment_folder, paste0("scenarios/", expName, "_7.xml")), schema = NULL, scenarios = scens)
## Run Open Malaria
## SLURM
slurmPrepareSimulations(expName = expName, scenarios = scens, nCPU = 1, bSize = 1, qos = "6hours",
                        rModule = "R/4.4.1-gfbf-2023b", 
                        omModule = "OpenMalaria/44.0-intel-compilers-2023.1.0")

#slurmRunSimulation()
check_simulations_created(experiment_folder)

# failed=read.table(file.path(experiment_folder, "missing_om_simulations.txt"))
# failed$V2=gsub("_out.txt","", gsub(paste0(expName, "_"), "", failed$V1))
# #
# paste0(failed$V2,collapse=",")




loadExperiment(expName)
scens=readRDS(file.path(experiment_folder, "cache/scenarios.rds"))
## Post-processing
age_groups_of_interest <- c("0.5-14", "All")

## SLURM
slurmPrepareResults(

  expDir=experiment_folder,
  dbName="results_RCTprotopopoff_validation",
  dbDir = experiment_folder,
  resultsName = "results_RCT13",
  resultsCols = c("scenario_id", "date_aggregation", "date", "age_group",
                  "nUncomp", "nHost", "nSevere",
                  "incidenceRate", "prevalenceRate"),
  indexOn = list(c("results", "scenario_id")),
  ncoresDT = 1,
  strategy = "batch",
  appendResults = FALSE,
  aggrFun = CalcEpiOutputs,
  aggrFunArgs = list(aggregateByAgeGroup = age_groups_of_interest, aggregateByDate = "month"),
  splitBy = 30000,
  verbose = FALSE,
  ntasks = 1,
  mem = "32G",
  nCPU = 10,
  time = "06:00:00",
  qos = "6hours", 
  rModule = "R/4.4.1-gfbf-2023b"
)
#slurmRunResults()

####
#### Step 4: Fitting to historical data
####

####Fitting with Clara's method
dir.create(file.path(experiment_folder, "calibration/"))
FittingDir=file.path(experiment_folder, "calibration/")

# PREPARE THE DATA TO FIT

# open database connection
conn = DBI::dbConnect(RSQLite::SQLite(),file.path(experiment_folder,"results_RCTprotopopoff_validation.sqlite"))
# list the tables of the database
DBI::dbListTables(conn)
# extract the results table, this should have the same name as provided to slurmPrepareResults
all_simul_0 = DBI::dbReadTable(conn, "results_RCT13")
# always disconnect from the database
DBI::dbDisconnect(conn)

all_simul=merge(all_simul_0, scens %>%mutate(scenario_id=ID), by="scenario_id")%>%
  mutate(EIR=as.numeric(EIR),
         month=as.numeric(format(as.Date(date), "%m")),
         year=as.numeric(format(as.Date(date), "%y"))+2000)



####
#### Step 4: Fitting to historical data
####

# PREPARE THE DATA TO FIT
scens$EIR=as.numeric(scens$EIR)

all_simul_chirps=all_simul
#write.csv(all_simul_chirps, file.path(FittingDir, "all_simul_chirps.csv"), row.names = F)

#all_simul_chirps=read.csv(file.path(FittingDir, "all_simul_chirps.csv"))

agePR="0.5-14"

###############################
## fitdat
######################


#'Table 1 Protopopoff 2018 https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(18)30427-6/fulltext
#'timing: Figure 1
fitdat = data.frame(setting=c("Pyrethroid", "OlysetPlus_BIT055","OlysetPlus_BIT103","OlysetPlus_BIT059")
                    ,age="0.5-14"
                    ,year=2014
                    ,month=9
                    ,pos_children=c(600,606, 606, 606)
                    ,nb_children=c(885,991, 991, 991)) %>%
  rowwise() %>%
  mutate(PR_obs=pos_children/nb_children,
         sd=sqrt(PR_obs*(1-PR_obs)/nb_children),
         #sd=4*sqrt(PR_obs*(1-PR_obs)/nb_children),
         PR_obs_min=PR_obs-1.96*sd,
         PR_obs_max=PR_obs+1.96*sd) %>%
  separate(setting, c("trial","arm"),remove = F)


######################
## validat
######################

setting=c("Pyrethroid", "OlysetPlus_BIT055","OlysetPlus_BIT103", "OlysetPlus_BIT059")
#' Protopopoff et al. 2018 https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(18)30427-6/fulltext
#' Figure 1 for timing and Table 2 for numbers
validat = data.frame(setting=rep(c("Pyrethroid", "OlysetPlus_BIT055","OlysetPlus_BIT103","OlysetPlus_BIT059"),each=6)
                     ,age="0.5-14"
                     ,year=rep(c(2015,2015,2016,2016, 2017, 2017),4)
                     ,month=rep(c(6,11,6,11, 6, 11),4)
                     ,pos_children=c(553,515,548,710, 1643, 1045,
                                     445,275,342,440,1280,901,
                                     445,275,342,440,1280,901,
                                     445,275,342,440,1280,901
                     )
                     ,nb_children=c(997,932,1034,1044,2032,1784,
                                    971,883,984,958,1848,1807,
                                    971,883,984,958,1848,1807,
                                    971,883,984,958,1848,1807
                     )
)%>%
  rowwise() %>%
  mutate(PR_obs=pos_children/nb_children,
         sd=sqrt(PR_obs*(1-PR_obs)/nb_children),
         PR_obs_min=PR_obs-1.96*sd,
         PR_obs_max=PR_obs+1.96*sd)


##########################################################################
##############
##############################################
## with chirps data
all_simul_chirps$setting=gsub("PBO", "OlysetPlus",gsub("pyrethroid", "Pyrethroid",all_simul_chirps$setting))
calibrated_data=calibration(my_all_simul=all_simul_chirps, fitdat, select_subset =0.5)

#check calibration
plotQuadraticApproxMaxMinEIR(quadratic_l_approx = calibrated_data$quadratic_l_approx
                             , best_EIRs_l_ci_approx = calibrated_data$best_EIRs_l_ci_approx
                             ,list_setting = fitdat$setting)+
  xlim(0, 300)



calibrated_data_uncertainty=rbind(calibrated_data$calibration_reformat,
                                  calibrated_data$calibration_reformat %>% dplyr::filter(setting !="pyrethroid")%>% mutate(setting=paste0(setting,"_inf")),
                                  calibrated_data$calibration_reformat %>% dplyr::filter(setting !="pyrethroid")%>% mutate(setting=paste0(setting,"_sup")))

final_output=merge_calibration_simulation(calibration_reformat=calibrated_data_uncertainty, my_all_simul=all_simul_chirps )

# take the envelope per setting
final_output$seed_avg= final_output$simul_calib %>%
  mutate(sup=str_detect( setting, pattern = ("_sup")),
         inf=str_detect( setting, pattern = ("_inf")),
         setting=gsub("_inf", "", gsub("_sup", "", setting)))%>%
  group_by(setting, year, month, age )%>%
  summarise(PR_mini=min(PR_mini), PR_maxi=max(PR_maxi), PR_middle=mean(PR_middle[sup==FALSE & inf==FALSE]))

# move outlier into validation set
fitdat_plot=fitdat%>% full_join(validat %>% filter(setting=="Pyrethroid" & PR_obs<0.8))
validat_plot=validat%>% filter(setting!="Pyrethroid" | PR_obs>0.8)

# take the envelope per net type
seed_avg_pooledModels= final_output$simul_calib %>%
  mutate(sup=str_detect( setting, pattern = ("_sup")),
         inf=str_detect( setting, pattern = ("_inf")),
         setting=gsub("_inf", "", gsub("_sup", "", gsub("_BIT055", "", gsub("_BIT103", "", gsub("_BIT059", "", setting))))))%>%
  group_by(setting, year, month, age)%>%
  summarise(PR_mini=min(PR_mini), PR_maxi=max(PR_maxi), PR_middle=mean(PR_middle[sup==FALSE & inf==FALSE]))

fitdat_pooledAll= fitdat %>%
  mutate(setting=gsub("_briet", "", gsub("_denz", "", gsub("_BIT055", "", gsub("_BIT103", "", gsub("_BIT059", "", setting))))))%>%
  unique()

validat_pooledAll= validat %>%
  mutate(setting=gsub("_briet", "", gsub("_denz", "", gsub("_BIT055", "", gsub("_BIT103", "", gsub("_BIT059", "", setting))))))%>%
  unique()

fitdat_pooledAll=fitdat_pooledAll%>% full_join(validat_pooledAll %>% filter(setting=="Pyrethroid" & PR_obs<0.8))
validat_pooledAll=validat_pooledAll%>% filter(setting!="Pyrethroid" | PR_obs>0.8)


write.csv(fitdat_pooledAll %>% select(-arm)%>% unique(), file = file.path(experiment_folder, "fitdat_protopopoff.csv"), row.names = F)
write.csv(validat_pooledAll, file = file.path(experiment_folder, "validat_protopopoff.csv"), row.names = F)
write.csv(seed_avg_pooledModels, file = file.path(experiment_folder, "results_protopopoff.csv"), row.names = F)
write.csv(final_output$seed_avg, file = file.path(experiment_folder, "results_protopopoff_detail.csv"), row.names = F)

