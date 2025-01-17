rm(list=ls())

library(devtools)
#devtools::install_github("SwissTPH/r-openMalariaUtilities", ref = "master")
#devtools::install("/scicore/home/pothin/champa0000/omuaddons") #run this if you have pulled in your local directory
#devtools::install("/scicore/home/pothin/champa0000/omu-slurm") #run this if you have pulled in your local directory
#devtools::install("/scicore/home/pothin/champa0000/omu-compat") #run this if you have pulled in your local directory

library(openMalariaUtilities)
library(tidyverse)
library(OMAddons)
library(omuslurm)
library(omucompat)
#library(caliviz)

source("helpfunctions_rct.R")
source("helpfunction_omu_gvi.R")
main_folder=(".")

# SET EXPERIMENT ----------------------------------------------------------

COUNTRY     = "TZA" ## options are: Tza, Ben, Moz, Cmr
SIMSTART    = "1918-01-01"
versionnum  = 44L

expName  = "RCTMosha_calibration"
experiment_folder=file.path(main_folder, expName)


###
### Define vector bionomics
###

# which mosquitoes will be used ?
mosqs  = c("funestus_indoor","gambiaesl_indoor")
#' Mosha 2022 https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(21)02499-5/fulltext#sec1
#' page 1235

# contrib: the relative contribution of each 'mosquito' to EIR
contrib = c("@fin@", "@gin@")
cbind( mosqs, contrib)

###----------------------------------------------------
###--------begin base file-----------------------------
###----------------------------------------------------

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

## Create demoography
baseList <- write_demography_compat(
  baseList = baseList, pop = "@pop@", country = COUNTRY
)


## Create monitoring
baseList <- write_monitoring_compat(
  baseList = baseList, y1 = 2015, y2 = 2022, detect = 100, SIMSTART = SIMSTART,
  upperbounds = c(.5,10, 14,100)
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


histITN= create_vectorInterventionParameters(deterrency = 0,
                                             preprandial = "@preprandial@",
                                             postprandial = 0,
                                             kappa=1.6,
                                             L=1, myname = "histITN", decay = "step", mosqs = mosqs)

baseList= defineGVI_simple(baseList = baseList,
                           vectorInterventionParameters=histITN,
                           append = TRUE,
                           verbatim = TRUE,
                           hist = TRUE)

futITN= create_vectorInterventionParameters(deterrency = 0,
                                            preprandial = "@preprandial@",
                                            postprandial = 0,
                                            kappa=1.6,
                                            L=2.3, myname = "futITN", mosqs = mosqs)

baseList= defineGVI_simple(baseList = baseList,
                           vectorInterventionParameters=futITN,
                           append = TRUE,
                           verbatim = TRUE,
                           hist = FALSE)

####--------------------####
#### deployment section ####
####--------------------####
#- successive distributions after start of trial

## Historical ITN
baseList <- deploy_it_compat(
  baseList = baseList, component = "histITN", coverage = "@histITNcov@",
  ## Allowing for different hist coverage levels
  byyear = FALSE,
  ## Annual deployments
  y1 = 2015, y2 = 2018, every = 1, interval = "year",
  m1 = 1, m2 = 1, d1 = 25, d2 = 25, #to end just when the trial distribution starts
  SIMSTART = SIMSTART
)


## Future ITN
baseList <- deploy_it_compat(
  baseList = baseList, component = "futITN",
  coverage = "@futITNcov@",
  ## allowing for different hist coverage levels
  byyear = FALSE,
  ##  annual deployments
  y1 = 2019, y2 = 2019, every = 3, interval = "year",
  m1 = 1, m2 = 1, d1 = 27, d2 = 27,
  #' Mosha 2022 Abstract, Findings
  #' https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(21)02499-5/fulltext#sec1
  #'"LLINs were distributed among households between Jan 26 and Jan 28, 2019"
  SIMSTART = SIMSTART
)



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


#######################
### seasonality
season_dir="/scicore/home/pothin/champa0000/IG2modelling/"
seasonality_chirps=read.csv(file.path(season_dir, "seasonality_Misungwi.csv"))

seasonality_full=rbind(seasonality_chirps %>% mutate(season="chirps"))

seasonality_full_lag2m=seasonality_full
names(seasonality_full_lag2m)=c("m3", "m4", "m5", "m6","m7","m8","m9" ,"m10", "m11", "m12", "m1", "m2", "season")

seasonality_full_lag3m=seasonality_full
names(seasonality_full_lag3m)=c("m4", "m5", "m6","m7","m8","m9" ,"m10", "m11", "m12", "m1", "m2", "m3", "season")

seasonality_full_all=rbind(seasonality_full%>% mutate(seasonLag="1m"),
                           seasonality_full_lag2m%>% mutate(seasonLag="2m"),
                           seasonality_full_lag3m%>% mutate(seasonLag="3m"))
##-- creating the experiment
full             = list()
full $ seed      = 1:10
full $ setting   = c("pyrethroid")
full $ pop       = 10000
##' DHS 2017 for Mwanza province, effective coverage value
##' code by Katya and Rémi https://clintonhealth.box.com/s/02vk4o11ta7etymks0du8hbzsb925n69
full$ EffCovconv = convert_cm(0.6155460)
full $ EIR = seq(10,200,2)
full $ preprandial = seq(0,1,0.05)

full $ season = c("chirps")
full $ seasonLag = c("1m", "2m", "3m")
#' https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(21)02499-5/fulltext#sec1
#' Mosha 2022 page 1235
full$fin = 4700/4973*1 #90% indoor biting
full$gin = (4973-4700)/4973*1 #90% indoor biting

#### 'scens' will do all possible combinations of those factors
scens = expand.grid( full )
scens = scens %>% left_join( seasonality_full_all)

##' historical net coverage: https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(21)02499-5/fulltext#sec1
##' Table 1, Long-lasting insecticidal net use In household residents of all age groups
##' future coverage: computed above

NetCov <- data.frame(setting=full$setting,
                     histITNcov=c(2957/4962),
                     futITNcov=c(0.77)
)

scens <- left_join(scens,NetCov)

scens <- finalize_scenarios(scens)
storeScenarios(scens)

dim(scens)

length(unique(scens$setting))

## SLURM
setup_scenarios(scenarios=scens[1,])
slurmPrepareScenarios(expName = expName, scenarios = scens, nCPU = 8, qos = "30min", time = "00:10:00", 
                      rModule = "R/4.4.1-gfbf-2023b")
#slurmCreateScenarios()
check_scenarios_created(experiment_folder)


validateXML(xmlfile = file.path(experiment_folder, paste0("scenarios/", expName, "_7.xml")), schema = NULL, scenarios = scens)


## Run Open Malaria
## SLURM
## Run Open Malaria
## SLURM
slurmPrepareSimulations(expName = expName, scenarios = scens, nCPU = 1, bSize = 1,
                        rModule = "R/4.4.1-gfbf-2023b", 
                        omModule = "OpenMalaria/44.0-intel-compilers-2023.1.0")
#slurmRunSimulation()
check_simulations_created(experiment_folder)


# cluster_status(experiment_folder,  stime="15:10" )
#
# failed=read.table(file.path(experiment_folder, "missing_om_simulations.txt"))
# failed$V2=gsub("_out.txt","", gsub(paste0(expName, "_"), "", failed$V1))
# #
# paste0(failed$V2,collapse=",")

## Post-processing
age_groups_of_interest <- c("0.5-14")
## 6 months to 14 years old for prevalence
## 6 months to 10 years old for incidence

loadExperiment(expName)
## SLURM
slurmPrepareResults(

  expDir=experiment_folder,
  dbName="results_RCTMosha_controlArm",
  dbDir = experiment_folder,
  resultsName = "results",
  resultsCols = c("scenario_id", "date_aggregation", "date", "age_group",
                  "nUncomp", "nHost", "nSevere",
                  "incidenceRate", "prevalenceRate"),
  indexOn = list(c("results", "scenario_id")),
  ncoresDT = 1,
  strategy = "batch",
  appendResults = FALSE,
  aggrFun = CalcEpiOutputs,
  aggrFunArgs = list(aggregateByAgeGroup = age_groups_of_interest, aggregateByDate = "month"),
  splitBy = 20000,
  verbose = FALSE,
  ntasks = 1,
  mem = "64G",
  nCPU = 1,
  time = "06:00:00",
  qos = "6hours",
  rModule = "R/4.4.1-gfbf-2023b"
)
#

#slurmRunResults()

####
#### Step 4: Fitting to historical data
####

dir.create(file.path(experiment_folder, "calibration/"))
FittingDir=file.path(experiment_folder, "calibration/")


loadExperiment(expName)
scens=readRDS(file.path(experiment_folder, "cache/scenarios.rds"))

# PREPARE THE DATA TO FIT

# open database connection
conn = DBI::dbConnect(RSQLite::SQLite(),paste0(expName,"/results_RCTMosha_controlArm.sqlite"))
# list the tables of the database
DBI::dbListTables(conn)
# extract the results table, this should have the same name as provided to slurmPrepareResults
all_simul_0 = DBI::dbReadTable(conn, "results")
# always disconnect from the database
DBI::dbDisconnect(conn)

all_simul=merge(all_simul_0, scens %>%mutate(scenario_id=ID))%>%
  mutate(EIR=as.numeric(EIR),
         month=as.numeric(format(as.Date(date), "%m")),
         year=as.numeric(format(as.Date(date), "%y"))+2000)

all_simul_chirps=all_simul

# write.csv(all_simul_chirps, file.path(FittingDir, "all_simul_chirps.csv"), row.names = F)
#all_simul_chirps=read.csv(file.path(FittingDir, "all_simul_chirps.csv"))

agePR="0.5-14"

###############################
## fitdat
######################


#' Mosha 2022 Table 1
#' #'https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(21)02499-5/fulltext#sec1
fitdat = data.frame(setting=c("pyrethroid","pyriproxyfen","chlorfenapyr","PBO")
                    ,age=agePR
                    ,year=2018
                    ,month=9
                    ,pos_children=c(519,516,469,444)
                    ,nb_children=c(1130,1118,1099,1056)) %>%
  rowwise() %>%
  mutate(PR_obs=pos_children/nb_children,
         sd=sqrt(PR_obs*(1-PR_obs)/nb_children))%>%
  filter(setting %in% c("pyrethroid"))


CI_fitdat = as.data.frame(t(mapply(compute_CI, fitdat$pos_children/fitdat$nb_children, fitdat$nb_children)))

fitdat$PR_obs_min = unlist(CI_fitdat$LCI)
fitdat$PR_obs_max = unlist(CI_fitdat$UCI)

######################
## validat
######################
#' Mosha 2022 Table 2
#' #'https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(21)02499-5/fulltext#sec1
validat = data.frame(setting=rep(c("pyrethroid","pyriproxyfen","chlorfenapyr","PBO"),each=5)
                     ,age=agePR
                     ,year=rep(c(2020,2020,2021, 2021, 2022),4)
                     ,month=rep(c(1,8,1, 8, 2),4)
                     ,pos_children=c(350,642,549,507,407,
                                     232,583,472,426,302,
                                     176,509,326,436,261,
                                     206,502,512,488, 338)
                     ,nb_children=c(1123,1227,1199,956,1088,
                                    1069,1153,1258,1004,1050,
                                    1126,1246,1272,1045,1145,
                                    1071,1160,1259,992,1049)) %>%
  rowwise() %>%
  mutate(PR_obs=pos_children/nb_children,
         sd=sqrt(PR_obs*(1-PR_obs)/nb_children))%>%
  filter(setting %in% c("pyrethroid"))
CI_validat = as.data.frame(t(mapply(compute_CI, validat$pos_children/validat$nb_children, fitdat$nb_children)))

validat$PR_obs_min = unlist(CI_validat$LCI)
validat$PR_obs_max = unlist(CI_validat$UCI)

fitdat_all=rbind(fitdat, validat)
fitdat_all=fitdat_all[fitdat_all$setting=="pyrethroid",]



fitdat_all=rbind(fitdat, validat)
fitdat_all=fitdat_all[fitdat_all$setting=="pyrethroid",]
##############################################
# select best post prandial and best seasonality

simul_fit=all_simul_chirps %>% filter(age_group==agePR) %>%
  right_join(fitdat_all)  %>%
  dplyr::mutate(
    SE=(prevalenceRate-PR_obs)^2,
    llk_gaussian = stats::dnorm(PR_obs, mean = prevalenceRate,  sd = sd, log = TRUE),
    SE_2018=ifelse(year==2018, SE, 0),
    llk_gaussian_2018=ifelse(year==2018, llk_gaussian, 0))

simul_fit_seedavg=simul_fit%>%
  group_by(preprandial, EIR, seasonLag, seed) %>%
  summarise(llk=sum(llk_gaussian), llk_2019=sum(llk_gaussian_2018), SE=sum(SE), SE_2018=sum(SE_2018), n=n())%>%
  group_by(preprandial, EIR, seasonLag) %>%
  summarise(llk=mean(llk), llk_2018=mean(llk_2019), SE=mean(SE), SE_2018=mean(SE_2018), n=mean(n))

# best EIR for all
simul_fit_seedavg %>%
  ungroup()%>%
  group_by(seasonLag) %>%
  dplyr::slice(which.max(llk))

# best EIR for all
simul_fit_seedavg %>%
  ungroup()%>%
  dplyr::slice(which.max(llk))



# best EIR for all
simul_fit_seedavg %>%
  ungroup()%>%
  filter(preprandial>0)%>%
  dplyr::slice(which.max(llk))
