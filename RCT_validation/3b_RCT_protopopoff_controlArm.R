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
COUNTRY <- "TZA"
## Start date of monitoring. OM's recommendation is to start 20 to 50 years
## before the first measures are deployed.
SIMSTART <- "1918-01-01"
## openMalariaUtilities supports Open Malaria version 43 and up
versionnum <- 44L


source("helpfunctions_rct.R")
source("helpfunction_omu_gvi.R")
main_folder=(".")



expName <- "RCT_protopopoff_calibration"
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


histITN= create_vectorInterventionParameters(deterrency = 0,
                                                preprandial ="@preprandial@",
                                                postprandial = 0,
                                                kappa=3.1,
                                                L=1, myname = "histITN", decay = "step", mosqs = mosqs)

baseList= defineGVI_simple(baseList = baseList,
                           vectorInterventionParameters=histITN,
                           append = TRUE,
                           verbatim = TRUE,
                           hist = TRUE)

futITN= create_vectorInterventionParameters(deterrency = 0,
                                               preprandial =  "@preprandial@",
                                               postprandial = 0,
                                               kappa=3.1,
                                               L=2.7, myname = "futITN", mosqs = mosqs)

baseList= defineGVI_simple(baseList = baseList,
                           vectorInterventionParameters=futITN,
                           append = TRUE,
                           verbatim = TRUE,
                           hist = FALSE)


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
#c("pyrethroid","PBO_BIT055_briet", "PBO_BIT055_denz")
full $ pop       = 10000
full$ EffCovconv = convert_cm(0.3782992) #DHS 2016, effective coverage, code by Katya an Rémi
full $ EIR = seq(10,200,2)
#full $ preprandial = seq(0,1,0.05)
full $ preprandial = seq(0,1,0.05)
# species composition: personal communication
full$ain = .04 
#full$aout = .04*.1 #90% indoor biting
full$fin = .04 
#full$fout = .04*.1 #90% indoor biting
full$gin = .92 
#full$gout = .92*.1 #90% indoor biting
full$season = c("chirps")
full $ seasonLag = c("1m", "2m", "3m")

#### 'scens' will do all possible combinations of those factors
scens = expand.grid( full )
scens = scens %>% left_join(seasonality_full_all) #%>%

#### add vector control coverages
NetCov <- data.frame(setting=c("pyrethroid"),
                     histITNcov=c(902/2996),
                     futITNcov=c(0.77)
)
scens <- left_join(scens,NetCov)

scens <- data.frame(scens)



scens <- finalize_scenarios(scens)
storeScenarios(scens)

dim(scens)

length(unique(scens$setting))

# small check on the scenarios, to check if the scenario has the right placeholder names
# setup_scenarios(scenarios=scens[1,])

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
# bSize: nb of simulations per job, and they run in parallel on the nCPUs
# nCPUs if too high, hard to get the resources on scicore
# old behaviour: bSize=1 and nCPU=1
#

#slurmRunSimulation()
check_simulations_created(experiment_folder)
#
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
  dbName="results_protopopoff_calibration_chirps",
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
  splitBy = 80000,
  verbose = FALSE,
  ntasks = 1,
  mem = "64G",
  nCPU = 10,
  time = "06:00:00",
  qos = "6hours",
  rModule = "R/4.4.1-gfbf-2023b"
)
#slurmRunResults()

####
#### Step 4: Fitting to historical data
####



loadExperiment(expName)
scens=readRDS(file.path(experiment_folder, "cache/scenarios.rds"))

####Fitting with Clara's method
dir.create(file.path(experiment_folder, "calibration/"))
FittingDir=file.path(experiment_folder, "calibration/")

# PREPARE THE DATA TO FIT

# open database connection
conn = DBI::dbConnect(RSQLite::SQLite(),file.path(experimentDir,"results_protopopoff_calibration_chirps.sqlite"))
# list the tables of the database
DBI::dbListTables(conn)
# extract the results table, this should have the same name as provided to slurmPrepareResults
all_simul_0 = DBI::dbReadTable(conn, "results_RCT13")
# always disconnect from the database
DBI::dbDisconnect(conn)

all_simul=merge(all_simul_0  %>% filter(experiment_id==1), scens %>%mutate(scenario_id=ID), by="scenario_id")%>%
  mutate(EIR=as.numeric(EIR),
         month=as.numeric(format(as.Date(date), "%m")),
         year=as.numeric(format(as.Date(date), "%y"))+2000)


####
#### Step 4: Fitting to historical data
####

# PREPARE THE DATA TO FIT
all_simul_chirps=all_simul

#write.csv(all_simul_chirps, file.path(FittingDir, "all_simul_chirps.csv"), row.names = F)

all_simul_chirps=read.csv(file.path(FittingDir, "all_simul_chirps.csv"))

agePR="0.5-14"

###############################
## fitdat
######################


#'Table 1 Protopopoff 2018 https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(18)30427-6/fulltext
#'timing: Figure 1
fitdat = data.frame(setting=c("pyrethroid", "PBO_BIT055_denz","PBO_BIT055_briet")
                    ,age="0.5-14"
                    ,year=2014
                    ,month=9
                    ,pos_children=c(600,606, 606)
                    ,nb_children=c(885,991, 991)) %>%
  rowwise() %>%
  mutate(PR_obs=pos_children/nb_children,
         sd=sqrt(PR_obs*(1-PR_obs)/nb_children),
         #sd=4*sqrt(PR_obs*(1-PR_obs)/nb_children),
         PR_obs_min=PR_obs-1.96*sd,
         PR_obs_max=PR_obs+1.96*sd) %>%
  filter(setting %in% c("pyrethroid"))


#CI_fitdat = as.data.frame(t(mapply(compute_CI, fitdat$pos_children/fitdat$nb_children, fitdat$nb_children)))

#fitdat$PR_obs_min = unlist(CI_fitdat$LCI)
#fitdat$PR_obs_max = unlist(CI_fitdat$UCI)


######################
## validat
######################

setting=c("pyrethroid", "PBO_BIT055_denz","PBO_BIT055_briet")
#' Protopopoff et al. 2018 https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(18)30427-6/fulltext
#' Figure 1 for timing and Table 2 for numbers
validat = data.frame(setting=rep(c("pyrethroid", "PBO_BIT055_denz","PBO_BIT055_briet"),each=6)
                     ,age="0.5-14"
                     ,year=rep(c(2015,2015,2016,2016, 2017, 2017),3)
                     ,month=rep(c(6,11,6,11, 6, 11),3)
                     ,pos_children=c(553,515,548,710, 1643, 1045,
                                     445,275,342,440,1280,901,
                                     445,275,342,440,1280,901
                     )
                     ,nb_children=c(997,932,1034,1044,2032,1784,
                                    971,883,984,958,1848,1807,
                                    971,883,984,958,1848,1807
                     )
) %>%
  rowwise() %>%
  mutate(PR_obs=pos_children/nb_children,
         sd=sqrt(PR_obs*(1-PR_obs)/nb_children),
         PR_obs_min=PR_obs-1.96*sd,
         PR_obs_max=PR_obs+1.96*sd)%>%
  filter(setting %in% c("pyrethroid"))

fitdat_all=rbind(fitdat, validat)
fitdat_all=fitdat_all[fitdat_all$setting=="pyrethroid",]%>%
  filter(PR_obs<0.8)

##########################################################################
# FIND THE BEST PREPRANDIAL KILLING EFFECT, and THE BEST SEASONALITY LAG

simul_fit=all_simul_chirps %>% filter(age_group==agePR) %>%
  right_join(fitdat_all)  %>%
  dplyr::mutate(
    SE=(prevalenceRate-PR_obs)^2,
    llk_gaussian = stats::dnorm(PR_obs, mean = prevalenceRate,  sd = sd, log = TRUE))

simul_fit_seedavg=simul_fit%>%
  group_by(preprandial, EIR, seasonLag, seed) %>%
  summarise(llk=sum(llk_gaussian), SE=sum(SE), n=n())%>%
  group_by(preprandial, EIR, seasonLag) %>%
  summarise(llk=mean(llk), SE=mean(SE), n=mean(n))

# best EIR for all
simul_fit_seedavg %>%
  ungroup()%>%
  group_by(seasonLag) %>%
  dplyr::slice(which.max(llk))


# best EIR for all
simul_fit_seedavg %>%
  ungroup()%>%
  dplyr::slice(which.max(llk))
