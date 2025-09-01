library(devtools)
#devtools::install_github("SwissTPH/r-openMalariaUtilities",
#                         ref = "master")

library(openMalariaUtilities)
library(tidyverse)
library(ggplot2)
#library(omuslurm)
library(omucompat)

rm(list=ls())

## Demographic data for the following countries included:
COUNTRY <- "TZA"
## Start date of monitoring. OM's recommendation is to start 20 to 50 years
## before the first measures are deployed.
SIMSTART <- "1918-01-01"
## openMalariaUtilities supports Open Malaria version 43 and up
versionnum <- 44L


main_folder=(".")
script_folder=""
parameters_GVI=read.csv(file.path(script_folder, "EHT_fit/fitted_parameters_posteriormax.csv"))
source(file.path(script_folder,"RCT_validation/helpfunctions_rct.R"))
source(file.path(script_folder,"RCT_validation/helpfunction_omu_gvi.R"))
source(file.path(script_folder,"helpfunctions_convert_cm.R"))



expName <- "RCT_protopopoff_full"
setwd(main_folder)
experimentDir <- paste0(main_folder,expName)
experiment_folder=file.path(main_folder, expName)

## Define mosquitoes
mosqs <- c("arabiensis_indoor","funestus_indoor", "gambiae_indoor")#personal communication
contrib <- c("@ain@","@fin@", "@gin@" )
cbind( mosqs, contrib)

mosqs_exposures=c("arabiensis_indoor"=1, "funestus_indoor"=1, "gambiae_indoor"=1)
exposure_in_bed=0.69


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
kappa_PBO=1.8

halflife_pyrethroid=2.7
kappa_pyrethroid=3.0

# historical ITNs
hist_pyrethroid_BIT103= extract_GVI_params(EHT="Assenga",
                                           netType="Pyrethroid-only",
                                           halflife_functionalSurvival=1, kappa_functionalSurvival=kappa_pyrethroid,
                                           myname="histITN_BIT103", mosqs=mosqs,  parameters_GVI=parameters_GVI, mosqs_exposures=mosqs_exposures, insecticide_decay = F, decay = "step")

hist_pyrethroid_martin= extract_GVI_params(EHT="Martin",
                                           netType="Pyrethroid-only",
                                           halflife_functionalSurvival=1, kappa_functionalSurvival=kappa_pyrethroid,
                                           myname="histITN_martin", mosqs=mosqs,  parameters_GVI=parameters_GVI, mosqs_exposures=mosqs_exposures, insecticide_decay = F, decay = "step")

hist_pyrethroid_kibondo= extract_GVI_params(EHT="Kibondo",
                                            netType="Pyrethroid-only",
                                            halflife_functionalSurvival=1, kappa_functionalSurvival=kappa_pyrethroid,
                                            myname="histITN_kibondo", mosqs=mosqs,  parameters_GVI=parameters_GVI, mosqs_exposures=mosqs_exposures, insecticide_decay = F, decay = "step")


hist_pyrethroid_BIT055= extract_GVI_params(EHT="BIT055",
                                           netType="Pyrethroid-only",
                                           halflife_functionalSurvival=1, kappa_functionalSurvival=kappa_pyrethroid,
                                           myname="histITN_BIT055", mosqs=mosqs,  parameters_GVI=parameters_GVI, mosqs_exposures=mosqs_exposures, insecticide_decay = F, decay = "step")

hist_pyrethroid_odufuwa= extract_GVI_params(EHT="Odufuwa",
                                            netType="Pyrethroid-only",
                                            halflife_functionalSurvival=1, kappa_functionalSurvival=kappa_pyrethroid,
                                            myname="histITN_odufuwa", mosqs=mosqs,  parameters_GVI=parameters_GVI, mosqs_exposures=mosqs_exposures, insecticide_decay = F, decay = "step")


hist_pyrethroid_BIT080= extract_GVI_params(EHT="BIT080",
                                           netType="Pyrethroid-only",
                                           halflife_functionalSurvival=1, kappa_functionalSurvival=kappa_pyrethroid,
                                           myname="histITN_BIT080", mosqs=mosqs,  parameters_GVI=parameters_GVI, mosqs_exposures=mosqs_exposures, insecticide_decay = F, decay = "step")

hist_pyrethroid_martinf= extract_GVI_params(EHT="Martin, f",
                                            netType="Pyrethroid-only",
                                            halflife_functionalSurvival=1, kappa_functionalSurvival=kappa_pyrethroid,
                                            myname="histITN_martinf", mosqs=mosqs,  parameters_GVI=parameters_GVI, mosqs_exposures=mosqs_exposures, insecticide_decay = F, decay = "step")


# Future ITNs
GVI_pyrethroid_BIT103= extract_GVI_params(EHT="Assenga",
                                          netType="Pyrethroid-only",
                                          halflife_functionalSurvival=halflife_pyrethroid, kappa_functionalSurvival=kappa_pyrethroid,
                                          myname="futITN_BIT103", mosqs=mosqs,  parameters_GVI=parameters_GVI, mosqs_exposures=mosqs_exposures)

GVI_pyrethroid_martin= extract_GVI_params(EHT="Martin",
                                          netType="Pyrethroid-only",
                                          halflife_functionalSurvival=halflife_pyrethroid, kappa_functionalSurvival=kappa_pyrethroid,
                                          myname="futITN_martin", mosqs=mosqs,  parameters_GVI=parameters_GVI, mosqs_exposures=mosqs_exposures)

GVI_pyrethroid_kibondo= extract_GVI_params(EHT="Kibondo",
                                           netType="Pyrethroid-only",
                                           halflife_functionalSurvival=halflife_pyrethroid, kappa_functionalSurvival=kappa_pyrethroid,
                                           myname="futITN_kibondo", mosqs=mosqs,  parameters_GVI=parameters_GVI, mosqs_exposures=mosqs_exposures)

GVI_pyrethroid_BIT055= extract_GVI_params(EHT="BIT055",
                                          netType="Pyrethroid-only",
                                          halflife_functionalSurvival=halflife_pyrethroid, kappa_functionalSurvival=kappa_pyrethroid,
                                          myname="futITN_BIT055", mosqs=mosqs,  parameters_GVI=parameters_GVI, mosqs_exposures=mosqs_exposures)

GVI_pyrethroid_odufuwa= extract_GVI_params(EHT="Odufuwa",
                                           netType="Pyrethroid-only",
                                           halflife_functionalSurvival=halflife_pyrethroid, kappa_functionalSurvival=kappa_pyrethroid,
                                           myname="futITN_odufuwa", mosqs=mosqs,  parameters_GVI=parameters_GVI, mosqs_exposures=mosqs_exposures)

GVI_pyrethroid_BIT080= extract_GVI_params(EHT="BIT080",
                                          netType="Pyrethroid-only",
                                          halflife_functionalSurvival=halflife_pyrethroid, kappa_functionalSurvival=kappa_pyrethroid,
                                          myname="futITN_BIT080", mosqs=mosqs,  parameters_GVI=parameters_GVI, mosqs_exposures=mosqs_exposures)

GVI_pyrethroid_martinf= extract_GVI_params(EHT="Martin, f",
                                           netType="Pyrethroid-only",
                                           halflife_functionalSurvival=halflife_pyrethroid, kappa_functionalSurvival=kappa_pyrethroid,
                                           myname="futITN_martinf", mosqs=mosqs,  parameters_GVI=parameters_GVI, mosqs_exposures=mosqs_exposures)

# future PBOs

GVI_PBO_BIT055= extract_GVI_params(EHT="BIT055",
                                   netType="Olyset Plus",
                                   halflife_functionalSurvival=halflife_PBO, kappa_functionalSurvival=kappa_PBO,
                                   myname="futPBO_BIT055", mosqs=mosqs, parameters_GVI=parameters_GVI, mosqs_exposures=mosqs_exposures)

GVI_PBO_BIT059= extract_GVI_params(EHT="Odufuwa",
                                   netType="Olyset Plus",
                                   halflife_functionalSurvival=halflife_PBO, kappa_functionalSurvival=kappa_PBO,
                                   myname="futPBO_BIT059", mosqs=mosqs, parameters_GVI=parameters_GVI, mosqs_exposures=mosqs_exposures)

GVI_PBO_BIT103= extract_GVI_params(EHT="Assenga",
                                   netType="Olyset Plus",
                                   halflife_functionalSurvival=halflife_PBO, kappa_functionalSurvival=kappa_PBO,
                                   myname="futPBO_BIT103", mosqs=mosqs,  parameters_GVI=parameters_GVI, mosqs_exposures=mosqs_exposures)

GVI_PBO_martin= extract_GVI_params(EHT="Martin",
                                   netType="Olyset Plus",
                                   halflife_functionalSurvival=halflife_PBO, kappa_functionalSurvival=kappa_PBO,
                                   myname="futPBO_martin", mosqs=mosqs,  parameters_GVI=parameters_GVI, mosqs_exposures=mosqs_exposures)

GVI_PBO_martinf= extract_GVI_params(EHT="Martin, f",
                                    netType="Olyset Plus",
                                    halflife_functionalSurvival=halflife_PBO, kappa_functionalSurvival=kappa_PBO,
                                    myname="futPBO_martinf", mosqs=mosqs,  parameters_GVI=parameters_GVI, mosqs_exposures=mosqs_exposures)


list_GVI_snippets=list(hist_pyrethroid_BIT103,hist_pyrethroid_martin,hist_pyrethroid_kibondo,hist_pyrethroid_BIT055,hist_pyrethroid_odufuwa,hist_pyrethroid_BIT080,hist_pyrethroid_martinf)
list_GVI_snippets=append(list_GVI_snippets, GVI_pyrethroid_BIT103)
list_GVI_snippets=append(list_GVI_snippets, GVI_pyrethroid_martin)
list_GVI_snippets=append(list_GVI_snippets, GVI_pyrethroid_kibondo)
list_GVI_snippets=append(list_GVI_snippets, GVI_pyrethroid_BIT055)
list_GVI_snippets=append(list_GVI_snippets, GVI_pyrethroid_odufuwa)
list_GVI_snippets=append(list_GVI_snippets, GVI_pyrethroid_BIT080)
list_GVI_snippets=append(list_GVI_snippets, GVI_pyrethroid_martinf)

list_GVI_snippets=append(list_GVI_snippets, GVI_PBO_BIT103)
list_GVI_snippets=append(list_GVI_snippets, GVI_PBO_BIT055)
list_GVI_snippets=append(list_GVI_snippets, GVI_PBO_BIT059)
list_GVI_snippets=append(list_GVI_snippets, GVI_PBO_martin)
list_GVI_snippets=append(list_GVI_snippets, GVI_PBO_martinf)

for(GVI_param in list_GVI_snippets){
  baseList= defineGVI_simple(baseList = baseList,
                             vectorInterventionParameters=GVI_param,
                             append = TRUE,
                             verbatim = TRUE,
                             hist = FALSE)
}

## Deployment section


## historical nets

list_histITN_snippets_id=c("histITN_BIT103","histITN_martin","histITN_kibondo",
                           "histITN_BIT055","histITN_odufuwa","histITN_BIT080","histITN_martinf"
)


for(GVI_snippet_id in list_histITN_snippets_id){
  baseList <- deploy_it_compat(
    baseList = baseList, component = paste0(GVI_snippet_id),
    coverage = paste0("@", GVI_snippet_id,"_cov@"),
    ## allowing for different hist coverage levels
    byyear = FALSE,
    ##  annual deployments
    y1 = 2012, y2 = 2014, every = 1, interval = "year",
    m1 = 2, m2 = 2, d1 = 25, d2 = 25,
    #' Mosha 2022 Abstract, Findings
    #' https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(21)02499-5/fulltext#sec1
    #'"LLINs were distributed among households between Jan 26 and Jan 28, 2019"
    SIMSTART = SIMSTART
  )
  
}


## Future IG2
list_GVI_snippets_id=c("futPBO_BIT055",
                       "futPBO_BIT059",
                       "futPBO_BIT103",
                       "futPBO_martin",
                       "futPBO_martinf",
                       "futITN_martin",
                       "futITN_BIT103",
                       "futITN_kibondo",
                       "futITN_BIT055",
                       "futITN_odufuwa",
                       "futITN_BIT080",
                       "futITN_martinf"
)


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
createBaseXml(baseList, replace = TRUE)


## Copy necessary Open Malaria files. Currently needs to be called after
## createBaseXml, otherwise the cache is not set up.
setupOM(version = versionnum, dir=experiment_folder)


#########################################
### seasonality
season_dir="/scicore/home/pothin/champa0000/IG2modelling/"
seasonality_chirps=read.csv(file.path(season_dir, "seasonality_Muleba.csv"))

seasonality_full=rbind(seasonality_chirps %>% mutate(season="chirps"))

seasonality_full_lag1m=rbind(seasonality_chirps %>% mutate(season="chirps"))

seasonality_full_lag2m=seasonality_full
names(seasonality_full_lag2m)=c("m3", "m4", "m5", "m6","m7","m8","m9" ,"m10", "m11", "m12", "m1", "m2", "season")

seasonality_full_lag3m=seasonality_full
names(seasonality_full_lag3m)=c("m4", "m5", "m6","m7","m8","m9" ,"m10", "m11", "m12", "m1", "m2", "m3", "season")

seasonality_full_all=rbind(seasonality_full_lag1m%>% mutate(seasonLag="1m"),
                           seasonality_full_lag2m%>% mutate(seasonLag="2m"),
                           seasonality_full_lag3m%>% mutate(seasonLag="3m"))

seasonality_full=seasonality_full_lag3m%>% mutate(seasonLag="3m")
seasonality_full %>%
  pivot_longer(cols=starts_with("m"),names_to = "month",values_to="intensity",names_prefix = "m") %>%
  mutate(month=as.numeric(month), dataset="CHIRPS +1m") %>%
  ggplot()+
  geom_line(aes(x=month,y=intensity, color=dataset))+
  theme_minimal()+
  scale_color_manual(values=c("red", "orange", "dodgerblue"), name="")


##-- creating the experiment
full             = list()
full $ seed      = 1:10
full $ setting   = c("pyrethroid",
                     "PBO_BIT055",
                     "PBO_BIT059",
                     "PBO_BIT103",
                     "PBO_martin",
                     "PBO_martinf"
)
full$StandardITN=c("BIT103","martin", "kibondo", "BIT055", "odufuwa", "BIT080", "martinf")
full $ pop       = 10000
full$ EffCovconv = convert_cm(0.3782992) #DHS 2016, effective coverage, code by Katya an Rémi
full $ EIR =c(1:9,seq(10,200,2))
# species composition: personal communication
full$ain = .04 
full$fin = .04 
full$gin = .92 
full$season = c("chirps")
full $ seasonLag = c("3m")

#### 'scens' will do all possible combinations of those factors
scens = expand.grid( full )
scens = scens %>% left_join(seasonality_full) #%>%



futPBO_cov=0.80*exposure_in_bed
NetCov <- data.frame(setting=c("pyrethroid","PBO"),
                     histITNcov=c(902/2996,810/3078)*exposure_in_bed,
                     futITNcov=c(0.77 , 0)*exposure_in_bed
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


PyrNetCov <- scens %>% select(setting, StandardITN)%>% unique %>%
  left_join(NetCov_full)
PyrNetCov_full=PyrNetCov
for(pyr_net in full$StandardITN){
  myPyrNetCov=PyrNetCov %>%
    filter(StandardITN ==pyr_net)
  names(myPyrNetCov)[3]=paste0("histITN_", pyr_net, "_cov")
  names(myPyrNetCov)[4]=paste0("futITN_", pyr_net, "_cov")
  PyrNetCov_full=merge(PyrNetCov_full, myPyrNetCov, all=TRUE)
}
PyrNetCov_full[is.na(PyrNetCov_full)]=0



scens <- left_join(scens,PyrNetCov_full)

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
scens=readRDS(file.path(experiment_folder, "cache/scenarios.rds"))

# PREPARE THE DATA TO FIT

# open database connection
conn = DBI::dbConnect(RSQLite::SQLite(),file.path(experiment_folder,"results_RCTprotopopoff_validation.sqlite"))
# list the tables of the database
DBI::dbListTables(conn)
# extract the results table, this should have the same name as provided to slurmPrepareResults
all_simul_0 = DBI::dbReadTable(conn, "results_RCT13")
# always disconnect from the database
DBI::dbDisconnect(conn)

####
#### Step 4: Fitting to historical data
####

agePR="0.5-14"

###############################
## fitdat
######################


#'Table 1 Protopopoff 2018 https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(18)30427-6/fulltext
#'timing: Figure 1
fitdat = data.frame(setting=c("pyrethroid","PBO")
                    ,age="0.5-14"
                    ,year=2014
                    ,month=9
                    ,pos_children=c(600,606)
                    ,nb_children=c(885,991)) %>%
  rowwise() %>%
  mutate(PR_obs=pos_children/nb_children,
         sd=sqrt(PR_obs*(1-PR_obs)/nb_children),
         #sd=4*sqrt(PR_obs*(1-PR_obs)/nb_children),
         PR_obs_min=PR_obs-1.96*sd,
         PR_obs_max=PR_obs+1.96*sd) %>%
  separate(setting, c("trial","arm"),remove = F)

fitdat=rbind(fitdat,
             fitdat %>% filter(setting=="PBO") %>% mutate(setting="PBO_BIT055"),
             fitdat %>% filter(setting=="PBO") %>% mutate(setting="PBO_BIT103"),
             fitdat %>% filter(setting=="PBO") %>% mutate(setting="PBO_martin"),
             fitdat %>% filter(setting=="PBO") %>% mutate(setting="PBO_martinf"),
             fitdat %>% filter(setting=="PBO") %>% mutate(setting="PBO_BIT059")
)
fitdat=fitdat[fitdat$setting !="PBO",]

######################
## validat
######################

#' Protopopoff et al. 2018 https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(18)30427-6/fulltext
#' Figure 1 for timing and Table 2 for numbers
validat0 = data.frame(setting=rep(c("pyrethroid", "PBO"),each=6)
                      ,age="0.5-14"
                      ,year=rep(c(2015,2015,2016,2016, 2017, 2017),2)
                      ,month=rep(c(6,11,6,11, 6, 11),2)
                      ,pos_children=c(553,515,548,710, 1643, 1045,
                                      445,275,342,440,1280,901
                      )
                      ,nb_children=c(997,932,1034,1044,2032,1784,
                                     971,883,984,958,1848,1807
                      )
)%>%
  rowwise() %>%
  mutate(PR_obs=pos_children/nb_children,
         sd=sqrt(PR_obs*(1-PR_obs)/nb_children),
         PR_obs_min=PR_obs-1.96*sd,
         PR_obs_max=PR_obs+1.96*sd)


validat=rbind(validat0,
              validat0 %>% filter(setting=="PBO") %>% mutate(setting="PBO_BIT055"),
              validat0 %>% filter(setting=="PBO") %>% mutate(setting="PBO_BIT103"),
              validat0 %>% filter(setting=="PBO") %>% mutate(setting="PBO_BIT059"),
              validat0 %>% filter(setting=="PBO") %>% mutate(setting="PBO_martin"),
              validat0 %>% filter(setting=="PBO") %>% mutate(setting="PBO_martinf"))

validat=validat[validat$setting !="chlorfenapyr" & validat$setting !="PBO",]


fitdat_pooledAll=fitdat%>%
  mutate(setting0=setting,
         setting=ifelse(str_detect( setting0, pattern = ("IG2")), "IG2",
                        ifelse(str_detect( setting, pattern = ("PBO")), "OlysetPlus", "Pyrethroid")))

validat_pooledAll=validat%>%
  mutate(setting0=setting,
         setting=ifelse(str_detect( setting0, pattern = ("IG2")), "IG2",
                        ifelse(str_detect( setting, pattern = ("PBO")), "OlysetPlus", "Pyrethroid")))

##########################################################################
##############
##############################################
## with chirps data

outputs_assenga=calibration_wrapper(my_all_simul=all_simul_0, scens=scens, fitdat=fitdat, select_subset =.5, mystandardITN="BIT103")
outputs_martin=calibration_wrapper(my_all_simul=all_simul_0, scens=scens, fitdat=fitdat, select_subset =.5, mystandardITN="martin")
outputs_kibondo=calibration_wrapper(my_all_simul=all_simul_0, scens=scens, fitdat=fitdat, select_subset =.5, mystandardITN="kibondo")
outputs_BIT055=calibration_wrapper(my_all_simul=all_simul_0, scens=scens, fitdat=fitdat, select_subset =.5, mystandardITN="BIT055")
outputs_odufuwa=calibration_wrapper(my_all_simul=all_simul_0, scens=scens, fitdat=fitdat, select_subset =.5, mystandardITN="odufuwa")
outputs_BIT080=calibration_wrapper(my_all_simul=all_simul_0, scens=scens, fitdat=fitdat, select_subset =.5, mystandardITN="BIT080")
outputs_martinf=calibration_wrapper(my_all_simul=all_simul_0, scens=scens, fitdat=fitdat, select_subset =.5, mystandardITN="martinf")

results_protopopoff_pooled=rbind(outputs_assenga$seed_avg_pooled %>% mutate(pyrethroid="assenga"),
                             outputs_martin$seed_avg_pooled %>% mutate(pyrethroid="martin"),
                             outputs_kibondo$seed_avg_pooled %>% mutate(pyrethroid="kibondo"),
                             outputs_BIT055$seed_avg_pooled %>% mutate(pyrethroid="BIT055"),
                             outputs_odufuwa$seed_avg_pooled %>% mutate(pyrethroid="odufuwa"),
                             outputs_BIT080$seed_avg_pooled %>% mutate(pyrethroid="BIT080"),
                             outputs_martinf$seed_avg_pooled %>% mutate(pyrethroid="martinf"))

results_protopopoff_details=rbind(outputs_assenga$final_output$seed_avg %>% mutate(pyrethroid="assenga"),
                              outputs_martin$final_output$seed_avg  %>% mutate(pyrethroid="martin"),
                              outputs_kibondo$final_output$seed_avg %>% mutate(pyrethroid="kibondo"),
                              outputs_BIT055$final_output$seed_avg %>% mutate(pyrethroid="BIT055"),
                              outputs_odufuwa$final_output$seed_avg %>% mutate(pyrethroid="odufuwa"),
                              outputs_BIT080$final_output$seed_avg %>% mutate(pyrethroid="BIT080"),
                              outputs_martinf$final_output$seed_avg %>% mutate(pyrethroid="martinf"))


write.csv(fitdat_pooledAll %>% select(-setting0)%>% unique(), file = file.path(experiment_folder, "fitdat_protopopoff.csv"), row.names = F)
write.csv(validat_pooledAll %>% select(-setting0)%>% ungroup()%>%unique(), file = file.path(experiment_folder, "validat_protopopoff.csv"), row.names = F)
write.csv(results_protopopoff_pooled, file = file.path(experiment_folder, "results_protopopoff.csv"), row.names = F)
write.csv(results_protopopoff_details, file = file.path(experiment_folder, "results_protopopoff_detail.csv"), row.names = F)


EIR_range=rbind(outputs_assenga$calibrated_data$best_EIRs_l_ci_approx %>% mutate(pyrethroid="assenga"),
                outputs_martin$calibrated_data$best_EIRs_l_ci_approx %>% mutate(pyrethroid="martin"),
                outputs_kibondo$calibrated_data$best_EIRs_l_ci_approx %>% mutate(pyrethroid="kibondo"),
                outputs_BIT055$calibrated_data$best_EIRs_l_ci_approx %>% mutate(pyrethroid="BIT055"),
                outputs_odufuwa$calibrated_data$best_EIRs_l_ci_approx %>% mutate(pyrethroid="odufuwa"),
                outputs_BIT080$calibrated_data$best_EIRs_l_ci_approx %>% mutate(pyrethroid="BIT080"),
                outputs_martinf$calibrated_data$best_EIRs_l_ci_approx %>% mutate(pyrethroid="martinf")) 

EIR_range%>%
  mutate(arm=str_detect("pyrethroid", sub))%>%
  group_by(arm)%>%
  summarise(EIR_min=min(EIR_lci),EIR_max=max(EIR_uci))


plot_validation(my_seed_avg=outputs_martin$seed_avg_pooled %>% filter(year>2013), fitdat_pooledAll, validat_pooledAll, myAgePR = agePR,
                mylimits=c(2014,2018), mybreaks=2014:2018,distribution_time = 2015+27/365,
                mylabel="Prevalence in 6 months\nto 14 years old")+ylim(0,1)
ggsave(filename=file.path(experiment_folder,"Protopopoff_validation_Martin_pooled.png")
       ,width = 39,height=15)

plot_validation(my_seed_avg=outputs_assenga$seed_avg_pooled %>% filter(year>2013), fitdat_pooledAll, validat_pooledAll, myAgePR = agePR,
                mylimits=c(2014,2018), mybreaks=2014:2018,distribution_time = 2015+27/365,
                mylabel="Prevalence in 6 months\nto 14 years old")+ylim(0,1)
ggsave(filename=file.path(experiment_folder,"Protopopoff_validation_Assenga_pooled.png")
       ,width = 39,height=15)

