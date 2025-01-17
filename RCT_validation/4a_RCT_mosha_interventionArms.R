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

source("helpfunctions_rct.R")
source("helpfunction_omu_gvi.R")
main_folder=(".")
parameters_GVI=read.csv(file.path(main_folder, "fitted_parameters_all_new_VCred.csv"))


# SET EXPERIMENT ----------------------------------------------------------

COUNTRY     = "TZA"
SIMSTART    = "1918-01-01"
versionnum  = 44L

expName  = "RCTMosha_validation"
experiment_folder=file.path(main_folder, expName)


###
### Define vector bionomics
###

# which mosquitoes will be used ?
mosqs  = c("funestus_indoor","gambiaesl_indoor")
#' Mosha 2022 https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(21)02499-5/fulltext#sec1
#' page 1235

# contrib: the relative contribution of each 'mosquito' to EIR
contrib = c("@fin@", "@gin@" )
cbind( mosqs, contrib)


distrib_year=2019
distrib_month=1
#' Mosha 2022 https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(21)02499-5/fulltext#sec1
#' Abtract, Findings


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


futITN= create_vectorInterventionParameters(deterrency = 0,
                                            preprandial = "@preprandial@",
                                            postprandial = 0,
                                            kappa=1.6,
                                            L=2.3, myname = "futITN", mosqs = mosqs)


exposure=0.69
halflife_PBO=1.7
kappa_PBO=1.9

halflife_IG2=2.4
kappa_IG2=2.4

GVI_IG2_kibondo= extract_GVI_params(EHT="Kibondo",
                                    netType="Interceptor G2",
                                    halflife_functionalSurvival=halflife_IG2, kappa_functionalSurvival=kappa_IG2,
                                    exposure_correction=exposure, myname="futIG2_kibondo", mosqs=mosqs,  parameters_GVI=parameters_GVI)


GVI_IG2_BIT103= extract_GVI_params(EHT="BIT103",
                                   netType="Interceptor G2",
                                   halflife_functionalSurvival=halflife_IG2, kappa_functionalSurvival=kappa_IG2,
                                   exposure_correction=exposure, myname="futIG2_BIT103", mosqs=mosqs,  parameters_GVI=parameters_GVI)


GVI_IG2_BIT080= extract_GVI_params(EHT="BIT080",
                                   netType="Interceptor G2",
                                   halflife_functionalSurvival=halflife_IG2, kappa_functionalSurvival=kappa_IG2,
                                   exposure_correction=exposure, myname="futIG2_BIT080", mosqs=mosqs,  parameters_GVI=parameters_GVI)


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
list_GVI_snippets=append(list_GVI_snippets, GVI_IG2_kibondo)
list_GVI_snippets=append(list_GVI_snippets, GVI_IG2_BIT103)
list_GVI_snippets=append(list_GVI_snippets, GVI_IG2_BIT080)
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
  y1 = 2015, y2 = 2018 , every = 1, interval = "year",
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

## Future IG2
list_GVI_snippets_id=c("futIG2_kibondo","futIG2_kibondo_inf","futIG2_kibondo_sup",
                       "futIG2_BIT103","futIG2_BIT103_sup","futIG2_BIT103_inf",
                       "futIG2_BIT080","futIG2_BIT080_sup","futIG2_BIT080_inf",
                       "futPBO_BIT059","futPBO_BIT059_sup","futPBO_BIT059_inf",
                       "futPBO_BIT055","futPBO_BIT055_sup","futPBO_BIT055_inf",
                       "futPBO_BIT103","futPBO_BIT103_sup","futPBO_BIT103_inf")


for(GVI_snippet_id in list_GVI_snippets_id){
  baseList <- deploy_it_compat(
    baseList = baseList, component = paste0(GVI_snippet_id, "_deterrency"),
    coverage = paste0("@", GVI_snippet_id,"_cov@"),
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

#######################
### seasonality
#' Mosha 2022 Methods, Mwanza
#' https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(21)02499-5/fulltext#sec1

season_dir="/scicore/home/pothin/champa0000/IG2modelling/"
seasonality_chirps=read.csv(file.path(season_dir, "seasonality_Misungwi.csv"))

seasonality_full=rbind(seasonality_chirps %>% mutate(season="chirps"))
names(seasonality_full)=c("m4", "m5", "m6","m7","m8","m9" ,"m10", "m11", "m12", "m1", "m2", "m3", "season")

seasonality_full %>%
  pivot_longer(cols=starts_with("m"),names_to = "month",values_to="intensity",names_prefix = "m") %>%
  mutate(month=as.numeric(month), dataset="CHIRPS +3m") %>%
  ggplot()+
  geom_line(aes(x=month,y=intensity, color=dataset))+
  theme_minimal()+
  scale_color_manual(values=c("red", "orange", "dodgerblue"), name="")

##-- creating the experiment
full             = list()
full $ seed      = 1:10
full $ preprandial = c(0, 0.05)
full $ setting   = c("pyrethroid",
                     "IG2_kibondo","IG2_kibondo_inf","IG2_kibondo_sup",
                     "PBO_BIT055", "PBO_BIT055_sup", "PBO_BIT055_inf",
                     "PBO_BIT059", "PBO_BIT059_sup", "PBO_BIT059_inf",
                     "PBO_BIT103", "PBO_BIT103_sup", "PBO_BIT103_inf",
                     "IG2_BIT103", "IG2_BIT103_sup", "IG2_BIT103_inf",
                     "IG2_BIT080", "IG2_BIT080_sup", "IG2_BIT080_inf")
full $ pop       = 10000
##' DHS 2017 for Mwanza province, effective coverage value
##' code by Katya and Rémi https://clintonhealth.box.com/s/02vk4o11ta7etymks0du8hbzsb925n69
full$ EffCovconv = convert_cm(0.6155460)
#full $ EIR_type = c("middle", "lower", "upper")
full $ EIR = seq(10,200,2)

full $ season = c("chirps")
#' https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(21)02499-5/fulltext#sec1
#' Mosha 2022 page 1235
full$fin = 4700/4973 
full$gin = (4973-4700)/4973 

#### 'scens' will do all possible combinations of those factors
scens = expand.grid( full )
scens = scens %>% left_join( seasonality_full)

##' historical net coverage: https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(21)02499-5/fulltext#sec1
##' Table 1, Long-lasting insecticidal net use In household residents of all age groups
##' future coverage: computed above


futIG2_cov=0.69
futPBO_cov=0.75
NetCov <- data.frame(setting=c("pyrethroid","chlorfenapyr","PBO"),
                     histITNcov=c(2957/4962,2849/4803,2695/4369),
                     futITNcov=c(0.77 , 0,0)
)

NetCov_full=NetCov
for(mysetting in full$setting[-1]){
  print(mysetting)
  is.pbo=str_detect( mysetting,"PBO")
  print(is.pbo)
  myNetCov=NetCov %>%
    filter(setting %in% ifelse(is.pbo, "PBO","chlorfenapyr"))%>%
    mutate(futcov=ifelse(is.pbo, futPBO_cov,futIG2_cov), setting=mysetting)
  names(myNetCov)[4]=paste0("fut", mysetting, "_cov")
  NetCov_full=merge(NetCov_full, myNetCov, all=TRUE)
}
NetCov_full[is.na(NetCov_full)]=0

scens <- left_join(scens,NetCov_full)

scens <- data.frame(scens)


####
#### writing the scenarios.csv
####
scens <- finalize_scenarios(scens)
storeScenarios(scens)

dim(scens)

length(unique(scens$setting))


## Generate scenarios

setup_scenarios(scenarios=scens[1,])
## SLURM
slurmPrepareScenarios(expName = expName, scenarios = scens, nCPU = 1, qos = "30min", time = "00:10:00", 
                      rModule = "R/4.4.1-gfbf-2023b")
#slurmCreateScenarios()
check_scenarios_created(experiment_folder)

validateXML(xmlfile = file.path(experiment_folder, paste0("scenarios/", expName, "_7.xml")), schema = NULL, scenarios = scens)

## Run Open Malaria
## SLURM
slurmPrepareSimulations(expName = expName, scenarios = scens, nCPU = 1, bSize = 1,
                        rModule = "R/4.4.1-gfbf-2023b", 
                        omModule = "OpenMalaria/44.0-intel-compilers-2023.1.0")
#slurmRunSimulation()
check_simulations_created(experiment_folder)

# failed=read.table(file.path(experiment_folder, "missing_om_simulations.txt"))
# failed$V2=gsub("_out.txt","", gsub(paste0(expName, "_"), "", failed$V1))
# #
# paste0(failed$V2,collapse=",")


####
#### Step 3: preparation for doing the post-processing
####

## Post-processing
age_groups_of_interest <- c("0.5-14","0.5-10", "All")
## 6 months to 14 years old for prevalence
## 6 months to 10 years old for incidence

loadExperiment(expName)
## SLURM
slurmPrepareResults(

  expDir=experiment_folder,
  dbName="results_RCTMosha_interventionArms_decay",
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
  splitBy = 10000,
  verbose = FALSE,
  ntasks = 1,
  mem = "64G",
  nCPU = 1,
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


loadExperiment(expName)
scens=readRDS(file.path(experiment_folder, "cache/scenarios.rds"))

# PREPARE THE DATA TO FIT

# open database connection
conn = DBI::dbConnect(RSQLite::SQLite(),file.path(expName,"results_RCTMosha_interventionArms_decay.sqlite"))
# list the tables of the database
DBI::dbListTables(conn)
# extract the results table, this should have the same name as provided to slurmPrepareResults
all_simul_0 = DBI::dbReadTable(conn, "results")
# always disconnect from the database
DBI::dbDisconnect(conn)

all_simul_chirps=merge(all_simul_0 , scens %>%mutate(scenario_id=ID))%>%
  mutate(EIR=as.numeric(EIR),
         month=as.numeric(format(as.Date(date), "%m")),
         year=as.numeric(format(as.Date(date), "%y"))+2000)

agePR="0.5-14"

###############################
## fitdat
######################

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
         sd=sqrt(PR_obs*(1-PR_obs)/nb_children),
         PR_obs_min=PR_obs-1.96*sd,
         PR_obs_max=PR_obs+1.96*sd)%>%
  filter(setting %in% c("pyrethroid","chlorfenapyr", "PBO"))

fitdat=rbind(fitdat,
             fitdat %>% filter(setting=="chlorfenapyr") %>% mutate(setting="IG2_kibondo"),
             fitdat %>% filter(setting=="chlorfenapyr") %>% mutate(setting="IG2_BIT103"),
             fitdat %>% filter(setting=="PBO") %>% mutate(setting="PBO_BIT055"),
             fitdat %>% filter(setting=="PBO") %>% mutate(setting="PBO_BIT103")
)
fitdat=fitdat[fitdat$setting !="chlorfenapyr" & fitdat$setting !="PBO",]

######################
## validat
######################
#' Mosha 2022 Table 2
#' #'https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(21)02499-5/fulltext#sec1
validat0 =  data.frame(setting=rep(c("pyrethroid","pyriproxyfen","chlorfenapyr","PBO"),each=5)
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
                                      1071,1160,1259,992,1049))  %>%
  rowwise() %>%
  mutate(PR_obs=pos_children/nb_children,
         sd=sqrt(PR_obs*(1-PR_obs)/nb_children),
         PR_obs_min=PR_obs-1.96*sd,
         PR_obs_max=PR_obs+1.96*sd)%>%
  filter(setting %in% c("pyrethroid","chlorfenapyr", "PBO"))
validat=rbind(validat0,
              validat0 %>% filter(setting=="chlorfenapyr") %>% mutate(setting="IG2_kibondo"),
              validat0 %>% filter(setting=="chlorfenapyr") %>% mutate(setting="IG2_BIT103"),
              validat0 %>% filter(setting=="PBO") %>% mutate(setting="PBO_BIT055"),
              validat0 %>% filter(setting=="PBO") %>% mutate(setting="PBO_BIT103"))

validat=validat[validat$setting !="chlorfenapyr" & validat$setting !="PBO",]

##############################################
##extract only month in fitdat and age

calibrated_data=calibration(all_simul_chirps %>% filter(preprandial==0.05), fitdat = fitdat)
calibrated_data0=calibration(all_simul_chirps %>% filter(preprandial==0), fitdat = fitdat)

# check calibration
plotQuadraticApproxMaxMinEIR(quadratic_l_approx = calibrated_data$quadratic_l_approx
                             , best_EIRs_l_ci_approx = calibrated_data$best_EIRs_l_ci_approx
                             ,list_setting = fitdat$setting)

final_output=merge_calibration_simulation(calibration_reformat=calibrated_data$calibration_reformat, my_all_simul=all_simul_chirps %>%filter(preprandial==0.05) )
final_output0=merge_calibration_simulation(calibration_reformat=calibrated_data0$calibration_reformat, my_all_simul=all_simul_chirps %>%filter(preprandial==0) )


# take the envelope per setting
final_output$seed_avg= final_output$simul_calib %>%
  mutate(sup=str_detect( setting, pattern = ("_sup")), inf=str_detect( setting, pattern = ("_inf")), setting=gsub("_inf", "", gsub("_sup", "", setting)))%>%
  group_by(setting, year, month, age )%>%
  summarise(PR_mini=min(PR_mini), PR_maxi=max(PR_maxi), PR_middle=mean(PR_middle[sup==FALSE & inf==FALSE]))

final_output0$seed_avg= final_output0$simul_calib %>%
  mutate(sup=str_detect( setting, pattern = ("_sup")), inf=str_detect( setting, pattern = ("_inf")), setting=gsub("_inf", "", gsub("_sup", "", setting)))%>%
  group_by(setting, year, month, age )%>%
  summarise(PR_mini=min(PR_mini), PR_maxi=max(PR_maxi), PR_middle=mean(PR_middle[sup==FALSE & inf==FALSE]))


# take the envelope per net type
fitdat_pooledAll=fitdat%>%
  mutate(setting0=setting,
         setting=ifelse(str_detect( setting0, pattern = ("IG2")), "IG2",
                        ifelse(str_detect( setting, pattern = ("PBO")), "OlysetPlus", "Pyrethroid")))

validat_pooledAll=validat%>%
  mutate(setting0=setting,
         setting=ifelse(str_detect( setting0, pattern = ("IG2")), "IG2",
                        ifelse(str_detect( setting, pattern = ("PBO")), "OlysetPlus", "Pyrethroid")))

seed_avg_pooledAll_TZ= final_output$simul_calib %>%
  filter(setting !="IG2_nguessan")%>%
  mutate(IG2=str_detect( setting, pattern = ("IG2")), PBO=str_detect( setting, pattern = ("PBO")),
         sup=str_detect( setting, pattern = ("_sup")), inf=str_detect( setting, pattern = ("_inf")),
         setting0=gsub("_inf", "", gsub("_sup", "", setting)),
         setting=ifelse(IG2, "IG2", ifelse(PBO, "OlysetPlus", "Pyrethroid")),
  )%>%
  group_by(setting, year, month, age )%>%
  summarise(PR_mini=min(PR_mini), PR_maxi=max(PR_maxi), PR_middle=mean(PR_middle[sup==FALSE & inf==FALSE]))

seed_avg_pooledAll0_TZ= final_output0$simul_calib %>%
  filter(setting !="IG2_nguessan")%>%
  mutate(IG2=str_detect( setting, pattern = ("IG2")), PBO=str_detect( setting, pattern = ("PBO")),
         sup=str_detect( setting, pattern = ("_sup")), inf=str_detect( setting, pattern = ("_inf")),
         setting0=gsub("_inf", "", gsub("_sup", "", setting)),
         setting=ifelse(IG2, "IG2", ifelse(PBO, "OlysetPlus", "Pyrethroid")),
  )%>%
  group_by(setting, year, month, age )%>%
  summarise(PR_mini=min(PR_mini), PR_maxi=max(PR_maxi), PR_middle=mean(PR_middle[sup==FALSE & inf==FALSE]))


fitdat_pooledAll=fitdat_pooledAll%>% full_join(validat_pooledAll %>% filter(setting=="Pyrethroid"))
validat_pooledAll=validat_pooledAll%>% filter(setting!="Pyrethroid")

fitdat_plot=fitdat%>% full_join(validat %>% filter(setting=="pyrethroid"))
validat_plot=validat%>% filter(setting!="pyrethroid")


write.csv(fitdat_pooledAll %>% select(-setting0)%>% unique(), file = file.path(experiment_folder, "fitdat_mosha.csv"), row.names = F)
write.csv(validat_pooledAll %>% select(-setting0)%>% ungroup()%>%unique(), file = file.path(experiment_folder, "validat_mosha.csv"), row.names = F)
write.csv(seed_avg_pooledAll_TZ, file = file.path(experiment_folder, "results_mosha.csv"), row.names = F)
write.csv(seed_avg_pooledAll0_TZ, file = file.path(experiment_folder, "results_mosha0.csv"), row.names = F)
write.csv(final_output$seed_avg, file = file.path(experiment_folder, "results_mosha_detail.csv"), row.names = F)
