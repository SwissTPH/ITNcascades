rm(list=ls())

library(devtools)
#devtools::install_github("SwissTPH/r-openMalariaUtilities", ref = "master")
#devtools::install("/scicore/home/pothin/champa0000/omu-slurm") #run this if you have pulled in your local directory

library(openMalariaUtilities)
library(tidyverse)
#library(omuslurm)

main_folder=(".")
script_folder=""
parameters_GVI=read.csv(file.path(script_folder, "EHT_fit/fitted_parameters_posteriormax.csv"))
source(file.path(script_folder,"RCT_validation/helpfunctions_rct.R"))
source(file.path(script_folder,"RCT_validation/helpfunction_omu_gvi.R"))
source(file.path(script_folder,"helpfunctions_convert_cm.R"))



# SET EXPERIMENT ----------------------------------------------------------

COUNTRY     = "BEN"
SIMSTART    = "1918-01-01"
versionnum  = 44L


expName  = "RCTAccrombessi_full"
experiment_folder=file.path(main_folder, expName)


###
### Define vector bionomics
###

# which mosquitoes will be used ?
mosqs  = c("gambiaesl_indoor")
#' Mosha 2022 https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(21)02499-5/fulltext#sec1
#' page 1235

mosqs_exposures=c("gambiaesl_indoor"=1)


# contrib: the relative contribution of each 'mosquito' to EIR
contrib = c("@gin@" )
cbind( mosqs, contrib)


distrib_year=2020
distrib_month=3


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
  baseList = baseList, y1 = 2015, y2 = 2023, detect = 100, SIMSTART = SIMSTART,
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





halflife_pyrethroid=3
kappa_pyrethroid=0.73

halflife_IG2=3
kappa_IG2=0.93




hist_pyrethroid_nguessan= extract_GVI_params(EHT="Nguessan",
                                             netType="Pyrethroid-only",
                                             halflife_functionalSurvival=1, kappa_functionalSurvival=kappa_pyrethroid,
                                             myname="histITN_nguessan", mosqs=mosqs,  parameters_GVI=parameters_GVI, mosqs_exposures=mosqs_exposures, insecticide_decay = F, decay = "step")


hist_pyrethroid_sovegnon= extract_GVI_params(EHT="Sovegnon",
                                             netType="Pyrethroid-only",
                                             halflife_functionalSurvival=1, kappa_functionalSurvival=kappa_pyrethroid,
                                             myname="histITN_sovegnon", mosqs=mosqs,  parameters_GVI=parameters_GVI, mosqs_exposures=mosqs_exposures, insecticide_decay = F, decay = "step")


hist_pyrethroid_assenga= extract_GVI_params(EHT="Assenga, CI",
                                            netType="Pyrethroid-only",
                                            halflife_functionalSurvival=1, kappa_functionalSurvival=kappa_pyrethroid,
                                            myname="histITN_assenga", mosqs=mosqs,  parameters_GVI=parameters_GVI, mosqs_exposures=mosqs_exposures, insecticide_decay = F, decay = "step")


GVI_pyrethroid_nguessan= extract_GVI_params(EHT="Nguessan",
                                            netType="Pyrethroid-only",
                                            halflife_functionalSurvival=halflife_pyrethroid, kappa_functionalSurvival=kappa_pyrethroid,
                                            myname="futITN_nguessan", mosqs=mosqs,  parameters_GVI=parameters_GVI, mosqs_exposures=mosqs_exposures)


GVI_pyrethroid_sovegnon= extract_GVI_params(EHT="Sovegnon",
                                            netType="Pyrethroid-only",
                                            halflife_functionalSurvival=halflife_pyrethroid, kappa_functionalSurvival=kappa_pyrethroid,
                                            myname="futITN_sovegnon", mosqs=mosqs,  parameters_GVI=parameters_GVI, mosqs_exposures=mosqs_exposures)

GVI_pyrethroid_assenga= extract_GVI_params(EHT="Assenga, CI",
                                           netType="Pyrethroid-only",
                                           halflife_functionalSurvival=halflife_pyrethroid, kappa_functionalSurvival=kappa_pyrethroid,
                                           myname="futITN_assenga", mosqs=mosqs,  parameters_GVI=parameters_GVI, mosqs_exposures=mosqs_exposures)


GVI_IG2_nguessan= extract_GVI_params(EHT="Nguessan",
                                     netType="Interceptor G2",
                                     halflife_functionalSurvival=halflife_IG2, kappa_functionalSurvival=kappa_IG2,
                                     mosqs_exposures=mosqs_exposures, myname="futIG2_nguessan", mosqs=mosqs,  parameters_GVI=parameters_GVI)


GVI_IG2_sovegnon= extract_GVI_params(EHT="Sovegnon",
                                     netType="Interceptor G2",
                                     halflife_functionalSurvival=halflife_IG2, kappa_functionalSurvival=kappa_IG2,
                                     mosqs_exposures=mosqs_exposures, myname="futIG2_sovegnon", mosqs=mosqs,  parameters_GVI=parameters_GVI)


GVI_IG2_assenga= extract_GVI_params(EHT="Assenga, CI",
                                    netType="Interceptor G2",
                                    halflife_functionalSurvival=halflife_IG2, kappa_functionalSurvival=kappa_IG2,
                                    mosqs_exposures=mosqs_exposures, myname="futIG2_assenga", mosqs=mosqs,  parameters_GVI=parameters_GVI)


list_GVI_snippets=list(hist_pyrethroid_nguessan,hist_pyrethroid_sovegnon, hist_pyrethroid_assenga)

list_GVI_snippets=append(list_GVI_snippets, GVI_pyrethroid_nguessan)
list_GVI_snippets=append(list_GVI_snippets, GVI_pyrethroid_sovegnon)
list_GVI_snippets=append(list_GVI_snippets, GVI_pyrethroid_assenga)


list_GVI_snippets=append(list_GVI_snippets, GVI_IG2_nguessan)
list_GVI_snippets=append(list_GVI_snippets, GVI_IG2_sovegnon)
list_GVI_snippets=append(list_GVI_snippets, GVI_IG2_assenga)

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


## historical nets
list_histITN_snippets_id=c("histITN_nguessan","histITN_assenga",
                           "histITN_sovegnon"
)


for(GVI_snippet_id in list_histITN_snippets_id){
  baseList <- deploy_it_compat(
    baseList = baseList, component = paste0(GVI_snippet_id),
    coverage = paste0("@", GVI_snippet_id,"_cov@"),
    ## allowing for different hist coverage levels
    byyear = FALSE,
    ##  annual deployments
    y1 = 2015, y2 = 2019, every = 1, interval = "year",
    m1 = 3, m2 = 3, d1 = 20, d2 = 20, #to end just when the trial distribution starts
    SIMSTART = SIMSTART
  )
  
}



## Future IG2
list_GVI_snippets_id=c("futIG2_nguessan",
                       "futIG2_sovegnon",
                       "futIG2_assenga",
                       "futITN_assenga",
                       "futITN_nguessan",
                       "futITN_sovegnon"
)


for(GVI_snippet_id in list_GVI_snippets_id){
  baseList <- deploy_it_compat(
    baseList = baseList, component = paste0(GVI_snippet_id, "_deterrency"),
    coverage = paste0("@", GVI_snippet_id,"_cov@"),
    ## allowing for different hist coverage levels
    byyear = FALSE,
    ##  annual deployments
    y1 = 2020, y2 = 2020, every = 3, interval = "year",
    m1 = 3, m2 = 3, d1 = 20, d2 = 20, 
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
           # list(component = list(id = "futITN")),
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

#######################
### seasonality
season_dir="/scicore/home/pothin/champa0000/IG2modelling/"
seasonality_chirps=read.csv(file.path(season_dir, "seasonality_Zou.csv")) %>% mutate(season="chirps")

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
                     "IG2_nguessan",#"IG2_nguessan_inf","IG2_nguessan_sup",
                     "IG2_sovegnon",# "IG2_sovegnon_sup", "IG2_sovegnon_inf"
                     "IG2_assenga"
)
full$StandardITN=c("nguessan","sovegnon","assenga")
full $ pop       = 10000
full$ EffCovconv = convert_cm(0.401)
full $ EIR = seq(10,200,2)

full $ season = c("chirps")
full$gin = 1

#### 'scens' will do all possible combinations of those factors
scens = expand.grid( full )
scens = scens %>% left_join( seasonality_full)

exposure_in_bed=0.77
futIG2_cov=1*exposure_in_bed


NetCov <- data.frame(setting=c("pyrethroid","chlorfenapyr"),
                     histITNcov=c(0.80, 0.80)*exposure_in_bed,
                     futITNcov=c(1 , 0)*exposure_in_bed
)
NetCov_full=NetCov
for(mysetting in full$setting[-1]){
  #print(mysetting)
  is.pbo=str_detect( mysetting,"PBO")
  #print(is.pbo)
  myNetCov=NetCov %>%
    filter(setting %in% ifelse(is.pbo, "PBO","chlorfenapyr"))%>%
    mutate(futcov=ifelse(is.pbo, futPBO_cov,futIG2_cov), setting=mysetting)
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
age_groups_of_interest <- c("All")
## 6 months to 14 years old for prevalence
## 6 months to 10 years old for incidence

loadExperiment(expName)
## SLURM
slurmPrepareResults(
  
  expDir=experiment_folder,
  dbName="results_RCTAccrombessi",
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
  mem = "20G",
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
conn = DBI::dbConnect(RSQLite::SQLite(),file.path(expName,"results_RCTAccrombessi.sqlite"))
# list the tables of the database
DBI::dbListTables(conn)
# extract the results table, this should have the same name as provided to slurmPrepareResults
all_simul_0 = DBI::dbReadTable(conn, "results")
# always disconnect from the database
DBI::dbDisconnect(conn)

agePR="All"


all_simul_0=all_simul_0 %>% filter(age_group==agePR)



###############################
## fitdat
######################

fitdat = data.frame(setting=c("pyrethroid","pyriproxyfen","chlorfenapyr")
                    ,age=agePR
                    ,year=2019
                    ,month=10
                    ,pos_children=c(690,636,598)
                    ,nb_children=c(1485,1475,1468)) %>%
  rowwise() %>%
  mutate(PR_obs=pos_children/nb_children,
         sd=sqrt(PR_obs*(1-PR_obs)/nb_children),
         PR_obs_min=PR_obs-1.96*sd,
         PR_obs_max=PR_obs+1.96*sd)%>%
  filter(setting %in% c("pyrethroid","chlorfenapyr", "PBO"))

fitdat=rbind(fitdat,
             fitdat %>% filter(setting=="chlorfenapyr") %>% mutate(setting="IG2_nguessan"),
             fitdat %>% filter(setting=="chlorfenapyr") %>% mutate(setting="IG2_sovegnon"),
             fitdat %>% filter(setting=="chlorfenapyr") %>% mutate(setting="IG2_assenga")
)
fitdat=fitdat[fitdat$setting !="chlorfenapyr" & fitdat$setting !="PBO",]

######################
## validat
######################
validat0 = data.frame(setting=rep(c("pyrethroid","pyriproxyfen","chlorfenapyr"),3)
                      ,age=agePR
                      ,year=rep(c(2020,2021, 2022),each=3)
                      ,month=10
                      ,pos_children=c(412,394,231,
                                      576,564,414,
                                      386,422,331)
                      ,nb_children=c(1471,1463,1475,
                                     1489,1478,1483,
                                     1474, 1485, 1469)) %>%
  rowwise() %>%
  mutate(PR_obs=pos_children/nb_children,
         sd=sqrt(PR_obs*(1-PR_obs)/nb_children),
         PR_obs_min=PR_obs-1.96*sd,
         PR_obs_max=PR_obs+1.96*sd)%>%
  filter(setting %in% c("pyrethroid","chlorfenapyr"))

validat=rbind(validat0,
              validat0 %>% filter(setting=="chlorfenapyr") %>% mutate(setting="IG2_nguessan"),
              validat0 %>% filter(setting=="chlorfenapyr") %>% mutate(setting="IG2_sovegnon"),
              validat0 %>% filter(setting=="chlorfenapyr") %>% mutate(setting="IG2_assenga")
)

validat=validat[validat$setting !="chlorfenapyr" & validat$setting !="PBO",]



fitdat_pooledAll=fitdat%>%
  mutate(setting0=setting,
         setting=ifelse(str_detect( setting0, pattern = ("IG2")), "IG2",
                        ifelse(str_detect( setting, pattern = ("PBO")), "OlysetPlus", "Pyrethroid")))

validat_pooledAll=validat%>%
  mutate(setting0=setting,
         setting=ifelse(str_detect( setting0, pattern = ("IG2")), "IG2",
                        ifelse(str_detect( setting, pattern = ("PBO")), "OlysetPlus", "Pyrethroid")))


##############################################
##extract only month in fitdat and age
outputs_nguessan=calibration_wrapper(my_all_simul=all_simul_0, scens=scens, fitdat=fitdat, select_subset =.5, mystandardITN="nguessan")
outputs_sovegnon=calibration_wrapper(my_all_simul=all_simul_0, scens=scens, fitdat=fitdat, select_subset =.5, mystandardITN="sovegnon")
outputs_assenga=calibration_wrapper(my_all_simul=all_simul_0, scens=scens, fitdat=fitdat, select_subset =.5, mystandardITN="assenga")


results_accrombessi_pooled=rbind(outputs_nguessan$seed_avg_pooled %>% mutate(pyrethroid="nguessan"),
                                 outputs_assenga$seed_avg_pooled %>% mutate(pyrethroid="assenga"),
                                 outputs_sovegnon$seed_avg_pooled %>% mutate(pyrethroid="sovegnon"))

results_accrombessi_details=rbind(outputs_nguessan$final_output$seed_avg %>% mutate(pyrethroid="nguessan"),
                                  outputs_assenga$final_output$seed_avg %>% mutate(pyrethroid="assenga"),
                                  outputs_sovegnon$final_output$seed_avg  %>% mutate(pyrethroid="sovegnon"))



write.csv(fitdat_pooledAll %>% select(-setting0)%>% unique(), file = file.path(experiment_folder, "fitdat_accrombessi.csv"), row.names = F)
write.csv(validat_pooledAll %>% select(-setting0)%>% ungroup()%>%unique(), file = file.path(experiment_folder, "validat_accrombessi.csv"), row.names = F)
write.csv(results_accrombessi_pooled, file = file.path(experiment_folder, "results_accrombessi.csv"), row.names = F)
write.csv(results_accrombessi_details, file = file.path(experiment_folder, "results_accrombessi_detail.csv"), row.names = F)


EIR_range=rbind(outputs_assenga$calibrated_data$best_EIRs_l_ci_approx %>% mutate(pyrethroid="assenga"),
                outputs_nguessan$calibrated_data$best_EIRs_l_ci_approx %>% mutate(pyrethroid="nguessan"),
                outputs_sovegnon$calibrated_data$best_EIRs_l_ci_approx %>% mutate(pyrethroid="sovegnon")) 

EIR_range%>%
  mutate(arm=str_detect("pyrethroid", sub))%>%
  group_by(arm)%>%
  summarise(EIR_min=min(EIR_lci),EIR_max=max(EIR_uci))


plot_validation(my_seed_avg=outputs_nguessan$seed_avg_pooled %>% filter(year<2023, year>2017), fitdat_pooledAll, validat_pooledAll, myAgePR = agePR,
                mylimits=c(2019,2023), mybreaks=2019:2022,distribution_time = 2020+3/12,
                mylabel="Prevalence in all ages"
)
ggsave(filename=file.path(experiment_folder,"Accrombessi_validation_pooled_nguessan.png")
       ,width = 49,height=12)

plot_validation(my_seed_avg=outputs_sovegnon$seed_avg_pooled %>% filter(year<2023, year>2017), fitdat_pooledAll, validat_pooledAll, myAgePR = agePR,
                mylimits=c(2019,2023), mybreaks=2019:2022,distribution_time = 2020+3/12,
                mylabel="Prevalence in all ages"
)
ggsave(filename=file.path(experiment_folder,"Accrombessi_validation_pooled_sovegnon.png")
       ,width = 49,height=12)

plot_validation(my_seed_avg=outputs_assenga$seed_avg_pooled %>% filter(year<2023, year>2017), fitdat_pooledAll, validat_pooledAll, myAgePR = agePR,
                mylimits=c(2019,2023), mybreaks=2019:2022,distribution_time = 2020+3/12,
                mylabel="Prevalence in all ages"
)
ggsave(filename=file.path(experiment_folder,"Accrombessi_validation_pooled_assenga.png")
       ,width = 49,height=12)

