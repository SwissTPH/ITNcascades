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

COUNTRY     = "TZA"
SIMSTART    = "1918-01-01"
versionnum  = 44L


expName  = "RCTMosha_full"
experiment_folder=file.path(main_folder, expName)


###
### Define vector bionomics
###

# which mosquitoes will be used ?
mosqs  = c("funestus_indoor","gambiaesl_indoor")
#' Mosha 2022 https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(21)02499-5/fulltext#sec1
#' page 1235

exposure_in_bed=0.69

mosqs_exposures=c("funestus_indoor"=1, "gambiaesl_indoor"=1) # because exposures are entered in the coverage with this approach


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





halflife_pyrethroid=2.3
kappa_pyrethroid=1.6

halflife_PBO=1.7
kappa_PBO=1.9

halflife_IG2=2.4
kappa_IG2=2.4

# historical ITNs
hist_pyrethroid_martin= extract_GVI_params(EHT="Martin",
                                           netType="Pyrethroid-only",
                                           halflife_functionalSurvival=1, kappa_functionalSurvival=kappa_pyrethroid,
                                           myname="histITN_martin", mosqs=mosqs,  parameters_GVI=parameters_GVI, mosqs_exposures=mosqs_exposures, insecticide_decay = F, decay = "step")

hist_pyrethroid_BIT103= extract_GVI_params(EHT="Assenga",
                                           netType="Pyrethroid-only",
                                           halflife_functionalSurvival=1, kappa_functionalSurvival=kappa_pyrethroid,
                                           myname="histITN_BIT103", mosqs=mosqs,  parameters_GVI=parameters_GVI, mosqs_exposures=mosqs_exposures, insecticide_decay = F, decay = "step")

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

# future ITNs
GVI_pyrethroid_martin= extract_GVI_params(EHT="Martin",
                                          netType="Pyrethroid-only",
                                          halflife_functionalSurvival=halflife_pyrethroid, kappa_functionalSurvival=kappa_pyrethroid,
                                          myname="futITN_martin", mosqs=mosqs,  parameters_GVI=parameters_GVI, mosqs_exposures=mosqs_exposures)


GVI_pyrethroid_BIT103= extract_GVI_params(EHT="Assenga",
                                          netType="Pyrethroid-only",
                                          halflife_functionalSurvival=halflife_pyrethroid, kappa_functionalSurvival=kappa_pyrethroid,
                                          myname="futITN_BIT103", mosqs=mosqs,  parameters_GVI=parameters_GVI, mosqs_exposures=mosqs_exposures)

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

# future IG2
GVI_IG2_kibondo= extract_GVI_params(EHT="Kibondo",
                                    netType="Interceptor G2",
                                    halflife_functionalSurvival=halflife_IG2, kappa_functionalSurvival=kappa_IG2,
                                    mosqs_exposures=mosqs_exposures, myname="futIG2_kibondo", mosqs=mosqs,  parameters_GVI=parameters_GVI)

GVI_IG2_martin= extract_GVI_params(EHT="Martin",
                                   netType="Interceptor G2",
                                   halflife_functionalSurvival=halflife_IG2, kappa_functionalSurvival=kappa_IG2,
                                   mosqs_exposures=mosqs_exposures, myname="futIG2_martin", mosqs=mosqs,  parameters_GVI=parameters_GVI)

GVI_IG2_martinf= extract_GVI_params(EHT="Martin, f",
                                    netType="Interceptor G2",
                                    halflife_functionalSurvival=halflife_IG2, kappa_functionalSurvival=kappa_IG2,
                                    mosqs_exposures=mosqs_exposures, myname="futIG2_martinf", mosqs=mosqs,  parameters_GVI=parameters_GVI)

GVI_IG2_BIT103= extract_GVI_params(EHT="Assenga",
                                   netType="Interceptor G2",
                                   halflife_functionalSurvival=halflife_IG2, kappa_functionalSurvival=kappa_IG2,
                                   mosqs_exposures=mosqs_exposures, myname="futIG2_BIT103", mosqs=mosqs,  parameters_GVI=parameters_GVI)


GVI_IG2_BIT080= extract_GVI_params(EHT="BIT080",
                                   netType="Interceptor G2",
                                   halflife_functionalSurvival=halflife_IG2, kappa_functionalSurvival=kappa_IG2,
                                   mosqs_exposures=mosqs_exposures, myname="futIG2_BIT080", mosqs=mosqs,  parameters_GVI=parameters_GVI)

# future PBO
GVI_PBO_BIT055= extract_GVI_params(EHT="BIT055",
                                   netType="Olyset Plus",
                                   halflife_functionalSurvival=halflife_PBO, kappa_functionalSurvival=kappa_PBO,
                                   mosqs_exposures=mosqs_exposures, myname="futPBO_BIT055", mosqs=mosqs, parameters_GVI=parameters_GVI)

GVI_PBO_BIT059= extract_GVI_params(EHT="Odufuwa",
                                   netType="Olyset Plus",
                                   halflife_functionalSurvival=halflife_PBO, kappa_functionalSurvival=kappa_PBO,
                                   mosqs_exposures=mosqs_exposures, myname="futPBO_BIT059", mosqs=mosqs, parameters_GVI=parameters_GVI)

GVI_PBO_BIT103= extract_GVI_params(EHT="Assenga",
                                   netType="Olyset Plus",
                                   halflife_functionalSurvival=halflife_PBO, kappa_functionalSurvival=kappa_PBO,
                                   mosqs_exposures=mosqs_exposures, myname="futPBO_BIT103", mosqs=mosqs,  parameters_GVI=parameters_GVI)

GVI_PBO_martin= extract_GVI_params(EHT="Martin",
                                   netType="Olyset Plus",
                                   halflife_functionalSurvival=halflife_PBO, kappa_functionalSurvival=kappa_PBO,
                                   mosqs_exposures=mosqs_exposures, myname="futPBO_martin", mosqs=mosqs,  parameters_GVI=parameters_GVI)

GVI_PBO_martinf= extract_GVI_params(EHT="Martin, f",
                                    netType="Olyset Plus",
                                    halflife_functionalSurvival=halflife_PBO, kappa_functionalSurvival=kappa_PBO,
                                    mosqs_exposures=mosqs_exposures, myname="futPBO_martinf", mosqs=mosqs,  parameters_GVI=parameters_GVI)

list_GVI_snippets=list(hist_pyrethroid_BIT103,hist_pyrethroid_martin,hist_pyrethroid_kibondo,hist_pyrethroid_BIT055,hist_pyrethroid_odufuwa,hist_pyrethroid_BIT080,hist_pyrethroid_martinf)
list_GVI_snippets=append(list_GVI_snippets, GVI_pyrethroid_BIT103)
list_GVI_snippets=append(list_GVI_snippets, GVI_pyrethroid_martin)
list_GVI_snippets=append(list_GVI_snippets, GVI_pyrethroid_kibondo)
list_GVI_snippets=append(list_GVI_snippets, GVI_pyrethroid_BIT055)
list_GVI_snippets=append(list_GVI_snippets, GVI_pyrethroid_odufuwa)
list_GVI_snippets=append(list_GVI_snippets, GVI_pyrethroid_BIT080)
list_GVI_snippets=append(list_GVI_snippets, GVI_pyrethroid_martinf)


list_GVI_snippets=append(list_GVI_snippets, GVI_IG2_kibondo)
list_GVI_snippets=append(list_GVI_snippets, GVI_IG2_martin)
list_GVI_snippets=append(list_GVI_snippets, GVI_IG2_martinf)
list_GVI_snippets=append(list_GVI_snippets, GVI_IG2_BIT103)
list_GVI_snippets=append(list_GVI_snippets, GVI_IG2_BIT080)
list_GVI_snippets=append(list_GVI_snippets, GVI_PBO_BIT103)
list_GVI_snippets=append(list_GVI_snippets, GVI_PBO_martin)
list_GVI_snippets=append(list_GVI_snippets, GVI_PBO_martinf)
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
    y1 = 2015, y2 = 2018 , every = 1, interval = "year",
    m1 = 1, m2 = 1, d1 = 25, d2 = 25,
    #' Mosha 2022 Abstract, Findings
    #' https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(21)02499-5/fulltext#sec1
    #'"LLINs were distributed among households between Jan 26 and Jan 28, 2019"
    SIMSTART = SIMSTART
  )
}



## Future IG2
list_GVI_snippets_id=c("futIG2_kibondo",
                       "futIG2_BIT103",
                       "futIG2_BIT080",
                       "futIG2_martin",
                       "futPBO_martin",
                       "futIG2_martinf",
                       "futPBO_martinf",
                       "futPBO_BIT059",
                       "futPBO_BIT055",
                       "futPBO_BIT103",
                       "futITN_BIT103",
                       "futITN_martin",
                       "futITN_kibondo",
                       "futITN_BIT055",
                       "futITN_odufuwa",
                       "futITN_BIT080",
                       "futITN_martinf"
)


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
#' Mosha 2022 Methods, Mwanza
#' https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(21)02499-5/fulltext#sec1

season_dir="/scicore/home/pothin/champa0000/IG2modelling/"
seasonality_chirps=read.csv(file.path(season_dir, "seasonality_Misungwi.csv"))

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
                     "IG2_kibondo",
                     "PBO_BIT055", 
                     "PBO_BIT059", 
                     "PBO_BIT103", 
                     "IG2_martin",
                     "PBO_martin", 
                     "IG2_martinf",
                     "PBO_martinf",
                     "IG2_BIT103", 
                     "IG2_BIT080"
)
full$StandardITN=c("BIT103","martin", "kibondo", "BIT055", "odufuwa", "BIT080", "martinf")
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


futIG2_cov=0.69*exposure_in_bed
futPBO_cov=0.75*exposure_in_bed

NetCov <- data.frame(setting=c("pyrethroid","chlorfenapyr","PBO"),
                     histITNcov=c(2957/4962,2849/4803,2695/4369)*exposure_in_bed,
                     futITNcov=c(0.77 , 0,0)*exposure_in_bed
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

#scens=scens%>% filter(ID !=100000)
scens <- finalize_scenarios(scens)
storeScenarios(scens)

dim(scens)

length(unique(scens$setting))


## Generate scenarios

setup_scenarios(scenarios=scens[1,])
## SLURM
slurmPrepareScenarios(expName = expName, scenarios = scens, nCPU = 1, qos = "6hours", time = "06:00:00", 
                      rModule = "R/4.4.1-gfbf-2023b")
#slurmCreateScenarios()
check_scenarios_created(experiment_folder)

validateXML(xmlfile = file.path(experiment_folder, paste0("scenarios/", expName, "_7.xml")), schema = NULL, scenarios = scens)


validateXML(xmlfile = file.path(experiment_folder, paste0("scenarios/", expName, "_37314.xml")), schema = NULL, scenarios = scens)
validateXML(xmlfile = file.path(experiment_folder, paste0("scenarios/", expName, "_40310.xml")), schema = NULL, scenarios = scens)

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
age_groups_of_interest <- c("0.5-14")

loadExperiment(expName)
## SLURM
slurmPrepareResults(
  
  expDir=experiment_folder,
  dbName="results_RCTMosha",
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
conn = DBI::dbConnect(RSQLite::SQLite(),file.path(expName,"results_RCTMosha.sqlite"))
# list the tables of the database
DBI::dbListTables(conn)
# extract the results table, this should have the same name as provided to slurmPrepareResults
all_simul_0 = DBI::dbReadTable(conn, "results")
# always disconnect from the database
DBI::dbDisconnect(conn)

agePR="0.5-14"

all_simul_0=all_simul_0 %>% filter(age_group==agePR)


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
             fitdat %>% filter(setting=="chlorfenapyr") %>% mutate(setting="IG2_BIT080"),
             fitdat %>% filter(setting=="chlorfenapyr") %>% mutate(setting="IG2_martin"),
             fitdat %>% filter(setting=="chlorfenapyr") %>% mutate(setting="IG2_martinf"),
             fitdat %>% filter(setting=="PBO") %>% mutate(setting="PBO_BIT055"),
             fitdat %>% filter(setting=="PBO") %>% mutate(setting="PBO_BIT059"),
             fitdat %>% filter(setting=="PBO") %>% mutate(setting="PBO_martin"),
             fitdat %>% filter(setting=="PBO") %>% mutate(setting="PBO_martinf"),
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
              validat0 %>% filter(setting=="chlorfenapyr") %>% mutate(setting="IG2_BIT080"),
              validat0 %>% filter(setting=="chlorfenapyr") %>% mutate(setting="IG2_martin"),
              validat0 %>% filter(setting=="chlorfenapyr") %>% mutate(setting="IG2_martinf"),
              validat0 %>% filter(setting=="PBO") %>% mutate(setting="PBO_BIT055"),
              validat0 %>% filter(setting=="PBO") %>% mutate(setting="PBO_BIT059"),
              validat0 %>% filter(setting=="PBO") %>% mutate(setting="PBO_martin"),
              validat0 %>% filter(setting=="PBO") %>% mutate(setting="PBO_martinf"),
              validat0 %>% filter(setting=="PBO") %>% mutate(setting="PBO_BIT103"))

validat=validat[validat$setting !="chlorfenapyr" & validat$setting !="PBO",]

##############################################
##extract only month in fitdat and age
outputs_assenga=calibration_wrapper(my_all_simul=all_simul_0, scens=scens, fitdat=fitdat, select_subset =.5, mystandardITN="BIT103")
outputs_martin=calibration_wrapper(my_all_simul=all_simul_0, scens=scens, fitdat=fitdat, select_subset =.5, mystandardITN="martin")
outputs_kibondo=calibration_wrapper(my_all_simul=all_simul_0, scens=scens, fitdat=fitdat, select_subset =.5, mystandardITN="kibondo")
outputs_BIT055=calibration_wrapper(my_all_simul=all_simul_0, scens=scens, fitdat=fitdat, select_subset =.5, mystandardITN="BIT055")
outputs_odufuwa=calibration_wrapper(my_all_simul=all_simul_0, scens=scens, fitdat=fitdat, select_subset =.5, mystandardITN="odufuwa")
outputs_BIT080=calibration_wrapper(my_all_simul=all_simul_0, scens=scens, fitdat=fitdat, select_subset =.5, mystandardITN="BIT080")
outputs_martinf=calibration_wrapper(my_all_simul=all_simul_0, scens=scens, fitdat=fitdat, select_subset =.5, mystandardITN="martinf")

results_mosha_pooled=rbind(outputs_assenga$seed_avg_pooled %>% mutate(pyrethroid="assenga"),
                                 outputs_martin$seed_avg_pooled %>% mutate(pyrethroid="martin"),
                                 outputs_kibondo$seed_avg_pooled %>% mutate(pyrethroid="kibondo"),
                                 outputs_BIT055$seed_avg_pooled %>% mutate(pyrethroid="BIT055"),
                                 outputs_odufuwa$seed_avg_pooled %>% mutate(pyrethroid="odufuwa"),
                                 outputs_BIT080$seed_avg_pooled %>% mutate(pyrethroid="BIT080"),
                                 outputs_martinf$seed_avg_pooled %>% mutate(pyrethroid="martinf"))

results_mosha_details=rbind(outputs_assenga$final_output$seed_avg %>% mutate(pyrethroid="assenga"),
                                  outputs_martin$final_output$seed_avg  %>% mutate(pyrethroid="martin"),
                                  outputs_kibondo$final_output$seed_avg %>% mutate(pyrethroid="kibondo"),
                                  outputs_BIT055$final_output$seed_avg %>% mutate(pyrethroid="BIT055"),
                                  outputs_odufuwa$final_output$seed_avg %>% mutate(pyrethroid="odufuwa"),
                                  outputs_BIT080$final_output$seed_avg %>% mutate(pyrethroid="BIT080"),
                                  outputs_martinf$final_output$seed_avg %>% mutate(pyrethroid="martinf"))


# take the envelope per net type
fitdat_pooledAll=fitdat%>%
  mutate(setting0=setting,
         setting=ifelse(str_detect( setting0, pattern = ("IG2")), "IG2",
                        ifelse(str_detect( setting, pattern = ("PBO")), "OlysetPlus", "Pyrethroid")))

validat_pooledAll=validat%>%
  mutate(setting0=setting,
         setting=ifelse(str_detect( setting0, pattern = ("IG2")), "IG2",
                        ifelse(str_detect( setting, pattern = ("PBO")), "OlysetPlus", "Pyrethroid")))


write.csv(fitdat_pooledAll %>% select(-setting0)%>% unique(), file = file.path(experiment_folder, "fitdat_mosha.csv"), row.names = F)
write.csv(validat_pooledAll %>% select(-setting0)%>% ungroup()%>%unique(), file = file.path(experiment_folder, "validat_mosha.csv"), row.names = F)
write.csv(results_mosha_pooled, file = file.path(experiment_folder, "results_mosha.csv"), row.names = F)
write.csv(results_mosha_details, file = file.path(experiment_folder, "results_mosha_detail.csv"), row.names = F)


EIR_range=rbind(outputs_assenga$calibrated_data$best_EIRs_l_ci_approx %>% mutate(pyrethroid="assenga"),
                outputs_martin$calibrated_data$best_EIRs_l_ci_approx %>% mutate(pyrethroid="martin"),
                outputs_kibondo$calibrated_data$best_EIRs_l_ci_approx %>% mutate(pyrethroid="kibondo"),
                outputs_BIT055$calibrated_data$best_EIRs_l_ci_approx %>% mutate(pyrethroid="BIT055"),
                outputs_odufuwa$calibrated_data$best_EIRs_l_ci_approx %>% mutate(pyrethroid="odufuwa"),
                outputs_BIT080$calibrated_data$best_EIRs_l_ci_approx %>% mutate(pyrethroid="BIT080"),
                outputs_martinf$calibrated_data$best_EIRs_l_ci_approx %>% mutate(pyrethroid="martinf")) 

EIR_range%>%
  group_by(sub)%>%
  summarise(EIR_min=min(EIR_lci),EIR_max=max(EIR_uci))

plot_validation(my_seed_avg=outputs_assenga$seed_avg_pooled %>% filter(year<2023, year>2017), fitdat_pooledAll, validat_pooledAll, myAgePR = agePR,
                mylimits=c(2018.5,2023), mybreaks=2018:2022,distribution_time = 2019+27/365,
                mylabel="Prevalence in 6 months\nto 14 years old"
)
ggsave(filename=file.path(experiment_folder,"Mosha_validation_assenga.png")
       ,width = 49,height=12)

plot_validation(my_seed_avg=outputs_martin$seed_avg_pooled %>% filter(year<2023, year>2017), fitdat_pooledAll, validat_pooledAll, myAgePR = agePR,
                mylimits=c(2018.5,2023), mybreaks=2018:2022,distribution_time = 2019+27/365,
                mylabel="Prevalence in 6 months\nto 14 years old"
)
ggsave(filename=file.path(experiment_folder,"Mosha_validation_martin.png")
       ,width = 49,height=12)

