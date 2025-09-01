
<br/>
<div align="center">
<a href="https://github.com/ShaanCoding/ReadME-Generator">
<img src="https://avatars.githubusercontent.com/u/5869420?s=200&v=4" alt="Logo" width="80" height="80">
</a>
<h3 align="center">Cascades of effectiveness of new-generation insecticide-treated nets against malaria, from entomological trials to real-life conditions</h3>
<p align="center">
<br/>

</p>
</div>

## About The Project

This code is associated with the following manuscript: 

__Cascades of effectiveness of new-generation insecticide-treated nets against malaria, from entomological trials to real-life conditions__  
Clara Champagne, Jeanne Lemant, Alphonce Assenga, Ummi A. Kibondo, Ruth G. Lekundayo, Emmanuel Mbuba, Jason Moore, Joseph B. Muganga, Watson S. Ntabaliba, Olukayode G. Odufuwa, Johnson Kyeba Swai, Maria Alexa, Roland Goers, Monica Golumbeanu, Nakul Chitnis, Amanda Ross, Raphael N'Guessan, Sarah Moore, Emilie Pothin


The user should refer to the [manuscript](https://www.medrxiv.org/content/10.1101/2025.02.07.25321565v1) for methodological details. 
Additionally, a dashboard to explore the results is available [here](https://swisstph.shinyapps.io/ITNcascadesdashboard/).



### Built With


- [AnophelesModel](https://github.com/SwissTPH/AnophelesModel/tree/main), version 1.1.0 and later
- [OpenMalaria](https://github.com/SwissTPH/openmalaria), version 44
- [OpenMalariaUtilities](https://github.com/SwissTPH/r-openMalariaUtilities), version 23.02
- [rstan](https://cran.r-project.org/web/packages/rstan/index.html), version 2.37.7

## Description of the code

The code is structured in two folders:

- __EHT_fit__: includes the code for the statistical analysis of Experimental Hut Trial data, as well as for building effectiveness cascades
- __RCT_validation__: includes the code for mathematical modelling of Randomized Controlled Trials with OpenMalaria
- __ITNcascadesdashboard__: includes the code for the online dashboard

### EHT_fit

The main scripts to be executed in that order are the following:

- __1_Data_formatting__: contains the code used to process the input EHT data, convert it into the proper format for inference, and compute descriptive statistics
- __2_EHT_fitting__: contains the code to perform the statistical inference as well as provide the effectiveness cascades
- __2b_convergence_diagnostics__: contains the convergence diagnostics for MCMC outputs with stan
- __3_compute_vectorial_capacity__: computes the vectorial capacity and summarises the of the EHT analysis outcomes (Figure 1)
- __3b_compute_vectorial_capacity_24__: computes the vectorial capacity and summarises the of the EHT analysis outcomes, suing 24h holding time for all EHTs
- __4_EffectivenessCascades__: calculates and visualise the effectiveness cascades

Other scripts are functions to be sourced while running the analysis.

The folder __processed_data__ contains the data used for the EHT fitting.
The file __fitted_parameters_posteriormax.csv__ contains the fitted parameters (posterior maximum). 
The file __fitted_parameters_pyr.csv__ contains summary statistics on the posterior distribution for each fitted parameter.

### RCT_validation

The main scripts to be executed in that order are the following:

- __1_seasonality__: contains the code used to export seasonality patterns from CHIRPS
- __2_ITNfunctionalSurvival__: contains the code to compute functional survival and net durability parameters for each RCT
- __3_get_exposure_all_trials__: contains the code to calculate in-bed exposure coefficients in each trial
- __4a_RCT_mosha__ , __4b_RCT_protopopoff__, __4c_RCT_accrombessi__ and __4b_RCT_staedke__: contains the code required to run OpenMalaria for the intervention arms of each RCT (validation)
- __5_validation_summary_figure__: contains the code for final visualisation of the results as well as goodness of fit diagnostics

Other scripts are functions to be sourced while running the analysis.
The folder __csv_inputs__ contains the data and simulation outputs required to reproduce the manuscript's figures.


## Authors contact

- Clara CHAMPAGNE - [@GitHub](https://github.com/clchampag) - clara.champagne@swisstph.ch


## Credits
This readme has been built with Bilal BENHANA thanks to a template
- [makeread.me](https://github.com/ShaanCoding/ReadME-Generator)
- [othneildrew](https://github.com/othneildrew/Best-README-Template)