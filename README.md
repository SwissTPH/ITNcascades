
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
Clara Champagne, Jeanne Lemant, Alphonce Assenga, Ummi A. Kibondo, Ruth G. Lekundayo, Emmanuel Mbuba, Jason Moore, Joseph B. Muganga, Watson S. Ntabaliba, Olukayode G. Odufuwa, Johnson Kyeba Swai, Maria Alexa, Roland Goers, Monica Golumbeanu, Nakul Chitnis, Amanda Ross, Sarah Moore, Emilie Pothin


The user should refer to the manuscript for methodological details. 
Additionally, a dashboard to explore the results is available [here](https://aimswisstph.shinyapps.io/ITNcascadesdashboard/).



### Built With


- [AnophelesModel](https://github.com/SwissTPH/AnophelesModel/tree/main)
- [OpenMalaria](https://github.com/SwissTPH/openmalaria)
- [OpenMalariaUtilities](https://github.com/SwissTPH/r-openMalariaUtilities)
- [rstan](https://cran.r-project.org/web/packages/rstan/index.html)

## Description of the code

The code is structured in two folders:

- __EHT_fit__: includes the code for the statistical analysis of Experimental Hut Trial data, as well as for building effectiveness cascades
- __RCT_validation__: includes the code for mathematical modelling of Randomized Controlled Trials with OpenMalaria

### EHT_fit

The main scripts to be executed in that order are the following:

- __1_Data_formatting__: contains the code used to process the input EHT data, convert it into the proper format for inference, and compute descriptive statistics
- __2_EHT_fitting_multinomial__: contains the code to perform the statistical inference as well as provide the effectiveness cascades
- __3_convergence_diagnostics__: contains the convergence diagnostics for MCMC outputs with stan
- __4_create_inputs_dasboards__: formats and summarises the outputs for the dashboard

Other scripts are functions to be sourced while running the analysis.

### RCT_validation

The main scripts to be executed in that order are the following:

- __1_seasonality__: contains the code used to export seasonality patterns from CHIRPS
- __2a_fit_ITNcov_mosha__ and __2b_fit_ITNcov_protopopoff__: contains the code to compute functional survival and net durability parameters for each RCT
- __3a_RCT_mosha_controlArm__ and __3b_RCT_protopopoff_controlArm__: contains the code required to run OpenMalaria for the control arms of each RCT
- __4a_RCT_mosha_interventionArms__ and __4b_RCT_protopopoff_interventionArms__: contains the code required to run OpenMalaria for the intervention arms of each RCT (validation)
- __5a_validation_final_plot__: contains the code for final visualisation of the results as well as goodness of fit diagnostics

Other scripts are functions to be sourced while running the analysis.

## Authors contact

- Clara CHAMPAGNE - [@GitHub](https://github.com/clchampag) - clara.champagne@swisstph.ch


## Credits
This readme has been built with Bilal BENHANA thanks to a template
- [makeread.me](https://github.com/ShaanCoding/ReadME-Generator)
- [othneildrew](https://github.com/othneildrew/Best-README-Template)