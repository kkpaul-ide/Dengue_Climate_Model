# The implications of climate change on dengue virus transmission in Bangladesh

This repository contains all the code for climate-driven dengue transmission model, extracted data for parameter fitting, processed climate model output, model parameters, model simulation, and analysis of the transmission model output.

[![DOI](https://zenodo.org/badge/649274053.svg)](https://zenodo.org/badge/latestdoi/649274053)

All codes were created and run using R statistical software version 4.2.2 and Rstudio version 2022.02.3+492 "Prairie Trillium" release with associated packages (tidyverse 1.3.2; tidymodels 1.0.0; reshape2 1.4.4; minpack.lm 1.2-3; readr 2.1.4; lubridate 1.9.2; R.utils 2.12.2; cowplot 1.1.1). Some of the data files required to run the scrits are not publicly available and have been excluded from the repository (but can be provided on request). Please report an issue if you have any problems running the code or contact either Kpaul@kirby.unsw.edu.au or Rgray@kirby.unsw.edu.au.

# Aims

1. To develop and validate a mechanistic model of dengue transmission that can make long-term projections of dengue incidence and sero-prevalence over time in endemic regions, under different climate change scenarios.

2. To utilize the mechanistic model and project the changes in magnitude and seasonality of dengue outbreaks in Dhaka, the capital city of Bangladesh, throughout the 21st Century, under different climate change scenarios.


# Contributors
The main developers for the project code are Dr Kishor K. Paul and Dr Richard T. Gray who are in the Surveillance and Evaluation Research Program (SERP) at the The Kirby Institute, UNSW Sydney, Sydney, NSW, Australia. Kishor Paul was the primarily developer of the model, data analysis and calculation scripts under the supervision of Dr Gray. Both maintain the repository.

# Project organization

*The project is organized with materials in main directory and additional code and data in sub-directories*

* 0-mosquito_paramater_fitting.Rmd Script used to (re)fit a Briere function to published mosquito data to generate best estimates and uncertainty in key variables that are dependent on temperature.

* 0-mosquito_paramater_mortality_fitting.Rmd Script used to (re)fit a Quadratic function to published mosquito data to generate best estimates and uncertainty in key variables that are dependent on temperature, humidity, and rainfall.

* 0-setup_project.Rmd Script used to set-up a dengue modelling project.

* 1-setup_model.Rmd Script used to set-up specific dengue model parameters.

* 2-test_calibrate_model.Rmd Script used to simulate the specified dengue model with the saved input parameters and data and then calibrate parameters.

* 3-setup_sims.Rmd Script used to set-up simulations of a specific dengue model.

* 4-run_scenarios.Rmd Script used to run simulations of a specific dengue model for available scenarios.

* 5-observed_analysis.Rmd Script used to analyse model output for the period 1995-2014.

* 6-projections_analysis.Rmd Script used to analyse model output for 2030-2049 and 2080-2099.

## code
Contains the dengue transmission model coded in R and specific R functions and scripts of used for data processing, reshaping ISIMIP netCDF data files to run dengue transmission model and analysis of the model simulation output.

## data
Contains some raw data used in the analysis. Temperature output of ISIMIP models are available as netCDF files at https://www.isimip.org/outputdata/isimip-data-on-the-esgf-server/ . Data for weather station locations need to be extraced first with the help of Climate Data Operator (CDO) (https://doi.org/10.5281/zenodo.3539275), a collection of command line operators, to be able to use user defined functions and analysis scripts.

## outputs
Contains R data pertaining to project specifications and calibrated model inputs.

## templates
Contains parameter values, ranges and type available from literature that has been used to run the model. 

