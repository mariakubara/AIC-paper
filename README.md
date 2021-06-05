# AIC-paper

This repository contains codes for the article 'Akaike Information Criteria in testing optimal spatial neighbourhood' by Maria Kubara and Katarzyna Kopczewska.

## Structure of files

### data

Contains shapefiles (used for irregular bounding box specification) and empirical point data structure.

- "geoloc data.\*" - point data regarding companies in Lubelskie region
- "powiaty.\*" - shapefiles regarding NUTS2 region structure in Poland
- "wojewodztwa.\*" - shapefiles regarding NUTS4 region structure in Poland

### codes_main_paper

Contains R scripts used for generating analysis and simulations presented in the paper. 

- "data preparation and setting up the environment.R" - main initialising script necessary for running all simulation and visualisation codes
- "SimulationScriptA.R" - script for general simulation scheme used in the paper
- "semiVarianceKnn.R" - script with self-defined function allowing for semiVariance calculation based on incremental knn change
- "calculating semivariance for the simulation A.R" - script in which semivariance results are calculated to match the simulation of type A
- "AIC behaviour on empirical sample.R" - R script in which subsamples of the empirical dataset are created and AIC behaviour is analysed

### codes_additional_simulation

Contains R scripts used for generating additional simulations discussed in the appendix

- "SimulationScriptB.R" - script for simulation of type B
- "SimulationScriptC.R" - script for simulation of type C
- "SimulationScriptD.R" - script for simulation of type D
- "SimulationScriptE.R" - script for simulation of type E
- "SimulationScriptF.R" - script for simulation of type F

### codes_visualisation

Contains general R scripts used to prepare graphics for the paper

- "Main visualisation script.R" - general visualisation script which allows to prepare most of the plots used in the paper. Script is prepared in a flexible way, which allows to generate visualisations for different simulation calls, depending on the result data loaded.
- "StructureX_Plot with general sample structure in three sizes \*.R" - group of scripts which are necessary to create sample overview figures for each simulation type, additionally to the figure creation Clark-Evans statistics are calculated

### results_simulation

Contains csv files with results from the main simulation

- "simulationCluster.csv" - results from simulation A (main simulation presented in the paper)
- "semiVarianceCluster.csv" - values of semiVariance matching the structure of simulation A 
- "resultsEmpiricalSimulation.csv" - results from analysing AIC behaviour in subsamples of the empirical point dataset

### results_additional

Contains csv files with results from additional simulations discussed in the appendix

- "simulationUniform.csv" - results from simulation B
- "simulationUniformCluster25.csv" - results from simulation C
- "simulationTwoUniformCluster25.csv" - results from simulation D
- "simulationUniformTwoClaster25.csv" - results from simulation E
- "TwoCluster25.csv" - results from simulation F
