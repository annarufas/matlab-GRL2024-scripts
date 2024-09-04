# Can we Constrain Geographical Variability in the Biological Carbon Pump's Transfer Efficiency from Observations?

This repository contains the MATLAB scripts and datasets used for the data analysis and figure generation in the study: 

**"Can We Constrain Geographical Variability in the Biological Carbon Pump's Transfer Efficiency from Observations?"**

A. Rufas<sup>1</sup>, S. Khatiwala<sup>1</sup>, K. M. Bisson<sup>2</sup>, A. P. Martin<sup>3</sup>, H. A. Bouman<sup>1</sup>

<sup>1</sup>Department of Earth Sciences, University of Oxford, Oxford, UK  
<sup>2</sup>NASA Headquarters Science Mission Directorate, Santa Barbara, CA, US  
<sup>3</sup>National Oceanography Centre, Southampton, UK

Find the pre-print in the [Earth and Space Science Open Archive](https://essopenarchive.org/users/806280/articles/1197117-can-we-constrain-geographical-variability-in-the-biological-carbon-pump-s-transfer-efficiency-from-observations).

## Requirements

To use the content of this repository, ensure you have the following.

- MATLAB version R2021a or later installed. 
- Third-party functions downloaded from [MATLAB'S File Exchange](https://mathworks.com/matlabcentral/fileexchange/): `worstcase`, `MCErrorPropagation`, `swtest.m`, `FMINSEARCHBND`, `m_map`, `brewermap` and `longhurst_v4_2010`. Once downloaded, please place the functions in the `./resources/external/` directory.
- MATLAB toolboxes: the Optimization toolbox, necessary to run `worstcase`, and the Statistics and Machine Learning toolbox, necessary to run `MCErrorPropagation` and `swtest.m`.

## Data used in this repository

This project uses various raw oceanographic datasets downloaded from online repositories, stored in the `./data/raw/` folder.

- Particle concentration measurements from the Underwater Vision Profiler 5 (UVP5) downloaded from the Ecotaxa repository hosted by [IFREMER](https://ecopart.obs-vlfr.fr). This dataset is stored in the subfolder `./data/raw/UVP5`.
- Particulate organic carbon (POC) flux measurements from sediment traps and radionuclides, a compilation that we made for this study. The file, `dataset_s0_trap_and_radionuclide_compilation.xlsx`, is not available within this repository as it contains data owned by other authors that are not in a preservation repository. References for constructing this dataset are provided in the Supporting Information of our paper, Tables S1-S6.

Due to their large size, the following datasets are not included in `./data/raw/` but instructions for acquiring them are provided. Place the corresponding `.mat` files in `./data/raw/`. 

- The World Ocean Atlas 2018 monthly climatology for **temperature**, which you can download from the [NOAA website](https://accession.nodc.noaa.gov/NCEI-WOA18). Once downloaded, read all the .nc files and store them in a single array and save it as `woa18_global.mat`. This file is necessary to run the Marsay et al. (2015) algorithm.
- The AVHRR Pathfinder v.5.0 global 4 km monthly climatology (1985--2001) for **sea surface temperature**, which you can download from the [NOAA website](https://www.ncei.noaa.gov/archive/accession/AVHRR_Pathfinder-NODC-v5.0-climatologies). Once downloaded, read all the .nc files and store them in a single array and save it as `sst_global.mat`. This file is necessary to run the Henson et al. (2012) algorithm.
- The SeaWiFS 9 km monthly climatology (1997--2010) for **chlorophyll concentration** and **photosynthetically available radiation**, which you can downlaod from the [NASA website](https://oceancolor.gsfc.nasa.gov/l3/). Once downloaded, read all the .nc files and store them in a single array and save it as `chla_seawifs_global.mat` and `par0_seawifs_global.mat` respectively. These files are necessary to run the Henson et al. (2012) algorithm.

## Available MATLAB scripts in the folder `./code`

The following scripts have been run in this order to analyse the data and reproduce the figures in our paper.

| Num| Script name                                  | Script action                                    |
|----|----------------------------------------------|--------------------------------------------------
| 1  | plotGlobalMapWithTimeseriesStations.m        | Creates Figure 1                                 |
| 2  | processPocFluxFromTrapAndRadCompilation.m    | Analyses Dataset S0                              |
| 3  | plotPocFluxFromTrapAndRadCompilation.m       | Creates Figure 2, S1 and S2                      |
| 4  | processPocFluxFromUvp.m                      | Anayses the UVP5 dataset downloaded from Ecotaxa |
| 5  | plotUvpDataset.m                             | Creates Figure S3                                | 
| 6  | findAndPlotUvpVsObsMatchups.m                | Creates Figure 3                                 |
| 7  | processBcpMetrics.m                          | Calculates b, z* and Teff                        |
| 8  | plotBcpMetrics.m                             | Creates Figure 4                                 |
| 9  | calculateBcpMetricsFromTrapAndRadionuclide.m | Called by script num. 7                          |
| 10 | calculateBcpMetricsFromUvp.m                 | Called by script num. 7                          |
| 11 | calculateBcpMetricsFromHenson2012.m          | Called by script num. 7                          |
| 12 | calculateBcpMetricsFromMarsay2015.m          | Called by script num. 7                          |
| 13 | extractHensonForcingDataForChosenLocations.m | Called by script num. 11                         |
| 14 | extractMarsayForcingDataForChosenLocations.m | Called by script num. 12                         |
| 15 | propagateErrorWithMCforMartinbAndZstar.m     | Called by script num. 9 and 10                   |
| 16 | propagateErrorWithMCforPEeffAndTeff.m        | Called by script num. 9 and 10                   |
| 17 | defineFittypesForPocFluxAttenuationCurves.m  | Called by script num. 15                         |



## Reproducibility

Our scripts showcase the application of Monte Carlo-based techniques for error propagation across diverse datasets of biological carbon pump (BCP) mesopelagic transfer efficiency metrics. Specifically, the scripts are hard-wired to handle data from our six designated study sites:
- the Hawaii Ocean Time-series (HOT) station ALOHA (HOT/ALOHA), in the subtropical NE Pacific (22.45ºN, 158ºW);
- the Bermuda Atlantic Time-Series/Oceanic Flux Program joint site (BATS/OFP), in the subtropical NW Atlantic (31.6ºN, 64.2ºW);
- the US JGOFS Equatorial Pacific process study experimental site (EqPac), in the central equatorial Pacific upwelling system (–2 to 2ºN, 140ºW);
- the Porcupine Abyssal Plain time-Series Observatory (PAP-SO), in the subpolar NE Atlantic (49.0ºN, 16.5ºW);
- Ocean Station Papa (OSP), in the HNLC region of the subpolar NE Pacific (50ºN, 145ºW), and
- the Long-Term Ecological Research (LTER) observatory HAUSGARTEN, in the polar Atlantic-Arctic boundary (79ºN, 4.0ºE).

## Author

* [**Anna Rufas**](mailto://Anna.RufasBlanco@earth.ox.ac.uk)

## Acknowledgments

This work was completed as part of my PhD project at the University of Oxford under the NERC large grant COMICS (Controls over Ocean Mesopelagic Interior Carbon Storage, NE/M020835/2). I also acknowledge funding from the University of Oxford's Covid-19 Scholarship Extension Fund and Oxford's Wolfson College Covid-19 Hardship Fund.