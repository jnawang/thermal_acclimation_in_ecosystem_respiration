# Strength and drivers of thermal responses in ecosystem respiration

## Project objectives:

**1.** Develop an algorithm to estimate site-level direct and total thermal response strength in ecosystem respiration (ER) using global eddy-covariance data.

**2.** Estimate direct and total thermal response strength at 110 terrestrial ecosystems.

**3.** Identify drivers of direct and total thermal response strength.

**4.** Predict how thermal responses will influence future warming-induced increases in ER.

## Structure of Rscripts:

01_data_preparation_for_thermal_response_estimate

-   01_01_estimate_soil_temperature_at_some_sites
-   01_02a_filter_high_quality_night_respiration_EuroFlux
-   01_02b_filter_high_quality_night_respiration_AmeriFlux

02_calculate_thermal_response_strength(TAS)\_for_all_sites

-   02_01a_estimate_total_TAS_use_moving_window
-   02_01b_estimate_direct_TAS_use_moving_window
-   02_02_compare_different_TAS

03_identify_drivers_of_thermal_response_strength

-   03_01_prepare_data_for_driver_analysis
    -   data include climate, primary productivity, and soil organic carbon
-   03_02_identify_thermal_response_driver
    -   use random forests

04_calculate_thermal_response_effects_on_future_respiration

-   04_01_get_future_night_soil_temperature_change
-   04_02_thermal_response_effect_on_future_ecosystem_respiration

## How to reproduce the workflow:

1.1 Download raw eddy covariance data from FLUXNET2015, FLUXNET2020, FLUXNET07202025, and ICOS (Ecosystem final quality (L2) product in ETC-Archive format - release 2025-1)

1.2 Download ERA5 daily soil water content data in 1990-2014 for all sites from <https://cds.climate.copernicus.eu/>

1.3 Download global soil carbon stock map from soilwise: <https://repository.soilwise-he.eu/cat/collections/metadata:main/items/7730e747-eb73-49c9-bfe6-84ebae718743>

1.4 Download remotely sensed vegetation index, gpp, lai, and fpar for all sites from NASA earth data: <https://www.earthdata.nasa.gov/data/tools/appeears>

1.5 Download BADM Data Product of all AmeriFlux sites: <https://ameriflux.lbl.gov/data/badm/>

-   all the data are now stored on Yale data server.

2.  Use the *site_info.csv* file in *data* folder to get started.

-   This file includes basic information for each site, including lat, long, IGBP, climate class, availability of measured SWC, key variable names in each flux file
-   It is pre-required to run the workflow.

3.  Scripts in the *workflows* folder should be executed in numerical order. The scrips that share the same number but have different letters (e.g., 01_02a and 01_02b) can be run in parallel.

4.  In the beginning of each script, it states the approximate time for running that script. The time is tested on a Mac mini, so it may vary on different computers.

## Output of the workflow:

1.  18 csv files containing predicted soil temperature for 18 sites with incomplete soil temperature measurements.

2.  110 csv files containing gap-filled subhourly meteorological and NEE time series (one file per site).

3.  110 csv files containing measured subhourly nightlight NEE and corresponding meteorological data (one file per site).

4.  2 csv files containing estimated thermal response strength (TAS): one for *direct* TAS and one for *total* TAS.

5.  1 csv file containing site-level climate, soil, and productivity characteristics used in the TAS driver analysis.

6.  2 csv files containing the drivers of estimated TAS: one for *direct* TAS and one for *total* TAS.

7.  1 csv file containing projected future changes in monthly nighttime soil temperature at each site.

8.  1 csv file containing projected future mean nighttime ecosystem respiration, with and without considering thermal response strength.
