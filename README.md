# Strength and drivers of thermal responses in ecosystem respiration

## Project objectives:

**1.** Develop an algorithm to estimate site-level total and direct thermal response strength in ecosystem respiration (ER) using global eddy-covariance data.

**2.** Estimate total and direct thermal response strength at 110 terrestrial ecosystems.

**3.** Identify drivers of direct and total thermal response strength.

**4.** Predict how thermal responses will influence future warming-induced increases in ER.

## Structure of R Scripts in *workflows*:

1.  Data Preparation for Thermal Response Estimates

-   `01_01_estimate_soil_temperature_at_some_sites.R`:
    -   Estimates top soil (\<10 cm) temperature at sites with lots of missing soil temperature data.
    -   Exports 16 csv files containing predicted soil temperature for 16 sites with incomplete soil temperature measurements. csv files are saved as SITE_TS_rfp.csv, where SITE is the name of the site.
-   `01_02a_filter_high_quality_night_respiration_EuroFlux.R`:
    -   Prepares data for fitting temperature-respiration curves for 48 EuroFlux sites.
    -   Exports 48 csv files containing gap-filled subhourly meteorological and NEE time series (one file per site).
    -   Exports 48 csv files containing measured high-quality subhourly nightlight NEE and corresponding meteorological data (one file per site).
    -   Exports 1 csv file containing the start and end date of growing seasons at all EuroFlux sites.
-   `01_02b_filter_high_quality_night_respiration_AmeriFlux`:
    -   Prepares data for fitting temperature-respiration curves for 62 AmeriFlux sites.
    -   Exports 62 csv files containing gap-filled subhourly meteorological and NEE time series (one file per site).
    -   Exports 62 csv files containing measured high-quality subhourly nightlight NEE and corresponding meteorological data (one file per site).
    -   Exports 1 csv file containing the start and end date of growing seasons at all AmeriFlux sites.

2.  Calculate Thermal Response Strength (TAS) for All Sites

-   `02_01a_estimate_total_TAS_use_moving_window.R`:
    -   Estimates *total* thermal response strength in ecosystem respiration using moving window methods.
    -   Exports 2 csv files (outcome_siteyear_temp.csv and outcome_temp.csv) containing estimated *total* thermal response strength (TAS) for all sites.
-   `02_01b_estimate_direct_TAS_use_moving_window.R`:
    -   Estimates *direct* thermal response strength in ecosystem respiration using moving window methods.
    -   Exports 2 csv files (outcome_siteyear_temp_water_gpp.csv and outcome_temp_water_gpp.csv) containing estimated *direct* thermal response strength (TAS) for all sites.
-   `02_02_compare_different_TAS.R`:
    -   Explores how *total* and *direct* thermal response strengths vary with climate class and vegetation class.
    -   No export.

3.  Identify Drivers of Total and Direct Thermal Response Strength

-   `03_01_prepare_data_for_driver_analysis.R`:
    -   Obtains environmental and biological conditions potentially affecting thermal response strength
    -   The information includes climate, soil organic carbon, and remotely sensed GPP, LAI, and vegetation index.
    -   Exports 1 csv file (acclimation_data.csv) containing climate, soil, and vegetation conditions at all sites.
-   `03_02_identify_thermal_response_driver.R`:
    -   Identifies drivers of direct and total thermal response strength using random forests.
    -   Exports 2 csv files. One records relative importance of predictors (varImp_plot.csv). The other contains data for making partial dependence plot (partial_plot.csv).

4.  Calculate Effects of Total Thermal Response on Future Ecosystem Respiration

-   `04_01_get_future_night_soil_temperature_change.R`:
    -   Estimates monthly nighttime air temperature change from the current (2000-2020) to the future (2041-2060) under medium emission scenario.
    -   Exports 1 csv file containing projected future changes in monthly nighttime air temperature at all sites (Tmin_month_ssp245_wc.csv).
-   `04_02_thermal_response_effect_on_future_ecosystem_respiration.R`:
    -   Estimates the effects of *total* thermal responses on future (2041-2060) ecosystem respiration at all sites.
    -   Exports 1 csv file containing projected future nighttime ecosystem respiration, with and without considering thermal response strength (acclimation_data_future_ssp245.csv).

## How to reproduce the workflow:

1.  Download raw data required for this project.

-   1.1 Download raw eddy covariance data from FLUXNET2015, FLUXNET2020, FLUXNET07202025, and ICOS (Ecosystem final quality (L2) product in ETC-Archive format - release 2025-1)

-   1.2 Download ERA5 daily soil water content data in 1990-2014 for all sites from <https://cds.climate.copernicus.eu/>

-   1.3 Download global soil carbon stock map from soilwise: <https://repository.soilwise-he.eu/cat/collections/metadata:main/items/7730e747-eb73-49c9-bfe6-84ebae718743>

-   1.4 Download remotely sensed vegetation index (NDVI and EVI), gpp, lai, and fpar for all sites using site coordinates from NASA earth data: <https://www.earthdata.nasa.gov/data/tools/appeears>

-   1.5 Download BADM Data Product of all AmeriFlux sites: <https://ameriflux.lbl.gov/data/badm/>

    -   All the data are now stored on Yale data server.

2.  Use the *site_info.csv* file in *data* folder to get started.

-   This file includes basic information for each site, including lat, long, IGBP, climate class, availability of measured SWC, key variable names in each flux file
-   It is pre-required to run the workflow.

3.  Scripts in the *workflows* folder should be executed in numerical order. The scrips that share the same number but have different letters (e.g., 01_02a and 01_02b) can be run in parallel.

4.  In the beginning of each script, it states the approximate time for running that script. The time is tested on a Mac mini, so it may vary on different computers.
