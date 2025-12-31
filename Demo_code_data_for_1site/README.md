This **README** includes how to run the scripts in this demonstration folder and how to reproduce the results for the AmeriFlux site: US-Kon in a manuscript on thermal responses in ecosystem respiration. The scripts have been tested using R version 4.3.2 on the platform aarch64-apple-darwin20

**Step 1**: Download the entire project in this github repo. If you only download this folder, you need to change the location of "dir_rawdata" to where you places the folder.

**Step 2**: To run these R scripts, some R packages such as amerifluxr, REddyProc, brms, and cmdstanr may not be installed automatically. In this case, users need to install them manually using devtools.

-   How to install amerifluxr? use devtools::install_github("chuhousen/amerifluxr")

-   How to install cmdstanr and CmdStan? please refer to this page: <https://mc-stan.org/cmdstanr/articles/cmdstanr.html>

**Step 3**: Please run the five R scripts starting with D0 in the numerical order. For each script, the time needed to run is listed in the beginning.

Here is the **introduction** of the 5 scripts:

**D01**: Prepares data for fitting temperature-respiration curves for this site.

-   Exports 1 csv files containing gap-filled subhourly meteorological and NEE time series.
-   Exports 1 csv files containing measured high-quality subhourly nightlight NEE and corresponding meteorological data.
-   Exports 1 csv file containing the start and end date of growing seasons.

**D02a**: Estimates *total* thermal response strength in ecosystem respiration using moving window methods.

-   Exports 2 csv files (outcome_siteyear_temp_1site.csv and outcome_temp_1site.csv) containing estimated *total* thermal response strength (TAS) for this site.

**D02b**: Estimates *direct* thermal response strength in ecosystem respiration using moving window methods.

-   Exports 2 csv files (outcome_siteyear_temp_water_gpp_1site.csv and outcome_temp_water_gpp_1site.csv) containing estimated *direct* thermal response strength (TAS) for this site.

**D03**: Identifies drivers of direct and total thermal response strength using random forests.

-   Exports 2 csv files. One records relative importance of predictors (varImp_plot.csv). The other contains data for making partial dependence plot (partial_plot.csv).

**D04**: Estimates the effects of *total* thermal responses on future (2041-2060) ecosystem respiration this site.

-   Exports 1 csv file containing projected future nighttime ecosystem respiration, with and without considering thermal response strength (acclimation_data_future_ssp245_1site.csv).

**Note**: you may want to compare your results for this site with author's complete results for all sites in acclimation_data.csv and acclimation_data_future_ssp245.csv.

It takes \~1 hour to finish running the five scripts. Thank you!
