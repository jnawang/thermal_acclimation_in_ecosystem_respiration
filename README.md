# Strength and drivers of thermal acclimation in ecosystem respiration

## Project objectives:

**1.** Develop an algorithm to estimate thermal acclimation strength in ecosystem respiration (ER) using eddy-covariance data.

**2.** Estimate thermal acclimation strength at 93 global terrestrial ecosystems.

**3.** Identify drivers of thermal acclimation strength.

**4.** Predict consequences of thermal acclimation on mediating future warming-induced ER increases.

## Structure of Rscripts:

01_data_preparation_for_thermal_acclimation_estimate

-   01_01_estimate_soil_temperature_at_some_sites
-   01_02_filter_high_quality_night_respiration_measurement
-   01_03_estimate_gross_primary_productivity (GPP) for each day at each study sites
    -   hourly daytime NEE, hourly GPP
    -   biweekly vegetation index: EVI and NDVI
    -   daily SIF for some sites
-   01_04_download_measured_SIF(solar-induced chlorophyll fluorescence)\_LAI_data_from_NASA
-   01_05_develop_different_respiration_models
-   01_06_get_annual_temperature_respiration_curves_for_acclimation_strength_estimate

02_calculate_thermal_acclimation_for_all_sites

-   02_01_get_annual_TS_ER_curves_for_each_site
-   02_02_estimate_thermal_acclimation_strenght_for_each_site

03_identify_drivers_thermal_acclimation_strength

04_calculate_thermal_acclimation_effects_on_future_respiration

05_supplementary_scripts_address_reviewer_comments

-   05_01_why not using air temperature instead of soil temperature?
-   05_02_performance difference of different ER models?
    -   only temperature
    -   temperature + moisture
    -   temperature + moisture + GPP
-   05_03_how about use an average year or the year whose temperature close to MAT as control conditions?
