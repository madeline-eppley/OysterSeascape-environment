# Eastern oyster seascape environmental data
##### last updated by MGE on 5/11/2026


### Methods
Temperature and salinity data for each sampling location (Supp. Table S2) was downloaded from sources with routine monitoring (e.g., National Estuarine Research Reserve System, National Parks Service) (Helmuth et al., 2006; Nadeau et al., 2017). Raw environmental data were converted into a standardized format using a custom R script (Data Accessibility).

Sampling locations were retained in genetic-environment association analysis if their data were high-quality by meeting the following minimum requirements: 1) at least 1 year of measurements, 2) a minimum of 6 months of data per year, with at least 1 month capturing seasonal highs and lows, 3) at least 1 measurement per day, and 4) recent data collected within the 5 years prior to sampling. Out of the 40 sampling locations, 32 locations met these quality control measures. Extreme values were filtered out to remove data errors (valid ranges: 0-40 ppt for salinity, 0-42℃ for temperature). For each sampling location, we compiled the mean, maximum, and minimum annual values, standard deviation, and quantiles (10th - low and 90th - high percentiles) for salinity and temperature, hereafter referred to as SalinityQ90, SalinityQ10, TemperatureQ90, and TemperatureQ10. Variables with significant pairwise correlations (r > 0.8) were removed. Quantiles and standard deviations were used in downstream analysis.

