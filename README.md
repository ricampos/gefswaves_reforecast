# NOAA Wave Ensemble Reforecast, Week2 Probabilistic Forecast, and Validation
#

AWS data storage noaa-nws-gefswaves-reforecast-pds

This is a project in cooperation with the National Weather Service’s Ocean Prediction Center and Environmental Modeling Center, along with NOAA’s Atlantic Oceanographic and Meteorological Laboratory (AOML) and the University of Miami’s Cooperative Institute for Marine and Atmospheric Studies (CIMAS).

# Description:

This repository is designed to help users to download and process the 20-year global wave reforecast generated 
 by WAVEWATCH III model (https://github.com/NOAA-EMC/WW3) forced by GEFSv12 winds (https://noaa-gefs-retrospective.s3.amazonaws.com/index.html).
The wave ensemble has been run with one cycle per day (at 03Z), spatial resolution of 0.25°X0.25° and temporal resolution of 3 hours. 
The reforecast ensemble comprises 5 members with one cycle per day, providing a forecast range of 16 days. 
Once a week (on Wednesdays), it extends to 35 days with 11 members.
More information about the wave modeling, wave grids and calibration can be found in the WAVEWATCH III regtest ww3_ufs1.3:
https://github.com/NOAA-EMC/WW3/tree/develop/regtests/ww3_ufs1.3
The 20 years of reforecast results have been analyzed and quality-controlled. Three output types are available:
1) Global wave fields, in grib2 format, with several variables including significant wave height, period, direction, and partitions;
2) Point output tables, in netcdf format, containing time-series of significant wave height, period and direction, for 658 points (latitude/longitude informed) at the positions of wave buoys;
3) For the same positions, spectral outputs are available, in netcdf format, containing the full spectra (2D directional spectrum).
Each file refers to one forecast cycle with date (year,month,day) written in the file name.

Details on the development and validation of this product can be found in the publication: \
Campos, R.M., Abdolali, A., Alves, J.H., Masarik, M., Meixner, J., Banihashemi, S., Mehra, A., Figurskey, D., Sienkiewicz, J., Lumpkin, R., 2024. Development and Validation of NOAA’s 20-year global wave ensemble reforecast. Weather & Forecasting. Under review.

This repository also contains scripts for the Week 2 Probabilistic Wave Forecast:
https://ocean.weather.gov/week2/
And validation and Fuzzy Verification of the product. The methodology is explained in details at:
https://journals.ametsoc.org/view/journals/wefo/aop/WAF-D-24-0154.1/WAF-D-24-0154.1.xml

# Update Frequency:

Quarterly

# License:

Open Data. There are no restrictions on the use of this data.

# Contact:

For questions related to wave modeling and the ensemble reforecast available, please contact Ricardo Martins Campos, email ricardo.campos@noaa.gov.

# Usage Examples:
