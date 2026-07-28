#!/bin/bash

########################################################################
# probmaps_gefs.sh and probmaps_gefs.py
#
# VERSION AND LAST UPDATE:
#   v1.0  06/09/2023
#   v1.1  01/04/2024
#
# PURPOSE:
#  Script to generate Probability Maps (Global Hazard Outlooks) of significant
#   wave height (Hs) and 10-m wind speed (WS10) for week-2 forecast based on
#   NOAA GEFSv12 wave ensemble forecast.
#  This shell script cheks forecast files exist, and runs the python script
#   probmaps_gefs.py reading the path where it is located in the configuration
#   file probmaps_gefs.yaml
#
# USAGE:
#  This code expects the GEFv12 files are downloaded (see download_GEFSwaves.sh),
#    the path where .grib2 files are located must be saved in the probmaps_gefs.yaml
#    file, variable gefspath
#  It reads probmaps_gefs.yaml, where the path+name of the python script
#   probmaps_gefs.py is located (users must edit variable pyscript).
#  The last configuration edit is the variable outpath in probmaps_gefs.yaml, which
#    is the location where the final plots will be saved. Once the .yaml
#    is configured, there is no need to change in the daily basis, unless
#    you want to modify the destination path or any other configuration.
#  Before running the python script, python must be loaded and activated. Please 
#    customize this part at the end (lines 102 and 103).
#
#  Example:
#    bash probmaps_gefs.sh /media/name/test/probmaps_gefs.yaml
#
# OUTPUT:
#  Figures containing the probability maps, saved in the directory
#    outpath informed in the probmaps_gefs.yaml file.
#
# DEPENDENCIES:
#  The python code probmaps_gefs.py contains the module dependencies.
#
# AUTHOR and DATE:
#  06/09/2023: Ricardo M. Campos, first version 
#  01/04/2024: Ricardo M. Campos, the download of GEFS was removed from here,
#    which is now download_GEFSwaves.sh
#
# PERSON OF CONTACT:
#  Ricardo M. Campos: ricardo.campos@noaa.gov
#
########################################################################

# INPUT ARGUMENT
# .yaml configuration file containing paths and information for this
#   shell script as well as for the python code.
PYCYAML="$1"
# PYCYAML="/media/ricardo/ssdrmc/analysis/products/probmaps/probmaps_gefs.yaml"

# Read the YAML as a text file:
#  GEFS data path
gefspath_line=$(grep 'gefspath' "${PYCYAML}")
GEFSMDIR=$(echo "$gefspath_line" | awk -F': ' '{print $2}')
#  Python script (probability maps)
pyscript_line=$(grep 'pyscript' "${PYCYAML}")
PYSCRIPT=$(echo "$pyscript_line" | awk -F': ' '{print $2}')
#  Variable names, for the python processing (probability maps)
mvars_line=$(grep 'mvars' "${PYCYAML}")
MVARS=$(echo "$mvars_line" | awk -F': ' '{gsub(/"/, "", $2); print $2}')
#  Output path
outpath_line=$(grep 'outpath' "${PYCYAML}")
OUTPATH=$(echo "$outpath_line" | awk -F': ' '{print $2}')

# Forecast cycle
# YEAR=`date +%Y`
# MONTH=`date +%m`
# DAY=`date +%d`

pa=0 #  days into the past. pa=1 runs using yesterday's cycle etc
YEAR=`date --date=-$pa' day' '+%Y'`
MONTH=`date --date=-$pa' day' '+%m'`
DAY=`date --date=-$pa' day' '+%d'`
HOUR="00" # first cycle 00Z

# Check GEFv12 is complete and ready.
# If not, it waits for 5 min and then try again (max 12 hours)
FSIZE=0
TRIES=1

while [ "$FSIZE" -lt 1000000 ] && [ "$TRIES" -le 144 ]; do

  # wait 5 minutes until next try
  if [ ${TRIES} -gt 5 ]; then
    sleep 300
  fi
  # Check if the last file (member 30, lead time 384h) is complete
  test -f $GEFSMDIR/GEFSv12Waves_$YEAR$MONTH$DAY$HOUR/gefs.wave.$YEAR$MONTH$DAY.30.global.0p25.f384.grib2
  TE=$?
  if [ ${TE} -eq 1 ]; then
    FSIZE=0
  else
    FSIZE=$(du -sb "$GEFSMDIR/GEFSv12Waves_$YEAR$MONTH$DAY$HOUR/gefs.wave.$YEAR$MONTH$DAY.30.global.0p25.f384.grib2" | awk '{print $1}')
  fi

  TRIES=`expr $TRIES + 1`

done

# Module load python and activate environment when necessary.
source /home/ricardo/work/python/anaconda3/setanaconda3.sh

echo "  "
echo " PYTHON PROCESSING: GLOBAL HAZARDS OUTLOOK - PROBABILITY MAPS, $YEAR$MONTH$DAY$HOUR "
echo "  "
# loop through variables
for WW3VAR in ${MVARS[*]}; do
  # 7 14 is the time intervall (days) for week 2
  python3 ${PYSCRIPT} ${PYCYAML} $YEAR$MONTH$DAY$HOUR 7 14 ${WW3VAR}
  echo " Probability maps for ${WW3VAR} Ok." 

done

echo "  "
echo " PYTHON PROCESSING COMPLETE."

# ----
cd ${OUTPATH}
mkdir $YEAR$MONTH$DAY
mkdir $YEAR$MONTH$DAY/Hs
mkdir $YEAR$MONTH$DAY/WS10
mv *Hs* $YEAR$MONTH$DAY/Hs/
mv *WS10* $YEAR$MONTH$DAY/WS10/

