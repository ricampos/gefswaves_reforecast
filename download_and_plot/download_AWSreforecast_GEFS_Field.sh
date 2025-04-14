#!/bin/bash

########################################################################
# download_AWSreforecast_GEFS_Field.sh
#
# VERSION AND LAST UPDATE:
#   v1.0  08/25/2023
#
# PURPOSE:
#  Script to download the Wave Reforecast from AWS.
#   NOAA Global Ensemble Forecast System (GEFSv12) using WAVEWATCH III.
#   It downloads the control and all perturbed members of the ensemble reforecast.
#   Global wind and wave fields.
#  https://noaa-nws-gefswaves-reforecast-pds.s3.amazonaws.com/index.html
#
# USAGE:
#  Two input arguments, date and path, must be entered.
#  Example:
#    bash download_AWSreforecast_GEFS_Field.sh 20150823 /home/ricardo/data/gefs
#
# OUTPUT:
#  Multiple grib2 files, for each ensemble member.
#
# DEPENDENCIES:
#  wget
#
# AUTHOR and DATE:
#  08/25/2023: Ricardo M. Campos, first version 
#
# PERSON OF CONTACT:
#  Ricardo M. Campos: ricardo.campos@noaa.gov
#
#  If you are interested in operational forecasts from NOAA ftp, see:
#  https://github.com/NOAA-EMC/WW3-tools/tree/develop/opforecast
#
########################################################################

# Two input arguments
# date
CTIME="$1"
# destination path
DIRW="$2"
# ensemble members
ensblm="`seq -f "%02g" 0 1 4`" # 4 or 10 (on Wednesdays)
# server address
SERVER=https://noaa-nws-gefswaves-reforecast-pds.s3.amazonaws.com/index.html#
# https://noaa-nws-gefswaves-reforecast-pds.s3.amazonaws.com/index.html#GEFSv12/reforecast/2000/20000108/gridded/
# gefs.wave.20000108.p01.global.0p25.grib2

echo $ensblm
cd ${DIRW}

echo " ======== GEFS Forecast, AWS archive: ${CTIME} 00Z ========"
for e in $ensblm;do
  echo $e
  # size TAM and tries TRIES will control the process
  TAM=0
  TRIES=1

  while [ $TAM -lt 1000000000 ] && [ $TRIES -le 130 ]; do
    # sleep 5 minutes between attemps
    if [ ${TRIES} -gt 5 ]; then
      sleep 30
    fi

    if [ ${TAM} -lt 1000000000 ]; then
        # Main line, download
        if [ ${e} == "00" ]; then
          wget -l1 -H -t1 -nd -N -np -erobots=off --tries=3 ${SERVER}GEFSv12/reforecast/${CTIME:0:4}/${CTIME}/gridded/gefs.wave.${CTIME}.c${e}.global.0p25.grib2 -O $DIRW/gefs.wave.${CTIME}.c${e}.global.0p25.grib2 2>&1
          wait $!
        else
          wget -l1 -H -t1 -nd -N -np -erobots=off --tries=3 ${SERVER}GEFSv12/reforecast/${CTIME:0:4}/${CTIME}/gridded/gefs.wave.${CTIME}.p${e}.global.0p25.grib2 -O $DIRW/gefs.wave.${CTIME}.p${e}.global.0p25.grib2 2>&1
          wait $!
        fi

        # test if the downloaded file exists
        if [ ${e} == "00" ]; then
          test -f $DIRW/gefs.wave.${CTIME}.c${e}.global.0p25.grib2
        else
          test -f $DIRW/gefs.wave.${CTIME}.p${e}.global.0p25.grib2
        fi

        TE=$?
        if [ ${TE} -eq 1 ]; then
          TAM=0
        else
          # check size of each file
          if [ ${e} == "00" ]; then
            TAM=`du -sb $DIRW/gefs.wave.${CTIME}.c${e}.global.0p25.grib2 | awk '{ print $1 }'` 2>&1
          else
            TAM=`du -sb $DIRW/gefs.wave.${CTIME}.p${e}.global.0p25.grib2 | awk '{ print $1 }'` 2>&1
          fi
        fi

        if [ ${e} == "00" ]; then
          echo $DIRW/gefs.wave.${CTIME}.c${e}.global.0p25.grib2
        else
          echo $DIRW/gefs.wave.${CTIME}.p${e}.global.0p25.grib2
        fi

    fi

    TRIES=`expr $TRIES + 1`
  done
done
# remove empty data
find $DIR -empty -type f -delete
echo " "
echo " Done ${CTIME}."

