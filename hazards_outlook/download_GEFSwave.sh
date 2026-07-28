#!/bin/bash

########################################################################
# download_GEFSwave.sh
#
# VERSION AND LAST UPDATE:
#   v1.0  02/15/2023
#   v1.1  04/30/2025
#   v1.2  07/22/2026 - Fixed critical bug: under `set -e`, a failed wget
#                       (404/timeout/reset) or a failed `test -f` check
#                       (file missing) would silently terminate the ENTIRE
#                       script mid-loop, skipping all remaining forecast
#                       hours and members. Both are now shielded so the
#                       script keeps going and simply retries/reports
#                       instead of dying partway through.
#                       Also reconciled the two mismatched size thresholds
#                       (10MB outer skip-check vs 1.2MB inner retry-check)
#                       to avoid redundant re-downloads, and added a
#                       missing-files summary at the end.
#   v1.3  07/28/2026 - Fixed silent pileup bug: a 404 (data aged off
#                       NOMADS, which only keeps ~2 days online) was being
#                       retried up to 130 times with sleeps, so a single
#                       stale file could burn over an hour. Combined with
#                       daily cron and no locking, stuck runs stacked on
#                       top of each other across multiple days, garbling
#                       the log and starving bandwidth. wget stderr is now
#                       captured; a 404 is treated as PERMANENT and skipped
#                       immediately instead of retried. Non-404 failures
#                       (timeouts, resets, connection errors) still get
#                       the full retry treatment since those are transient.
#
# PURPOSE:
#  Script to download NOAA Global Ensemble Forecast System (GEFS), Wave 
#   Forecast from WAVEWATCH III operational. Download from AWS archive
#   and save the grib2 files without any conversion or processing. It
#   includes the control and all perturbed members of the ensemble.
#   Global wind and wave fields.
#
# USAGE:
#  Two input arguments, date and path, must be entered.
#  Example:
#    bash download_GEFSwave.sh 20220823 00 /home/ricardo/data/gefs
#
# OUTPUT:
#  Multiple grib2 files, for each time step and ensemble member.
#
# DEPENDENCIES:
#  wget
#
# AUTHOR and DATE:
#  02/15/2023: Ricardo M. Campos, first version 
#  04/30/2025: Ricardo M. Campos, flexible cycle time
#  07/22/2026: Ricardo M. Campos, fixed set -e early-exit bug on failed downloads
#  07/28/2026: Ricardo M. Campos, fixed 404 retry pileup / added error logging
#
# PERSON OF CONTACT:
#  Ricardo M. Campos: ricardo.campos@noaa.gov
#
#  If you are interested in operational forecasts from NOAA ftp, see:
#  https://github.com/NOAA-EMC/WW3-tools/tree/develop/opforecast
#
########################################################################

set -euo pipefail
export USER_IS_ROOT=0
export MODULEPATH=/etc/scl/modulefiles:/apps/lmod/lmod/modulefiles/Core:/apps/modules/modulefiles/Linux:/apps/modules/modulefiles
source /apps/lmod/lmod/init/bash

# Two input arguments
# date
CTIME="$1"
# cycle 00,06,12,18
# HCYCLE="00"
HCYCLE="$2"
# destination path
DIRW="$3"
# server address
# SERVER=https://noaa-gefs-pds.s3.amazonaws.com/
# SERVER=ftp://ftpprd.ncep.noaa.gov/pub/data/nccf/com/gens/prod/
SERVER=https://nomads.ncep.noaa.gov/pub/data/nccf/com/gens/prod/
# ensemble members
ensblm="`seq -f "%02g" 0 1 30`"
# Forecast lead time (hours) to download
fleads="`seq -f "%03g" 0 6 384`"
# minimum acceptable size (bytes) for a completed grib2 file - used both
# for the outer "already downloaded" skip check and the inner retry loop,
# so the two are now consistent
MINSIZE=10000000
# scratch file for capturing wget's stderr so we can tell a permanent
# 404 (data aged off NOMADS) apart from a transient network hiccup
WGETERR=$(mktemp)
trap 'rm -f "$WGETERR"' EXIT
# track failures so we can report them at the end instead of finding out later
FAILED_FILES=()
# separately track permanent 404s (expected/benign - data no longer on
# NOMADS) vs other failures (worth investigating)
NOTFOUND_FILES=()
cd ${DIRW}
for h in $fleads;do
  echo " ======== GEFS Forecast, AWS archive: ${CTIME} ${HCYCLE}Z $h ========"
  for e in $ensblm;do
    echo $e
    FILE=$DIRW/gefs.wave.${CTIME}.${e}.global.0p25.f$(printf "%03.f" $h).grib2
    # Skip if file exists and is large enough
    if [ -f "$FILE" ]; then
      TAM=$(du -sb "$FILE" | awk '{ print $1 }')
      if [ "$TAM" -ge "$MINSIZE" ]; then
        echo "File $FILE already exists and is large enough. Skipping download."
        continue
      fi
    fi
    # size TAM and tries TRIES will control the process
    TAM=0
    TRIES=1
    while [ $TAM -lt $MINSIZE ] && [ $TRIES -le 130 ]; do
      # sleep 5 minutes between attemps
      if [ ${TRIES} -gt 5 ]; then
        sleep 30
      fi
      if [ ${TAM} -lt $MINSIZE ]; then
          # Main line, download
          # `|| true` prevents a non-zero wget exit (404/timeout/reset)
          # from killing the whole script under `set -e`
          # stderr is captured to $WGETERR so we can inspect it for a
          # permanent 404 right after
          if [ ${e} == "00" ]; then
            wget --wait=1 --random-wait --limit-rate=5m --continue -l1 -H -t1 -nd -N -np -erobots=off --tries=3 ${SERVER}gefs.${CTIME}/${HCYCLE}/wave/gridded/gefs.wave.t${HCYCLE}z.c${e}.global.0p25.f"$(printf "%03.f" $h)".grib2 -O $DIRW/gefs.wave.${CTIME}.${e}.global.0p25.f"$(printf "%03.f" $h)".grib2 > "$WGETERR" 2>&1 || true
          else
            wget --wait=1 --random-wait --limit-rate=5m --continue -l1 -H -t1 -nd -N -np -erobots=off --tries=3 ${SERVER}gefs.${CTIME}/${HCYCLE}/wave/gridded/gefs.wave.t${HCYCLE}z.p${e}.global.0p25.f"$(printf "%03.f" $h)".grib2 -O $DIRW/gefs.wave.${CTIME}.${e}.global.0p25.f"$(printf "%03.f" $h)".grib2 > "$WGETERR" 2>&1 || true
          fi
          cat "$WGETERR"
          # test if the downloaded file exists
          # rewritten so a missing file (test -f returns 1) does NOT
          # trigger `set -e` and kill the script
          if [ -f $DIRW/gefs.wave.${CTIME}.${e}.global.0p25.f"$(printf "%03.f" $h)".grib2 ]; then
            TE=0
          else
            TE=1
          fi
          if [ ${TE} -eq 1 ]; then
            TAM=0
          else
            # check size of each file
            TAM=`du -sb $DIRW/gefs.wave.${CTIME}.${e}.global.0p25.f"$(printf "%03.f" $h)".grib2 | awk '{ print $1 }'`
          fi
          echo $DIRW/gefs.wave.${CTIME}.${e}.global.0p25.f"$(printf "%03.f" $h)".grib2
          # If NOMADS returned a clean 404, the data has aged off the
          # server (NOMADS only keeps ~2 days online). Retrying up to
          # 130 times against a permanent 404 is what caused runs to
          # stack up under cron - short-circuit out of the retry loop
          # immediately instead.
          if [ $TAM -lt $MINSIZE ] && grep -q "404 Not Found" "$WGETERR"; then
            echo "PERMANENT 404 (data likely aged off NOMADS) for $FILE - not retrying."
            NOTFOUND_FILES+=("$FILE")
            TRIES=131
            TAM=0
            break
          fi
	  # sleep 2
      fi
      TRIES=`expr $TRIES + 1`
      sleep 1
    done
    # If we exhausted all tries (or hit a permanent 404) and still don't
    # have a valid file, log it and move on instead of failing the whole run
    if [ $TAM -lt $MINSIZE ]; then
      echo "WARNING: failed to download $FILE after ${TRIES} tries."
      FAILED_FILES+=("$FILE")
    fi
  done
done
echo " Done ${CTIME}."
if [ ${#NOTFOUND_FILES[@]} -gt 0 ]; then
  echo " "
  echo "==== SUMMARY: ${#NOTFOUND_FILES[@]} file(s) returned 404 (likely aged off NOMADS) ===="
  for f in "${NOTFOUND_FILES[@]}"; do
    echo "  $f"
  done
fi
if [ ${#FAILED_FILES[@]} -gt 0 ]; then
  echo " "
  echo "==== SUMMARY: ${#FAILED_FILES[@]} file(s) failed to download ===="
  for f in "${FAILED_FILES[@]}"; do
    echo "  $f"
  done
fi
