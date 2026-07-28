#!/bin/bash
# bash rdownload_GEFSwave.sh 00 1
#
# 07/28/2026: Ricardo M. Campos - added a lock file so a new cron-triggered
#   run cannot start while a previous one is still going. Without this,
#   a run stuck retrying stale/404 data could still be alive when the
#   next day's cron fired, and the two (or more) instances would run
#   concurrently, interleaving output in cron_log.txt and competing for
#   bandwidth. Also wrapped the actual download in `timeout` as a backstop
#   so a run can't linger indefinitely even if something unexpected stalls.

set -euo pipefail

# ---- prevent overlapping runs -------------------------------------------
LOCKFILE="/tmp/rdownload_GEFSwave.lock"
exec 200>"$LOCKFILE"
if ! flock -n 200; then
  echo "$(date): another instance of rdownload_GEFSwave.sh is already running - exiting."
  exit 1
fi
# --------------------------------------------------------------------------

export USER_IS_ROOT=0
export MODULEPATH=/etc/scl/modulefiles:/apps/lmod/lmod/modulefiles/Core:/apps/modules/modulefiles/Linux:/apps/modules/modulefiles
source /apps/lmod/lmod/init/bash
# cycle 00,06,12,18
# HCYCLE="00"
HCYCLE="$1"
# time lag in days (past days)
pa="$2"
# Date cycle for the download
YEAR=`date --date=-$pa' day' '+%Y'`
MONTH=`date --date=-$pa' day' '+%m'`
DAY=`date --date=-$pa' day' '+%d'`
# Path
DIRS="/scratch4/AOML/aoml-phod/Ricardo.Campos/data/archives/GEFS"
cd ${DIRS}
WTIME=${YEAR}`printf %2.2d $MONTH``printf %2.2d $DAY`
# Flat structure: all cycles live directly under DIRS, cycle hour folded
# into the directory name, e.g. GEFS/GEFSv12Waves_2026072100
DIRW=${DIRS}/GEFSv12Waves_${WTIME}${HCYCLE}
if [ ! -d "${DIRW}" ]; then
   mkdir -p "${DIRW}"
fi
# Wrapped in `timeout` as a backstop: if a run somehow can't finish a
# full cycle in 3 hours, kill it cleanly instead of letting it linger
# into the next scheduled cron trigger.
timeout 3h bash ${DIRS}/download_GEFSwave.sh ${WTIME} ${HCYCLE} ${DIRW}
DL_RC=$?
if [ $DL_RC -eq 124 ]; then
  echo "$(date): WARNING - download_GEFSwave.sh for ${WTIME}${HCYCLE} hit the 3h timeout and was killed."
fi

# clean old data
# Get the cutoff date in YYYYMMDD format
CUTOFFDAYS=16
CUTOFF=$(date -d "${CUTOFFDAYS} days ago" +%Y%m%d)
# Loop through directories matching pattern (now directly under DIRS,
# named GEFSv12Waves_YYYYMMDDHH)
for dir in "$DIRS"/GEFSv12Waves_*; do
    # Extract the date+cycle part from directory name
    BASENAME=$(basename "$dir")
    DIR_DATETIME=${BASENAME#GEFSv12Waves_}
    # Check it's a valid 10-digit YYYYMMDDHH
    if [[ $DIR_DATETIME =~ ^[0-9]{10}$ ]]; then
        # Compare just the YYYYMMDD portion (first 8 digits) against cutoff
        DIR_DATE=${DIR_DATETIME:0:8}
        if [[ $DIR_DATE -lt $CUTOFF ]]; then
            echo "Deleting $dir (older than ${CUTOFFDAYS} days)"
            rm -rf "$dir"
        fi
    fi
done

