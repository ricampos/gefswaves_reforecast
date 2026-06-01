#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=128G
#SBATCH -q batch
#SBATCH -t 08:00:00
#SBATCH -A marine-cpu
#SBATCH -p orion

# This job run the code buildfuzzydataset.py on Orion for one month
#
# DIRJOUT="/work/noaa/marine/ricardo.campos/work/analysis/3assessments/fuzzy_verification/jobs"
# export YEAR=2023
# export MONTH=1
# sbatch --output=${DIRJOUT}/jbuildfuzzydataset_${YEAR}${MONTH}.out jbuildfuzzydataset.sh

echo " Starting at "$(date +"%T")

ulimit -s unlimited
ulimit -c 0

DIRSCRIPTS="/work/noaa/marine/ricardo.campos/work/analysis/3assessments/fuzzy_verification"

export YEAR=${YEAR}
export MONTH=${MONTH}

# python env
source /work/noaa/marine/ricardo.campos/progs/python/setanaconda3.sh
sh /work/noaa/marine/ricardo.campos/progs/python/setanaconda3.sh

echo "  "
echo " BUILD FUZZY DATASET"
echo "  "

# work dir
cd ${DIRSCRIPTS}

for DAY in `seq 1 31`; do

  export DAY=${DAY}
  DATE=${YEAR}$(printf "%02.f" $MONTH)$(printf "%02.f" $DAY)
  export DATE=${DATE}
  echo "  "
  echo " "${DATE}
  echo "  "
  # run
  python3 ${DIRSCRIPTS}/buildfuzzydataset.py ${DATE}
  wait $!
  echo " buildfuzzydataset.py for "${DATE}" OK at "$(date +"%T")
  echo "  "
  sleep 2

done

echo " Complete at "$(date +"%T")

