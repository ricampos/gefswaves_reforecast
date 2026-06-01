#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=128G
#SBATCH -q batch
#SBATCH -t 08:00:00
#SBATCH -A marine-cpu
#SBATCH -p orion

# This job script runs fuzzy_verification_GEFS.py 

echo " Starting at "$(date +"%T")

ulimit -s unlimited
ulimit -c 0

DIRSCRIPTS="/work/noaa/marine/ricardo.campos/work/analysis/3assessments/fuzzy_verification"

export NOCEAN=${NOCEAN}
export GRPID=${GRPID}
export WVAR=${WVAR}

# python env
source /work/noaa/marine/ricardo.campos/progs/python/setanaconda3.sh
sh /work/noaa/marine/ricardo.campos/progs/python/setanaconda3.sh

echo "  "
echo " Fuzzy Verification ${NOCEAN} ${GRPID} ${WVAR}"
echo "  "

# work dir
cd ${DIRSCRIPTS}
# Run
python3 ${DIRSCRIPTS}/fuzzy_verification_GEFS.py ${NOCEAN} ${GRPID} ${WVAR}
wait $!
echo " fuzzy_verification_GEFS.py for ${NOCEAN} ${GRPID} ${WVAR} OK at $(date +"%T")"
echo " "
sleep 1
echo " Complete at "$(date +"%T")

