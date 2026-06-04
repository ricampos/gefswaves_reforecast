#!/bin/bash

DIRJOUT="/work/noaa/marine/ricardo.campos/work/analysis/3assessments/fuzzy_verification/jobs"
DIRSCRIPTS="/work/noaa/marine/ricardo.campos/work/analysis/3assessments/fuzzy_verification"

cd ${DIRSCRIPTS}

for NOCEAN in Atlantic Pacific; do
  # total_lines=$(wc -l < "groups_"${NOCEAN}".txt")
  total_lines=$(wc -l < "groups_seasonal_"${NOCEAN}".txt")
  for ((GRPID=0; GRPID<total_lines; GRPID++)); do
    for WVAR in hs u10; do
      export NOCEAN=${NOCEAN}
      export GRPID=${GRPID}
      export WVAR=${WVAR}
      sbatch --output=${DIRJOUT}"/jfuzzy_verification_GEFS_"${NOCEAN}"_"${GRPID}"_"${WVAR}".out" ${DIRSCRIPTS}"/jfuzzy_verification_GEFS.sh"
      echo " job jfuzzy_verification_GEFS_"${NOCEAN}"_"${GRPID}"_"${WVAR}" submitted OK at "$(date +"%T")
      sleep 1
    done
  done
done

