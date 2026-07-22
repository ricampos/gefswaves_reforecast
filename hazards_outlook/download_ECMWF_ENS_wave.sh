#!/bin/bash

########################################################################
# download_ECMWF_ENS_wave.sh
#
# VERSION AND LAST UPDATE:
#   v1.0  06/13/2025
#   v1.1  07/22/2026 - Updated for IFS Cycle 50r1 (effective 2026-05-12):
#                       stream=waef/type=cf is deprecated. The wave
#                       ensemble control forecast (member 0) is now
#                       retrieved via stream=wave/type=fc instead.
#                       Perturbed members (1-50) are unchanged:
#                       stream=waef/type=pf.
#                       Script intended for use with 00Z and 12Z cycles only.
#
# PURPOSE:
#  Script to download ECMWF Ensemble Forecast Data from OpenData
#  https://www.ecmwf.int/en/forecasts/datasets/open-data
#  using python https://pypi.org/project/ecmwf-opendata/
#  https://www.ecmwf.int/sites/default/files/elibrary/2021/81274-ifs-documentation-cy47r3-part-vii-ecmwf-wave-model_1.pdf
#
# USAGE:
#  Two input arguments must be entered: cycle time (00 or 12) and the output path.
#  Example:
#    bash download_ECMWF_ENS_wave.sh 00 /home/ricardo/data/ECMWF_ENS
#  it will download the current day (00Z cycle)
#
# OUTPUT:
#  Multiple netcdf files containing the forecast steps, one file per ensemble member.
#
# DEPENDENCIES:
#  wget, python (ecmwf.opendata), CDO, NCO
#  user must activate python/environment (see source below)
#
# AUTHOR and DATE:
#  06/13/2025: Ricardo M. Campos, first version
#  07/22/2026: Updated for IFS Cycle 50r1 open-data changes
#
# PERSON OF CONTACT:
#  Ricardo M. Campos: ricardo.campos@noaa.gov
#
########################################################################

set -euo pipefail
export USER_IS_ROOT=0
export MODULEPATH=/etc/scl/modulefiles:/apps/lmod/lmod/modulefiles/Core:/apps/modules/modulefiles/Linux:/apps/modules/modulefiles
source /apps/lmod/lmod/init/bash
module load cdo
module load nco

# Directory, cycle time, and python environment setup
CHOUR="$1"
CHOUR=$(printf "%02.f" $CHOUR)
DIR="$2"
# Forecast date
# DATE=$(date '+%Y%m%d')
pa=1
DATE=$(date --date="-${pa} day" '+%Y%m%d')
exec > >(tee -a "$DIR/download_ECMWF_ENS_wave_${DATE}${CHOUR}.log") 2>&1
cd "$DIR"
if [ ! -d "$DIR/work_${DATE}${CHOUR}" ]; then
  mkdir "$DIR/work_${DATE}${CHOUR}"
fi
cd "$DIR/work_${DATE}${CHOUR}"
source /home/Ricardo.Campos/python/envs/intelpy_env/bin/activate
# Function to convert and compress GRIB2 to NetCDF
function ccompress() {
  local arqn=$1
  local DIR=$2
  cd "$DIR/work_${DATE}${CHOUR}"
  cdo -f nc4 copy "${arqn}.grib2" "${arqn}.saux1.nc"
  ncks -4 -L 1 -d lat,${latmin},${latmax} "${arqn}.saux1.nc" "${arqn}.saux2.nc"
  ncks --ppc default=.${dp} "${arqn}.saux2.nc" "${arqn}.nc"
  rm -f "${arqn}".saux* "${arqn}".*idx* "${arqn}".*ncks* "${arqn}".*tmp "${arqn}.grib2" "${sname}"
  echo "File ${arqn} converted to NetCDF and compressed."
}
# Variables
wparam_wave="swh/mwd/pp1d/mwp"
# restricting domain
latmin=-82.
latmax=89.
# decimals/resolution for compression
dp=2
# Lead times: 3h to 144h, then 6h to 360h (valid for 00Z/12Z cycles)
fleads=$(echo $(seq -f "%g" 0 3 144) $(seq -f "%g" 150 6 360) | tr ' ' '/')
# Ensemble members
ensblm=$(seq -f "%02g" 0 1 50)
# Loop over members
for e in $ensblm; do
  member_num=$((10#$e))
  arqn="ECMWF_ENS_wave_${DATE}${CHOUR}.${e}"
  FILE="${DIR}/${arqn}.nc"
  # Skip if existing file is large enough
  if [[ -f "$FILE" ]]; then
    size=$(du -sb "$FILE" | awk '{print $1}')
    if (( size >= 210000000 )); then
      echo " $FILE exists and is large enough. Skipping."
      continue
    else
      rm -f "$FILE"
    fi
  fi
  sname="download_ECMWF_ENS_wave_m${e}.py"
  echo " Running ${sname} ..."
  # IFS Cycle 50r1 (2026-05-12): stream=waef/type=cf is deprecated.
  # Member 0 (control) now comes from stream=wave/type=fc.
  # Members 1-50 (perturbed) remain stream=waef/type=pf.
  if (( member_num == 0 )); then
    estream="wave"
    etype="fc"
    number_line=""
  else
    estream="waef"
    etype="pf"
    number_line="        \"number\": \"${member_num}\","
  fi
  # Create download script
  cat > "${sname}" << EOF
#!/usr/bin/env python3
from ecmwf.opendata import Client
client = Client()
client.retrieve({
    "class": "od",
    "date": "${DATE}",
    "time": "${CHOUR}",
    "stream": "${estream}",
    "type": "${etype}",
${number_line}
    "step": "${fleads}",
    "levtype": "sfc",
    "param": "${wparam_wave}",
    "target": "${arqn}.grib2"
})
EOF
  rm -f "${arqn}".saux* "${arqn}".*idx* "${arqn}".*ncks* "${arqn}".*tmp "${arqn}.grib2"
  chmod +x "${sname}"
  python3 "${sname}"
  ccompress "${arqn}" "${DIR}"
  chmod 775 "${arqn}.nc"
  mv "${arqn}.nc" ${DIR}
done
sleep 1
cd $DIR
rm -rf "work_${DATE}${CHOUR}"
echo " "
echo " Done download_ECMWF_ENS_wave.sh ${DATE} ${CHOUR}Z "

