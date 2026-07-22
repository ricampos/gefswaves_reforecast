#!/bin/bash

########################################################################
# download_GEWPS_wave.sh
#
# VERSION AND LAST UPDATE:
#   v1.0  06/18/2025
#
# PURPOSE:
#  Script to download Environment Canada's wave ensemble forecast
#   https://dd.meteo.gc.ca/
#   https://weather.gc.ca/ensemble/index_e.html
#   https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_table4-2-10-0.shtml
#
# USAGE:
#  Two input arguments must be entered: cycle time (0 or 12) and the output path.
#  Example:
#    bash download_GEWPS_wave.sh 00 /home/ricardo/data/EnvCanada/wave
#  it will download the current day (00Z cycle)
#
# OUTPUT:
#  Multiple netcdf files containing the forecast steps, one file per ensemble member.
#
# DEPENDENCIES:
#  wget, CDO, NCO
#
# AUTHOR and DATE:
#  06/18/2025: Ricardo M. Campos, first version 
#
# PERSON OF CONTACT:
#  Ricardo M. Campos: ricardo.campos@noaa.gov
#
########################################################################

set -euo pipefail

# When working on the cluster
# export USER_IS_ROOT=0
# export MODULEPATH=/etc/scl/modulefiles:/apps/lmod/lmod/modulefiles/Core:/apps/modules/modulefiles/Linux:/apps/modules/modulefiles
# source /apps/lmod/lmod/init/bash
# module load cdo
# module load nco

# Function to check if file has been previously downloaded
checkfile() {
  local FILE="$1"

  if [[ -f "$FILE" ]]; then
    local TAM
    TAM=$(stat -c%s "$FILE" 2>/dev/null)

    if [[ "$TAM" -ge 230000000 ]]; then
      echo "0"  # File exists and is big enough
      return
    fi
  fi

  echo "1"  # File missing or too small
}

# Directory and cycle time
CHOUR="$1"
CHOUR=$(printf "%02.f" $CHOUR)
WDIR="$2"
# Forecast cycle date
# DATE=$(date '+%Y%m%d')
pa=1
DATE=$(date --date="-${pa} day" '+%Y%m%d')
exec > >(tee -a "$WDIR/download_GEWPS_wave_${DATE}${CHOUR}.log") 2>&1

# Wave variables: Hs, Dp, Tp, Tm
WVARS=("HTSGW" "PWAVEDIR" "PWPER" "MZWPER")
# corresponding variable names after netcdf conversion
wvarname=("swh" "param46.0.10" "pp1d" "mp2")

# Ensemble members
ensbl=($(seq 1 21))
ensblm=($(seq -f "%02g" 0 1 20))

cd $WDIR
if [ ! -d "$WDIR/work_${DATE}${CHOUR}" ]; then
  mkdir "$WDIR/work_${DATE}${CHOUR}"
fi
cd "$WDIR/work_${DATE}${CHOUR}"

# EnvCanada server address
SERVER=https://dd.meteo.gc.ca
# Forecast lead time
fleads=$(echo $(seq -f "%03g" 0 3 168) $(seq -f "%03g" 174 6 384))

# restricting domain
latmin=-82.
latmax=89.
# decimals/resolution for compression
dp=2

# Dowload loop
for h in $fleads;do

  for i in "${!WVARS[@]}"; do

    wv="${WVARS[$i]}"
    wvarn="${wvarname[$i]}"

    echo " ======== GEWPS Forecast: ${DATE} ${CHOUR}Z $h ${wv} ========"

    arqn="gewps.wave.${DATE}T${CHOUR}Z_"${h}"H_${wv}"

    FILE="$WDIR/GEWPS_wave_${DATE}${CHOUR}.00.nc"
    cfile=$(checkfile "$FILE")
    if [ "$cfile" -eq 0 ]; then
      echo "File $FILE already exists and is large enough. Skipping download 1."
      continue
    fi

    wfile="${DATE}T${CHOUR}Z_MSC_GEWPS_${wv}_Sfc_LatLon0.25_PT"${h}"H.grib2"  
    wget -l1 -H -t1 -nd -N -np -erobots=off --tries=3 "${SERVER}/${DATE}/WXO-DD/model_gewps/25km/${CHOUR}/${wfile}" -O "$WDIR/work_${DATE}${CHOUR}/${arqn}.grib2"
    echo "wget $WDIR/work_${DATE}${CHOUR}/${arqn}.grib2"

    test -f "$WDIR/work_${DATE}${CHOUR}/${arqn}.grib2"
    TE=$?
    if [ ${TE} -eq 0 ]; then

      # wgrib2 "${arqn}.grib2" -netcdf ${arqn}.saux.nc
      cdo -f nc4 copy "${arqn}.grib2" "${arqn}.saux.nc"

      for j in "${!ensbl[@]}"; do
        e="${ensbl[$j]}"
        ne="${ensblm[$j]}"

        if [ "$j" -eq 0 ]; then
          ncks -v ${wvarn} ${arqn}.saux.nc -o ${arqn}_${ne}.saux.nc
          cdo sellonlatbox,-180,180,-90,90 ${arqn}_${ne}.saux.nc ${arqn}_${ne}.saux1.nc
          ncrename -v ${wvarn},${wv} ${arqn}_${ne}.saux1.nc ${arqn}_${ne}.nc
        else
          ncks -v "${wvarn}_$e" ${arqn}.saux.nc -o ${arqn}_${ne}.saux.nc
          cdo sellonlatbox,-180,180,-90,90 ${arqn}_${ne}.saux.nc ${arqn}_${ne}.saux1.nc
          ncrename -v ${wvarn}_${e},${wv} ${arqn}_${ne}.saux1.nc ${arqn}_${ne}.nc
        fi

        rm -f ${arqn}_${ne}.saux*.nc

      done

      rm -f $arqn.grib2
      rm -f $arqn.saux.nc

    else
      exit 1
    fi

  done

  for i in "${!ensbl[@]}"; do
    ne="${ensblm[$i]}"

    FILE="$WDIR/GEWPS_wave_${DATE}${CHOUR}.${ne}.nc"
    cfile=$(checkfile "$FILE")
    if [ "$cfile" -eq 0 ]; then
      echo "File $FILE already exists and is large enough. Skipping download 2."
      continue
    fi

    for j in "${!WVARS[@]}"; do
      wv="${WVARS[$j]}"
      if [ "$j" -eq 0 ]; then
        iname="gewps.wave.${DATE}T${CHOUR}Z_${h}H_${ne}.nc"
        mv "gewps.wave.${DATE}T${CHOUR}Z_${h}H_${wv}_${ne}.nc" "$iname"
      else
        src="gewps.wave.${DATE}T${CHOUR}Z_${h}H_${wv}_${ne}.nc"
        ncks -A -v $wv "$src" "$iname"
        rm -f $src
      fi
    done

  done

done


for i in "${!ensbl[@]}"; do

  ne="${ensblm[$i]}"

  FILE="$WDIR/GEWPS_wave_${DATE}${CHOUR}.${ne}.nc"
  cfile=$(checkfile "$FILE")
  if [ "$cfile" -eq 0 ]; then
    echo "File $FILE already exists and is large enough. Skipping download 3."
    continue
  fi

  arqn="GEWPS_wave_"${DATE}${CHOUR}.${ne}
  ncrcat gewps.wave.${DATE}T${CHOUR}Z_*H_${ne}.nc ${arqn}.saux.nc
  rm -f gewps.wave.${DATE}T${CHOUR}Z_*H_${ne}.nc

  ncks -4 -L 1 -d lat,${latmin},${latmax} ${arqn}".saux.nc" ${arqn}".saux1.nc"
  ncks --ppc default=.$dp ${arqn}".saux1.nc" "${arqn}.nc"
  ncatted -a _FillValue,,o,f,NaN "${arqn}.nc"
  chmod 775 "${arqn}.nc"
  mv "${arqn}.nc" $WDIR

  rm -f ${arqn}.saux*
  rm -f ${arqn}.*idx*
  rm -f *ncks* *tmp
  echo " File ${arqn} converted to netcdf and compressed with success. "
  sleep 1

done

sleep 1
cd $WDIR
rm -rf "work_${DATE}${CHOUR}"

echo " "
echo " Done download_GEWPS_wave.sh ${DATE} ${CHOUR}Z "

