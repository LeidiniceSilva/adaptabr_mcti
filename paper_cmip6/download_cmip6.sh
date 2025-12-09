#!/bin/bash

#SBATCH -N 1
#SBATCH -t 24:00:00
#SBATCH -J CMIP6
#SBATCH -p esp
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=mda_silv@ictp.it

model_list=( 'ACCESS-CM2' 'BCC-CSM2-MR' 'CanESM5' 'CMCC-ESM2' 'CNRM-CM6-1' 'CNRM-ESM2-1' 'GFDL-ESM4' 'INM-CM4-8' 'INM-CM5-0' 'KIOST-ESM' 'MIROC6' 'MIROC-ES2L' 'MPI-ESM1-2-HR' 'MPI-ESM1-2-LR' 'MRI-ESM2-0' 'NESM3' 'NorESM2-MM' ) 

var_list=( 'zg')     
exp='historical'

echo "Starting download"
for mdl in ${model_list[@]}; do
    for var in ${var_list[@]}; do

        if [ ${mdl} == 'ACCESS-CM2' ]; then
	    mdl_family='CSIRO-ARCCSS'
	    type='gn'
	    member='r1i1p1f1'
	    version='v20191108'
            if [ ${var} == 'hus' ] || [ ${var} == 'ta' ] || [ ${var} == 'ua' ] || [ ${var} == 'va' ]; then
	        declare -a YEARS=('195001-201412')
	    else
	        declare -a YEARS=('185001-201412')
	    fi
	    
        elif [ ${mdl} == 'BCC-CSM2-MR' ]; then
	     mdl_family='BCC'
	     type='gn'
	     member='r1i1p1f1'
	     version='v20181126'
            if [ ${var} == 'hus' ] || [ ${var} == 'ta' ] || [ ${var} == 'ua' ] || [ ${var} == 'va' ]; then
	        declare -a YEARS=('197001-200912' '201001-201412')
	    else
	        declare -a YEARS=('185001-201412')
	    fi
	    
        elif [ ${mdl} == 'CanESM5' ]; then
	     mdl_family='CCCma'
	     type='gn'
	     member='r1i1p1f1'
	     version='v20190429'
	        declare -a YEARS=('185001-201412')
	    
        elif [ ${mdl} == 'CMCC-ESM2' ]; then
	     mdl_family='CMCC'
	     type='gn'
	     member='r1i1p1f1'
	     version='v20210114'
            if [ ${var} == 'hus' ] || [ ${var} == 'ta' ] || [ ${var} == 'ua' ] || [ ${var} == 'va' ]; then
	        declare -a YEARS=('197501-199912' '200001-201412')
	    else
	        declare -a YEARS=('185001-201412')
	    fi
	    
        elif [ ${mdl} == 'CNRM-CM6-1' ]; then
	     mdl_family='CNRM-CERFACS'
	     type='gr'
	     member='r1i1p1f2'
	     version='v20180917'
	        declare -a YEARS=('185001-201412')
	    	     
        elif [ ${mdl} == 'CNRM-ESM2-1' ]; then
	     mdl_family='CNRM-CERFACS'
	     type='gr'
	     member='r1i1p1f2'
	     version='v20181206'	
	     declare -a YEARS=('185001-201412')
           
        elif [ ${mdl} == 'GFDL-ESM4' ]; then
	     mdl_family='NOAA-GFDL'
	     type='gr1'
	     member='r1i1p1f1'
	     version='v20190726'
	     declare -a YEARS=('195001-201412')
	     
        elif [ ${mdl} == 'INM-CM4-8' ]; then
	     mdl_family='INM'
	     type='gr1'
	     member='r1i1p1f1'
            if [ ${var} == 'hus' ] || [ ${var} == 'ta' ] || [ ${var} == 'ua' ] || [ ${var} == 'va' ]; then
                version='v20190605'
	        declare -a YEARS=('197001-197912' '198001-198912' '199001-199912' '200001-200912' '201001-201412')
	    else     
	    	version='v20190530'
	        declare -a YEARS=('195001-201412')
	    fi

        elif [ ${mdl} == 'INM-CM5-0' ]; then
	     mdl_family='INM'
	     type='gr1'
	     member='r1i1p1f1'
	     version='v20190610'
            if [ ${var} == 'hus' ] || [ ${var} == 'ta' ] || [ ${var} == 'ua' ] || [ ${var} == 'va' ]; then
	        declare -a YEARS=('197001-197912' '198001-198912' '199001-199912' '200001-200912' '201001-201412')
	    else
	        declare -a YEARS=('195001-201412')
	    fi
	    
        elif [ ${mdl} == 'KIOST-ESM' ]; then
	     mdl_family='KIOST'
	     type='gr1'
	     member='r1i1p1f1'
            if [ ${var} == 'hus' ] || [ ${var} == 'ta' ] || [ ${var} == 'ua' ] || [ ${var} == 'va' ]; then
                version='v20191104'
	    else     
	    	version='v20210928'
	    fi
	     
	elif [ ${mdl} == 'MIROC6' ]; then
	     mdl_family='MIROC'
	     type='gn'
	     member='r1i1p1f1'
             if [ ${var} == 'hus' ] || [ ${var} == 'ta' ] || [ ${var} == 'ua' ] || [ ${var} == 'va' ]; then
                 version='v20190311'
	         declare -a YEARS=('197001-197912' '198001-198912' '199001-199912' '200001-200912' '201001-201412')
	     else
                 version='v20181212'
	         declare -a YEARS=('195001-201412')
	     fi
	    
	elif [ ${mdl} == 'MIROC-ES2L' ]; then
	     mdl_family='MIROC'
	     type='gn'
	     member='r1i1p1f2'
	     version='v20190823'
	     declare -a YEARS=('185001-201412')
	      
	elif [ ${mdl} == 'MPI-ESM1-2-HR' ]; then
	     mdl_family='MPI-M'
	     type='gn'
	     member='r1i1p1f1'
	     version='v20190710'
	     declare -a YEARS=('197001-197412' '197501-197912' '198001-198412' '198501-198912' '199001-199412' '199501-199912' '200001-200412' '200501-200912' '201001-201412')
	     
	elif [ ${mdl} == 'MPI-ESM1-2-LR' ]; then
	     mdl_family='MPI-M'
	     type='gn'
	     member='r1i1p1f1'
	     version='v20190710'
	     declare -a YEARS=('197001-198912' '199001-200912' '201001-201412')
	  	     
	elif [ ${mdl} == 'MRI-ESM2-0' ]; then
	     mdl_family='MRI'
	     type='gn'
	     member='r1i1p1f1'
	     version='v20190308'
	     declare -a YEARS=('195001-199912' '199001-200912' '201001-201412')
             if [ ${var} == 'hus' ] || [ ${var} == 'ta' ] || [ ${var} == 'ua' ] || [ ${var} == 'va' ]; then
                 version='v20190308'
	         declare -a YEARS=('195001-199912')
	         #declare -a YEARS=('195001-199912' '200001-201412')
	     else
                 version='v20190222'
	         declare -a YEARS=('185001-201412')
	     fi

	elif [ ${mdl} == 'NESM3' ]; then
	     mdl_family='NUIST'
	     type='gn'
	     member='r1i1p1f1'
	     version='v20190630'
	     declare -a YEARS=('185001-201412')
	   
	else
	     mdl_family='NCC'
	     type='gn'
	     member='r1i1p1f1'
	     version='v20191108'
	     declare -a YEARS=('197001-197912' '198001-198912' '199001-199912' '200001-200912' '201001-201412')

	fi
	
	base_url="https://data.ceda.ac.uk/badc/cmip6/data/CMIP6/CMIP/${mdl_family}/${mdl}/${exp}/${member}/Amon/${var}/${type}/${version}"
	dir="/home/mda_silv/users/AdaptaBr_MCTI/database/paper_cmip6/cmip6/${mdl}"
	
	for year in "${YEARS[@]}"; do

	    filename="${var}_Amon_${mdl}_${exp}_${member}_${type}_${year}.nc"
	    url="${base_url}/${filename}"
	    file_dir="/home/mda_silv/users/AdaptaBr_MCTI/database/paper_cmip6/cmip6/${mdl}/${filename}"

	    if [ -f "$file_dir" ]; then
		echo "File exists: ${filename}"
	    else
		echo "Downloading: $filename"
		wget -P "$dir" ${url}
	    fi
	    
	done
    done
done

echo "Finished downloads."
