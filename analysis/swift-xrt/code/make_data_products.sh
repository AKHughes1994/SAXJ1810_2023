#!/bin/bash

################################################################################################
# For this to work you need both HEASOFT and the Swift-XRT Python API installed and functional #
################################################################################################

# Basic arrays
segments=("05" "06" "07" "08" "09" "10" "11" "12" "13" "14" "15" "17" "18" "19" "20" "21" "22")
modes=('wt' 'wt' 'wt' 'wt' 'wt' 'wt' 'wt' 'wt' 'wt' 'wt' 'wt' 'pc' 'wt' 'wt' 'wt' 'wt' 'wt')
len=${#segments[@]}
directory="$(pwd)"

# First get your data products
python3 get_data_products.py --segments 5 6 7 8 9 10 11 12 13 14 15 17 18 19 20 21 22 --no-centroid

# Next run grphha on the spectra that correspond to detections -- spectra with chi2 fitting
for ((i=0; i<$len; i++));do
    name="${directory/code/"spectra"}"/Obs_000324590${segments[$i]}${modes[$i]}
    echo ${name}source.pi
    grppha ${name}source.pi !${name}group.pi comm="chkey backfile ${name}back.pi & chkey respfile ${name}.rmf & chkey ancrfile ${name}.arf & exit"
    python new_grppha.py -i ${name}group.pi -o ${name}final.pi -l 500 -u 10000 -c 25 -e 0.0
done

# Next run grphha on the spectra that correspond to detections -- spectra with cashstat fitting
segments=("23" "24")
modes=("wt" "pc")
len=${#segments[@]}

# First get your data products
python3 get_data_products.py --segments 23 24 --no-centroid

for ((i=0; i<$len; i++));do
    name="${directory/code/"spectra"}"/Obs_000324590${segments[$i]}${modes[$i]}
    echo ${name}source.pi
    grppha ${name}source.pi !${name}group.pi comm="chkey backfile ${name}back.pi & chkey respfile ${name}.rmf & chkey ancrfile ${name}.arf & exit"
    python new_grppha.py -i ${name}group.pi -o ${name}final.pi -l 500 -u 10000 -c 1 -e 0.0
done
