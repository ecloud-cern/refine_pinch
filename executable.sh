#!/bin/bash


PINCH=$1
parent=root://eosproject.cern.ch//eos/project/e/ecloud-simulations/kparasch/LHC_Triplets
optics=run3_bet30cm_160urad_1.2e11ppb_2.0um/
eos_folder=/eos/project/e/ecloud-simulations/...

SEY=`echo $PINCH | grep -o -P '_sey.{4}_' | cut -c 5-8`

source_file=${parent}/Pinches/${optics}/SEY${SEY}/${PINCH}
xrdcp ${source_file} .

wget https://github.com/ecloud-cern/refine_pinch/raw/main/ecloud_xsuite_filemanager.py
wget https://github.com/ecloud-cern/refine_pinch/raw/main/refine_pinch.py
wget https://github.com/ecloud-cern/refine_pinch/raw/main/refinement_helpers.py
wget https://github.com/ecloud-cern/refine_pinch/raw/main/reorder_slices.py
wget https://github.com/ecloud-cern/refine_pinch/raw/main/volume_helpers.py

python reorder_slices.py $PINCH
python refine_pinch.py $PINCH --MTI $2 --MLI $2

destination_folder=${parent}/Refined_pinches/${optics}/SEY${SEY}
xrdcp ${PINCH} ${destination_folder}