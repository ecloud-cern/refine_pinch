#!/bin/bash


export OMP_NUM_THREADS=24

collider=$1
PINCH=$2
parent=root://eosproject.cern.ch//eos/project/e/ecloud-simulations/kparasch/LHC_Triplets

SEY=`echo $PINCH | grep -o -P '_sey.{4}_' | cut -c 5-8`

source_file=${parent}/Pinches/${collider}/${PINCH}
xrdcp ${source_file} .

wget https://github.com/ecloud-cern/refine_pinch/raw/main/ecloud_xsuite_filemanager.py
wget https://github.com/ecloud-cern/refine_pinch/raw/main/refine_pinch.py
wget https://github.com/ecloud-cern/refine_pinch/raw/main/refinement_helpers.py
wget https://github.com/ecloud-cern/refine_pinch/raw/main/reorder_slices.py
wget https://github.com/ecloud-cern/refine_pinch/raw/main/volume_helpers.py
wget https://github.com/ecloud-cern/refine_pinch/raw/main/generate_rho_slices.py

python generate_rho_slices.py $PINCH
python reorder_slices.py $PINCH
# python refine_pinch.py $PINCH --MTI $3 --MLI $3
python refine_pinch.py $PINCH --MTI $3 --MLI $3 --DTO 2

destination_folder=${parent}/Refined_pinches/${collider}/
xrdcp -f refined_*.h5 ${destination_folder}
