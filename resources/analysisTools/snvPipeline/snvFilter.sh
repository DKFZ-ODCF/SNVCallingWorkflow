#!/bin/bash

#PBS -l walltime=5:00:00
#PBS -l nodes=1:ppn=2
#PBS -l mem=2g
#PBS -m a

set -o pipefail

[[ -f ${FILENAME_CHECKPOINT} ]] && rm ${FILENAME_CHECKPOINT}

wait 60

touch ${FILENAME_CHECKPOINT}