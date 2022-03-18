#!/bin/bash

# @Author: Tobias Jakobi <tjakobi>
# @Email:  tobias.jakobi@med.uni-heidelberg.de
# @Project: University Hospital Heidelberg, Section of Bioinformatics and Systems Cardiology
# @License: CC BY-NC-SA

#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 2
#SBATCH --mem=10G
#SBATCH -J "multiqc"

module load multiqc

export LC_ALL=C.UTF-8
export LANG=C.UTF-8

# check if we have 4 arguments
if [ ! $# == 4 ]; then
  echo "Usage: $0 [STAR] [FLEXBAR] [BOWTIE] [FASTQC]"
  exit
fi


multiqc $1 $2 $3 $4 -f
