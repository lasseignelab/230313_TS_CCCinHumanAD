#!/bin/bash

#SBATCH --job-name=shortest_path
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tsoelter@uab.edu
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=65000
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --partition=short
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --array=0-15

########################################
### PUT YOUR COMMANDS BELOW THIS BOX ###
########################################

# load modules
module load Singularity/3.5.2-GCC-5.4.0-2.26

# set variables
wd="/data/user/tsoelter/projects/230313_TS_CCCinHumanAD"
src="/data/user/tsoelter/projects/230313_TS_CCCinHumanAD/src/shortest_path"

export SINGULARITYENV_PASSWORD='pass'
export SINGULARITYENV_USER='tsoelter'

# sample list of igraph objects
SAMPLE_LIST="${wd}/results/intermediate_outputs/inputs/shortest_path_inputs.txt"
SAMPLE_ARRAY=(`cat ${SAMPLE_LIST}`)
INPUT=`echo ${SAMPLE_ARRAY[$SLURM_ARRAY_TASK_ID]}`

# execute docker
singularity exec --cleanenv --containall -B ${wd} ${wd}/bin/docker/rstudio_ccc_ad_1.0.1.sif Rscript --vanilla ${src}/05_shortest_path.R ${INPUT}