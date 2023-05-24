#!/bin/bash

#SBATCH --job-name=ppi_contruction
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tsoelter@uab.edu
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=65000
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --partition=short
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --array=0-2

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

# sample list of seurat objects
SAMPLE_LIST="${wd}/results/intermediate_outputs/inputs/ppi_construction_gex_inputs.txt"
SAMPLE_ARRAY=(`cat ${SAMPLE_LIST}`)
INPUT=`echo ${SAMPLE_ARRAY[$SLURM_ARRAY_TASK_ID]}`

# sample list of mapped gene inputs
SAMPLE_LIST2="${wd}/results/intermediate_outputs/inputs/ppi_construction_gene_inputs.txt"
SAMPLE_ARRAY2=(`cat ${SAMPLE_LIST2}`)
INPUT2=`echo ${SAMPLE_ARRAY2[$SLURM_ARRAY_TASK_ID]}`

# sample list of tmp ppi objects
SAMPLE_LIST3="${wd}/results/intermediate_outputs/inputs/ppi_construction_ppi_inputs.txt"
SAMPLE_ARRAY3=(`cat ${SAMPLE_LIST3}`)
INPUT3=`echo ${SAMPLE_ARRAY3[$SLURM_ARRAY_TASK_ID]}`

# execute docker
singularity exec --cleanenv --containall -B ${wd} ${wd}/bin/docker/rstudio_ccc_ad_1.0.1.sif Rscript --vanilla ${src}/03_inhibitory_ppi_construction.R ${INPUT} ${INPUT2} ${INPUT3} 