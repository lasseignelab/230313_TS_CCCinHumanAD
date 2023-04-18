#!/bin/bash
#
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tsoelter@uab.edu
#SBATCH --job-name=CellRanger
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu=64000
#SBATCH --nodes=1
#SABTCH --cpus-per-task=4
#SBATCH --time=2-02:00:00
#SBATCH --share
#SBATCH --partition=medium
#SBATCH --error=%A_%a.err.txt
#SBATCH --output=S%A_%a.out.txt
#SBATCH --array=0-7

########################################
### PUT YOUR COMMANDS BELOW THIS BOX ###
########################################
SAMPLE_LIST=($(<sample_sheet.csv))
SAMPLE=${SAMPLE_LIST[${SLURM_ARRAY_TASK_ID}]}

cd /data/project/lasseigne_lab/TabeaSoelter/tms_ccc/VINE_cortex/CellRanger/outputs/

module load CellRanger/6.1.1

cellranger count --id=$SAMPLE \
                 --transcriptome=/data/user/tsoelter/apps/refdata-gex-GRCh38-2020-A \
                 --include-introns \
                 --fastqs=/data/project/lasseigne_lab/TabeaSoelter/tms_ccc/VINE_cortex/rawData/ \
                 --sample=$SAMPLE \