#!/bin/bash
#
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tsoelter@uab.edu
#SBATCH --job-name=CellRanger
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu=8000
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --share
#SBATCH --partition=short
#SBATCH --error=%j.%N.err.txt
#SBATCH --output=%j.%N.out.txt

########################################
### PUT YOUR COMMANDS BELOW THIS BOX ###
########################################
cd /data/project/lasseigne_lab/TabeaSoelter/tms_CaSpER_data/CellRangerOutputs/outputs/HumanHIPIntrons/

cellranger count --id=SRR13278449_S13_AD \
                 --transcriptome=/data/user/tsoelter/apps/refdata-gex-GRCh38-2020-A \
                 --fastqs=/data/project/lasseigne_lab/TabeaSoelter/tms_CaSpER_data/rawData \
                 --sample=SRR13278449 \
                 --expect-cells=6000