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
#SBATCH --error=SAMN19128602_S10_AD.err.txt
#SBATCH --output=SAMN19128602_S10_AD.out.txt

########################################
### PUT YOUR COMMANDS BELOW THIS BOX ###
########################################
cd /data/project/lasseigne_lab/TabeaSoelter/tms_ccc/GSE174367/CellRanger/outputs/

module load CellRanger/6.1.1

cellranger count --id=SAMN19128602_S10_AD \
                 --transcriptome=/data/user/tsoelter/apps/refdata-gex-GRCh38-2020-A \
                 --include-introns \
                 --fastqs=/data/project/lasseigne_lab/TabeaSoelter/tms_ccc/GSE174367/rawData/ \
                 --sample=SRR14514049,SRR14514050,SRR14514051,SRR14514052,SRR14514053,SRR14514054,SRR14514055,SRR14514056 \