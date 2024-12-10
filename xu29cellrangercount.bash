#!/bin/bash
######################
##SLURM CONFIGURATION#
######################
##Job name
#SBATCH -J xu29cellrangercount.bash
# Standard out and Standard Error output files with the job number in the name.
#SBATCH -o logs/xu29cellrangercount.bash_res.txt
#SBATCH -e logs/xu29cellrangercount.bash_err.txt
#Number of cores/tasks
#SBATCH -n 8
#Number of Nodes
#SBATCH -N 1
#SBATCH --mem=64G
#SBATCH -t 64:00:00
#SBATCH --mail-type=ALL # other options are ALL, NONE, BEGIN, FAIL
#SBATCH --mail-user=jhruan@ucdavis.edu

module load fastx_toolkit
module load samtools
module load hisat2
module load gatk
module load picard-tools/2.6.0
module load star
module load cellranger

cellranger count --id=xu29CountOut \
 --fastqs=/share/xulab/seanr/strAnalyz/xu29/ \
 --sample=xx29 \
 --transcriptome=/share/xulab/seanr/cellranger/analysis/yard/run_cellranger_mkref/MaizeV3_genome \

