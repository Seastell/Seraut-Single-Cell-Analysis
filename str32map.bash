#!/bin/bash
######################
##SLURM CONFIGURATION#
######################
##Job name
#SBATCH -J str32map.bash
# Standard out and Standard Error output files with the job number in the name.
#SBATCH -o logs/str32map.bash_res.txt
#SBATCH -e logs/str32map.bash_err.txt
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

STAR --runThreadN 8 \
 --genomeDir /share/xulab/seanr/strAnalyz/maizev3Index \
 --readFilesIn /share/xulab/seanr/strAnalyz/xu32/xx32_S1_L001_R2_001.fastq.gz /share/xulab/seanr/strAnalyz/xu32/xx32_S1_L001_R1_001.fastq.gz \
 --readFilesCommand zcat \
 --soloType CB_UMI_Simple \
 --soloCBwhitelist /share/xulab/software/cellranger-8.0.1/lib/python/cellranger/barcodes/3M-february-2018.txt\
 --soloUMIlen 12 \
 --soloBarcodeReadLength 0 \

