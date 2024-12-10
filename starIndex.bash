#!/bin/bash
######################
##SLURM CONFIGURATION#
######################
##Job name
#SBATCH -J starIndex.bash
# Standard out and Standard Error output files with the job number in the name.
#SBATCH -o logs/starIndex.bash_res.txt
#SBATCH -e logs/starIndex.bash_err.txt
#Number of cores/tasks
#SBATCH -n 8
#Number of Nodes
#SBATCH -N 1
#SBATCH --mem=128G
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
 --runMode genomeGenerate \
 --genomeDir /share/xulab/seanr/strAnalyz/maizev3Index \
 --genomeFastaFiles /share/xulab/seanr/strAnalyz/maize.v3.genome.gtf/Zea_mays.AGPv3.dna.toplevel.fa\
 --sjdbGTFfile /share/xulab/seanr/strAnalyz/maize.v3.genome.gtf/genes_apr9.gtf \
 --sjdbOverhang 150 \

