#!/bin/bash
#$ -l mem_free=6G
#$ -l h_vmem=6G
#$ -l h_rt=24:00:00
#$ -cwd
#$ -j y
#$ -R y
#$ -t 1-8000
module load conda_R
Rscript phmrc_adult_5_causes.R $SGE_TASK_ID