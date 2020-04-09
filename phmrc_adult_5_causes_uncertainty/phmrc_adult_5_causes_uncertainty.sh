#!/bin/bash
#$ -l mem_free=2G
#$ -l h_vmem=2G
#$ -l h_rt=24:00:00
#$ -cwd
#$ -j y
#$ -R y
#$ -t 1-8000
module load conda_R
Rscript phmrc_adult_5_causes_uncertainty.R $SGE_TASK_ID