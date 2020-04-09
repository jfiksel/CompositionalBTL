#!/bin/bash
#$ -l mem_free=2G
#$ -l h_vmem=2G
#$ -l h_rt=24:00:00
#$ -cwd
#$ -j y
#$ -R y
#$ -t 1-1500
module load conda_R
Rscript CompBTL_Dirichlet.R $SGE_TASK_ID