#!/bin/bash
#$ -l mem_free=2G
#$ -l h_vmem=2G
#$ -l h_rt=24:00:00
#$ -cwd
#$ -pe local 3
#$ -j y
#$ -R y
#$ -t 1-4000
module load conda_R
Rscript CompBTL_vs_Dirichlet_Single_Cause.R $SGE_TASK_ID