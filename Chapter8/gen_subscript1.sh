#$ -V -cwd

#$ -l h_rt=08:00:00

#$ -l h_vmem=2.5G

#$ -tc 100


R CMD BATCH singlesim_1_$SGE_TASK_ID.R singlesim_1_$SGE_TASK_ID.Rout #run r in batch command mode, run averages.R script, direct output to averages.Rout
