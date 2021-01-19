#$ -V -cwd

#$ -l h_rt=01:00:00

#$ -l h_vmem=3G

#$ -tc 100



R CMD BATCH singlesim_3_$SGE_TASK_ID.R singlesim_3_$SGE_TASK_ID.Rout #run r in batch command mode, run averages.R script, direct output to averages.Rout
