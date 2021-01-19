#$ -V -cwd

#$ -l h_rt=02:00:00

#$ -l h_vmem=5G

#$ -tc 100

R CMD BATCH singlesim_0_$SGE_TASK_ID.R singlesim_0_$SGE_TASK_ID.Rout #run r in batch command mode, run averages.R script, direct output to averages.Rout
  
