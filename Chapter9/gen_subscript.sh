#$ -V -cwd

#$ -l h_rt=10:00:00

#$ -l h_vmem=5G

#$ -tc 50

R CMD BATCH bigsim_pattern_$SGE_TASK_ID.R bigsim_pattern_$SGE_TASK_ID.Rout #run r in batch command mode, run averages.R script, direct output to averages.Rout
