#!/bin/bash
#SBATCH -n 48
#SBATCH -N 1     # Number of cores
#SBATCH --job-name=jags-model
#SBATCH --output=out-%j.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --time=48:00:00
#SBATCH --mail-type=END  	  # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=jdeleeuw@vassar.edu

Rscript ~/shape-of-statistical-learning/experiment-1/junior/test-cores.R