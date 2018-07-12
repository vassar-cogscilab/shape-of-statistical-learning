#!/bin/bash
#SBATCH -n 64
#SBATCH -N 1
#SBATCH --partition=emc
#SBATCH --job-name=jags-model
#SBATCH --output=out-%j.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=48:00:00
#SBATCH --mail-type=END  	  # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=jdeleeuw@vassar.edu

Rscript ~/shape-of-statistical-learning/experiment-2/junior/junior-jags-fit-exp-2.R