#!/bin/bash
#SBATCH --job-name=fullWorkflow
#SBATCH --partition=sixhour
#SBATCH --mail-type=ALL
#SBATCH --mail-user=e378m007@ku.edu
#SBATCH --time=0-06:00:00
#SBATCH --output=fullWorkflow_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32gb

echo "loading programs"

module load python/3.7
module load repeatmasker/4.1.2
module load bwa/0.7.17 
module load samtools/1.14 
module load bedtools/2.30.0 

echo "running workflow"

python command_line.py