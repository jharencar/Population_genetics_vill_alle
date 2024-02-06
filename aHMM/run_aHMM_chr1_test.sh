#!/bin/bash
#SBATCH --partition 128x24
#SBATCH --job-name=QTL_aHMM_1scaf
#SBATCH --output=QTL_aHMM_1scaf.out
#SBATCH -e QTL_aHMM_1scaf.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jharenca@ucsc.edu
#SBATCH --mem=124G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20

# load conda env containing many dependencies
module load miniconda3.9
source activate JGH_conda_env
cd /hb/groups/kay_lab/Julia_popgen/aHMM/old_py_test

# Run ancestry_hmm: 
# This will fit a model of a single pulse at time 10 with 50% of individuals from the first population, and 50% from the second. 
# Genotypes are used with a uniform error rate of 1e-3. 
# All individuals are ploidy 2. 

/hb/home/jharenca/Software/Ancestry_HMM/src/ancestry_hmm -i chr1_old_py.panel -s fixed_ahmm.ploidy -a 2 0.5 0.5 -p 0 100000 0.5 -p 1 10 0.5 -e 1e-3


## Example from github page: 
# ./ancestry_hmm -i example.panel -s example.sample -a 2 0.8 0.2 -p 0 100000 0.8 -p 1 10 0.2 -g -e 1e-3
#This will fit a model of a single pulse at time 10 with 80% of individuals from the first population, and 20% from the second. 
#Genotypes are used with a uniform error rate of 1e-3. 
#All individuals are ploidy 2. 
