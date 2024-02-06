#!/bin/bash
#SBATCH --partition 128x24
#SBATCH --job-name=VCF_to_aHMM
#SBATCH --output=VCF_to_aHMM.out
#SBATCH -e VCF_to_aHMM.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jharenca@ucsc.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24

# Converting VCF tyo aHMM input file

# --minGT : for parentals, min ave # of genotyped called for a site (prevents including site with no info)
# --minDP : min depth at the site (0.5)
# --rate : recombination rate
# --dist : distance between sites 1000 a pretty good default/min
# --minDif : allele freq difference between parentals (approximates Fst)

# load python
module load python-3.6.2

# Change directory 
cd /hb/groups/kay_lab/Julia_popgen/aHMM/old_py_test

# Convert!
python3 vcf2aHMM.py --vcf /hb/groups/kay_lab/Julia_popgen/aHMM/AHVF1_filtered_scaf_1.vcf --pop aHMM_pop_parent_subset_old_style.pop --rate 2 --minGT 1 --minDif .8 --dist 1000 --minDP 0.5 1>chr1_old_py.panel
