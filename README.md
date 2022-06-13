# EnhancerGene
This is a repository for scripts and programs developped for the identification of enhancer/gene relationships in animal genomes

Example of script that was sent to our slurm cluster

#!/bin/sh

cd ~/dynagen/sdjebali/enhancer.gene/exploratory.analysis/methylation.rnaseq/chicken

module load system/Python-3.6.3 

python /home/sdjebali/dynagen/sdjebali/enhancer.gene/software/EnhancerGene/compute_correlations.py -d 50000 -t 4 cpgid.methylpcent.7samples.bed > cpgid.methylpcent.7samples.allcorr500kb.tsv 2> cpgid.methylpcent.7samples.allcorr500kb.err
