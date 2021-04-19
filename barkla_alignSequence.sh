#!/bin/sh

# Request GPU partitions
#SBATCH -p gpu,gpuc
# Request the number of GPUs to be used (if more than 1 GPU is required, change 1 into Ngpu, where Ngpu=2,3,4)
#SBATCH --gres=gpu:1
# Request the number of nodes
#SBATCH -N 1
# Request the number of CPU cores (There are 24 cores on the GPU node). 1 core is enough for us.
#SBATCH -n 1 -t 30


# Usual invocation, eg.: ./alignSequence -g data/dna/dna_01.txt data/dna/dna_02.txt
./alignSequence -g <TEXT_FILE> <PATTERN_FILE>
