#!/bin/bash

module load R

Rscript biowulf.slingshot.tradeSeq.R > biowulf.slingshot.tradeSeq.R.out

# run:
# sbatch --mem=96g --cpus-per-task=8 --gres=lscratch:100 --time=10:00:00 biowulf.slingshot.tradeSeq.sbatch.sh
