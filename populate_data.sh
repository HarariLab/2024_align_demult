#!/usr/bin/bash

# Link to genotypes. I'm assuming all genotypes are in the same directory and are named as following: <pool>.vcf.gz.
incoming_dir="/fs/ess/PAS2694/Data/Genotypes/2024_knight_phase4_multiregion/step10_demuxlet_ready"

mkdir -p data/genotypes
ln -s ${incoming_dir}/*.vcf.gz data/genotypes/

