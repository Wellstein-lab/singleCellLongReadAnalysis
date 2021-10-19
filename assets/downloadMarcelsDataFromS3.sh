#!/bin/bash
#----------------------------------------------------------
#
# Obtained from Marcel's Dropbox Dropbox\Single_cell_analysis\PacBio\single_cell_pieline\collapse_merged_abund
#
#----------------------------------------------------------
# mv ~/Downloads/collapse_merged_abund.zip .
# unzip collapse_merged_abund.zip
#----------------------------------------------------------
#
# Obtained from WellsteinLab Amazon Bucket Locations
#
#----------------------------------------------------------
#
# merged fasta from all
#
#----------------------------------------------------------
aws s3 cp s3://marcel-scrna/PacBio/single_cell_pipeline/trim_3pass_BC/BM_merged.fasta .
#----------------------------------------------------------
#
# Each of the samples separately
#
#----------------------------------------------------------
#
# Sample A Lin negative
#
#----------------------------------------------------------
aws s3 cp s3://marcel-scrna/PacBio/single_cell_pipeline/trim_3pass_BC/Alin_neg.flnc_BC.bam .
aws s3 cp s3://marcel-scrna/PacBio/single_cell_pipeline/trim_3pass_BC/Alin_neg.flnc_BC.bam.bpi .
aws s3 cp s3://marcel-scrna/PacBio/single_cell_pipeline/trim_3pass_BC/Alin_neg.flnc_BC.bam.pbi .
aws s3 cp s3://marcel-scrna/PacBio/single_cell_pipeline/trim_3pass_BC/Alin_neg.flnc_BC.consensusreadset.xml .
aws s3 cp s3://marcel-scrna/PacBio/single_cell_pipeline/trim_3pass_BC/Alin_neg.flnc_BC.fasta .
aws s3 cp s3://marcel-scrna/PacBio/single_cell_pipeline/trim_3pass_BC/Alin_neg.flnc_BC.filter_summary.json .
aws s3 cp s3://marcel-scrna/PacBio/single_cell_pipeline/trim_3pass_BC/Alin_neg.flnc_BC.report.csv .
aws s3 cp s3://marcel-scrna/PacBio/single_cell_pipeline/trim_3pass_BC/Alin_neg.19x_BC.trimmed.bam .
aws s3 cp s3://marcel-scrna/PacBio/single_cell_pipeline/trim_3pass_BC/Alin_neg.10x_BC.trimmed.bam .
aws s3 cp s3://marcel-scrna/PacBio/single_cell_pipeline/trim_3pass_BC/Alin_neg.10x_BC.trimmed.csv .
aws s3 cp s3://marcel-scrna/PacBio/single_cell_pipeline/trim_3pass_BC/Alin_neg_10x_BC.trimmed.bam .
aws s3 cp s3://marcel-scrna/PacBio/single_cell_pipeline/trim_3pass_BC/Alin_neg_10x_BC.trimmed.csv .
#----------------------------------------------------------
#
# Sample B BM total
#
#----------------------------------------------------------
aws s3 cp s3://marcel-scrna/PacBio/single_cell_pipeline/trim_3pass_BC/B_BM_tot.flnc_BC.bam .
aws s3 cp s3://marcel-scrna/PacBio/single_cell_pipeline/trim_3pass_BC/B_BM_tot.flnc_BC.bam.bpi .
aws s3 cp s3://marcel-scrna/PacBio/single_cell_pipeline/trim_3pass_BC/B_BM_tot.flnc_BC.bam.pbi .
aws s3 cp s3://marcel-scrna/PacBio/single_cell_pipeline/trim_3pass_BC/B_BM_tot.flnc_BC.consensusreadset.xml .
aws s3 cp s3://marcel-scrna/PacBio/single_cell_pipeline/trim_3pass_BC/B_BM_tot.flnc_BC.fasta .
aws s3 cp s3://marcel-scrna/PacBio/single_cell_pipeline/trim_3pass_BC/B_BM_tot.flnc_BC.filter_summary.json .
aws s3 cp s3://marcel-scrna/PacBio/single_cell_pipeline/trim_3pass_BC/B_BM_tot.flnc_BC.report.csv .
aws s3 cp s3://marcel-scrna/PacBio/single_cell_pipeline/trim_3pass_BC/B_BM_tot.19x_BC.trimmed.bam .
aws s3 cp s3://marcel-scrna/PacBio/single_cell_pipeline/trim_3pass_BC/B_BM_tot.10x_BC.trimmed.bam .
aws s3 cp s3://marcel-scrna/PacBio/single_cell_pipeline/trim_3pass_BC/B_BM_tot.10x_BC.trimmed.csv .
aws s3 cp s3://marcel-scrna/PacBio/single_cell_pipeline/trim_3pass_BC/B_BM_tot_10x_BC.trimmed.bam .
aws s3 cp s3://marcel-scrna/PacBio/single_cell_pipeline/trim_3pass_BC/B_BM_tot_10x_BC.trimmed.csv .
#----------------------------------------------------------
#
# Sample B Lin negative
#
#----------------------------------------------------------
aws s3 cp s3://marcel-scrna/PacBio/single_cell_pipeline/trim_3pass_BC/Blin_neg.flnc_BC.bam .
aws s3 cp s3://marcel-scrna/PacBio/single_cell_pipeline/trim_3pass_BC/Blin_neg.flnc_BC.bam.bpi .
aws s3 cp s3://marcel-scrna/PacBio/single_cell_pipeline/trim_3pass_BC/Blin_neg.flnc_BC.bam.pbi .
aws s3 cp s3://marcel-scrna/PacBio/single_cell_pipeline/trim_3pass_BC/Blin_neg.flnc_BC.consensusreadset.xml .
aws s3 cp s3://marcel-scrna/PacBio/single_cell_pipeline/trim_3pass_BC/Blin_neg.flnc_BC.fasta .
aws s3 cp s3://marcel-scrna/PacBio/single_cell_pipeline/trim_3pass_BC/Blin_neg.flnc_BC.filter_summary.json .
aws s3 cp s3://marcel-scrna/PacBio/single_cell_pipeline/trim_3pass_BC/Blin_neg.flnc_BC.report.csv .
aws s3 cp s3://marcel-scrna/PacBio/single_cell_pipeline/trim_3pass_BC/Blin_neg.19x_BC.trimmed.bam .
aws s3 cp s3://marcel-scrna/PacBio/single_cell_pipeline/trim_3pass_BC/Blin_neg.10x_BC.trimmed.bam .
aws s3 cp s3://marcel-scrna/PacBio/single_cell_pipeline/trim_3pass_BC/Blin_neg.10x_BC.trimmed.csv .
aws s3 cp s3://marcel-scrna/PacBio/single_cell_pipeline/trim_3pass_BC/Blin_neg_10x_BC.trimmed.bam .
aws s3 cp s3://marcel-scrna/PacBio/single_cell_pipeline/trim_3pass_BC/Blin_neg_10x_BC.trimmed.csv 
#----------------------------------------------------------
#
# Sample KPC (unrelated non-bone marrow sample - used as a negative control ?)
#
#----------------------------------------------------------
aws s3 cp s3://marcel-scrna/PacBio/single_cell_pipeline/trim_3pass_BC/KPC.flnc_BC.bam .
aws s3 cp s3://marcel-scrna/PacBio/single_cell_pipeline/trim_3pass_BC/KPC.flnc_BC.bam.bpi .
aws s3 cp s3://marcel-scrna/PacBio/single_cell_pipeline/trim_3pass_BC/KPC.flnc_BC.bam.pbi .
aws s3 cp s3://marcel-scrna/PacBio/single_cell_pipeline/trim_3pass_BC/KPC.flnc_BC.consensusreadset.xml .
aws s3 cp s3://marcel-scrna/PacBio/single_cell_pipeline/trim_3pass_BC/KPC.flnc_BC.fasta .
aws s3 cp s3://marcel-scrna/PacBio/single_cell_pipeline/trim_3pass_BC/KPC.flnc_BC.filter_summary.json .
aws s3 cp s3://marcel-scrna/PacBio/single_cell_pipeline/trim_3pass_BC/KPC.flnc_BC.report.csv .
aws s3 cp s3://marcel-scrna/PacBio/single_cell_pipeline/trim_3pass_BC/KPC.19x_BC.trimmed.bam .
aws s3 cp s3://marcel-scrna/PacBio/single_cell_pipeline/trim_3pass_BC/KPC.10x_BC.trimmed.bam .
aws s3 cp s3://marcel-scrna/PacBio/single_cell_pipeline/trim_3pass_BC/KPC.10x_BC.trimmed.csv .
aws s3 cp s3://marcel-scrna/PacBio/single_cell_pipeline/trim_3pass_BC/KPC_10x_BC.trimmed.bam .
aws s3 cp s3://marcel-scrna/PacBio/single_cell_pipeline/trim_3pass_BC/KPC_10x_BC.trimmed.csv .

