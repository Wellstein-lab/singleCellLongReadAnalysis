#!/bin/bash

#
#
#
cd ../data/BC_ranked_isoforms
PWD=$(pwd)

#python3 ../../IsoAnnotLite/IsoAnnotLite_SQ3.py \
#	filtered_Blin_neg_filt_ranked_BC_clust0_B_ccsids.merge5.collapsed_corrected.gtf \
#	Blin_neg_filt_ranked_BC_clust0_B_ccsids.merge5.collapsed_classification.txt \
#	Blin_neg_filt_ranked_BC_clust0_B_ccsids.merge5.collapsed_junctions.txt \
#	-gff3 ../../IsoAnnotLite/Homo_sapiens_GRCh38_Ensembl_86.gff3 \
#	-o  Blin_neg_clust0_B_new.gff3 \
#	-stdout Blin_neg_clust0_B_statisticalResults \
#	-novel -nointronic -saveTranscriptIDs

python3 ../../IsoAnnotLite/IsoAnnotLite_SQ3.py \
	filtered_Blin_neg_filt_ranked_BC_clust1_B_ccsids.merge5.collapsed_corrected.gtf \
	Blin_neg_filt_ranked_BC_clust1_B_ccsids.merge5.collapsed_classification.txt \
	Blin_neg_filt_ranked_BC_clust1_B_ccsids.merge5.collapsed_junctions.txt \
	-gff3 ../../IsoAnnotLite/Homo_sapiens_GRCh38_Ensembl_86.gff3 \
	-o  Blin_neg_clust1_B_new.gff3 \
	-stdout Blin_neg_clust1_B_statisticalResults \
	-novel -nointronic -saveTranscriptIDs

