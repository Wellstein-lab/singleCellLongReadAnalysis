#!/bin/bash

#
#
#
cd ../data/BC_ranked_isoforms
PWD=$(pwd)

new_gff3="_new.gff3"
statistical_results="_statisticalResults"
gtf="_B_ccsids.merge5.collapsed_corrected.gtf"
class="_B_ccsids.merge5.collapsed_classification.txt"
junctions="_B_ccsids.merge5.collapsed_junctions.txt"
Blin_neg="Blin_neg_filt_ranked_BC_clust"
B_BM_tot="B_BM_tot_filt_ranked_BC_clust"

for ((clust=0;clust<=9;clust++)); do
    blin_clust_gtf=$Blin_neg$clust$gtf
    blin_clust_class=$Blin_neg$clust$class
    blin_clust_junctions=$Blin_neg$clust$junctions
    blin_clust_new_gff3=$Blin_neg$clust$new_gff3
    blin_clust_statistics=$Blin_neg$clust$statistical_results

    btot_clust_gtf=$B_BM_tot$clust$gtf
    btot_clust_class=$B_BM_tot$clust$class
    btot_clust_junctions=$B_BM_tot$clust$junctions
    btot_clust_new_gff3=$B_BM_tot$clust$new_gff3
    btot_clust_statistics=$B_BM_tot$clust$statistical_results

    echo "blin_clust_gtf        = " $blin_clust_gtf
    echo "blin_clust_class      = " $blin_clust_class
    echo "blin_clust_junc       = " $blin_clust_junctions
    echo "blin_clust_new_gff3   = " $blin_clust_new_gff3
    echo "blin_clust_statistics = " $blin_clust_statistics

    echo "btot_clust_gtf        = " $btot_clust_gtf
    echo "btot_clust_class      = " $btot_clust_class
    echo "btot_clust_junc       = " $btot_clust_junctions
    echo "btot_clust_new_gff3   = " $btot_clust_new_gff3
    echo "btot_clust_statistics = " $btot_clust_statistics

    python3 ../../IsoAnnotLite/IsoAnnotLite_SQ3.py \
	    $blin_clust_gtf \
	    $blin_clust_class \
	    $blin_clust_junctions \
	    -gff3 ../../IsoAnnotLite/Homo_sapiens_GRCh38_Ensembl_86.gff3 \
	    -o  $blin_clust_new_gff3 \
	    -stdout $blin_clust_statistics \
	    -novel -nointronic -saveTranscriptIDs

    python3 ../../IsoAnnotLite/IsoAnnotLite_SQ3.py \
	    $btot_clust_gtf \
	    $btot_clust_class \
	    $btot_clust_junctions \
	    -gff3 ../../IsoAnnotLite/Homo_sapiens_GRCh38_Ensembl_86.gff3 \
	    -o  $btot_clust_new_gff3 \
	    -stdout $btot_clust_statistics \
	    -novel -nointronic -saveTranscriptIDs

done
