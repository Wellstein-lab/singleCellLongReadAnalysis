#!/bin/bash
#
#/*--------------------------------------------------
#Transcriptome Summary
# * Compares the abundance (CPM) based on long-read sequencing
# * to the abundances (TPM) inferred from short-read sequencing,
# * as computed by Kallisto (analyzed outside of this pipeline).
# * Additionally produces a pacbio-gene reference table
#
#
# THE MOST IMPORTANT THING HERE IS THE PRODUCTION OF THE PACBIO-GENE REFERENCE TABLE
# In our case here we will have one pacbio-gene reference per cluster
#
#  input:
#    file(sqanti_classification) from ch_sample_classification_transcriptome
#    file(tpm) from ch_sample_kallisto
#    file(ribo) from ch_normalized_ribo_kallisto
#    file(ensg_to_gene) from ch_ensg_gene
#    file(enst_to_isoname) from ch_enst_isoname
#    file(len_stats) from ch_gene_lens_transcriptome
# output: 
#  file("gene_level_tab.tsv") into ch_gene_level
#  file("sqanti_isoform_info.tsv") into ch_sqanti_isoform_info
#  file("pb_gene.tsv") into ch_pb_gene
#---------------------------------------------------*/

allccsidsfasta="*ccsids.fasta"
reference_genome="GRCh38.primary_assembly.genome.fa"
ensg_to_gene="ensg_gene.tsv"
enst_to_isoname="enst_isoname.tsv"
lens_stat="isoname_lens.tsv"

merge5=".merge5.collapsed"
filtered="filtered_"
corrected="_corrected"
degfilter=".5degfilter"
classification=".5degfilter_classification.5degfilter"
pb_gene=".pb_gene.tsv"
pb_gene_tsv="pb_gene.tsv"


cd ../data/BC_ranked_isoforms
PWD=$(pwd)
allsqanticlassification="filtered*5degfilter*classification*tsv"*
echo "Current Working Directory is = " $PWD

#
# an error occurs but we get our pb_gene file!
#
for file in $allsqanticlassification; do

    name="${file%%.*}"
    name_merge5_corrected_5degfilter_classification_pb_gene=$name$merge5$corrected$classification$pb_gene

    echo "file                                                    = " $file
    echo "name                                                    = " $name
    echo "name_merge5_corrected_5degfilter_classification_pb_gene = " $name_merge5_corrected_5degfilter_classification_pb_gene
    
#    docker run --rm -v $PWD:$PWD -w $PWD -it ghcr.io/adeslatt/transcriptome-summary-docker transcriptome_summary.py \
#	   --sq_out $file \
#	   --ensg_to_gene ensg_gene.tsv \
#	   --enst_to_isoname isoname_lens.tsv \
#	   --len_stats gene_lens.tsv \
#	   --odir junk
    
#    mv $pb_gene_tsv $name_merge5_corrected_5degfilter_classification_pb_gene
    
done
