#!/bin/bash
#
#  Filename:  run_filter_sqanti.sh
#
# 
#
#Filter SQANTI
#* Filter SQANTI results based on several criteria
#* - protein coding only
#*      PB transcript aligns to a GENCODE-annotated protein coding gene.
#* - percent A downstream
#*      perc_A_downstreamTTS : percent of genomic "A"s in the downstream 20 bp window.
# *      If this number if high (> 80%), the 3' end have arisen from intra-priming during the RT step 
# * - RTS stage
# *      RTS_stage: TRUE if one of the junctions could be an RT template switching artifact.
# * - Structural Category
# *      keep only transcripts that have a isoform structural category of:
# *        -novel_not_in_catalog
# *        -novel_in_catalog
# *        -incomplete-splice_match
# *        -full-splice_match
#
#
# From https://github.com/sheynkman-lab/Long-Read-Proteogenomics/blob/main/main.nf
# Specifically with permalink
# https://github.com/sheynkman-lab/Long-Read-Proteogenomics/blob/23e345dafb0ef90e479cac94a29e3d702472e37#0/main.nf#L655-L731
#
#  input:
#    file(classification) from ch_sample_unfiltered_classification
#    file(sample_fasta) from ch_sample_unfiltered_fasta
#    file(sample_gtf) from ch_sample_unfiltered_gtf
#    file(protein_coding_genes) from ch_protein_coding_genes_filter_sqanti
#    file(ensg_gene) from ch_ensg_gene_filter
#
#  output:
#    file("${params.name}_classification.5degfilter.tsv") into ch_sample_classification
#    file("${params.name}_corrected.5degfilter.fasta") into ch_sample_fasta
#    file("${params.name}_corrected.5degfilter.gff") into ch_sample_gtf
#    file("*")
#
#  script:
#    """
#    filter_sqanti.py \
#    --sqanti_classification $classification \
#    --sqanti_corrected_fasta $sample_fasta \
#    --sqanti_corrected_gtf $sample_gtf \
#    --protein_coding_genes $protein_coding_genes \
#    --ensg_gene $ensg_gene \
#    --filter_protein_coding yes \
#    --filter_intra_polyA yes \
#    --filter_template_switching yes \
#    --percent_A_downstream_threshold 95 \
#    --structural_categories_level strict \
#    --minimum_illumina_coverage 3 \
#
#  collapse_isoforms.py \
#    --name ${params.name} \
#    --sqanti_gtf filtered_${params.name}_corrected.gtf \
#    --sqanti_fasta filtered_${params.name}_corrected.fasta
#
#    collapse_classification.py \
#    --name ${params.name} \
#    --collapsed_fasta ${params.name}_corrected.5degfilter.fasta \
#    --classification filtered_${params.name}_classification.tsv
#
#---------------------------------------------------*/
#
#


gtf=".gtf"
fasta=".fasta"
collapsed=".collapsed.filtered.rep"
corrected="_corrected"
classification="_classification.txt"
allclassified="*classification.txt"

#cd ../data/BC_ranked_isoforms
cd $1

PWD=$(pwd)
echo "Current Working Directory is = " $PWD
#echo "allclassified                = " $allclassified

# loop through 
#for file in $allclassified; do
#    name="${file%%.*}"
    name_classification=$2
    name_corrected_fasta=$3
    name_corrected_gtf=$4

    echo "name_classification  = " $name_classification
    echo "name_corrected_fasta = " $name_corrected_fasta
    echo "name_corrected_gtf   = " $name_corrected_gtf
    
    docker run --rm -v $PWD:$PWD -w $PWD -it gsheynkmanlab/proteogenomics-base filter_sqanti.py \
	   --sqanti_classification  $name_classification \
	   --sqanti_corrected_fasta $name_corrected_fasta \
	   --sqanti_corrected_gtf   $name_corrected_gtf \
	   --protein_coding_genes   protein_coding_genes.txt \
	   --ensg_gene              ensg_gene.tsv \
	   --filter_protein_coding  yes \
           --filter_intra_polyA     yes \
	   --filter_template_switching yes \
	   --percent_A_downstream_threshold 95 \
	   --structural_categories_level strict \
	   --minimum_illumina_coverage 3

#done

