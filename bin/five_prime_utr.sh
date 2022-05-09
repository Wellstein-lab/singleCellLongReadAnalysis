#!/bin/bash
#/*--------------------------------------------------
# repurposing five_prime_utr process from the
#  Long Read Proteogenomics Workflow
#  https://github.com/sheynkman-lab/Long-Read-Proteogenomics/main.nf
#
# 2022 MAY 9 - step 1 ran but not step 2
# 
#  5' UTR Status
# * Intermediate step for protein classification
# * Dtermines the 5' UTR status of the protein in order 
# * to classify protein category in latter step
#
#---------------------------------------------------*/

gencode_primary_assembly_annotation="gencode.v32.primary_assembly.annotation.gtf"
reference_genome="GRCh38.primary_assembly.genome.fa"
pb_gene="pb_gene.tsv"

echo "gencode_primary_assembly_annotation = " $gencode_primary_assembly_annotation
echo "reference_genome                    = " $reference_genome

gtf=".gtf"
fasta=".fasta"
merge5=".merge5.collapsed"
filtered="filtered_"
corrected="_corrected"
degfilter=".5degfilter"
refined_orf="_refined_orf_refined.fasta"
refined_tsv="_refined_orf_refined.tsv"
orf_prob_best=".ORF_prob.best.tsv"
orf_prob=".ORF_prob.tsv"
orf_fasta=".ORF_seqs.fa"
no_orf=".no_ORF.txt"
best_orf=".best_orf.tsv"
classification=".5degfilter_classification.5degfilter"
pb_gene=".pb_gene.tsv"
with_transcript_with_cds="_with_transcript_with_cds"
cds_renamed_exon=".cds_renamed_exon"
transcript_exons_only=".transcript_exons_only"
sqanti_w_5utr="_sqanti_protein_classification_w_5utr_info.tsv"

cd ../data/BC_ranked_isoforms
PWD=$(pwd)

allsqantiproteinclasstsv="*sqanti_protein_classification.tsv"

echo "Current Working Directory is = " $PWD
echo "all sqantiproteinclasstsv    = " $allsqantiproteinclasstsv


#
# this is run just once - if the version of Gencode changes
# -- TODO - move this in the data preparation step -- the reference set up step
#  this is the command as run on the gencode.v32.primary_assembly.annotation.gtf

#docker run -v $PWD:$PWD -w $PWD gsheynkmanlab/proteogenomics-base:v1.0 1_get_gc_exon_and_5utr_info.py --gencode_gtf gencode.v32.primary_assembly.annotation.gtf --odir ./


# loop through cds scripts
# step 2
for file in $allwithcds; do
    name="${file%%.*}"
    echo "name = " $name
    echo "file = " $file
    docker run -v $PWD:$PWD -w $PWD gsheynkmanlab/proteogenomics-base:v1.0 \
         2_classify_5utr_status.py \
         --gencode_exons_bed gencode_exons_for_cds_containing_ensts.bed \
         --gencode_exons_chain gc_exon_chain_strings_for_cds_containing_transcripts.tsv \
         --sample_cds_gtf $file \
         --odir ./ 
done

    

# loop through
# step 3
#for file in $allsqantiproteinclasstsv; do
#    name="${file%%.*}"
#    name_sqanti_w_5utr=$name$sqanti_w_5utr
#
#    echo "name              = "$name
#    echo "file              = "$file
#    echo "name_sqanti_w5utr = "$name_sqanti_w_5utr

#    3_merge_5utr_info_to_pclass_table.py \
#    --name ${params.name} \
#    --utr_info pb_5utr_categories.tsv \
#    --sqanti_protein_classification $sqanti_protein_classification \
#    --odir ./

done


    
