#!/bin/bash
#/*--------------------------------------------------
#Then we do some clean up with the CDS GTF
#/*--------------------------------------------------
#Rename CDS to Exon
# * Preprocessing step to SQANTI Protein
# * CDS is renamed to exon and transcript stop and start
# * locations are updated to reflect CDS start and stop
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

cd ../data/BC_ranked_isoforms
PWD=$(pwd)

allwithgtf="*cds.transcript_exons_only.gtf"

echo "Current Working Directory is = " $PWD
echo "all with transcript gtf      = " $allwithgtf

#
# this is run just once - if the version of Gencode changes
# -- TODO - move this in the data preparation step -- the reference set up step
#  this is the command as run on the gencode.v32.primary_assembly.annotation.gtf
#
#  docker run -v $PWD:$PWD -w $PWD gsheynkmanlab/proteogenomics-base:v1.0 1_get_gc_exon_and_5utr_info.py --gencode_gtf gencode.v32.primary_assembly.annotation.gtf --odir ./

# loop through
#
for file in $allwithgtf; do
    name="${file%%.*}"
    name_merge5_corrected_5degfilter_with_transcript_with_cds_exons_only_gtf=$name$merge5$corrected$degfilter$with_transcript_with_cds$transcript_exons_only$gtf
    name_merge5_corrected_5degfilter_with_transcript_with_cds_renamed_exon_gtf=$name$merge5$corrected$degfilter$with_transcript_with_cds$cds_renamed_exon$gtf
    name_merge5_corrected_5degfilter_best_orf=$name$merge5$corrected$degfilter$best_orf
    param_name=$name$merge5$corrected$degfilter

    echo "name                                                                       = " $name
    echo "file                                                                       = " $file
    echo "name_merge5_corrected_5degfilter_with_transcript_with_cds_exon_only_gtf    = " $name_merge5_corrected_5degfilter_with_transcript_with_cds_exons_only_gtf
    echo "name_merge5_corrected_5degfilter_with_transcript_with_cds_renamed_exon_gtf = " $name_merge5_corrected_5degfilter_with_transcript_with_cds_renamed_exon_gtf
    echo "name_merge5_corrected_5degfilter_best_orf                                  = " $name_merge5_corrected_5degfilter_best_orf
    
    docker run -v $PWD:$PWD -w $PWD gsheynkmanlab/proteogenomics-base:v1.0 2_classify_5utr_status.py \
	   --gencode_exons_bed gencode_exons_for_cds_containing_ensts.bed --gencode_exons_chain gc_exon_chain_strings_for_cds_containing_transcripts.tsv --sample_cds_gtf filtered_Blin_neg_filt_ranked_BC_clust0_B_ccsids.merge5.collapsed_corrected.5degfilter_with_transcript_with_cds.gtf --odir ./
    # variables need to be provided in a specific order
    docker run -v $PWD:$PWD -w $PWD gsheynkmanlab/sqanti_protein:sing sqanti3_protein.py \
	   $name_merge5_corrected_5degfilter_with_transcript_with_cds_exons_only_gtf \
	   $name_merge5_corrected_5degfilter_with_transcript_with_cds_renamed_exon_gtf \
	   $name_merge5_corrected_5degfilter_best_orf \
	   gencode.transcript_exon_only.gtf \
	   gencode.cds_renamed_exon.gtf \
	   -d . \
           -p $param_name
   
done

