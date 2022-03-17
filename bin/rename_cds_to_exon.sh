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
with_transcript="_with_transcript"
no_transcript="_no_transcript"

cd ../data/BC_ranked_isoforms
PWD=$(pwd)

allwithgtf="*with_transcript_with_cds.gtf"

echo "Current Working Directory is = " $PWD
echo "all with transcript gtf      = " $allwithgtf

# loop through 
for file in $allwithgtf; do
    name="${file%%.*}"
    name_merge5_corrected_5degfilter_with_transcript=$name$merge5$degfilter$with_transcript
    name_merge5_corrected_5degfilter_with_transcript_gtf=$name$merge5$degfilter$with_transcript$gtf

    echo "file                                                 = " $file
    echo "name                                                 = " $name
    echo "name_merge5_corrected_5degfilter_with_transcript     = " $name_merge5_corrected_5degfilter_with_transcript
    echo "name_merge5_corrected_5degfilter_with_transcript_gtf = " $name_merge5_corrected_5degfilter_with_transcript_gtf

            rename_cds_to_exon.py \

    docker run -v $PWD:$PWD -w $PWD -it gsheynkmanlab/proteogenomics-base:v1.0 rename_cds_to_exon \
	   --sample_gtf             $name_merge5_corrected_5degfilter_with_transcript_gtf \
           --sample_name            $name_merge5_corrected_5degfilter_with_transcript \
           --reference_gtf          $reference_gtf \
           --reference_name         gencode \
           --num_cores 8

done

