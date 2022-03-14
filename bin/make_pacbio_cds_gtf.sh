#!/bin/bash
#/*--------------------------------------------------
#Then we do some clean up with the CDS GTF
#/*--------------------------------------------------
#PacBio CDS GTF 
# * derive a GTF file that includes the ORF regions (as CDS features)
#https://github.com/sheynkman-lab/Long-Read-Proteogenomics/blob/23e345dafb0ef90e479cac94a29e3d702472e370/main.nf#L966-L1008

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
refined="_refined"
orf_prob_best=".ORF_prob.best.tsv"
orf_prob=".ORF_prob.tsv"
orf_fasta=".ORF_seqs.fa"
no_orf=".no_ORF.txt"
best_orf=".best_orf.tsv"
classification=".5degfilter_classification.5degfilter"
pb_gene=".pb_gene.tsv"

cd ../data/BC_ranked_isoforms
PWD=$(pwd)

allorfprobbest="*.ORF_prob.best.tsv"

echo "Current Working Directory is = " $PWD
echo "all orf prob best            = " $allorfprobbest

# loop through 
for file in $allorfprobbest; do
    name="${file%%.*}"

    name_merge5_corrected_5degfilter_with_transcript=$name$merge5$corrected$degfilter$with_transcript
    name_merge5_corrected_5degfilter_no_transcript=$name$merge5$corrected$degfilter$no_transcript
    name_merge5_corrected_5degfilter_orf_prob_best=$name$merge5$corrected$degfilter$orf_prob_best
    name_merge5_corrected_5degfilter_fasta=$name$merge5$corrected$degfilter$fasta
    name_merge5_corrected_5degfilter_classification_pb_gene=$name$merge5$corrected$classification$pb_gene
    name_merge5_corrected_gtf=$name$merge5$corrected$gtf
    
    echo "file                                            = " $file
    echo "name                                            = " $name
    echo "name_merge5_corrected_gtf                       = " $name_merge5_corrected_gtf
    echo "name_merge5_corrected_5degfilter_fasta          = " $name_merge5_corrected_5degfilter_fasta
    echo "name_merge5_corrected_5degfilter_with_transcript= " $name_merge5_corrected_5degfilter_with_transcript
    echo "name_merge5_corrected_5degfilter_no_transcript  = " $name_merge5_corrected_5degfilter_no_transcript
    echo "name_merge5_corrected_5degfilter_orf_prob_best  = " $name_merge5_corrected_5degfilter_orf_prob_best
    echo "name_merge5_corrected_5degfilter_classification_pb_gene = " $name_merge5_corrected_5degfilter_classification_pb_gene
    
#    docker run -v $PWD:$PWD -w $PWD -it gsheynkmanlab/proteogenomics-base make_pacbio_cds_gtf.py \
#	   --name $name_merge5_corrected_5degfilter_with_transcript \
#	   --sample_gtf $name_merge5_corrected_gtf \
#	   --refined_database $file \
#           --called_orfs $name_merge5_corrected_5degfilter_orf_prob_best \
#	   --pb_gene $name_merge5_corrected_5degfilter_classification_pb_gene \
#	   --include_transcript yes
#
#    docker run -v $PWD:$PWD -w $PWD -it gsheynkmanlab/proteogenomics-base make_pacbio_cds_gtf.py \
#	   --name $name_merge5_corrected_5degfilter_no_transcript \
#	   --sample_gtf $name_merge5_corrected_gtf \
#	   --refined_database $file \
#           --called_orfs $name_merge5_corrected_5degfilter_orf_prob_best \
#	   --pb_gene $name_merge5_corrected_5degfilter_classification_pb_gene \
#	   --include_transcript no

done

