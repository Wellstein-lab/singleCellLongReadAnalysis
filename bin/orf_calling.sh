#!/bin/bash
#/*--------------------------------------------------
#ORF Calling 
# * Selects the most plausible ORF from each pacbio transcript,
# * using the following information
# *    comparison of ATG start to reference (GENCODE) 
# *        - selects ORF with ATG start matching the one in the reference, if it exists 
# *    coding probability score from CPAT
# *    number of upstream ATGs for the candidate ORF
# *        - decrease score as number of upstream ATGs increases 
# *         using sigmoid function
# *  Additionally provides calling confidence of each ORF called
# *      - Clear Best ORF  : best score and fewest upstream ATGs of all called ORFs
# *      - Plausible ORF   : not clear best, but decent CPAT coding_score (>0.364) 
# *      - Low Quality ORF : low CPAT coding_score (<0.364)       
#---------------------------------------------------*/

gencode_primary_assembly_annotation="gencode.v32.primary_assembly.annotation.gtf"
reference_genome="GRCh38.primary_assembly.genome.fa"
pb_gene="pb_gene.tsv"

echo "gencode_primary_assembly_annotation = " $gencode_primary_assembly_annotation
echo "reference_genome                    = " $reference_genome

tsv=".tsv"
gtf=".gtf"
fasta=".fasta"
merge5=".merge5.collapsed"
filtered="filtered_"
corrected="_corrected"
degfilter=".5degfilter"
orf_prob_best=".ORF_prob.best"
orf_prob=".ORF_prob"
orf_fasta=".ORF_seqs.fa"
no_orf=".no_ORF.txt"
best_orf=".best_orf"
classification="_classification.5degfilter"
pb_gene=".pb_gene"

cd ../data/BC_ranked_isoforms
PWD=$(pwd)

allorfs="*.ORF_prob.tsv"

echo "Current Working Directory is = " $PWD
echo "allorfs                      = " $allorfs

# loop through 
for file in $allorfs; do
    name="${file%%.*}"

    name_merge5_corrected_5degfilter=$name$merge5$corrected$degfilter
    name_merge5_corrected_5degfilter_orf_prob_best_tsv=$name$merge5$corrected$degfilter$orf_prob_best$tsv
    name_merge5_corrected_5degfilter_orf_prob_tsv=$name$merge5$corrected$degfilter$orf_prob$tsv
    name_merge5_corrected_5degfilter_orf_fasta=$name$merge5$corrected$degfilter$orf_fasta
    name_merge5_corrected_5degfilter_no_orf=$name$merge5$corrected$degfilter$no_orf
    name_merge5_corrected_5degfilter_best_orf_tsv=$name$merge5$corrected$degfilter$best_orf$tsv
    name_merge5_corrected_5degfilter_fasta=$name$merge5$corrected$degfilter$fasta
    name_merge5_corrected_5degfilter_classification_tsv=$name$merge5$corrected$degfilter$classification$tsv
    name_merge5_corrected_5degfilter_classification_pb_gene_tsv=$name$merge5$corrected$degfilter$classification$pb_gene$tsv
    name_merge5_corrected_gtf=$name$merge5$corrected$gtf
    
    echo "file                                                    = " $file
    echo "name                                                        = " $name
    echo "name_merge5_corrected_5degfilter                            = " $name_merge5_corrected_5degfilter
    echo "name_merge5_corrected_5degfilter_orf_prob_best_tsv          = " $name_merge5_corrected_5degfilter_orf_prob_best_tsv
    echo "name_merge5_corrected_5degfilter_orf_prob_tsv               = " $name_merge5_corrected_5degfilter_orf_prob_tsv
    echo "name_merge5_corrected_5degfilter_orf_fasta                  = " $name_merge5_corrected_5degfilter_orf_fasta
    echo "name_merge5_corrected_5degfilter_no_orf                     = " $name_merge5_corrected_5degfilter_no_orf
    echo "name_merge5_corrected_5degfilter_best_orf_tsv               = " $name_merge5_corrected_5degfilter_best_orf_tsv
    echo "name_merge5_corrected_5degfilter_fasta                      = " $name_merge5_corrected_5degfilter_fasta
    echo "name_merge5_corrected_5degfilter_classification_tsv         = " $name_merge5_corrected_5degfilter_classification_tsv
    echo "name_merge5_corrected_5degfilter_classification_pb_gene_tsv = " $name_merge5_corrected_5degfilter_classification_pb_gene_tsv
    echo "name_merge5_corrected_gtf                                   = " $name_merge5_corrected_gtf
    
    docker run -v $PWD:$PWD -w $PWD -it gsheynkmanlab/proteogenomics-base:v1.0 orf_calling.py \
	   --orf_coord $name_merge5_corrected_5degfilter_orf_prob_tsv \
	   --orf_fasta $name_merge5_corrected_5degfilter_orf_fasta \
	   --gencode_gtf gencode.v32.primary_assembly.annotation.gtf \
	   --sample_gtf $name_merge5_corrected_gtf \
	   --pb_gene $name_merge5_corrected_5degfilter_classification_pb_gene_tsv \
	   --classification $name_merge5_corrected_5degfilter_classification_tsv \
	   --sample_fasta $name_merge5_corrected_5degfilter_fasta \
	   --num_cores 8 \
	   --output $name_merge5_corrected_5degfilter_best_orf_tsv
    
done

