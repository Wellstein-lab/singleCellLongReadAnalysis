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
orf_prob_tsv=".ORF_prob.tsv"
orf_fasta=".ORF_seqs.fa"
no_orf=".no_ORF.txt"
best_orf_tsv=".best_orf.tsv"
classification_tsv="_classification.tsv"
pb_gene_tsv=".pb_gene.tsv"

#cd ../data/BC_ranked_isoforms
cd $1
PWD=$(pwd)

allorfs="*.ORF_prob.tsv"

echo "Current Working Directory is = " $PWD
echo "allorfs                      = " $allorfs

# loop through 
for file in $allorfs; do
    name="${file%%.ORF_prob.tsv}"
    name_orf_prob_best_tsv=$name$orf_prob_best
    name_orf_prob_tsv=$name$orf_prob_tsv
    name_orf_fasta=$name$orf_fasta
    name_no_orf=$name$no_orf
    name_best_orf_tsv=$name$best_orf_tsv
    name_fasta=$name$fasta
    name_pb_gene_tsv=$name$pb_gene_tsv
    name_classification_tsv=$name$classification_tsv
    name_gtf=$name$gtf
    
    echo "file                            = " $file
    echo "name                            = " $name
    echo "name_orf_prob_best_tsv          = " $name_orf_prob_best_tsv
    echo "name_orf_prob_tsv               = " $name_orf_prob_tsv
    echo "name_orf_fasta                  = " $name_orf_fasta
    echo "name_no_orf                     = " $name_no_orf
    echo "name_best_orf_tsv               = " $name_best_orf_tsv
    echo "name_fasta                      = " $name_fasta
    echo "name_classification_tsv         = " $name_classification_tsv
    echo "name_pb_gene_tsv                = " $name_pb_gene_tsv
    echo "name_gtf                        = " $name_gtf
    
    docker run -v $PWD:$PWD -w $PWD -it gsheynkmanlab/proteogenomics-base:v1.0 orf_calling.py \
	   --orf_coord $name_orf_prob_tsv \
	   --orf_fasta $name_orf_fasta \
	   --gencode_gtf gencode.v32.primary_assembly.annotation.gtf \
	   --sample_gtf $name_gtf \
	   --pb_gene $name_pb_gene_tsv \
	   --classification $name_classification_tsv \
	   --sample_fasta $name_fasta \
	   --num_cores 8 \
	   --output $name_best_orf_tsv
    
done

