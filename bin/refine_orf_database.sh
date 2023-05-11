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
classification="5degfilter_classification.5degfilter.tsv"
coding_score_cutoff=0.364

#cd ../data/BC_ranked_isoforms
cd $1
PWD=$(pwd)

#
# the input to refine is the output of orf_calling
#  confusing to call it best_orf - this is the output of orf_calling.sh is the *.best_orf.tsv
# 
allbestorfs="*.best_orf.tsv"

echo "Current Working Directory is = " $PWD
echo "allbestorfsorfs              = " $allbestorfs
echo "coding_score_cutoff          = " $coding_score_cutoff

# loop through 
for file in $allbestorfs; do
    name="${file%%.best_orf.tsv}"

    name_merge5_corrected_5degfilter_refined=$name$merge5$corrected$degfilter$refined
    name_merge5_corrected_5degfilter_fasta=$name$merge5$corrected$degfilter$fasta
    
    echo "file                                            = " $file
    echo "name                                            = " $name
    echo "name_merge5_corrected_5degfilter_refined        = " $name_merge5_corrected_5degfilter_refined
    echo "name_merge5_corrected_5degfilter_fasta          = " $name_merge5_corrected_5degfilter_fasta
    
    docker run -v $PWD:$PWD -w $PWD -it gsheynkmanlab/proteogenomics-base:v1.0 refine_orf_database.py \
	   --name         $name_merge5_corrected_5degfilter_refined \
	   -io            $file \
	   -if            $name_merge5_corrected_5degfilter_fasta \
           -cut           0.364

done

