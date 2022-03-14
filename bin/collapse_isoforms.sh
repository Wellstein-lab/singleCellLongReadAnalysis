#!/bin/bash
#
# collapse_isoforms
#

gtf=".gtf"
fasta=".fasta"
merge5=".merge5.collapsed"
filtered="filtered_"
corrected="_corrected"
classification="_classification.txt"
allclassified="*classification.txt"
allfiltered="filtered_*.tsv"
degfilter=".5degfilter"
classification_tsv="_classification.tsv"

cd ../data/BC_ranked_isoforms

PWD=$(pwd)
echo "Current Working Directory is = " $PWD
echo "allfiltered                  = " $allfiltered

# loop through 
for file in $allfiltered; do
    name="${file%%.*}"
    filtered_name_merge5_corrected_5degfilter=$name$merge5$corrected$degfilter
    filtered_name_merge5_corrected_gtf=$name$merge5$corrected$gtf
    filtered_name_merge5_corrected_fasta=$name$merge5$corrected$fasta
    filtered_name_merge5_corrected_5degfilter_fasta=$name$merge5$corrected$degfilter$fasta
    filtered_name_merge5_classification_tsv=$name$merge5$classification_tsv
    

    echo "filtered_name_merge5_corrected_5degfilter       = "$filtered_name_merge5_corrected_5degfilter
    echo "filtered_name_merge5_corrected_fasta            = "$filtered_name_merge5_corrected_fasta
    echo "filtered_name_merge5_corrected_gtf              = "$filtered_name_merge5_corrected_gtf
    echo "filtered_name_merge5_corrected_5degfilter_fasta = "$filtered_name_merge5_corrected_5degfilter_fasta
    echo "filtered_name_merge5_classification_tsv         = "$filtered_name_merge5_classification_tsv

#    docker run --rm -v $PWD:$PWD -w $PWD -it gsheynkmanlab/proteogenomics-base collapse_isoforms.py --name $name_merge5 --sqanti_gtf $filtered_name_merge5_corrected_gtf --sqanti_fasta $filtered_name_merge5_corrected_fasta
    docker run --rm -v $PWD:$PWD -w $PWD -it gsheynkmanlab/proteogenomics-base collapse_classification.py --name $filtered_name_merge5_corrected_5degfilter --collapsed_fasta $filtered_name_merge5_corrected_5degfilter_fasta --classification $filtered_name_merge5_classification_tsv 
    
done

