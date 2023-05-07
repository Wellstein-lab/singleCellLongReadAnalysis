#!/bin/bash
#/*--------------------------------------------------
#Six-Frame Translation
# * Generates a fasta file of all possible protein sequences
# * derivable from each PacBio transcript, by translating the
# * fasta file in all six frames (3+, 3-). This is used to examine
# * what peptides could theoretically match the peptides found via
# * a mass spectrometry search against GENCODE. 
#---------------------------------------------------*/

gtf=".gtf"
fasta=".fasta"
merge5=".merge5.collapsed"
filtered="filtered_"
corrected="_corrected"
degfilter=".5degfilter"
sixframe=".6frame"
classification_tsv=".5degfilter_classification.5degfilter.tsv"

#cd ../data/BC_ranked_isoforms
cd $1

PWD=$(pwd)
echo "Current Working Directory is = " $PWD
echo "allfiltered                  = " $allfiltered

# pull the docker image
docker pull ghcr.io/adeslatt/six_frame_translation:sha256-ae3c04adb8bdeb613ededc737b555af5740a2329293a3549014d1c8873390d43.sig

#allfiltered="*_corrected.5degfilter.fasta"

# loop through 
#for file in $allfiltered; do
    name="${file%%.*}"

#    name_merge5_corrected_classification=$name$merge5$corrected$classification_tsv
    name_classification_tsv=$2
    name_fasta=$3
    name_sixframe_fasta=$name".6frame.fasta"
    
    echo "file                                            = " $file
    echo "name                                            = " $name
    echo "name_merge5_corrected_classification            = " $name_merge5_corrected_classification
    echo "name_merge5_corrected_5degfilter_fasta          = " $name_merge5_corrected_5degfilter_fasta
    echo "name_merge5_corrected_5degfilter_sixframe_fasta = " $name_merge5_corrected_5degfilter_sixframe_fasta
    
    docker run --rm -v $PWD:$PWD -w $PWD -it sixframe six_frame_translation.py \
	   --iso_annot $name$merge5$corrected$classification_tsv \
	   --ensg_gene ensg_gene.tsv \
	   --sample_fasta $name_merge5_corrected_5degfilter_fasta \
	   --output_fasta $name_merge5_corrected_5degfilter_sixframe_fasta
        
done

