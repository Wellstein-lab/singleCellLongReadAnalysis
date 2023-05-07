#!/bin/bash
#Then we run CPAT
#/*--------------------------------------------------
#CPAT
# * CPAT is a bioinformatics tool to predict an RNA’s coding probability 
# * based on the RNA sequence characteristics. 
# * To achieve this goal, CPAT calculates scores of sequence-based features 
# * from a set of known protein-coding genes and background set of non-coding genes.
# *     ORF size
# *     ORF coverage
# *     Fickett score
# *     Hexamer usage bias
# * 
# * CPAT will then builds a logistic regression model using these 4 features as 
# * predictor variables and the “protein-coding status” as the response variable. 
# * After evaluating the performance and determining the probability cutoff, 
# * the model can be used to predict new RNA sequences.
# *
# * https://cpat.readthedocs.io/en/latest/
#---------------------------------------------------*/
#https://github.com/sheynkman-lab/Long-Read-Proteogenomics/blob/23e345dafb0ef90e479cac94a29e3d702472e370/main.nf#L810-L867
#
# get the Human hexamer and logit models needed by cpat
#
#--------------------------------
#  name: run_cpat.sh
#  parameters: $1 - the directory (can be absolute or relative)
#              containing the data desired.
#  output: cpat


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

# Get the files needed by cpat
#wget https://zenodo.org/record/5703754/files/Human_Hexamer.tsv
#wget https://zenodo.org/record/5703754/files/Human_logitModel.RData.gz

human_hexamer="Human_Hexamer.tsv"
human_logitmodel="Human_logitModel.RData"

cpat_output="_cpat.out"
cpat_error="_cpat.error"

#allfiltered="*_corrected.5degfilter.fasta"
allfiltered="*corrected.fasta"

# loop through 
for file in $allfiltered; do
    name="${file%%.fasta}"

#    name_merge5_corrected_5degfilter=$name$merge5$corrected$degfilter
#    name_merge5_corrected_5degfilter_cpat_output=$name$merge5$corrected$degfilter$cpat_output
#    name_merge5_corrected_5degfilter_cpat_error=$name$merge5$corrected$degfilter$cpat_error

    name_cpat_output=$name$cpat_output
    name_cpat_error=$name$cpat_error
    
    echo "file                = " $file
    echo "name                = " $name"_cpat"
    echo "name_cpat_output    = " $name_cpat_output
    echo "name_cpat_error     = " $name_cpat_error
    
    docker run -v $PWD:$PWD -w $PWD -it gsheynkmanlab/cpat:addr cpat.py \
       -x $human_hexamer \
       -d $human_logitmodel \
       -g $file \
       --min-orf=50 \
       --top-orf=50 \
       -o $name \
       1> $name_cpat_output \
       2> $name_cpat_error

done

