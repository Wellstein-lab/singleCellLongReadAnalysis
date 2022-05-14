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
# assumes you are in the directory with the data.
# assumes you have transeq
#
PWD=$(pwd)
allprecollapse="*_ccsids.fasta"

for file in $allprecollapse; do
    name="${file%%.*}"

    cpat_out="_cpat.out"
    cpat_error="_cpat.error"
    cpat="_cpat"
    cpat_orf_dna="_cpat.ORF_seqs.fa"
    cpat_aa="cpat.ORF_seqs_aa.fa"

    # get the best ORFs
    docker run -v $PWD:$PWD -it gsheynkmanlab/cpat:addr cpat.py \
	   -x $PWD/Human_Hexamer.tsv \
	   -d $PWD/Human_logitModel.RData \
	   -g $PWD/$file \
	   --min-orf=50 \
	   --top-orf=50 \
	   -o $PWD/$name$cpat \
	   1>$PWD/$name$cpat_out \
	   2>$PWD/$name$cpat_error 


    # translate to amino acid sequence
    transeq -sequence $PWD/$name$cpat_orf_dna \
	    --outseq $PWD/$name$cpat_aa \
	    -frame 1
    
done

