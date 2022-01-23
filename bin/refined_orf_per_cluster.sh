#!/bin/bash
#Then we make a refined database (even though we do not have mass spec data - because it does some collapsing of the transcripts into the shared proteins (or ORFs).
#/*--------------------------------------------------
#Refined DB Generation
# * - Filteres ORF database to only include accessions 
# *   with a CPAT coding score above a threshold (default 0.0)
# * - Filters ORFs to only include ORFs that have a stop codon 
# * - Collapses transcripts that produce the same protein
# *   into one entry, keeping a base accession (first alphanumeric).
# *   Abundances of transcripts (CPM) are collapsed during this process.
#---------------------------------------------------*/
#https://github.com/sheynkman-lab/Long-Read-Proteogenomics/blob/23e345dafb0ef90e479cac94a29e3d702472e370/main.nf#L923-L962

for file in "*ORF_prob.best.tsv"; do
    name="${file%%.*}"
    fasta=".ORF_seqs.fa"
    name_fasta="$name$fasta"
    echo "running " $file
    echo "refining to " $name_fasta
    docker run --rm -v $PWD:$PWD -w $PWD -it gsheynkmanlab/proteogenomics-base refine_orf_database.py --name $name -io $file -if $name_fasta -cut 0.364

done
