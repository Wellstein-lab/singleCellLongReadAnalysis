#!/bin/bash
cd ../data/BC_ranked_isoforms

Alin_neg_fasta="Alin_neg.flnc_BC.fasta"
echo "Alin_neg_fasta    = " $Alin_neg_fasta
for file in Alin_neg*ccsids.csv; do
    
    name="${file%%.*}"
    fasta=".fasta"
    file_fasta="$name$fasta"
    
    echo "name=" $name
    echo "file=" $file
    echo "file_fasta=" $file_fasta
    
    echo "extracting Alin clusters and making " $file_fasta 
    seqkit grep -n -f $file $Alin_neg_fasta > $file_fasta
done

Blin_neg_fasta="Blin_neg.flnc_BC.fasta"
echo "Blin_neg_fasta    = " $Blin_neg_fasta
for file in Blin_neg*ccsids.csv; do
    
    name="${file%%.*}"
    fasta=".fasta"
    file_fasta="$name$fasta"
    
    echo "name=" $name
    echo "file=" $file
    echo "file_fasta=" $file_fasta
    
    echo "extracting Blin clusters and making " $file_fasta 
    seqkit grep -n -f $file $Blin_neg_fasta > $file_fasta
done

B_BM_tot_fasta="B_BM_tot.flnc_BC.fasta"
echo "B_BM_tot_fasta    = " $B_BM_tot_fasta
for file in B_BM*ccsids.csv; do
    
    name="${file%%.*}"
    fasta=".fasta"
    file_fasta="$name$fasta"
    
    echo "name=" $name
    echo "file=" $file
    echo "file_fasta=" $file_fasta
    
    echo "extracting B_BM_tot clusters and making " $file_fasta 
    seqkit grep -n -f $file $B_BM_tot_fasta > $file_fasta
done

