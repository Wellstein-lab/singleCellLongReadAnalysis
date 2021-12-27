#!/bin/bash
cd ../data/BC_ranked_isoforms

allcsv="*.csv"

for file in $allcsv; do
    
    name="${file%%.*}"
    
    echo "name=" $name
    echo "file=" $file
    
    pbids="_pbids.csv"
    ccsids="_ccsids.csv"
    ccsfasta="_ccsids.fasta"
    file_pbids="$name$pbids"
    file_ccsids="$name$ccsids"
    file_ccsids_fasta="$name$ccsfasta"
    Alin_neg_fasta="Alin_neg.flnc_BC.fasta"
    Blin_neg_fasta="Blin_neg.flnc_BC.fasta"
    B_BM_tot_fasta="B_BM_tot.flnc.BC.fasta"
    
    echo "file_pbids        = " $file_pbids
    echo "file_ccsids       = " $file_ccsids
    echo "Alin_neg_fasta    = " $Alin_neg_fasta
    echo "Blin_neg_fasta    = " $Blin_neg_fasta
    echo "B_BM_tot_fasta    = " $B_BM_tot_fasta
    
    cut -f 2 $file > $file_pbids
    cut -f 1 $file > $file_ccsids
done

for file in "Alin*"; do
    echo "extracting Alin clusters and making " $file_ccsids_fasta 
    seqkit grep -n -f $file_ccsids $Alin_neg_fasta > $file_ccsids_fasta
done

for file in "Blin*"; do
    echo "extracting Blin clusters and making " $file_ccsids_fasta
    seqkit grep -n -f $file_ccsids $Blin_neg_fasta > $file_ccsids_fasta
done

for file in "B_BM_tot*"; do
    echo "extracting B_BM_tot clusters and making " $file_ccsids_fasta
    seqkit grep -n -f $file_ccsids $B_BM_tot_fasta > $file_ccsids_fasta
done

    
	  
