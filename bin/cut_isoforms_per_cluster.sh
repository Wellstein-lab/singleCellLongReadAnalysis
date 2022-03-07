#!/bin/bash
#
#  Filename:  cut_isoforms_per_cluster.sh
#
#  Setup:     Make a symbolic link to the directory containing the csv files for the ranks of the isoforms
#             Put the files (copy or symbolic link)
#                Alin_neg.flnc_BC.fasta
#                Blin_neg.flnc_BC.fasta
#                B_BM_tot.flnc.BC.fasta
#
#  Pseudo code:
#             The rank files contain 
#
#
#----------------------------------
# Get the data
mkdir -p ../data
mkdir -p ../data/BC_ranked_isoforms

cd ../data/BC_ranked_isoforms

#
# Get the Full-length PacBio Consensus Fasta Files for the samples
#
cp /Users/annedeslattesmays/Scitechcon\ Dropbox/Anne\ DeslattesMays/Single_cell_analysis/scRNA\ paper/single\ cell\ Long\ Read\ Protein\ Analysis/Alin_neg.flnc_BC.fasta .
cp /Users/annedeslattesmays/Scitechcon\ Dropbox/Anne\ DeslattesMays/Single_cell_analysis/scRNA\ paper/single\ cell\ Long\ Read\ Protein\ Analysis/Blin_neg.flnc_BC.fasta .
cp /Users/annedeslattesmays/Scitechcon\ Dropbox/Anne\ DeslattesMays/Single_cell_analysis/scRNA\ paper/single\ cell\ Long\ Read\ Protein\ Analysis/B_BM_tot.flnc_BC.fasta .

#
# Get the barcoded ranked clusters
#
cp /Users/annedeslattesmays/Scitechcon\ Dropbox/Anne\ DeslattesMays/Single_cell_analysis/PacBio/single_cell_pipeline/BC-ranked_isoforms/BC_ranked_isoforms/*.csv .

Alin_neg_fasta="Alin_neg.flnc_BC.fasta"
Blin_neg_fasta="Blin_neg.flnc_BC.fasta"
B_BM_tot_fasta="B_BM_tot.flnc_BC.fasta"

echo "Alin_neg_fasta    = " $Alin_neg_fasta
echo "Blin_neg_fasta    = " $Blin_neg_fasta
echo "B_BM_tot_fasta    = " $B_BM_tot_fasta

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
    
    echo "file_pbids        = " $file_pbids
    echo "file_ccsids       = " $file_ccsids
    echo "file_ccsids_fasta = " $file_ccsids_fasta
    
    cut -f 2 $file > $file_pbids
    cut -f 1 $file > $file_ccsids


    # in Bash - spaces are important!  without a space around the string, it interrupts as a command
    if [[ $file_ccsids == "Alin"* ]]; then 
	echo "extracting Alin clusters and making " $file_ccsids_fasta
        seqkit grep -n -f $file_ccsids $Alin_neg_fasta > $file_ccsids_fasta
    elif [[ $file_ccsids == "^Blin"* ]]; then
	echo "extracting Blin clusters and making " $file_ccsids_fasta
	seqkit grep -n -f $file_ccsids $Blin_neg_fasta > $file_ccsids_fasta
    elif [[ $file_ccsids == "^B_BM_tot"* ]]; then
	echo "extracting B_BM_tot clusters and making " $file_ccsids_fasta
	seqkit grep -n -f $file_ccsids $B_BM_tot_fasta > $file_ccsids_fasta
    else
	echo "no match"
    fi

done

    
	  
