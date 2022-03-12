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

    # simple to understand - but need to know the ending
    name2="$(basename "$file" .csv)"
    # more difficult to parse -- no need to know the ending -- so what is happening
    # What we want is the prefix of the file so we can name all of the files
    # Bash has special letters used to do that
    # $() means evaluate this and then the rest of the line
    # ${} means expand the name first and then the line -- since file is a variable
    #    We need this.
    name="${file%%.*}"
    
    bam=".bam"
    sam=".sam"
    sortedbam=".sorted.bam"
    sortedsam=".sorted.sam"
    merge5collapsedabundance="_merge5.collapsed.abundance.txt"
    merge5collapsedgff="_merge5.collapsed.gff"
    merge5collapsedgffunfuzzy="_merge5.collapsed.gff.unfuzzy"
    merge5collapsedgrouptxt="_merge5.collapsed.group.txt"
    merge5collapsedgrouptxtunfuzzy="_merge5.collapsed.group.txt.unfuzzy"
    merge5collapsedrepfa="_merge5.collapsed.rep.fa"
    merge5ignoredids="_merge5.ignored_ids.txt"
    pbids="_pbids.csv"
    ccsids="_ccsids.csv"
    ccsfasta="_ccsids.fasta"
    file_pbids="$name$pbids"
    file_ccsids="$name$ccsids"
    file_ccsids_fasta="$name$ccsfasta"
    file_ccsids_bam="$name$bam"
    file_ccsids_sorted_bam="$name$sortedbam"
    file_ccsids_sam="$name$sam"
    file_ccsids_sorted_sam="$name$sortedsam"
    file_ccsids_merge5="$name$merge5"
    file_ccids_merge5collapsedabundance="$name$merge5collapsedabundance"
    file_ccids_merge5collapsedgff="$name$merge5collapsedgff"
    file_ccids_merge5collapsedgffunfuzzy="$name$merge5collapsedgffunfuzzy"
    file_ccids_merge5collapsedgrouptxt="$name$merge5collapsedgrouptxt"
    file_ccids_merge5collapsedgrouptxtunfuzzy="$name$merge5collapsedgrouptxtunfuzzy"
    file_ccids_merge5collapsedrepfa="$name$merge5collapsedrepfa"
    file_merge5ignoredids="$name$merge5ignoredids"
    
    echo "name=" $name
    echo "file=" $file
    echo "file_pbids        = " $file_pbids
    echo "file_ccsids       = " $file_ccsids
    echo "file_ccsids_fasta = " $file_ccsids_fasta
    
    cut -f 2 $file > $file_pbids
    cut -f 1 $file > $file_ccsids


    # in Bash - spaces are important!  without a space around the string, it interrupts as a command
    if [[ $file_ccsids == "Alin"* ]]; then 
	echo "extracting Alin clusters and making " $file_ccsids_fasta
        seqkit grep -n -f $file_ccsids $Alin_neg_fasta > $file_ccsids_fasta
    else
	echo "no match to Alin"
    fi
    
    if [[ $file_ccsids == "Blin"* ]]; then
	echo "extracting Blin clusters and making " $file_ccsids_fasta
	seqkit grep -n -f $file_ccsids $Blin_neg_fasta > $file_ccsids_fasta
    else
	echo "no match to Blin"
    fi
    
    if [[ $file_ccsids == "B_BM_tot"* ]]; then
	echo "extracting B_BM_tot clusters and making " $file_ccsids_fasta
	seqkit grep -n -f $file_ccsids $B_BM_tot_fasta > $file_ccsids_fasta
    else
	echo "no match to B_BM_tot"
    fi

done

    
	  
