#!/bin/bash
#
#
#----------------------------------

allseqsaa="*.ORF_seqs_aa.fa"
seqsaa="_seqs_aa.fa"

for file in $allseqsaa; do

    # simple to understand - but need to know the ending
    name2="$(basename "$file" .csv)"
    # more difficult to parse -- no need to know the ending -- so what is happening
    # What we want is the prefix of the file so we can name all of the files
    # Bash has special letters used to do that
    # $() means evaluate this and then the rest of the line
    # ${} means expand the name first and then the line -- since file is a variable
    #    We need this.
    name="${file%%.*}"
    outname=$name$seqsaa
    
    echo "file    = " $file
    echo "name    = " $name
    echo "outname = " $outname

    seqkit grep -n -f DDX5.essential.evidence.ids.csv     $file >"DDX5."$outname
    seqkit grep -n -f HNRNPA1.essential.evidence.ids.csv  $file >"HNRNPA1."$outname
    seqkit grep -n -f HNRNPF.essential.evidence.ids.csv   $file >"HNRNPF."$outname
    seqkit grep -n -f PNF1.essential.evidence.ids.csv     $file >"PNF1."$outname
    seqkit grep -n -f SMG1.essential.evidence.ids.csv     $file >"SMG1."$outname
    seqkit grep -n -f SRSF5.essential.evidence.ids.csv    $file >"SRSF5."$outname
    seqkit grep -n -f SRSF7.essential.evidence.ids.csv    $file >"SRSF7."$outname

done

    
	  
