#!/bin/bash
#
#
#----------------------------------

allseqsaa="*.ORF_seqs_aa.fa"
seqsaa="_linear_aa.fa"

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

    while read line ; do
	if [ "${line:0:1}" == ">" ];
 	   then echo -e "\n"$line ;
	else echo $line | tr -d '\n' ;
	fi ;
    done < $file | tail -n+2 > $outname
    
done

    
	  
