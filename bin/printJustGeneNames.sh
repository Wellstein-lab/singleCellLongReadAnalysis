#!/bin/bash
#/*--------------------------------------------------
# print Just Gene Names
#---------------------------------------------------*/


cd ../data/BC_ranked_isoforms
PWD=$(pwd)

allfilteredgenepred="filtered*transcript_exons_only.genePred"
end="filtered_genePred_names.txt"
echo "Current Working Directory is = " $PWD

# loop through all filtered genePred files and print the gene name
for file in $allfilteredgenepred; do
    name="${file%%.*}"
    name_end=$name$end
    
    echo "file     = " $file
    echo "name     = " $name
    echo "name_end = " $name_end

    # column 12 contains the gene name
    awk -F "\t" '{print $12}' $file > $name_end
    
done
