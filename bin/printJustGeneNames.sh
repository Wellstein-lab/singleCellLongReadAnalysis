#!/bin/bash
#/*--------------------------------------------------
# print Just Gene Names
#---------------------------------------------------*/


cd ../data/BC_ranked_isoforms
PWD=$(pwd)

allfilteredgenepred="filtered*transcript_exons_only.genePred"
end="_filtered_genePred_names.txt"
sorted="_filtered_genePred_names_sorted_unique.csv"
echo "Current Working Directory is = " $PWD

# loop through all filtered genePred files and print the gene name
for file in $allfilteredgenepred; do
    name="${file%%.*}"
    name_end=$name$end
    name_sorted=$name$sorted
    
    echo "file        = " $file
    echo "name        = " $name
    echo "name_end    = " $name_end
    echo "name_sorted = "$name_sorted

    # column 12 contains the gene name
    awk -F "\t" '{print $12}' $file > $name_end
    sort -u $name_end > $name_sorted
    
done
