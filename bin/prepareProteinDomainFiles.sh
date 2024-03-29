#!/bin/bash
#---------------------------------------------------------------------
#
# name:  prepareProteinDomainFiles.sh
#
# philosophy - always an elements of style approach
#
#  input - process - output
#
#
# input:   the directory of protein files
#
# process:
#
#  step 1 - countDomainPerSequence.sh
#           input: Experiment Directory
#                  protein domain info from Uniprot
#
#                  arranged as:
#
#                  Protein:Domain_Name:Amino Acid Range: Amino Acid Sequence
#
#          process
#
#              Using grep -- count the number of reads that match the experiment_file
#                  focus only on exact matches
#
#          output
#
#             Protein:Domain_Name:Amino Acid Range: Amino Acid Sequence and read count
#
#   step 5 - add an ID
#
#   step 6 - sort the individual files
#
#   step 7 - normalize the individual files
#
#   step 8 - break normalized files into IJC, SJC files
#          - coordinate file master will be ID and coordinates as it is
#          - norm IJC file will be ID IJC
#          - norm SJC file will be ID SJC
#
#   step 9 - create SE.coordinates.matrix.txt
#                   SE.IJC.w.coordinates.matrix.txt
#                   SE.SJC.w.coordinates.matrix.txt
#
#            join all the files together by ID
#            as each file is joined add to the header the name.
#            finish with adding the header to the final matrix
#            all.SE.IJC.txt
#            all.SE.SJC.txt
#
#   step 10 - create the SE.coordinates.bed file for display
#            in a genome browser, such as the UCSC genome browser
#            
#
# output:  the coordinates and gene identifiers, strand for SE
#---------------------------------------------------------------------

# get to the directory where output will be put
cd $1

# what directory holds the experiment files
experiment_dir=$2

# what directory holds the protein files
protein_dir=$3

tmp="tmp.txt"
PWD=$(pwd)
echo "Current Working Directory is = " $PWD

#
#   step 1 - count the Domains per sequence for each protein in the file
#
for file in "$protein_dir/*.txt"; do
    # make the name of the output file contain the protein
    name="${file%-*.txt}"
    
    echo "file                = " $file
    echo "name                = " $name
    echo "name_se             = " $name_se
    
    cut -f 2-11,13-14 $file > $name_se

    # remove the header
    tail -n +2 $name_se > $tmp && mv $tmp $name_se
done

#
#   step 2 - create the master union file by cat'ing these files together
#
allSEend="*.SE.txt"
allSE="SE.all.txt"

cat $allSEend > $allSE

#
#   step 3 - unique sort master union file (all.SE.txt)
#
# sort unique
#   our column structure has changed
#          col 1 - GeneID
#          col 2 - geneSymbol
#          col 3 - chr
#          col 4 - strand
#          col 5 - exonStart_0base
#          col 6 - exonEnd
#          col 7 - upstreamES
#          col 8 - upstreamEE
#          col 9 - downstreamES
#          col 10 - downstreamEE
#          col 11 - IJC_SAMPLE_1
#          col 12 - SJC_SAMPLE_2
#

allSorted="all.SE.sorted.txt"
sort -u -k3,3 -k4,4 -k5,5 -k6,6 -k7,7 -k8,8 -k9,9 -k10,10 $allSE > $allSorted


#
#   step 4 - cut the last two columns from sorted unique union file
#
# cut the last two columns - because they do not represent the sample specific data
allCut="all.SE.sorted.cut.txt"

cut -f 1-10 $allSorted > $allCut

#
#   step 5 - add an ID
#
# add an ID - counting each row
# afterwards we have the following order in the allUnionFile
#
#          col 1 - ID new -- adding now
#          col 2 - GeneID
#          col 3 - geneSymbol
#          col 4 - chr
#          col 5 - strand
#          col 6 - exonStart_0base
#          col 7 - exonEnd
#          col 8 - upstreamES
#          col 9 - upstreamEE
#          col 10 - downstreamES
#          col 11 - downstreamEE

seCoordinatesFile="SE.coordinates.matrix.txt"
nl $allCut > $seCoordinatesFile

#
# step 6 - sort the individual files
#
# sort the individual files (here we keep the sample count)
sortedSEend=".sorted.SE.txt"

for file in $allSEend; do
    name="${file%.SE.txt}"
    name_sorted_se=$name$sortedSEend

    echo "file                = " $file
    echo "name                = " $name
    echo "name_sorted_se      = " $name_sorted_se

    sort -u -k3,3 -k4,4 -k5,5 -k6,6 -k7,7 -k8,8 -k9,9 -k10,10 $file > $name_sorted_se

done

#
# step 7 - normalize the individual files
#
allSortedSE="*.sorted.SE.txt"
normSEend=".sorted.norm.SE.txt"

for file in $allSortedSE; do
    name="${file%.sorted.SE.txt}"
    name_norm_se=$name$normSEend
    
    echo "file                = " $file
    echo "name                = " $name
    echo "name_norm_se        = " $name_norm_se

    awk -f "../bin/match_se.awk" $file $seCoordinatesFile > $name$normSEend

done

#   step 8 - break normalized files into IJC, SJC files
#          - coordinate file master will be ID and coordinates as it is
#          - norm IJC file will be ID IJC
#          - norm SJC file will be ID SJC
#
allNormSE="*.sorted.norm.SE.txt"
ijc=".IJC.txt"
sjc=".SJC.txt"

for file in $allNormSE; do
    name="${file%.sorted.norm.SE.txt}"
    name_ijc=$name$ijc
    name_sjc=$name$sjc

    echo "file                = " $file
    echo "name                = " $name
    echo "name_ijc            = " $name_ijc
    echo "name_sjc            = " $name_sjc

    cut -f 1,12 $file > $name_ijc
    cut -f 1,13 $file > $name_sjc

done

    
#   step 9 - create final matrix
#            join all the files together by ID
#            as each file is joined add to the header the name.
#            finish with adding the header to the final matrix
#            all.SE.IJC.txt
#            all.SE.SJC.txt
#
seCoordinatesFile="SE.coordinates.matrix.txt"

IJC_matrix="SE.IJC.matrix.txt"
IJC_matrix_csv="SE.IJC.matrix.csv"
IJC_w_coordinates_matrix="SE.IJC.w.coordinates.matrix.txt"
IJC_w_coordinates_matrix_csv="SE.IJC.w.coordinates.matrix.csv"
SJC_matrix="SE.SJC.matrix.txt"
SJC_matrix_csv="SE.SJC.matrix.csv"
SJC_w_coordinates_matrix="SE.SJC.w.coordinates.matrix.txt"
SJC_w_coordinates_matrix_csv="SE.SJC.w.coordinates.matrix.csv"
allIJC="*.IJC.txt"
allSJC="*.SJC.txt"
IJCend=".IJC.txt"
SJCend=".SJC.txt"

#
# need tmp files for temporality
#
tmp_IJC="tmp_IJC.txt"
tmp_SJC="tmp_SJC.txt"
tmp_coordinates_IJC="tmp_coord_IJC.txt"
tmp_coordinates_SJC="tmp_coord_SJC.txt"


#
# headers
#
header_file="SE.header.txt"
coordinates_header_file="SE.coordinates.header.txt"

header="ID"
coordinates_header="ID	GeneID	geneSymbol	chr	strand	exonStart_0base	exonEnd	upstreamES	upstreamEE	downstreamES	downstreamEE"
tab="	"

#
# Populate the IJC and SJC matrices with the ID from the coordinate file
#
cut -f 1 $seCoordinatesFile > $IJC_matrix
cut -f 1 $seCoordinatesFile > $SJC_matrix

cp $seCoordinatesFile $IJC_w_coordinates_matrix
cp $seCoordinatesFile $SJC_w_coordinates_matrix


for file in $allIJC; do
    name="${file%.IJC.txt}"
    header=$header$tab$name
    coordinates_header=$coordinates_header$tab$name

    IJCfile=$name$IJCend
    SJCfile=$name$SJCend
    
    echo "IJCfile             = " $IJCfile
    echo "SJCfile             = " $SJCfile
    echo "name                = " $name
    echo "header              = " $header
    echo "coordinates_header  = " $coordinates_header
    
    join -1 1 -2 1 $IJC_matrix $IJCfile > $tmp_IJC
    join -1 1 -2 1 $SJC_matrix $SJCfile > $tmp_SJC

    join -1 1 -2 1 $IJC_w_coordinates_matrix $IJCfile > $tmp_coordinates_IJC
    join -1 1 -2 1 $SJC_w_coordinates_matrix $SJCfile > $tmp_coordinates_SJC
    
    cp $tmp_IJC $IJC_matrix
    cp $tmp_SJC $SJC_matrix

    cp $tmp_coordinates_IJC $IJC_w_coordinates_matrix
    cp $tmp_coordinates_SJC $SJC_w_coordinates_matrix
    
done

#
# Add header
#
echo $header > $header_file
echo $coordinates_header > $coordinates_header_file

cat $header_file $IJC_matrix > $tmp_IJC
cat $header_file $SJC_matrix > $tmp_SJC

cat $coordinates_header_file $IJC_w_coordinates_matrix > $tmp_coordinates_IJC
cat $coordinates_header_file $SJC_w_coordinates_matrix > $tmp_coordinates_SJC

cp $tmp_IJC $IJC_matrix
cp $tmp_SJC $SJC_matrix

cp $tmp_coordinates_IJC $IJC_w_coordinates_matrix
cp $tmp_coordinates_SJC $SJC_w_coordinates_matrix

#
# clean up
#
rm $tmp_IJC
rm $tmp_SJC

rm $tmp_coordinates_IJC
rm $tmp_coordinates_SJC

#
# add a csv version for files
#
sed  's/ /,/g' < $IJC_matrix > $IJC_matrix_csv
sed  's/ /,/g' < $SJC_matrix > $SJC_matrix_csv

sed 's/ /,/g' < $IJC_w_coordinates_matrix > $IJC_w_coordinates_matrix_csv
sed 's/ /,/g' < $SJC_w_coordinates_matrix > $SJC_w_coordinates_matrix_csv
#
#   step 10 - create the SE.coordinates.bed file for display
#            in a genome browser, such as the UCSC genome browser
#            

echo "track name=rMATS_SE description=\"rMATS SE Events DS-AML\"" > SE.coordinates.bed
awk -f ../bin/make_bed_se.awk SE.coordinates.matrix.txt >> SE.coordinates.bed
