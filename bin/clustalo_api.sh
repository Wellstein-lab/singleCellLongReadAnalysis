#!/bin/bash
# Running clustal Omega using REST API from embl (url at end)
# INPUT:  order - first argument - User can specify order [input|aligned]
#         file  - the input file containing sequences to align
# Documenting the runnning of the alignment - trying to sort out differences
#
order=$1
file=$2
#
# read the file into a single variable
#
quote=\'
#sequence_command=$'\'sequence='
#
# String needs to be between two quotes
# so begin with the quote
#
entireFileAsString=`cat $file`
#entireFileAsString=$entireFileAsString$quote
echo "$entireFileAsString"
curl --form 'email=adeslat@scitechcon.org' \
     --form 'dealign=false' \
     --form 'mbed=yes' \
     --form 'mbediterations=5' \
     --form 'iterations=5' \
     --form 'outfmt=clustal_num' \
     --form 'order='$order \
     --form 'stype=protein' \
     --form $entireFileAsString \
     https://www.ebi.ac.uk/Tools/services/rest/clustalo/run
