#
# Let the alignment be for the result
#
jobid=$1
curl -X GET --header 'Accept: text/x-clustalw-alignment' 'https://www.ebi.ac.uk/Tools/services/rest/clustalo/result/'$jobid'/aln-clustal_num' 
