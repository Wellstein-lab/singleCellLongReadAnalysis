#!/bin/bash

for file in *clust*fasta; do
    name="${file%%.*}"
    cpat="_cpat"
    cpat_output="cpat.output"
    cpat_error="cpat.error"
    namecpat=$name$cpat
    namecpatout=$name$cpat_output
    namecpaterror=$name$cpat_error

    echo "running " $file
    echo "        name = " $name
    echo "        namecpat = " $namecpat
    
    docker run --rm -v $PWD:$PWD -w $PWD -it gsheynkmanlab/cpat:addr cpat.py -x Human_Hexamer.tsv -d Human_logitModel.RData -g $file --min-orf=50 --top-orf=50 -o $namecpat 1>$namecpatout 2>$namecpaterror

done

