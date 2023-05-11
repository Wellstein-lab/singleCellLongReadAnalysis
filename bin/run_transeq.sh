#
# transeq is part of the EMBOSS suite
# Package available through Anaconda
#
# conda install -c bioconda emboss
#
# installs the package
#
# transeq is capable of translating all 6 frames
#
# given we are using the output of cpat -- it provides us
# the likely ORF frame specific so we will just use
# the provided frame details -- aka frame 1 is adequate
#
# but we can explore this in the future.
#
#!/bin/bash
cd $1
PWD=$(pwd)

cpat_aa="_cpat.ORF_seqs_aa.fa"

allorfseqsfa=$(ls *.ORF_seqs.fa)

for file in $allorfseqsfa; do
    name="${file%%.ORF_seqs.fa}"
    out=$name$cpat_aa
    echo "file is " $file
    echo "name is " $name
    echo "out  is " $out
    # translate to amino acid sequence
	       
    transeq -sequence $PWD/$file \
	    --outseq $PWD/$out \
	    -frame 1
    
done
