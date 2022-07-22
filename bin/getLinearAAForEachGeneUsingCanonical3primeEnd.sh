#!/bin/bash
#------------------------------------------
#
# Previous step linearized the aa sequence to prepare for this step,
# which is to use a regular expression to extract those molecule sequences
# which match the cannonical 3 prime end of our favorite genes.
#
# Using in this case a very liberal match of end 4 amino acides
# we will use this then to then do the alignment between clusters and
# samples - to get at functional differences
# 
#-----------------------------------------

allseqsaa="*_linear_aa.fa"
seqsaa="_linear_aa.fa"
outdir="../favorite_genes/"
testfile="B_BM_tot_filt_ranked_BC_clust1_B_ccsidscpat_linear_aa.fa"
file=$testfile
underscore="_"
dot="."

for file in $allseqsaa; do

    # simple to understand - but need to know the ending
    # using the function basename - we can strip off the end of the file
    # handy for cleaning up growing filenames
    # What we want is the prefix of the file so we can name all of the files
    # Bash has special letters used to do that
    # $() means evaluate this and then the rest of the line, i.e. use the function
    # ${} means expand the name first and then the line -- since file is a variable

    name="$(basename "$file" _ccsidscpat_linear_aa.fa)"

    outname=$outdir$name$seqsaa
 
    echo "file = " $file
    echo "name = " $name
    echo "outname = " $outname

    # ^ means match at the beginning
    # $ means match at the end
    # | pipe symbole
    # sed '/^--$/d' (piping the result of the grep through sed is required to remove the -- separator that grep adds)
    # always define a variable
    
    # DDX5 isoforms begin with "MSGY"
    DDX5_MSGY="DDX5.MSGY"
    grep -B1 "^MSGY" $file | sed '/^--$/d' | sed "s%^>\(.*\)%>$DDX5_MSGY$underscore$name$underscore\1%"    >$outdir$DDX5_MSGY$dot$name$seqsaa
    
    # DDX5 isoforms end in GYSQ
    DDX5_GYSQ="DDX5.GYSQ"
    grep -B1 "GYSQ$" $file | sed '/^--$/d' | sed "s%^>\(.*\)%>$DDX5_GYSQ$underscore$name$underscore\1%"    >$outdir$DDX5_GYSQ$dot$name$seqsaa
    DDX5_KRGG="DDX5.KRGG"
    grep -B1 "KRGG$" $file | sed '/^--$/d' | sed "s%^>\(.*\)%>$DDX5_KRGG$underscore$name$underscore\1%"    >$outdir$DDX5_KRGG$dot$name$seqsaa

    # MSKS is the beginning of both the HNRNPA1 isoforms
    HNRNPA1_MSKS="HNRNPA1.MSKS"
    grep -B1 "^MSKS" $file | sed '/^--$/d' | sed "s%^>\(.*\)%>$HNRNPA1_MSKS$underscore$name$underscore\1%" >$outdir$HNRNPA1_MSKS$dot$name$seqsaa
    
    # all HNRNPA1 isoforms end with GRRF
    HNRNPA1_GRRF="HNRNPA1.GRRF"
    grep -B1 "GRRF$" $file | sed '/^--$/d' | sed "s%^>\(.*\)%>$HNRNPA1_GRRF$underscore$name$underscore\1%" >$outdir$HNRNPA1_GRRF$dot$name$seqsaa
    
    # QGGY is the ending of one of the HNRNPA1 features
    HNRNPA1_QGGY="HNRNPA1.QGGY"
    grep -B1 "QGGY$" $file | sed '/^--$/d' | sed "s%^>\(.*\)%>$HNRNPA1_QGGY$underscore$name$underscore\1%" >$outdir$HNRNPA1_QGGY$dot$name$seqsaa

    # the HNRNPF isoform starts with MMLG
    HNRNPF_MMLG="HNRNPF.MMLG"
    grep -B1 "^MMLG" $file | sed '/^--$/d' | sed "s%^>\(.*\)%>$HNRNPF_MMLG$underscore$name$underscore\1%"  >$outdir$HNRNPF_MMLG$dot$name$seqsaa

    # the HNRNPF isoform ends with GGYD
    HNRNPF_GGYD="HNRNPF.GGYD"
    grep -B1 "GGYD$" $file | sed '/^--$/d' | sed "s%^>\(.*\)%>$HNRNPF_GGYD$underscore$name$underscore\1%"  >$outdir$HNRNPF_GGYD$dot$name$seqsaa

    # the PFN1 isoform starts with MAGW
    PFN1_MAGW="PFN1.MAGW"
    grep -B1 "^MAGW" $file | sed '/^--$/d' | sed "s%^>\(.*\)%>$PFN1_MAGW$underscore$name$underscore\1%"    >$outdir$PFN1_MAGW$dot$name$seqsaa

    # the PFN1 isoform ends with RSQY
    PFN1_RSQY="PRN1.RSQY"
    grep -B1 "RSQY$" $file | sed '/^--$/d' | sed "s%^>\(.*\)%>$PFN1_RSQY$underscore$name$underscore\1%"    >$outdir$PFN1_RSQY$dot$name$seqsaa

    # the SMG1 isoforms 1 starts with MSRR
    SMG1_MSRR="SMG1.MSRR"
    grep -B1 "^MSRR" $file | sed '/^--$/d' | sed "s%^>\(.*\)%>$SMG1_MSRR$underscore$name$underscore\1%"    >$outdir$SMG1_MSRR$dot$name$seqsaa

    # the SMG1 isoforms 2 starts with MSYS
    SMG1_MSYS="SMG1.MSYS"
    grep -B1 "^MSYS" $file | sed '/^--$/d' | sed "s%^>\(.*\)%>$SMG1_MSYS$underscore$name$underscore\1%"    >$outdir$SMG1_MSYS$dot$name$seqsaa

    # the SMG1 isoforms 3 starts with MWAL
    SMG1_MWAL="SMG1.MWAL"
    grep -B1 "^MWAL" $file | sed '/^--$/d' | sed "s%^>\(.*\)%>$SMG1_MWAL$underscore$name$underscore\1%"    >$outdir$SMG1_MWAL$dot$name$seqsaa

    # the SMG1 isoforms 4 starts with MKKL
    SMG1_MKKL="SMG1.MKKL"
    grep -B1 "^MKKL" $file | sed '/^--$/d' | sed "s%^>\(.*\)%>$SMG1_MKKL$underscore$name$underscore\1%"    >$outdir$SMG1_MKKL$dot$name$seqsaa
    
    # all the 4 SMG1 isoforms all end in TAWV only present in clust5 Blin_neg_B (clust4 Alin_neg_AB)
    SMG1_TAMV="SMG1.TAMV"
    grep -B1 "TAWV$" $file | sed '/^--$/d' | sed "s%^>\(.*\)%>$SMG1_TAMV$underscore$name$underscore\1%"    >$outdir$SMG1_TAMV$dot$name$seqsaa

    # the SRSF5 isoforms 1,2 and 4 start with MSGC
    SRSF5_MSGC="SRSF5.MSGC"
    grep -B1 "^MSGC" $file | sed '/^--$/d' | sed "s%^>\(.*\)%>$SRSF5_MSGC$underscore$name$underscore\1%"   >$outdir$SRSF5_MSGC$dot$name$seqsaa
    
    # there are 4 isoforms for SRSF5 with two different endings
    SRSF5_DSGN="SRSF5.DSGN"
    grep -B1 "DSGN$" $file | sed '/^--$/d' | sed "s%^>\(.*\)%>$SRSF5_DSGN$underscore$name$underscore\1%"   >$outdir$SRSF5_DSGN$dot$name$seqsaa
    SRSF5_GWLH="SRSF5.GWLH"
    grep -B1 "GWLH$" $file | sed '/^--$/d' | sed "s%^>\(.*\)%>$SRSF5_GWLH$underscore$name$underscore\1%"   >$outdir$SRSF5_GWLH$dot$name$seqsaa

    # all 4 SRSF7 isoforms start with MSRY
    SRSF7_MSRY="SRSF7.MSRY"
    grep -B1 "^MSRY" $file | sed '/^--$/d' | sed "s%^>\(.*\)%>$SRSF7_MSRY$underscore$name$underscore\1%"  >$outdir$SRSF7_MSRY$dot$name$seqsaa

    # there are 4 SRSF7 isoforms
    SRSF7_ERMD="SRSF7.ERMD"
    grep -B1 "ERMD$" $file | sed '/^--$/d' | sed "s%^>\(.*\)%>$SRSF7_ERMD$underscore$name$underscore\1%"  >$outdir$SRSF7_ERMD$dot$name$seqsaa
    SRSF7_NLRR="SRSF7.NLRR"
    grep -B1 "NLRR$" $file | sed '/^--$/d' | sed "s%^>\(.*\)%>$SRSF7_NLRR$underscore$name$underscore\1%"  >$outdir$SRSF7_NLRR$dot$name$seqsaa
    SRSF7_RYLF="SRSF7.RYLF"
    grep -B1 "RYLF$" $file | sed '/^--$/d' | sed "s%^>\(.*\)%>$SRSF7_RYLF$underscore$name$underscore\1%"  >$outdir$SRSF7_RYLF$dot$name$seqsaa
 
done
