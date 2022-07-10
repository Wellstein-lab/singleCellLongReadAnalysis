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
    
    # DDX5 isoforms begin with "MSGY"
    grep -B1 "^MSGY" $file | sed '/^--$/d' | sed "s%^>\(.*\)%>$name$underscore\1%" >$outdir"DDX5.MSGY."$name$seqsaa
    
    # DDX5 isoforms end in GYSQ
    grep -B1 "GYSQ$" $file | sed '/^--$/d' | sed "s%^>\(.*\)%>$name$underscore\1%"  >$outdir"DDX5.GYSQ."$name$seqsaa
    grep -B1 "KRGG$" $file | sed '/^--$/d' | sed "s%^>\(.*\)%>$name$underscore\1%"  >$outdir"DDX5.KRGG."$name$seqsaa

    # MSKS is the beginning of both the HNRNPA1 isoforms
    grep -B1 "^MSKS" $file | sed '/^--$/d' | sed "s%^>\(.*\)%>$name$underscore\1%"   > $outdir"HNRNPA1.MSKS."$name$seqsaa
    
    # all HNRNPA1 isoforms end with GRRF
    grep -B1 "GRRF$" $file | sed '/^--$/d' | sed "s%^>\(.*\)%>$name$underscore\1%"  >$outdir"HNRNPA1.GRRF."$name$seqsaa
    
    # QGGY is the ending of one of the HNRNPA1 features
    grep -B1 "QGGY$" $file | sed '/^--$/d' | sed "s%^>\(.*\)%>$name$underscore\1%"  >$outdir"HNRNPA1.QGGY."$name$seqsaa

    # the HNRNPF isoform starts with MMLG
    grep -B1 "^MMLG" $file | sed '/^--$/d' | sed "s%^>\(.*\)%>$name$underscore\1%"  >$outdir"HNRNPA1.MMLG."$name$seqsaa

    # the HNRNPF isoform ends with GGYD
    grep -B1 "GGYD$" $file | sed '/^--$/d' | sed "s%^>\(.*\)%>$name$underscore\1%"  >$outdir"HNRNPF.GGYD."$name$seqsaa

    # the PFN1 isoform starts with MAGW
    grep -B1 "^MAGW" $file | sed '/^--$/d' | sed "s%^>\(.*\)%>$name$underscore\1%"  >$outdir"PFN1.MAGW."$name$seqsaa

    # the PFN1 isoform ends with RSQY
    grep -B1 "RSQY$" $file | sed '/^--$/d' | sed "s%^>\(.*\)%>$name$underscore\1%"  >$outdir"PFN1.RSQY."$name$seqsaa

    # the SMG1 isoforms 1 starts with MSRR
    grep -B1 "^MSRR" $file | sed '/^--$/d' | sed "s%^>\(.*\)%>$name$underscore\1%"  >$outdir"SMG1.MSRR."$name$seqsaa

    # the SMG1 isoforms 2 starts with MSYS
    grep -B1 "^MSYS" $file | sed '/^--$/d' | sed "s%^>\(.*\)%>$name$underscore\1%"  >$outdir"SMG1.MSYS."$name$seqsaa

    # the SMG1 isoforms 3 starts with MWAL
    grep -B1 "^MWAL" $file | sed '/^--$/d' | sed "s%^>\(.*\)%>$name$underscore\1%"  >$outdir"SMG1.MWAL."$name$seqsaa

    # the SMG1 isoforms 4 starts with MKKL
    grep -B1 "^MKKL" $file | sed '/^--$/d' | sed "s%^>\(.*\)%>$name$underscore\1%"  >$outdir"SMG1.MKKL."$name$seqsaa
    
    # all the 4 SMG1 isoforms all end in TAWV only present in clust5 Blin_neg_B (clust4 Alin_neg_AB)
    grep -B1 "TAWV$" $file | sed '/^--$/d' | sed "s%^>\(.*\)%>$name$underscore\1%"  >$outdir"SMG1.TAWV."$name$seqsaa

    # the SRSF5 isoforms 1,2 and 4 start with MSGC
    grep -B1 "^MSGC" $file | sed '/^--$/d' | sed "s%^>\(.*\)%>$name$underscore\1%"  >$outdir"SRSF5.MSGC."$name$seqsaa
    
    # there are 4 isoforms for SRSF5 with two different endings
    grep -B1 "DSGN$" $file | sed '/^--$/d' | sed "s%^>\(.*\)%>$name$underscore\1%"  >$outdir"SRSF5.DSGN."$name$seqsaa
    grep -B1 "GWLH$" $file | sed '/^--$/d' | sed "s%^>\(.*\)%>$name$underscore\1%"  >$outdir"SRSF5.GWLH."$name$seqsaa

    # all 4 SRSF7 isoforms start with MSRY
    grep -B1 "^MSRY" $file | sed '/^--$/d' | sed "s%^>\(.*\)%>$name$underscore\1%"  >$outdir"SRSF7.MSRY."$name$seqsaa

    # there are 4 SRSF7 isoforms 
    grep -B1 "ERMD$" $file | sed '/^--$/d' | sed "s%^>\(.*\)%>$name$underscore\1%"  >$outdir"SRSF7.ERMD."$name$seqsaa
    grep -B1 "NLRR$" $file | sed '/^--$/d' | sed "s%^>\(.*\)%>$name$underscore\1%"  >$outdir"SRSF7.NLRR."$name$seqsaa
    grep -B1 "RYLF$" $file | sed '/^--$/d' | sed "s%^>\(.*\)%>$name$underscore\1%"  >$outdir"SRSF7.RYLF."$name$seqsaa
 
done
