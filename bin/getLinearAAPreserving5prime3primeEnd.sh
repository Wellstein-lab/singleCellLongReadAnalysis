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

allseqsaa="B*_linear_aa.fa"
seqsaa="_linear_aa.fa"
outdir="./"
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
    # DDX5 isoforms end in GYSQ
    DDX5_MSGY_GYSQ="DDX5.MSGY.GYSQ"
    DDX5_MSGY_KRGG="DDX5.MSGY.KRGG"
    grep -B1 "^MSGY" $file | sed '/^--$/d' | grep -B1 "GYSQ$" | sed '/^--$/d' | sed "s%^>\(.*\)%>$DDX5_MSGY_GYSQ$underscore$name$underscore\1%"    >$outdir$DDX5_MSGY_GYSQ$dot$name$seqsaa
    grep -B1 "^MSGY" $file | sed '/^--$/d' | grep -B1 "KRGG$" | sed '/^--$/d' | sed "s%^>\(.*\)%>$DDX5_MSGY_KRGG$underscore$name$underscore\1%"    >$outdir$DDX5_MSGY_KRGG$dot$name$seqsaa

    # HNRNPA1 isoforms begin with "MSKS"
    # HNRNPA1 isoforms end with "GRRF" 
    HNRNPA1_MSKS_GRRF="HNRNPA1.MSKS.GRRF"
    grep -B1 "^MSKS" $file | sed '/^--$/d' | grep -B1 "GRRF$" | sed '/^--$/d' | sed "s%^>\(.*\)%>$HNRNPA1_MSKS_GRRF$underscore$name$underscore\1%" >$outdir$HNRNPA1_MSKS_GRRF$dot$name$seqsaa
    
    # HNRNPA1 isoforms begin with "MSKS"
    # HNRNPA1 isoforms end with "QGGY" 
    HNRNPA1_MSKS_QGGY="HNRNPA1.MSKS.QGGY"
    grep -B1 "^MSKS" $file | sed '/^--$/d' | grep -B1 "QGGY$" | sed '/^--$/d' | sed "s%^>\(.*\)%>$HNRNPA1_MSKS_QGGY$underscore$name$underscore\1%" >$outdir$HNRNPA1_MSKS_QGGY$dot$name$seqsaa

    # HNRNPF isoforms starts with "MMLG"
    # HNRNPF isoforms end with "GGYD"
    HNRNPF_MMLG_GGYD="HNRNPF.MMLG.GGYD"
    grep -B1 "^MMLG" $file | sed '/^--$/d' | grep -B1 "GGYD$" | sed '/^--$/d' | sed "s%^>\(.*\)%>$HNRNPF_MMLG_GGYD$underscore$name$underscore\1%"  >$outdir$HNRNPF_MMLG_GGYD$dot$name$seqsaa

    # PFN1 isoforms start with MAGW
    # PFN1 isoforms end with RSQY
    PFN1_MAGW_RSQY="PFN1.MAGW_RSQY"
    grep -B1 "^MAGW" $file | sed '/^--$/d' | grep -B1 "RSQY$" | sed '/^--$/d' | sed "s%^>\(.*\)%>$PFN1_MAGW_RSQY$underscore$name$underscore\1%"   >$outdir$PFN1_MAGW_RSQY$dot$name$seqsaa

    # SMG1 isoform 1 starts with MSRR
    # SMG1 isoform 1 ends with TAMV
    SMG1_MSRR_TAMV="SMG1.MSRR.TAMV"
    grep -B1 "^MSRR" $file | sed '/^--$/d' | grep -B1 "TAWV$" | sed '/^--$/d' | sed "s%^>\(.*\)%>$SMG1_MSRR_TAMV$underscore$name$underscore\1%"    >$outdir$SMG1_MSRR_TAMV$dot$name$seqsaa

    # SMG1 isoform 2 starts with MSYS
    # SMG1 isoform 2 ends with TAMV
    SMG1_MSYS_TAMV="SMG1.MSYS.TAMV"
    grep -B1 "^MSYS" $file | sed '/^--$/d' | grep -B1 "TAWV$" | sed '/^--$/d' | sed "s%^>\(.*\)%>$SMG1_MSYS_TAMV$underscore$name$underscore\1%"    >$outdir$SMG1_MSYS_TAMV$dot$name$seqsaa

    # SMG1 isoform 3 starts with MWAL
    # SMG1 isoform 3 ends with TAMV
    SMG1_MWAL_TAMV="SMG1.MWAL.TAMV"
    grep -B1 "^MWAL" $file | sed '/^--$/d' | grep -B1 "TAWV$" | sed '/^--$/d' | sed "s%^>\(.*\)%>$SMG1_MWAL_TAMV$underscore$name$underscore\1%"    >$outdir$SMG1_MWAL_TAMV$dot$name$seqsaa

    # SMG1 isoform 4 starts with MKKL
    # SMG1 isoform 4 ends with TAMV
    SMG1_MKKL_TAMV="SMG1.MKKL.TAMV"
    grep -B1 "^MKKL" $file | sed '/^--$/d' | grep -B1 "TAWV$" | sed '/^--$/d' | sed "s%^>\(.*\)%>$SMG1_MKKL_TAMV$underscore$name$underscore\1%"    >$outdir$SMG1_MKKL_TAMV$dot$name$seqsaa

    # SRSF5 isoforms 1,2 and 4 start with MSGC
    # SRSF5 isoforms end with DSGN
    SRSF5_MSGC_DSGN="SRSF5.MSGC.DSGN"
    grep -B1 "^MSGC" $file | sed '/^--$/d' | grep -B1 "DSGN$" | sed '/^--$/d' | sed "s%^>\(.*\)%>$SRSF5_MSGC_DSGN$underscore$name$underscore\1%"   >$outdir$SRSF5_MSGC_DSGN$dot$name$seqsaa
   
    # SRSF5 isoforms 1,2 and 4 start with MSGC
    # SRSF5 isoforms end with GWLH
    SRSF5_MSGC_GWLH="SRSF5.MSGC.GWLH"
    grep -B1 "^MSGC" $file | sed '/^--$/d' | grep -B1 "GWLH$" | sed '/^--$/d' | sed "s%^>\(.*\)%>$SRSF5_MSGC_GWLH$underscore$name$underscore\1%"   >$outdir$SRSF5_MSGC_GWLH$dot$name$seqsaa

    # SRSF7 isoforms start with MSRY
    # SRSF7 isoforms 1, 4 end with ERMD
    SRSF7_MSRY_ERMD="SRSF7.MSRY.ERMD"
    grep -B1 "^MSRY" $file | sed '/^--$/d' | grep -B1 "ERMD$" | sed '/^--$/d' |sed "s%^>\(.*\)%>$SRSF7_MSRY_ERMD$underscore$name$underscore\1%"  >$outdir$SRSF7_MSRY_ERMD$dot$name$seqsaa

    # SRSF7 isoforms start with MSRY
    # SRSF7 isoforms 2 end with NLRR
    SRSF7_MSRY_NLRR="SRSF7.MSRY.NLRR"
    grep -B1 "^MSRY" $file | sed '/^--$/d' | grep -B1 "NLRR$" | sed '/^--$/d' |sed "s%^>\(.*\)%>$SRSF7_MSRY_NLRR$underscore$name$underscore\1%"  >$outdir$SRSF7_MSRY_NLRR$dot$name$seqsaa

    # SRSF7 isoforms start with MSRY
    # SRSF7 isoforms 3 end with RYLF
    SRSF7_MSRY_RYLF="SRSF7.MSRY.RYLF"
    grep -B1 "^MSRY" $file | sed '/^--$/d' | grep -B1 "RYLF$" | sed '/^--$/d' |sed "s%^>\(.*\)%>$SRSF7_MSRY_RYLF$underscore$name$underscore\1%"  >$outdir$SRSF7_MSRY_RYLF$dot$name$seqsaa
 
done
