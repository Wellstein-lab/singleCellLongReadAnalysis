#!/bin/bash

PWD=$(pwd)
cd ../data/BC_ranked_isoforms

allcollapsed="*.collapsed.merge5*"
wrong=".collapsed.merge5"
merge5=".merge5.collapsed"
alignedsam=".Aligned.out.sam"
logfinalout=".Log.final.out"
logout=".Log.out"
logprogressout=".Log.progress.out"
sjout=".SJ.out.tab"

echo "merge5= " $merge5

for file in $allcollapsed; do
    name="${file%%.*}"
    name_merge5=$name$merge5
    name_wrong=$name$wrong
    name_alignedsam=$name$merge5$alignedsam
    name_wrong_alignedsam=$name$wrong$alignedsam
    name_logfinalout=$name$merge5$logfinalout
    name_wrong_logfinalout=$name$wrong$logfinalout
    name_logout=$name$merge5$logout
    name_wrong_logout=$name$wrong$logout
    name_logprogressout=$name$merge5$logprogressout
    name_wrong_logprogressout=$name$wrong$logprogressout
    name_sjout=$name$merge5$sjout
    name_wrong_sjout=$name$wrong$sjout
    
    echo "file                      = " $file
    echo "name                      = " $name
    echo "name_merge5               = " $name_merge5
    echo "name_alignedsam           = " $name_alignedsam
    echo "name_wrong_alignedsam     = " $name_wrong_alignedsam
    echo "name_logfinalout          = " $name_logfinalout
    echo "name_wrong_logfinalout    = " $name_wrong_logfinalout
    echo "name_logout               = " $name_logout
    echo "name_wrong_logout         = " $name_wrong_logout
    echo "name_logprogressout       = " $name_logprogressout
    echo "name_wrong_logprogressout = " $name_wrong_logprogressout
    echo "name_sjout                = " $name_sjout
    echo "name_wrong_sjout          = " $name_wrong_sjout

    # rename the misnamed files
#    mv $name_wrong_alignedsam     $name_alignedsam
#    mv $name_wrong_logfinalout    $name_logfinalout
#    mv $name_wrong_logout         $name_logout
    #    mv $name_wrong_logprogressout $name_logprogressout
    mv $name_wrong_sjout $name_sjout
    
done
