#!/bin/bash
#
#  Filename:  run_sqanti3_long_read_clusters.sh
#
#And then sqanti3
#/*--------------------------------------------------
#SQANTI3
#* https://github.com/ConesaLab/SQANTI3
#* Corrects any errors in alignment from IsoSeq3 and
#* classifies each accession in relation to the reference
#* genome
#
#
# Container here Dockerfile
#  /Users/annedeslattesmays/Desktop/projects/Long-Read-Proteogenomics/modules/sqanti3:
#---------------------------------------------------*/
#https://github.com/sheynkman-lab/Long-Read-Proteogenomics/blob/23e345dafb0ef90e479cac94a29e3d702472e370/main.nf#L591-L653

cd ../data/BC_ranked_isoforms
PWD=$(pwd)
echo "Current Working Directory is = " $PWD

gencode_primary_assembly_annotation=gencode.v32.primary_assembly.annotation.gtf
reference_genome=GRCh38.primary_assembly.genome.fa

for file in $allccsidsfasta; do
    name="$(basename "$file" .fasta)"
    bam=".bam"
    bed=".bed"
    genepred=".genepred"
    gff=".gff"
    gtf=".gtf"
    file_ccsids_bam="$name$bam"
    file_ccsids_bed="$name$bed"
    file_ccsids_genepred="$name$genepred"
    file_ccsids_gtf="$name$gtf"
    file_ccsids_gff="$name$gff"

        
    echo "name=" $name
    echo "file=" $file
    echo "file_ccsids_bam      = " $file_ccsids_bam
    echo "file_ccsids_bed      = " $file_ccsids_bed
    echo "file_ccsids_genepred = " $file_ccsids_genepred
    echo "file_ccsids_gtf      = " $file_ccsids_gtf
    echo "file_ccsids_gff      = " $file_ccsids_gff

    docker run --rm -v $PWD:$PWD -w $PWD -it gsheynkmanlab/sqanti3 sqanti3_qc.py $file_ccsids_gff  $gencode_primary_assembly_annotation $reference_genome --skipORF -o $name --fl_count BM_merged3.sort.flnc_BC.merge5.collapsed.abundance.txt --gtf -c 100
    minimap2 -ax map-pb $reference_genome $file | samtools sort -O BAM - > $file_ccsids_bam
    bedtools bamtobed -bed12 -i $file_ccsids_bam > $file_ccsids_bed
    bedToGenePred $file_ccsids_bed $file_ccsids_genepred
    genePredToGtf "file" $file_ccsids_genepred $file_ccsids_gtf
    gffread -E $file_ccsids_gtf -o $file_ccsids_gff

done

    
	  
