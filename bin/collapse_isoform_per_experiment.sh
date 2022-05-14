#!/bin/bash
#
#  Filename:  collapse_isoform_per_experiment
#
#  Conda install
#     conda 
#     conda install -c bioconda minimap2 -y
#     conda install -c bioconda bedtools -y
#     conda install -c bioconda ucsc-bedtogenepred -y
#     conda install -c bioconda ucsc-genepredtogtf -y
#
#  Pseudo code:
#             Using minimap2 align each of the fasta sequences to the GRCh38 v32 genome
#             Output is a bam file - with that we can get a gff file
#
#  In fact here is a way to do that from https://github.com/lh3/minimap2/issues/455
#      and thanks to Guilherme Borges Dias @gbdias
#
#      For future googlers,
#            An option to convert spliced alignments from SAM/BAM to gff/gtf:
#
#            1. Align sequences and convert to BAM
#               minimap2 -ax splice --cs target.fa query.fa | samtools sort -O BAM - > alignments.bam
#            2. Convert to BED12 using BEDtools
#               bedtools bamtobed -bed12 -i alignments.bam > alignments.bed
#            3. Convert to genePred using UCSC tools
#               bedToGenePred alignments.bed alignments.genepred
#            4. Convert to GTF2 using UCSC tools
#               genePredToGtf "file" alignments.genepred alignments.gtf
#            5. genePredToGtf has additional options that might be useful in specific use cases.
#
#      All of these tools are available in Bioconda.
#
#             which can be used for sqanti annotation (and maybe tappas)
#
#----------------------------------

# assumes you are in the directory and all items are here

allexperiments="*_BC.fasta"
reference_genome="GRCh38.primary_assembly.genome.fa"

PWD=$(pwd)
echo "Current Working Directory is = " $PWD

for file in $allexeriments; do
    name="$(basename "$file" .fasta)"
    bam=".bam"
    sortedbam=".sorted.bam"
    sam=".sam"
    sortedsam=".sorted.sam"
    merge5=".merge5"
    bed=".bed"
    genepred=".genepred"
    gff=".gff"
    gtf=".gtf"
    file_bam="$name$bam"
    file_sorted_bam="$name$sortedbam"
    file_sam="$name$sam"
    file_sorted_sam="$name$sortedsam"
    file_merge5="$name$merge5"
    file_bed="$name$bed"
    file_genepred="$name$genepred"
    file_gtf="$name$gtf"
    file_gff="$name$gff"

        
    echo "name=" $name
    echo "file=" $file
    echo "file_merge5   = " $file_merge5
    echo "file_bam      = " $file_bam
    echo "file_sam      = " $file_sam
    echo "file_bed      = " $file_bed
    echo "file_genepred = " $file_genepred
    echo "file_gtf      = " $file_gtf
    echo "file_gff      = " $file_gff

    minimap2 -ax map-pb $reference_genome $file > $file_sam
    samtools view -bh $file_sam > $file_bam
    samtools sort $file_bam > $file_sorted_bam
    samtools view -h $file_sorted_bam > $file_sorted_sam
    collapse_isoforms_by_sam.py --input $file -s $file_sorted_sam -c 0.99 -i 0.95 --gen_mol_count -o  $file_merge5
#    bedtools bamtobed -bed12 -i $file_bam > $file_bed
#    bedToGenePred $file_bed $file_genepred
#    genePredToGtf "file" $file_genepred $file_gtf
#    gffread -E $file_gtf -o $file_gff

done
