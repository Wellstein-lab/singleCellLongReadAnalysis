#!/bin/bash
#
#  Filename:  align_long_read_clusters.sh
#
#  Conda install
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

#  Conda install
conda install -c bioconda minimap2 -y
conda install -c bioconda bedtools -y
conda install -c bioconda ucsc-bedtogenepred -y
conda install -c bioconda ucsc-genepredtogtf -y
#
# samtools is trying to find a dynamic library not loaded
# Using solution https://github.com/conda/conda/issues/8103 suggested by ReneKat
#
conda install -c bioconda samtools=1.9 -y
conda install -c bioconda gffread -y

allccsidsfasta="*ccsids.fasta"
reference_genome="../GRCh38.primary_assembly.genome.fa"

cd ../data/BC_ranked_isoforms
PWD=$(pwd)
echo "Current Working Directory is = " $PWD

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

    minimap2 -ax map-pb $reference_genome $file | samtools sort -O BAM - > $file_ccsids_bam
    bedtools bamtobed -bed12 -i $file_ccsids_bam > $file_ccsids_bed
    bedToGenePred $file_ccsids_bed $file_ccsids_genepred
    genePredToGtf "file" $file_ccsids_genepred $file_ccsids_gtf
    gffread -E $file_ccsids_gtf -o $file_ccsids_gff

done

    
	  
