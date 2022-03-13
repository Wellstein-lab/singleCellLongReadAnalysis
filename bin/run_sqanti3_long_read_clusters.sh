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

# install through bioconda channel agat (another gff/gtf analysis toolkit)
#conda install -c bioconda agat

cd ../data/BC_ranked_isoforms
PWD=$(pwd)
echo "Current Working Directory is = " $PWD

gencode_primary_assembly_annotation="../gencode.v32.primary_assembly.annotation.gtf"
reference_genome="../GRCh38.primary_assembly.genome.fa"
allccsidsfasta="*_ccsids.fasta"

echo "gencode_primary_assembly_annotation = " $gencode_primary_assembly_annotation
echo "reference_genome                    = " $reference_genome
echo "allccsidsfasta                      = " $allccsidsfasta

for file in $allccsidsfasta; do
    name="${file%%.*}"
    fasta=".fasta"
    sam=".sam"
    bam=".bam"
    sortedbam=".sorted.bam"
    sortedsam=".sorted.sam"
    bed=".bed"
    genepred=".genepred"
    merge5collapsedabundance=".merge5.collapsed.abundance.txt"
    merge5collapsedgff=".merge5.collapsed.gff"
    merge5collapsedgtf=".merge5.collapsed.gtf"
    merge5collapsedgffunfuzzy=".merge5.collapsed.gff.unfuzzy"
    merge5collapsedgrouptxt=".merge5.collapsed.group.txt"
    merge5collapsedgrouptxtunfuzzy=".merge5.collapsed.group.txt.unfuzzy"
    merge5collapsedrepfa=".merge5.collapsed.rep.fa"
    merge5ignoredids=".merge5.ignored_ids.txt"
    file_ccsids_fasta="$name$fasta"
    file_ccsids_bam="$name$bam"
    file_ccsids_sam="$name$sam"
    file_ccsids_sorted_bam="$name$sortedbam"
    file_ccsids_sorted_sam="$name$sortedsam"
    file_ccsids_merge5="$name$merge5"
    file_ccids_merge5collapsedabundance="$name$merge5collapsedabundance"
    file_ccids_merge5collapsedgff="$name$merge5collapsedgff"
    file_ccids_merge5collapsedgtf="$name$merge5collapsedgtf"
    file_ccids_merge5collapsedgffunfuzzy="$name$merge5collapsedgffunfuzzy"
    file_ccids_merge5collapsedgrouptxt="$name$merge5collapsedgrouptxt"
    file_ccids_merge5collapsedgrouptxtunfuzzy="$name$merge5collapsedgrouptxtunfuzzy"
    file_ccids_merge5collapsedrepfa="$name$merge5collapsedrepfa"
    file_ccids_merge5collapsedrepfa_sam="$name$merge5collapsedrepfa$sam"
    file_ccids_merge5collapsedrepfa_bam="$name$merge5collapsedrepfa$bam"
    file_ccids_merge5collapsedrepfa_sortedbam="$name$merge5collapsedrepfa$sortedbam"
    file_ccids_merge5collapsedrepfa_sortedsam="$name$merge5collapsedrepfa$sortedsam"
    file_ccids_merge5collapsedrepfa_bed="$name$merge5collapsedrepfa$sorted$bed"
    file_merge5ignoredids="$name$merge5ignoredids"

    echo "file                                      = " $file
    echo "name                                      = " $name
    echo "file_ccsids_fasta                         = " $file_ccsids_fasta
    echo "file_ccsids_bam                           = " $file_ccsids_bam
    echo "file_ccsids_sorted_bam                    = " $file_ccsids_sorted_bam
    echo "file_ccsids_sam                           = " $file_ccsids_sam
    echo "file_ccsids_sorted_sam                    = " $file_ccsids_sorted_sam
    echo "file_ccsids_merge5                        = " $file_ccsids_merge5
    echo "file_ccids_merge5collapsedabundance       = " $file_ccids_merge5collapsedabundance
    echo "file_ccids_merge5collapsedgff             = " $file_ccids_merge5collapsedgff
    echo "file_ccids_merge5collapsedgtf             = " $file_ccids_merge5collapsedgtf
    echo "file_ccids_merge5collapsedgffunfuzzy      = " $file_ccids_merge5collapsedgffunfuzzy
    echo "file_ccids_merge5collapsedgrouptxt        = " $file_ccids_merge5collapsedgrouptxt
    echo "file_ccids_merge5collapsedgrouptxtunfuzzy = " $file_ccids_merge5collapsedgrouptxtunfuzzy
    echo "file_ccids_merge5collapsedrepfa           = " $file_ccids_merge5collapsedrepfa
    echo "file_ccids_merge5collapsedrepfa_sam       = " $file_ccids_merge5collapsedrepfa_sam
    echo "file_ccids_merge5collapsedrepfa_bam       = " $file_ccids_merge5collapsedrepfa_bam
    echo "file_ccids_merge5collapsedrepfa_sortedbam = " $file_ccids_merge5collapsedrepfa_sortedbam
    echo "file_ccids_merge5collapsedrepfa_sortedsam = " $file_ccids_merge5collapsedrepfa_sortedsam
    echo "file_ccids_merge5collapsedrepfa_bed       = " $file_ccids_merge5collapsedrepfa_bed
    echo "file_merge5ignoredids                     = " $file_merge5ignoredids

 
#    agat_convert_sp_gff2gtf.pl --gff $file_ccids_merge5collapsedgff -o $file_ccids_merge5collapsedgtf
    docker run --rm -v $PWD:$PWD -w $PWD -it gsheynkmanlab/sqanti3 sqanti3_qc.py $file_ccsids_merge5collapsedgtf $gencode_primary_assembly_annotation $reference_genome --skipORF -o $name --fl_count $file_ccsids_merge5collapsedabundance --gtf -c 100
#    minimap2 -ax map-pb $reference_genome $file_ccids_merge5collapsedrepfa > 
#    bedtools bamtobed -bed12 -i $file_ccsids_bam > $file_ccsids_bed
#    bedToGenePred $file_ccsids_bed $file_ccsids_genepred
#    genePredToGtf "file" $file_ccsids_genepred $file_ccsids_gtf
#    gffread -E $file_ccsids_gtf -o $file_ccsids_gff

done

    
	  
