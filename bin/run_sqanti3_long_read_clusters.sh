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

gencode_primary_assembly_annotation="gencode.v32.primary_assembly.annotation.gtf"
reference_genome="GRCh38.primary_assembly.genome.fa"
allcollapsedrepfasta="*rep.fa"
gtf=".gtf"


echo "gencode_primary_assembly_annotation = " $gencode_primary_assembly_annotation
echo "reference_genome                    = " $reference_genome
echo "allcollapsedrepfasta                = " $allcollapsedrepfasta

for file in $allcollapsedrepfasta; do
    name="${file%%.*}"
    merge5=".merge5.collapsed"
    abundance=".abundance.txt"
    star_junction=".SJ.out.tab"
    name_merge5=$name$merge5
    name_merge5_abundance=$name$merge5$abundance
    name_merge5_star_junction=$name$merge5$star_junction
    name_merge5_gtf=$name$merge5$gtf

    docker run --rm -v $PWD:$PWD -w $PWD -it gsheynkmanlab/sqanti3 sqanti3_qc.py --aligner_choice minimap2 --skipORF -o $name_merge5 --fl_count $name_merge5_abundance --gtf $name_merge5_gtf $gencode_primary_assembly_annotation $reference_genome
#    minimap2 -ax map-pb $reference_genome $file_ccids_merge5collapsedrepfa > 
#    bedtools bamtobed -bed12 -i $file_ccsids_bam > $file_ccsids_bed
#    bedToGenePred $file_ccsids_bed $file_ccsids_genepred
#    genePredToGtf "file" $file_ccsids_genepred $file_ccsids_gtf
#    gffread -E $file_ccsids_gtf -o $file_ccsids_gff
done

    
	  
