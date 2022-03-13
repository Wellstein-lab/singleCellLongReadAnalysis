#!/bin/bash
#
#  Filename:  run_sqanti3_long_read_clusters.sh
#
#Filter SQANTI
#* Filter SQANTI results based on several criteria
#* - protein coding only
#*      PB transcript aligns to a GENCODE-annotated protein coding gene.
#* - percent A downstream
#*      perc_A_downstreamTTS : percent of genomic "A"s in the downstream 20 bp window.
# *      If this number if high (> 80%), the 3' end have arisen from intra-priming during the RT step 
# * - RTS stage
# *      RTS_stage: TRUE if one of the junctions could be an RT template switching artifact.
# * - Structural Category
# *      keep only transcripts that have a isoform structural category of:
# *        -novel_not_in_catalog
# *        -novel_in_catalog
# *        -incomplete-splice_match
# *        -full-splice_match 
#---------------------------------------------------*/
#https://github.com/sheynkman-lab/Long-Read-Proteogenomics/blob/23e345dafb0ef90e479cac94a29e3d702472e37#0/main.nf#L655-L731


cd ../data/BC_ranked_isoforms
PWD=$(pwd)
echo "Current Working Directory is = " $PWD

gtf=".gtf"
fasta=".fasta"
merge5=".merge5.collapsed"
corrected=".corrected"
classified=".classification.txt"
allclassfied="*classification.txt"

echo "gencode_primary_assembly_annotation = " $gencode_primary_assembly_annotation
echo "reference_genome                    = " $reference_genome
echo "allclassified                       = " $allclassified

cd ../data/BC_ranked_isoforms

PWD=$(pwd)

echo "current working directory = " $PWD

# loop through 
for file in $allcollapsedrepfasta; do
    name="${file%%.*}"
    name_merge5=$name$merge5
    name_merge5_classification=$name$merge5$classification
    name_merge5_corrected_fasta=$name$merge5$corrected$fasta
    name_merge5_corrected_gtf=$name$merge5$corrected$gtf

    echo "name_merge5_classification  = " $name_merge5_classification
    echo "name_merge5_corrected_fasta = " $name_merge5_corrected_fasta
    echo "name_merge5_corrected_gtf   = " $name_merge5_corrected_gtf
    
    docker run --rm -v $PWD:$PWD -w $PWD -it gsheynkmanlab/proteogenomics-base filter_sqanti.py --sqanti_classification $name_merge5_classification --sqanti_corrected_fasta $name_merge5_corrected_fasta --sqanti_corrected_gtf $name_merge5_corrected_gtf --protein_coding_genes protein_coding_genes.txt --ensg_gene ensg_gene.tsv --filter_protein_coding yes --filter_intra_polyA yes --filter_template_switching yes --percent_A_downstream_threshold 95 --structural_categories_level strict --minimum_illumina_coverage 3

done

