#!/bin/bash
#
#  Get the reference fasta and the gtf files
#
#----------------------------------------------------------
#wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/GRCh38.primary_assembly.genome.fa.gz;
#gunzip GRCh38.primary_assembly.genome.fa.gz
#wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.primary_assembly.annotation.gtf.gz
#gunzip gencode.v32.primary_assembly.annotation.gtf.gz
#wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.pc_transcripts.fa.gz
#gunzip gencode.v32.pc_transcripts.fa.gz
#wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.pc_translations.fa.gz
#gunzip gencode.v32.pc_translations.fa.gz
#
#----------------------------------------------------------
#
# Generate reference tables
#
# using two input files Protein-coding transcript sequences -
#  which are the nucleotide sequences of coding transcripts on the reference chromosomes.
#  Transcript biotypes: protein_coding, nonsense_mediated_decay, non_stop_decay, IG_*_gen,
#  TR*_gene, polymorphic_pseudogene
#
# First command - generate-reference-tables
#
# using container gsheynkmanlab/generate-reference-tables -- found here
# https://github.com/sheynkman-lab/Long-Read-Proteogenomics/blob/main/modules/generate_reference_tables/Dockerfile
#
#----------------------------------------------------------
PWD=$(pwd)
docker run --rm -v $PWD:$PWD -w $PWD -it gsheynkmanlab/generate-reference-tables prepare_reference_tables.py --gtf gencode.v32.primary_assembly.annotation.gtf --fa gencode.v32.pc_transcripts.fa --ensg_gene ensg_gene.tsv --enst_isoname enst_isoname.tsv --gene_ensp gene_ensp.tsv --gene_isoname gene_isoname.tsv --isoname_lens isoname_lens.tsv --gene_lens gene_lens.tsv --protein_coding_genes protein_coding_genes.txt

#----------------------------------------------------------
#
# Make_gencode_database.py
#
# using container gsheynkmanlab/proteogenomics-base -- found here
# https://github.com/sheynkman-lab/Long-Read-Proteogenomics/blob/main/modules/make_gencode_database/Dockerfile
#
#----------------------------------------------------------
PWD=$(pwd)
docker run --rm -v $PWD:$PWD -w $PWD -it gsheynkmanlab/proteogenomics-base make_gencode_database.py   --gencode_fasta gencode.v32.pc_translations.fa   --protein_coding_genes protein_coding_genes.txt   --output_fasta gencode_protein.fasta   --output_cluster gencode_isoname_clusters.tsv

