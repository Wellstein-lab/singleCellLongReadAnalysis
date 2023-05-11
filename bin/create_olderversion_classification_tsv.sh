#!/bin/bash
#
#  2023 May 7 -
#  Author:  Anne Deslattes Mays, PhD
#
#  file: create_olderversion_classification_tsv.sh
#
#  Justification:
#
#  Tewer output from Sqanti seems to have an extra column which contains ORF
#  to run the older script orf_calling from the Sheynkman Lab's proteogenomics workflow
#
#  a few lines to strip the last column of the file.
#
#  How was this discovered?
#  When I didn't find the file 'classification.tsv' I wondered what the difference was between
#  `classification.txt` and the older version file named `classification.tsv`
#
#  A nice way to sort this out is to turn column headers into rows and inspect their names side by side.
#
#  By "heading" the first line and then literally flip it on its head by "tr"anslating the tabs (\t) with newlines (\n).
#  As follows:
#
#   (eos) head -1 Lneg_7_2_hq_mmp_.collapsed.filtered.rep_classification.txt | tr "\t" "\n" 
#   isoform
#   chrom
#   strand
#   length
#   exons
#   structural_category
#   associated_gene
#   associated_transcript
#   ref_length
#   ref_exons
#   diff_to_TSS
#   diff_to_TTS
#   diff_to_gene_TSS
#   diff_to_gene_TTS
#   subcategory
#   RTS_stage
#   all_canonical
#   min_sample_cov
#   min_cov
#   min_cov_pos
#   sd_cov
#   FL
#   n_indels
#   n_indels_junc
#   bite
#   iso_exp
#   gene_exp
#   ratio_exp
#   FSM_class
#   coding
#   ORF_length
#   CDS_length
#   CDS_start
#   CDS_end
#   CDS_genomic_start
#   CDS_genomic_end
#   predicted_NMD
#   perc_A_downstream_TTS
#   seq_A_downstream_TTS
#   dist_to_cage_peak
#   within_cage_peak
#   dist_to_polya_site
#   within_polya_site
#   polyA_motif
#   polyA_dist
#   ORF_seq
#
# and then when we looked at an older file that was in the subdirectory "BC_ranked_isoforms"
#
# (eos) head -1  ../BC_ranked_isoforms/filtered_Blin_neg_filt_ranked_BC_clust9_B_ccsids.merge5.collapsed_classification.tsv| tr "\t" "\n" 
#   (eos) head -1 Lneg_7_2_hq_mmp_.collapsed.filtered.rep_classification.txt | tr "\t" "\n" 
#   isoform
#   chrom
#   strand
#   length
#   exons
#   structural_category
#   associated_gene
#   associated_transcript
#   ref_length
#   ref_exons
#   diff_to_TSS
#   diff_to_TTS
#   diff_to_gene_TSS
#   diff_to_gene_TTS
#   subcategory
#   RTS_stage
#   all_canonical
#   min_sample_cov
#   min_cov
#   min_cov_pos
#   sd_cov
#   FL
#   n_indels
#   n_indels_junc
#   bite
#   iso_exp
#   gene_exp
#   ratio_exp
#   FSM_class
#   coding
#   ORF_length
#   CDS_length
#   CDS_start
#   CDS_end
#   CDS_genomic_start
#   CDS_genomic_end
#   predicted_NMD
#   perc_A_downstream_TTS
#   seq_A_downstream_TTS
#   dist_to_cage_peak
#   within_cage_peak
#   dist_to_polya_site
#   within_polya_site
#   polyA_motif
#   polyA_dist
#
#  we see the last column "ORF_seq" is not there.  I will investigate later where it is generated and what the conditions are,
#  could be useful later.
#
#
# unfortunately awk defaults to printing spaces in stead of tabs.
#
cd $1
files=$(ls *rep_classification.txt)
for file in $files; do
    name="${file%%.txt}"
    in=$name".txt"
    out=$name".tsv"
    out2=$name".tsv2"
    # cut the last column
    awk 'NF{--NF};1' <$in > $out
    # change spaces to tabs
    cat $out | tr " " "\t" > $out2
    cp $out2 $out
    rm $out2
done

