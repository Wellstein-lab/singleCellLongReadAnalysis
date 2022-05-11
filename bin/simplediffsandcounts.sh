#!/bin/bash
#
# simple accounting differences
#


wc -l filtered_B_BM_tot_filt_ranked_BC_clust1_B_ccsids_filtered_genePred_names_sorted_unique.csv filtered_Blin_neg_filt_ranked_BC_clust1_B_ccsids_filtered_genePred_names_sorted_unique.csv 
wc -l filtered_B_BM_tot_filt_ranked_BC_clust2_B_ccsids_filtered_genePred_names_sorted_unique.csv filtered_Blin_neg_filt_ranked_BC_clust2_B_ccsids_filtered_genePred_names_sorted_unique.csv 
wc -l filtered_B_BM_tot_filt_ranked_BC_clust3_B_ccsids_filtered_genePred_names_sorted_unique.csv filtered_Blin_neg_filt_ranked_BC_clust3_B_ccsids_filtered_genePred_names_sorted_unique.csv 
wc -l filtered_B_BM_tot_filt_ranked_BC_clust4_B_ccsids_filtered_genePred_names_sorted_unique.csv filtered_Blin_neg_filt_ranked_BC_clust4_B_ccsids_filtered_genePred_names_sorted_unique.csv 
wc -l filtered_B_BM_tot_filt_ranked_BC_clust5_B_ccsids_filtered_genePred_names_sorted_unique.csv filtered_Blin_neg_filt_ranked_BC_clust5_B_ccsids_filtered_genePred_names_sorted_unique.csv 
wc -l filtered_B_BM_tot_filt_ranked_BC_clust6_B_ccsids_filtered_genePred_names_sorted_unique.csv filtered_Blin_neg_filt_ranked_BC_clust6_B_ccsids_filtered_genePred_names_sorted_unique.csv 
wc -l filtered_B_BM_tot_filt_ranked_BC_clust7_B_ccsids_filtered_genePred_names_sorted_unique.csv filtered_Blin_neg_filt_ranked_BC_clust7_B_ccsids_filtered_genePred_names_sorted_unique.csv 
wc -l filtered_B_BM_tot_filt_ranked_BC_clust8_B_ccsids_filtered_genePred_names_sorted_unique.csv filtered_Blin_neg_filt_ranked_BC_clust8_B_ccsids_filtered_genePred_names_sorted_unique.csv 
wc -l filtered_B_BM_tot_filt_ranked_BC_clust9_B_ccsids_filtered_genePred_names_sorted_unique.csv filtered_Blin_neg_filt_ranked_BC_clust9_B_ccsids_filtered_genePred_names_sorted_unique.csv 
diff filtered_B_BM_tot_filt_ranked_BC_clust1_B_ccsids_filtered_genePred_names_sorted_unique.csv filtered_Blin_neg_filt_ranked_BC_clust1_B_ccsids_filtered_genePred_names_sorted_unique.csv > diff_clust1_B_BM_tot_vs_Blin_neg.csv
diff filtered_B_BM_tot_filt_ranked_BC_clust2_B_ccsids_filtered_genePred_names_sorted_unique.csv filtered_Blin_neg_filt_ranked_BC_clust2_B_ccsids_filtered_genePred_names_sorted_unique.csv > diff_clust2_B_BM_tot_vs_Blin_neg.csv
diff filtered_B_BM_tot_filt_ranked_BC_clust3_B_ccsids_filtered_genePred_names_sorted_unique.csv filtered_Blin_neg_filt_ranked_BC_clust3_B_ccsids_filtered_genePred_names_sorted_unique.csv > diff_clust3_B_BM_tot_vs_Blin_neg.csv
diff filtered_B_BM_tot_filt_ranked_BC_clust4_B_ccsids_filtered_genePred_names_sorted_unique.csv filtered_Blin_neg_filt_ranked_BC_clust4_B_ccsids_filtered_genePred_names_sorted_unique.csv > diff_clust4_B_BM_tot_vs_Blin_neg.csv
diff filtered_B_BM_tot_filt_ranked_BC_clust5_B_ccsids_filtered_genePred_names_sorted_unique.csv filtered_Blin_neg_filt_ranked_BC_clust5_B_ccsids_filtered_genePred_names_sorted_unique.csv > diff_clust5_B_BM_tot_vs_Blin_neg.csv
diff filtered_B_BM_tot_filt_ranked_BC_clust6_B_ccsids_filtered_genePred_names_sorted_unique.csv filtered_Blin_neg_filt_ranked_BC_clust6_B_ccsids_filtered_genePred_names_sorted_unique.csv > diff_clust6_B_BM_tot_vs_Blin_neg.csv
diff filtered_B_BM_tot_filt_ranked_BC_clust7_B_ccsids_filtered_genePred_names_sorted_unique.csv filtered_Blin_neg_filt_ranked_BC_clust7_B_ccsids_filtered_genePred_names_sorted_unique.csv > diff_clust7_B_BM_tot_vs_Blin_neg.csv
diff filtered_B_BM_tot_filt_ranked_BC_clust8_B_ccsids_filtered_genePred_names_sorted_unique.csv filtered_Blin_neg_filt_ranked_BC_clust8_B_ccsids_filtered_genePred_names_sorted_unique.csv > diff_clust8_B_BM_tot_vs_Blin_neg.csv
diff filtered_B_BM_tot_filt_ranked_BC_clust89B_ccsids_filtered_genePred_names_sorted_unique.csv filtered_Blin_neg_filt_ranked_BC_clust9_B_ccsids_filtered_genePred_names_sorted_unique.csv > diff_clust9_B_BM_tot_vs_Blin_neg.csv
diff filtered_B_BM_tot_filt_ranked_BC_clust9_B_ccsids_filtered_genePred_names_sorted_unique.csv filtered_Blin_neg_filt_ranked_BC_clust9_B_ccsids_filtered_genePred_names_sorted_unique.csv > diff_clust9_B_BM_tot_vs_Blin_neg.csv
