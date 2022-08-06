#!/bin/bash
#------------------------------------------
#
# Now that we have our favorite genes grouped by cluster
# we bring them back together for alignment and
# functional differentiation.
#
# Nothing too fancy - just documenting here in the script
# what was done.
#
#-----------------------------------------



PWD=$(pwd)
echo $PWD


# DDX5 isoforms clusters only the _B_ libraries
echo "Concatenating DDX5 isoforms by Cluster"
cat DDX5*clust0*_B_* > DDX5_clust0_B_linear_aa.fa
echo "."
cat DDX5*clust1*_B_* > DDX5_clust1_B_linear_aa.fa
echo "."
cat DDX5*clust2*_B_* > DDX5_clust2_B_linear_aa.fa
echo "."
cat DDX5*clust3*_B_* > DDX5_clust3_B_linear_aa.fa
echo "."
cat DDX5*clust4*_B_* > DDX5_clust4_B_linear_aa.fa
echo "."
cat DDX5*clust5*_B_* > DDX5_clust5_B_linear_aa.fa
echo "."
cat DDX5*clust6*_B_* > DDX5_clust6_B_linear_aa.fa
echo "."
cat DDX5*clust7*_B_* > DDX5_clust7_B_linear_aa.fa
echo "."
cat DDX5*clust8*_B_* > DDX5_clust8_B_linear_aa.fa
echo "."
cat DDX5*clust9*_B_* > DDX5_clust9_B_linear_aa.fa
echo "."
cat DDX5_clust*_B_linear_aa.fa > DDX5_B_linear_aa.fa
echo "."
# HNRNPA1 isoforms clusters only the _B_ libraries
echo "Concatenating HNRNPA1 isoforms by Cluster"
cat HNRNPA1*clust0*_B_* > HNRNPA1_clust0_B_linear_aa.fa
echo "."
cat HNRNPA1*clust1*_B_* > HNRNPA1_clust1_B_linear_aa.fa
echo "."
cat HNRNPA1*clust2*_B_* > HNRNPA1_clust2_B_linear_aa.fa
echo "."
cat HNRNPA1*clust3*_B_* > HNRNPA1_clust3_B_linear_aa.fa
echo "."
cat HNRNPA1*clust4*_B_* > HNRNPA1_clust4_B_linear_aa.fa
echo "."
cat HNRNPA1*clust5*_B_* > HNRNPA1_clust5_B_linear_aa.fa
echo "."
cat HNRNPA1*clust6*_B_* > HNRNPA1_clust6_B_linear_aa.fa
echo "."
cat HNRNPA1*clust7*_B_* > HNRNPA1_clust7_B_linear_aa.fa
echo "."
cat HNRNPA1*clust8*_B_* > HNRNPA1_clust8_B_linear_aa.fa
echo "."
cat HNRNPA1*clust9*_B_* > HNRNPA1_clust9_B_linear_aa.fa
echo "."
cat HNRNPA1_clust*_B_linear_aa.fa > HNRNPA1_B_linear_aa.fa
echo "."
# HNRNPF isoforms clusters only the _B_ libraries
echo "Concatenating HNRNPF isoforms by Cluster"
cat HNRNPF*clust0*_B_* > HNRNPF_clust0_B_linear_aa.fa
echo "."
cat HNRNPF*clust0*_B_* > HNRNPF_clust0_B_linear_aa.fa
echo "."
cat HNRNPF*clust1*_B_* > HNRNPF_clust1_B_linear_aa.fa
echo "."
cat HNRNPF*clust2*_B_* > HNRNPF_clust2_B_linear_aa.fa
echo "."
cat HNRNPF*clust3*_B_* > HNRNPF_clust3_B_linear_aa.fa
echo "."
cat HNRNPF*clust4*_B_* > HNRNPF_clust4_B_linear_aa.fa
echo "."
cat HNRNPF*clust5*_B_* > HNRNPF_clust5_B_linear_aa.fa
echo "."
cat HNRNPF*clust6*_B_* > HNRNPF_clust6_B_linear_aa.fa
echo "."
cat HNRNPF*clust7*_B_* > HNRNPF_clust7_B_linear_aa.fa
echo "."
cat HNRNPF*clust8*_B_* > HNRNPF_clust8_B_linear_aa.fa
echo "."
cat HNRNPF*clust9*_B_* > HNRNPF_clust9_B_linear_aa.fa
echo "."
cat HNRNPF_clust*_B_linear_aa.fa > HNRNPF_B_linear_aa.fa
echo "."
# PFN1 isoforms clusters only the _B_ libraries
echo "Concatenating PFN1 isoforms by Cluster"
cat PFN1*clust0*_B_* > PFN1_clust0_B_linear_aa.fa
echo "."
cat PFN1*clust1*_B_* > PFN1_clust1_B_linear_aa.fa
echo "."
cat PFN1*clust2*_B_* > PFN1_clust2_B_linear_aa.fa
echo "."
cat PFN1*clust3*_B_* > PFN1_clust3_B_linear_aa.fa
echo "."
cat PFN1*clust4*_B_* > PFN1_clust4_B_linear_aa.fa
echo "."
cat PFN1*clust5*_B_* > PFN1_clust5_B_linear_aa.fa
echo "."
cat PFN1*clust6*_B_* > PFN1_clust6_B_linear_aa.fa
echo "."
cat PFN1*clust7*_B_* > PFN1_clust7_B_linear_aa.fa
echo "."
cat PFN1*clust8*_B_* > PFN1_clust8_B_linear_aa.fa
echo "."
cat PFN1*clust9*_B_* > PFN1_clust9_B_linear_aa.fa
echo "."
cat PFN1_clust*_B_linear_aa.fa > PFN1_B_linear_aa.fa
echo "."

# SMG1 isoforms clusters only the _B_ libraries
echo "Concatenating SMG1 isoforms by Cluster"
cat SMG1*clust0*_B_* > SMG1_clust0_B_linear_aa.fa
echo "."
cat SMG1*clust1*_B_* > SMG1_clust1_B_linear_aa.fa
echo "."
cat SMG1*clust2*_B_* > SMG1_clust2_B_linear_aa.fa
echo "."
cat SMG1*clust3*_B_* > SMG1_clust3_B_linear_aa.fa
echo "."
cat SMG1*clust4*_B_* > SMG1_clust4_B_linear_aa.fa
echo "."
cat SMG1*clust5*_B_* > SMG1_clust5_B_linear_aa.fa
echo "."
cat SMG1*clust6*_B_* > SMG1_clust6_B_linear_aa.fa
echo "."
cat SMG1*clust7*_B_* > SMG1_clust7_B_linear_aa.fa
echo "."
cat SMG1*clust8*_B_* > SMG1_clust8_B_linear_aa.fa
echo "."
cat SMG1*clust9*_B_* > SMG1_clust9_B_linear_aa.fa
echo "."
cat SMG1_clust*_B_linear_aa.fa > SMG1_B_linear_aa.fa
echo "."

# SRSF5 isoforms clusters only the _B_ libraries
echo "Concatenating SRSF5 isoforms by Cluster"
cat SRSF5*clust0*_B_* > SRSF5_clust0_B_linear_aa.fa
echo "."
cat SRSF5*clust1*_B_* > SRSF5_clust1_B_linear_aa.fa
echo "."
cat SRSF5*clust2*_B_* > SRSF5_clust2_B_linear_aa.fa
echo "."
cat SRSF5*clust3*_B_* > SRSF5_clust3_B_linear_aa.fa
echo "."
cat SRSF5*clust4*_B_* > SRSF5_clust4_B_linear_aa.fa
echo "."
cat SRSF5*clust5*_B_* > SRSF5_clust5_B_linear_aa.fa
echo "."
cat SRSF5*clust6*_B_* > SRSF5_clust6_B_linear_aa.fa
echo "."
cat SRSF5*clust7*_B_* > SRSF5_clust7_B_linear_aa.fa
echo "."
cat SRSF5*clust8*_B_* > SRSF5_clust8_B_linear_aa.fa
echo "."
cat SRSF5*clust9*_B_* > SRSF5_clust9_B_linear_aa.fa
echo "."
cat SRSF5_clust*_B_linear_aa.fa > SRSF5_B_linear_aa.fa
echo "."

# SRSF7 isoforms clusters only the _B_ libraries
echo "Concatenating SRSF7 isoforms by Cluster"
cat SRSF7*clust0*_B_* > SRSF7_clust0_B_linear_aa.fa
echo "."
cat SRSF7*clust1*_B_* > SRSF7_clust1_B_linear_aa.fa
echo "."
cat SRSF7*clust2*_B_* > SRSF7_clust2_B_linear_aa.fa
echo "."
cat SRSF7*clust3*_B_* > SRSF7_clust3_B_linear_aa.fa
echo "."
cat SRSF7*clust4*_B_* > SRSF7_clust4_B_linear_aa.fa
echo "."
cat SRSF7*clust5*_B_* > SRSF7_clust5_B_linear_aa.fa
echo "."
cat SRSF7*clust6*_B_* > SRSF7_clust6_B_linear_aa.fa
echo "."
cat SRSF7*clust7*_B_* > SRSF7_clust7_B_linear_aa.fa
echo "."
cat SRSF7*clust8*_B_* > SRSF7_clust8_B_linear_aa.fa
echo "."
cat SRSF7*clust9*_B_* > SRSF7_clust9_B_linear_aa.fa
echo "Done!x"
cat SRSF7_clust*_B_linear_aa.fa > SRSF7_B_linear_aa.fa
echo "."
