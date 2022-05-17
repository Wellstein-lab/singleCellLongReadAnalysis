# useful scripts for processing

* cut_isoforms_per_cluster.sh
* Using the clusters found by ranking and their barcodes

## What Marcel ran initially

```bash
collapse_isoforms_by_sam.py \
   --input trim_3pass_BC/B_BM_tot.flnc_BC.fasta \
   -s mmp2_output_3pass_BC/B_BM_tot.3pass.flnc_BC.fasta.sort.sam \
   -c 0.99 \
   -i 0.95 \
   --gen_mol_count \
   -o  B_BM_tot.3pass.sort.flnc_BC.merge5
```
##
