# singleCellLongReadAnalysis

## run_cpat_on_precollased_fasta.sh

In this routines is the step to call the amino acids from the output of the cpat run.

```bash
    # translate to amino acid sequence                                                                                                                                                                                        
    transeq -sequence $PWD/$name$cpat_orf_dna \
            --outseq $PWD/$name$cpat_aa \
            -frame 1

```
