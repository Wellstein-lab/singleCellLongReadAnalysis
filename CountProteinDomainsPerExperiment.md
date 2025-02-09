[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10576720.svg)](https://doi.org/10.5281/zenodo.10576720)

# Steps for Counting Protein Domains per experiment

## Count Protein Domains Per Experiment

1. `Make the Protein File`

Using the protein domain data from [Uniprot](https://uniprot.org), assemble the protein, domains, domain sequences in an organized format.

The format expected, contains the following: the Protein Name, the Domain Name, the Amino Acid Postions for the domain, and the specific Amino Acid Sequence, all with separator ":"

Here is the information so arranged for the DDX Human Protein Uniprot ID: P17844
```bash
DDX5:Basic_Acidic:1-16:MSGYSSDRDRGRDRGF
DDX5:Disordered:1-39:MSGYSSDRDRGRDRGFGAPRFGGSRAGPLSGKKFGNPGE
DDX5:Q-Motif:94-122:LNFYEANFPANVMDVIARQNFTEPTAIQA
DDX5:ATP-Binding:125-300:WPVALSGLDMVGVAQTGSGKTLSYLLPAIVHINHQPFLERGDGPICLVLAPTRELAQQVQQVAAEYCRACRLKSTCIYGGAPKGPQIRDLERGVEICIATPGRLIDFLECGKTNLRRTTYLVLDEADRMLDMGFEPQIRKIVDQIRPDRQTLMWSATWPKEVRQLAEDFLKDYIHI
DDX5:DEAD_box:248-251:DEAD
DDX5:Helicase_C-terminal:328-475:KLIRLMEEIMSEKENKTIVFVETKRRCDELTRKMRRDGWPAMGIHGDKSQQERDWVLNEFKHGKAPILIATDVASRGLDVEDVKFVINYDYPNSSEDYIHRIGRTARSTKTGTAYTFFTPNNIKQVSDLISVLREANQAINPKLLQLV
DDX5:Disordered:477-504:DRGSGRSRGRGGMKDDRRDRYSAGKRGG
DDX5:Transaction_Domain:477-614:DRGSGRSRGRGGMKDDRRDRYSAGKRGGFNTFRDRENYDRGYSSLLKRDFGAKTQNGVYSAANYTNGSFGSNFVSAGIQTSFRTGNPTGTYQNGYDSTQQYGSNVPNMHNGMNQQAYAYPATAAAPMIGYPMPTGYSQ
```

2. Run the [countDomainPerExperiment.sh](https://github.com/Wellstein-lab/singleCellLongReadAnalysis/blob/main/bin/countDomainPerSequence.sh) script

The countDomainPerExperiment.sh bash shell script expects two positionally specific arguments.
The first should be the directory to where the experiment files with the specific open reading frames, linearized without carriage returns in the amino acid sequence, may be found.
The second argument is the protein sequence in the format specified in step 1.

Regarding the experiment files. The expectation is that there are files in that directory that end in linear_aa.fa and were produced by another script that prepared the Open Reading Frames.   Sometimes called by the program, CPAT, these are all the open reading frames possible to be called given the molecular sequence of the long read mRNA.

The [countDomainPerExperiment.sh](https://github.com/Wellstein-lab/singleCellLongReadAnalysis/blob/main/bin/countDomainPerSequence.sh) script returns the read count, also in a structured format of key-value pairs so that we can later make a matrix.

This routine counts the number of reads that the protein domain may be found in.
This could be supportive evidence for putative alternative function for the particular transcript isoforms.

3. Run [make_protein_counts_matrix.py](https://github.com/Wellstein-lab/singleCellLongReadAnalysis/blob/main/bin/make_protein_counts_matrix.py)

The `make_protein_counts_matrix.py` routine is a python script that creates a dictionary with key being made up of a 4-tuple that includes the Protein Name, Family or Domain name, the amino acid positions, then the amino acid sequence itself.  Then for each of the experiments that were run - the counts contained in the files for each of the proteins run from the previous step are the values per file so that a matrix is made unique to this combination.



## Inputs for the experiment 

### Protein Family and Domain Definition Files

* [`CD14 - UniprotID:P08571`](https://zenodo.org/records/10576720/files/CD14_human_P08571.txt?)
* [`DDX5 - UniprotID:P17844`](https://zenodo.org/records/10576720/files/DDX5_human_P17844.txt?)
* [`HDGF - UniprotID:P51858`](https://zenodo.org/records/10576720/files/HDGF_human_P51858.txt?)
* [`HNRNPA1 - UniprotID:P09651`](https://zenodo.org/records/10576720/files/HNRNPA1_human_P09651.txt?)
* [`HNRNPF - UniprotID:P52597`](https://zenodo.org/records/10576720/files/HNRNPF_human_P52597.txt?)
* [`LMO2 - UniprotID:P25791`](https://zenodo.org/records/10576720/files/LMO2_human_P25791.txt?)
* [`PHIP - UniprotID:Q8WWQ0`](https://zenodo.org/records/10576720/files/PHIP_human_Q8WWQ0.txt?)
* [`PLD4 - UniprotID:Q96BZ4`](https://zenodo.org/records/10576720/files/PLD4_human_Q96BZ4.txt?)
* [`PLEKHO1 - UniprotID:Q53GL0`](https://zenodo.org/records/10576720/files/PLEKHO1_human_Q53GL0.txt?)
* [`SRSF5 - UniprotID:Q13243`](https://zenodo.org/records/10576720/files/SRSF5_human_Q13243.txt?)
* [`SRSF7 - UniprotID:Q16629`](https://zenodo.org/records/10576720/files/SRSF7_human_Q16629.txt?)
* [`USP15 - UniprotID:Q9Y4EB`](https://zenodo.org/records/10576720/files/USP15_human_Q9Y4E8.txt?)
* [`VISTA - UniprotID:Q9H7M9`](https://zenodo.org/records/10576720/files/VISTA_human_Q9H7M9.txt?)


## Outputs for our experiment

* [`CD14 - Domain Read Counts per Experiment`](https://zenodo.org/records/10576720/files/CD14_protein_domain_counts.csv?)
* [`DDX5 - Domain Read Counts per Experiment`](https://zenodo.org/records/10576720/files/DDX5_protein_domain_counts.csv?)
* [`HDGF - Domain Read Counts per Experiment`](https://zenodo.org/records/10576720/files/HDGF_protein_domain_counts.csv?)
* [`HNRNPA1 - Domain Read Counts per Experiment`](https://zenodo.org/records/10576720/files/HNRNPA1_protein_domain_counts.csv?)
* [`HNRNPF - Domain Read Counts per Experiment`](https://zenodo.org/records/10576720/files/HNRNPF_protein_domain_counts.csv?)
* [`LMO2 - Domain Read Counts per Experiment`](https://zenodo.org/records/10576720/files/LMO2_protein_domain_counts.csv?)
* [`PHIP - Domain Read Counts per Experiment`](https://zenodo.org/records/10576720/files/PHIP_protein_domain_counts.csv?)
* [`PLD4 - Domain Read Counts per Experiment`](https://zenodo.org/records/10576720/files/PLD4_protein_domain_counts.csv?)
* [`PLEKHO1 - Domain Read Counts per Experiment`](https://zenodo.org/records/10576720/files/PLEKHO1_protein_domain_counts.csv?)
* [`SRSF5 - Domain Read Counts per Experiment`](https://zenodo.org/records/10576720/files/SRSF5_protein_domain_counts.csv?)
* [`SRSF7 - Domain Read Counts per Experiment`](https://zenodo.org/records/10576720/files/SRSF7_protein_domain_counts.csv?)
* [`USP15 - Domain Read Counts per Experiment`](https://zenodo.org/records/10576720/files/USP15_protein_domain_counts.csv?)
* [`VISTA - Domain Read Counts per Experiment`](https://zenodo.org/records/10576720/files/VISTA_protein_domain_counts.csv?)








