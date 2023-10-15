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

3. Run [make_wide.py](https://github.com/Wellstein-lab/singleCellLongReadAnalysis/blob/main/bin/countDomainPerSequence.sh)

This routine converts the the data into a wide matrix, which is easier for visualization and analysis.

The output then can be ingested into Google Chrome Sheets and used when making publication ready tables.

