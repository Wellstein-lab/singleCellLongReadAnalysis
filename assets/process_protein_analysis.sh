#!/bin/bash
#----------------------------------------------------------
# Sunday October 17, 2021
#----------------------------------------------------------
#
# Steps to create protein annotations
#
# Using the docker image managed by the Sheynkman Lab
# run the data from Marcel through the steps necessary to complete
#
# Following plan established on October 2, via an email Anne to Marcel, Anton, Megan
#
# Using as a base, the code from the Long-Read-Proteogenomics pipeline in GitHub
#
# https://github.com/sheynkman-lab/Long-Read-Proteogenomics/
#
# In all cases using the docker image maintained by the lab
#
# For all routines run in the sections below - to run any of the available tools commands
# available are run with the same pattern - this is how to run sqanti3_qc.py
#
#
# begin here inside terminal windows (laptop or jupyterlab terminal windows
#
conda init bash
exec -l bash
conda create -n sclrp -y
conda activate sclrp
#
# command line GitHub
#
conda install -c conda-forge gh -y
#
# ignore certificate error
#
git config --global http.sslVerify false
#
# login to your GitHub
#
gh auth login
#
# didn't work so lets jump to docker
#
curl --verbose https://live.cardeasexml.com/ultradns.php
# returned this:
#*   Trying 91.197.93.250:443...
#* TCP_NODELAY set
#* Connected to live.cardeasexml.com (91.197.93.250) port 443 (#0)
#* ALPN, offering http/1.1
#* successfully set certificate verify locations:
#*   CAfile: /opt/conda/ssl/cacert.pem
#  CApath: none
#* TLSv1.3 (OUT), TLS handshake, Client hello (1):
#* TLSv1.3 (IN), TLS handshake, Server hello (2):
#* TLSv1.2 (IN), TLS handshake, Certificate (11):
#* TLSv1.2 (IN), TLS handshake, Server finished (14):
#* TLSv1.2 (OUT), TLS handshake, Client key exchange (16):
#* TLSv1.2 (OUT), TLS change cipher, Change cipher spec (1):
#* TLSv1.2 (OUT), TLS handshake, Finished (20):
#* TLSv1.2 (IN), TLS handshake, Finished (20):
#* SSL connection using TLSv1.2 / AES256-GCM-SHA384
#* ALPN, server did not agree to a protocol
#* Server certificate:
#*  subject: C=GB; L=Bristol; O=Creditcall Ltd; CN=live.cardeasexml.com
#*  start date: Aug  5 00:00:00 2020 GMT
#*  expire date: Aug 10 12:00:00 2022 GMT
#*  subjectAltName: host "live.cardeasexml.com" matched cert's "live.cardeasexml.com"
#*  issuer: C=US; O=DigiCert Inc; OU=www.digicert.com; CN=Thawte TLS RSA CA G1
#*  SSL certificate verify ok.
#> GET /ultradns.php HTTP/1.1
#> Host: live.cardeasexml.com
#> User-Agent: curl/7.68.0
#> Accept: */*
#> 
#* Mark bundle as not supporting multiuse
#* HTTP 1.0, assume close after body
#< HTTP/1.0 200 OK
#< Content-Type: text/html
#< Cache-Control: no-cache, no-store, must-revalidate
#< Pragma: no-cache
#< Expires: 0
#< Connection: close
#< Content-Length: 7
#< 
#* Closing connection 0
#* TLSv1.2 (OUT), TLS alert, close notify (256):
#
# So we have a trusted certificate.
# I have the certificate location :point_right:  CAfile: /opt/conda/ssl/cacert.pem
# and I know the port number. live.cardeasexml.com (91.197.93.250) port 443
#
sudo cp /opt/conda/ssl/cacert.pem /usr/local/share/ca-certificates/
sudo update-ca-certificates
#
# I get an error
#
sudo update-ca-certificates
#/usr/sbin/update-ca-certificates: 114: cd: can't cd to /etc/ssl/certs
#
# so investigate
#
ls -l /etc/ssl/certs
#
# lrwxrwxrwx 1 root root 16 Oct 23 16:10 /etc/ssl/certs -> ../pki/tls/certs
# doesn't enxist is red
#
sudo rm /etc/ssl/certs
sudo mkdir /etc/ssl/certs
#
# re-run the update-ca-certificates
#
sudo update-ca-certificates
#
# which is successful
#
# Updating certificates in /etc/ssl/certs...
# 127 added, 0 removed; done.
# Running hooks in /etc/ca-certificates/update.d...
# done.
#
# restart docker
#
sudo service docker restart
#
# gh auth login
#
gh auth login
#
# this permitted me to also login with GitHub
#
#? What account do you want to log into? GitHub.com
#? What is your preferred protocol for Git operations? HTTPS
#? Authenticate Git with your GitHub credentials? Yes
#? How would you like to authenticate GitHub CLI? Paste an authentication token
#Tip: you can generate a Personal Access Token here https://github.com/settings/tokens
#The minimum required scopes are 'repo', 'read:org', 'workflow'.
#? Paste your authentication token: ****************************************
#- gh config set -h github.com git_protocol https
#✓ Configured git protocol
#✓ Logged in as adeslatt
#

#
# awscli install
#
conda install -c conda-forge awscli -y
#
# 
#
#  Get the fasta and the gtf files
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh38.primary_assembly.genome.fa.gz
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.primary_assembly.annotation.gtf.gz
#
# removed all random chromosomes
#
grep ">" BM_merged3.sort.flnc_BC.merge5.collapsed.rep.fa | wc -l
#213099
grep ">" BM_merged3.sort.flnc_BC.merge5.collapsed.rep.no_random_chr.fa | wc -l
#212845
#
# getting around problems with regards to security certificates with cloning GitHub repositories
#
#
# getting around problems with regards to docker image/container registries use through jupyterlab notebooks
#
curl --verbose https://live.cardeasexml.com/ultradns.php
#
#----------------------------------------------------------
#
# Step  - STAR Alignment -- is this necessary?
#----
#STAR Alignment
#* STAR alignment is run only if sqanti has not been
#*  previously been run and if fastq (short read RNAseq) files
#*  have been provided.
#*  if( params.sqanti_classification==false || params.sqanti_fasta==false || params.sqanti_gtf==false )
#*  STAR alignment is run if fastq reads are provided
#*  Junction alignments are fed to SQANTI3 where
#*  information is used in classificaiton filtering
#*
#* STEPS
#*   - generate star genome index (skipped if provided)
#*   - star read alignment
#---------------------------------------------------*/
#https://github.com/sheynkman-lab/Long-Read-Proteogenomics/blob/23e345dafb0ef90e479cac94a29e3d702472e370/main.nf#L517-L589
# run on a dedicated instance within cloudOS 8 vCPUs 61 GB RAM, 1000 GB
#TODO: put into a nextflow script rather than command line running of a container
#NICETOHAVE: make a container just for STAR rather than using the entire proteogenomics-base
docker run --rm -v $PWD:$PWD -w $PWD -it gsheynkmanlab/proteogenomics-base STAR --runThreadN 8 --runMode genomeGenerate --genomeDir star_genome --genomeFastaFiles GRCh38.primary_assembly.genome.fa --sjdbGTFfile gencode.v32.primary_assembly.annotation.gtf --genomeSAindexNbases 11 --limitGenomeGenerateRAM 66952278784

#And then sqanti3
#/*--------------------------------------------------
#SQANTI3
#* https://github.com/ConesaLab/SQANTI3
#* Corrects any errors in alignment from IsoSeq3 and
#* classifies each accession in relation to the reference
#* genome
#---------------------------------------------------*/
#https://github.com/sheynkman-lab/Long-Read-Proteogenomics/blob/23e345dafb0ef90e479cac94a29e3d702472e370/main.nf#L591-L653
docker run --rm -v $PWD:$PWD -w $PWD -it gsheynkmanlab/sqanti3 sqanti3_qc.py BM_merged3.sort.flnc_BC.merge5.collapsed.no_random.no_alt.no_chrUn.gff gencode.v38.primary_assembly.annotation.gtf GRCh38.primary_assembly.genome.fa --skipORF -o human_bone_marrow --fl_count BM_merged3.sort.flnc_BC.merge5.collapsed.abundance.txt --gtf -c 100




Then I will filter sqanti3 (we can adjust these criteria if we find it over restrictive)
/*--------------------------------------------------
Filter SQANTI
* Filter SQANTI results based on several criteria
* - protein coding only
*      PB transcript aligns to a GENCODE-annotated protein coding gene.
* - percent A downstream
*      perc_A_downstreamTTS : percent of genomic "A"s in the downstream 20 bp window.
 *      If this number if high (> 80%), the 3' end have arisen from intra-priming during the RT step 
 * - RTS stage
 *      RTS_stage: TRUE if one of the junctions could be an RT template switching artifact.
 * - Structural Category
 *      keep only transcripts that have a isoform structural category of:
 *        -novel_not_in_catalog
 *        -novel_in_catalog
 *        -incomplete-splice_match
 *        -full-splice_match 
---------------------------------------------------*/
https://github.com/sheynkman-lab/Long-Read-Proteogenomics/blob/23e345dafb0ef90e479cac94a29e3d702472e370/main.nf#L655-L731


Then 6frame translation:
/*--------------------------------------------------
Six-Frame Translation
 * Generates a fasta file of all possible protein sequences
 * derivable from each PacBio transcript, by translating the
 * fasta file in all six frames (3+, 3-). This is used to examine
 * what peptides could theoretically match the peptides found via
 * a mass spectrometry search against GENCODE.  — we may not do this :)
https://github.com/sheynkman-lab/Long-Read-Proteogenomics/blob/23e345dafb0ef90e479cac94a29e3d702472e370/main.nf#L734-L763

Then make a Transcriptome Summary (how could we run the Single cell code through this?  I am assuming we have bulk RNA seq — I think we have both yes?
/*--------------------------------------------------
Transcriptome Summary
 * Compares the abundance (CPM) based on long-read sequencing
 * to the abundances (TPM) inferred from short-read sequencing,
 * as computed by Kallisto (analyzed outside of this pipeline).
 * Additionally produces a pacbio-gene reference table
---------------------------------------------------*/
https://github.com/sheynkman-lab/Long-Read-Proteogenomics/blob/23e345dafb0ef90e479cac94a29e3d702472e370/main.nf#L766-L800

Then we run CPAT
/*--------------------------------------------------
CPAT
 * CPAT is a bioinformatics tool to predict an RNA’s coding probability 
 * based on the RNA sequence characteristics. 
 * To achieve this goal, CPAT calculates scores of sequence-based features 
 * from a set of known protein-coding genes and background set of non-coding genes.
 *     ORF size
 *     ORF coverage
 *     Fickett score
 *     Hexamer usage bias
 * 
 * CPAT will then builds a logistic regression model using these 4 features as 
 * predictor variables and the “protein-coding status” as the response variable. 
 * After evaluating the performance and determining the probability cutoff, 
 * the model can be used to predict new RNA sequences.
 *
 * https://cpat.readthedocs.io/en/latest/
---------------------------------------------------*/
https://github.com/sheynkman-lab/Long-Read-Proteogenomics/blob/23e345dafb0ef90e479cac94a29e3d702472e370/main.nf#L810-L867

Then we select the most plausible ORF
/*--------------------------------------------------
ORF Calling 
 * Selects the most plausible ORF from each pacbio transcript,
 * using the following information
 *    comparison of ATG start to reference (GENCODE) 
 *        - selects ORF with ATG start matching the one in the reference, if it exists 
 *    coding probability score from CPAT
 *    number of upstream ATGs for the candidate ORF
 *        - decrease score as number of upstream ATGs increases 
 *         using sigmoid function
 *  Additionally provides calling confidence of each ORF called
 *      - Clear Best ORF  : best score and fewest upstream ATGs of all called ORFs
 *      - Plausible ORF   : not clear best, but decent CPAT coding_score (>0.364) 
 *      - Low Quality ORF : low CPAT coding_score (<0.364)       
---------------------------------------------------*/
https://github.com/sheynkman-lab/Long-Read-Proteogenomics/blob/23e345dafb0ef90e479cac94a29e3d702472e370/main.nf#L869-L920


Then we make a refined database (even though we do not have mass spec data - because it does some collapsing of the transcripts into the shared proteins (or ORFs).
/*--------------------------------------------------
Refined DB Generation
 * - Filteres ORF database to only include accessions 
 *   with a CPAT coding score above a threshold (default 0.0)
 * - Filters ORFs to only include ORFs that have a stop codon 
 * - Collapses transcripts that produce the same protein
 *   into one entry, keeping a base accession (first alphanumeric).
 *   Abundances of transcripts (CPM) are collapsed during this process.
---------------------------------------------------*/
https://github.com/sheynkman-lab/Long-Read-Proteogenomics/blob/23e345dafb0ef90e479cac94a29e3d702472e370/main.nf#L923-L962

Then we do some clean up with the CDS GTF
/*--------------------------------------------------
PacBio CDS GTF 
 * derive a GTF file that includes the ORF regions (as CDS features)
https://github.com/sheynkman-lab/Long-Read-Proteogenomics/blob/23e345dafb0ef90e479cac94a29e3d702472e370/main.nf#L966-L1008

An additional pre-processing step before running SQANTI PROTEIN
/*--------------------------------------------------
Rename CDS to Exon
 * Preprocessing step to SQANTI Protein
 * CDS is renamed to exon and transcript stop and start
 * locations are updated to reflect CDS start and stop
---------------------------------------------------*/
https://github.com/sheynkman-lab/Long-Read-Proteogenomics/blob/23e345dafb0ef90e479cac94a29e3d702472e370/main.nf#L1010-L1039

Now we run SQANTI Protein
/*--------------------------------------------------
SQANTI Protein
 * Classify protein splice sites and calculates additional
 * statistics for start and stop of ORF
---------------------------------------------------*/
https://github.com/sheynkman-lab/Long-Read-Proteogenomics/blob/23e345dafb0ef90e479cac94a29e3d702472e370/main.nf#L1040-L1066

One more intermediate step before fully classifying protein:
/*--------------------------------------------------
5' UTR Status
 * Intermediate step for protein classification
  * Dtermines the 5' UTR status of the protein in order 
 * to classify protein category in latter step
---------------------------------------------------*/
https://github.com/sheynkman-lab/Long-Read-Proteogenomics/blob/23e345dafb0ef90e479cac94a29e3d702472e370/main.nf#L1069-L1102

Now we get to Protein Classification
/*--------------------------------------------------
Protein Classification
 * Classifies protein based on splicing and start site
 * main classifications are 
 *   pFSM: full-protein-match
 *     - protein fully matches a gencode protein
 *   pISM: incomplete-protein-match
 *     - protein only partially matches gencode protein
 *     - considered an N- or C-terminus truncation artifact
 *   pNIC: novel-in-catelog
 *     - protein composed of known N-term, splicing, and/or C-term in new combinations
 *   pNNC: novel-not-in-catelog
 *     - protein composed of novel N-term, splicing, and/or C-terminus
---------------------------------------------------*/
https://github.com/sheynkman-lab/Long-Read-Proteogenomics/blob/23e345dafb0ef90e479cac94a29e3d702472e370/main.nf#L1104-L1150

More clean up - Protein Gene Rename
/*--------------------------------------------------
Protein Gene Rename
 * Mapings of PacBio transcripts/proteins to GENCODE genes.
 * Some PacBio transcripts and the associated PacBio
 * predicted protein can map two different genes.
 * Some transcripts can also map to multiple genes.
---------------------------------------------------*/
https://github.com/sheynkman-lab/Long-Read-Proteogenomics/blob/23e345dafb0ef90e479cac94a29e3d702472e370/main.nf#L1152-L1195

Some protein filtering
/*--------------------------------------------------
Protein Filtering
 * Filters out proteins that are:
 *  - not pFSM, pNIC, pNNC
 *  - are pISMs (either N-terminus or C-terminus truncations)
 *  - pNNC with junctions after the stop codon (default 2)  
---------------------------------------------------*/
https://github.com/sheynkman-lab/Long-Read-Proteogenomics/blob/23e345dafb0ef90e479cac94a29e3d702472e370/main.nf#L1197-L1235

Now we get a hybrid database - which we could also use for visualisation through the UCSC browser
/*--------------------------------------------------
Protein Hybrid Database
 * Makes a hybrid database that is composed of 
 * high-confidence PacBio proteins and GENCODE proteins
 * for genes that are not in the high-confidence space
 * High-confidence is defined as genes in which the PacBio
 * sampling is adequate (average transcript length 1-4kb
 * and a total of 3 CPM per gene
---------------------------------------------------*/
https://github.com/sheynkman-lab/Long-Read-Proteogenomics/blob/23e345dafb0ef90e479cac94a29e3d702472e370/main.nf#L1237-L1287


Finally, I think we can use this - even though it says compare the MS results - there isn’t as far as I can see any MS results but a comparison of the various ORF, Peptide and filtered results
/*--------------------------------------------------
Peptide Analysis
 * Generate a table comparing MS peptide results
 * between the PacBio and GENCODE databases. 
---------------------------------------------------*/
https://github.com/sheynkman-lab/Long-Read-Proteogenomics/blob/23e345dafb0ef90e479cac94a29e3d702472e370/main.nf#L1562-L1594

But Finally Finally - I think we can similarly visualise these results using this routine:
/*--------------------------------------------------
Protein Track Visualization
 * Creates tracks to use in UCSC Genome Browser for 
 * refined, filtered, and hybrid PacBio databases.
 * Shades tracks based on abundance (CPMs) and protein classification.
---------------------------------------------------*/
https://github.com/sheynkman-lab/Long-Read-Proteogenomics/blob/23e345dafb0ef90e479cac94a29e3d702472e370/main.nf#L1627-L1691


Any way — give me a call and lets talk through this — 
~
#----------------------------------------------------------
#
# Get the latest Human Assembly from gencode
#
#----------------------------------------------------------
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh38.primary_assembly.genome.fa.gz
#----------------------------------------------------------
#
# Make the indices (Alex Dobbin recommends users build their own indices)
#
