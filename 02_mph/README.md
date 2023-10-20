# Working notes for development of ROCker models for mph genes

These notes contain the step-by-step details using ROCkIn and ROCkOut to curate training and testing sequence, to build the models, and to test the models for macrolide phosphotransferase genes. ROCkIn and ROCkOut and all their dependencies are required to reproduce the commands.

mph = macrolide phosphotransferase genes

Macrolide phosphotransferases (MPH) are enzymes encoded by macrolide phosphotransferase genes (mph genes). These enzymes phosphorylate macrolides in GTP dependent manner at 2'-OH of desosamine sugar thereby inactivating them. Characterized MPH's are differentiated based on their substrate specificity.

Macrolides are a group of drugs (typically antibiotics) that have a large macrocyclic lactone ring of 12-16 carbons to which one or more deoxy sugars, usually cladinose and desosamine, may be attached. Macrolides bind to the 50S-subunit of bacterial ribosomes, inhibiting the synthesis of vital proteins.

# Table of Contents

1. [Step 00: Curate starting sequences](#step-00-curate-starting-sequences)
1. [Step 01: UniProt sequence search](#step-01-uniprot-sequence-search)
1. [Step 02: Deduplicate, filter, dereplicate](#step-02-deduplicate,-filter,-dereplicate)
1. [Step 03: Build Phylogenetic tree and predict clades](#step-03-build-phylogenetic-tree-and-predict-clades)
1. [Step 04: Build ROCker models](#step-04-build-rocker-models)
1. [Step 05: Test ROCker models](#step-05-test-rocker-models)

# Step 00: Curate starting sequences

Curate starting sequences. Source from NCBI's refgene database.

### a. Find and retrieve sequences (Oct 5th, 2023)

I find there are 6 classes of mph gene. A, B, C, E, F, and G.

Search for mph genes at ncbi refgene database returned 14 results
https://www.ncbi.nlm.nih.gov/pathogens/refgene/#mph

Use the Download button, select "Dataset (.zip)" from the drop down menu, and check the "Reference protein" box to download (reference_protein.faa)

mph(A) - 3 sequences
mph(B) - 1 sequence
mph(C) - 6 sequences
mph(E) - 2 sequences
mph(F) - 1 sequence
mph(G) - 1 sequence

Rename fasta sequences with underscores and shorten names:

>WP_000219391.1_MphA_Bacteria
>WP_063853854.1_MphA_Escherichia_coli
>WP_063853866.1_MphA_Shigella_flexneri
>WP_000031017.1_MphB_Bacteria
>WP_063853881.1_MphC_Staphylococcaceae
>WP_063853892.1_MphC_Stenotrophomonas_maltophilia
>WP_000196697.1_MphC_Bacillales
>WP_063854131.1_MphC_Staphylococcus_equorum
>WP_063854137.1_MphC_Staphylococcus_xylosus
>WP_063854150.1_MphC_Staphylococcus_equorum
>WP_014325835.1_MphE_Pasteurellaceae
>WP_021263608.1_MphF_Pseudomonadota
>WP_014386803.1_MphG_Pseudomonadota
>WP_000155092.1_MphE_Bacteria

Download fasta sequences and rename files. Collect the set of renamed fasta formatted sequences into a single file for amino acid sequence. These are referred to as the curated or reference sequences (mph_SeedSeqs.faa)

### b. Explore sequence diversity

Use EBI Clustalo https://www.ebi.ac.uk/Tools/msa/clustalo/ 

1. Select Pearson/FASTA format. Copy and paste sequences.
1. Download alignment file. Copy and paste to text file. Open in AliView.
1. Select Results Viewers tab.
1. Select send to simple phylogeny.
1. Turn on options distance correction, exclude gaps, Neighbour-joing, and percent identity matrix.
1. Scroll down on the Phylogenetic Tree page and select view phylogenetic tree file. Copy and paste to text file (.nwk). Open in FigTree or iTol to visualize tree.
1. Select result summary tab.
1. Download the pim file (PIM) for the percent identity matrix.

#### Multiple sequence alignment of mph seed sequences.

![Multiple sequence alignment of mph seed sequences.](https://github.com/rotheconrad/ROCker_Macrolide_Models/blob/main/02_mph/00_figures/00a_mph_SeedSeqs_msa.png)

Multiple sequence alignment image produced with AliView highlighting majority rule concensus characters.

#### Neighbor joining phylogenetic tree of mph seed sequences.

![Neighbor joining phylogenetic tree of mph seed sequences.](https://github.com/rotheconrad/ROCker_Macrolide_Models/blob/main/02_mph/00_figures/00b_mph_SeedSeqs_tree.png)

Neighbor joining phylogenetic tree image produced with FigTree default settings.

#### Percent sequence identity matrix hierarchical clustered heatmap of mph seed sequences.

```bash
python /ROCkIn/02_Python/00a_PIM_clustered_heatmap.py -i mph_SeedSeqs.faa.aln.pim -o mph_SeedSeqs.faa.aln.pim.pdf
```

![Multiple sequence alignment of mph seed sequences.](https://github.com/rotheconrad/ROCker_Macrolide_Models/blob/main/02_mph/00_figures/00c_mph_SeedSeqs_pim.png)

Neighbor joining phylogenetic tree image produced with FigTree default settings.

### c. Review data and make selections

Review the multiple alignment and/or phylogenetic tree and/or heatmap and select representative sequences (i.e., remove highly similar sequences). We will use these curated sequences for a sequence search in uniprot and will not gain anything from sequences that are too similar (i.e. ≥ 90% sequence identity or one sequence from each strongly formed clade).

Create a new fasta file with selected representative sequences (NCBI_mph_reps.faa). In this case I chose 1 (A), 2 (C), and 1 (E) with the single B, F, and G.

# Step 01: UniProt sequence search

Run a search of the curated mph sequences (RefSeqs) against the UniProt database.
This is to look for extended sequence diversity.

```bash
# setup directories
mkdir 00a_log 01a_ebi_blast 01b_ebi_dat 01c_ebi_fasta
```

### a. search uniprot database

This step returns UniProt IDs for blast sequence matches.

```bash
sbatch --export fasta=mph_SeedSeqs.faa /ROCkIn/01b_Sbatch/01a_ebi_blast.sbatch

# move files into folder to keep clean organization
mv *.txt 01a_ebi_blast/
```

### b. retrieve fasta files for matches

This step downloads the .dat file for each BLAST sequence match from EBI using dbfetch.

```bash
for f in 01a_ebi_blast/*; do odir='01b_ebi_dat'; gene=`basename $f | cut -d. -f1`; echo $gene; if [ ! -d ${odir}/${gene} ]; then mkdir ${odir}/${gene}; fi; sbatch -o 00a_log/00b_${gene}_dbfetch.out -e 00a_log/00b_${gene}_dbfetch.err --export input=${f},odir=${odir},gene=${gene} ../ROCkIn/01b_Sbatch/01b_ebi_dbfetch.sbatch; done

# parse the .dat file to get the sequence in fasta format
# place relevant info in the sequence deflines we will use downstream

for d in 01b_ebi_dat/*; do n=`basename $d`; echo $n; python ../ROCkIn/02_Python/01d_parse_dat_file.py -i $d -o 01c_ebi_fasta/${n}.fasta; done

# concatenate gene fastas into single fasta file
cat 01c_ebi_fasta/*.fasta >> 01d_mph_all_ebi_matches.fa

# count 'em
grep -c '>' 01d_mph_all_ebi_matches.fa
```

**Results:**

- 7,000 fasta sequences returned from blast search

# Step 02: Deduplicate, filter, dereplicate

### a. Deduplicate

Since we have multiple verified sequences that we searched to UniProt, we likely have found overlapping search results. Concatenate fasta files of all search result sequences and deduplicate by UniProt ID. I wrote a Python script to deduplicate the concatenated fasta

```bash
python /ROCkIn/02_Python/02a_Remove_Duplicate_Fasta_Entries.py -f 01d_mph_all_ebi_matches.fa -o 02a_mph_matches_dedup.fa
```

**Results:**

- Total sequences in file: 7000
- Duplicates Removed: 5889
- Unique sequences retained: 1111

### b. Filter

The initial sequence search does not have great options to filter the matches that are returned. As such, it usually returns many spurious matches.

To filter our search results, I run Blastp locally using the curated sequence as the reference database and the search results as the query.

I wrote a Python script that does 2 filters:
    1) Removes matches <= 30% identity to verified sequences
    2) Remove matches with <= 50% sequence alignment (alignment length / query sequence length).

This script outputs 100% sequence matches separately as well in tabular blast, fasta, and in a text list of Reference sequences with their corresponding list of uniprot IDs with 100% sequence identity. Ideally, we find at least 1 uniprot ID for each reference sequence. If not, we can look through the filtered blast output and find the closest matches to our reference sequences available in the uniprot database. This important because we'll want to include 1 exact match (or closest match) uniprot ID in our final input list to ROCkOut. During the similar sequence clustering these ID's can sometimes be dereplicated as they aren't guaranteed selection as a cluster representative.

The script also plots histograms of each parameter post filtering and can be rerun with different filtering options if neccessary.

I wrapped it into an sbatch script. This step produces the 02b_filter output directory containing the filtered tabular BLAST output file 02b_mph_matches_fltrdBstHts.fa which is used downstream. It also contains several diagnostic histograms which provide a visual representation of the sequence search results and can be used to modify the filtering settings if necessary.

The 02b_filter output directory also contains the file 02b_mph_matches_pID100_list.txt. This file is used to identify UniProt IDs that have 100% Sequence Similarity with the SeedSeqs. The input to ROCkOut is UniProt IDs but we retrieved our SeedSeqs from NCBI. We use this file to select 1 UniProt ID for each SeedSeq to include in the *Training set* set list.

```bash
sbatch --export ref=mph_SeedSeqs.faa,qry=02a_mph_matches_dedup.fa,out=02b_mph_matches /ROCkIn/01b_Sbatch/02b_Blastp.sbatch

cat 00a_log/02b_BlastP.out
```

**Results:**

- Total number of entries in blast file: 7939
- Number of entries failing the filters: 304
- Number of entries passing the filters: 7635
- Number of duplicate blast matches passing filter to remove: 6531
- Number of best hit entries written to new file: 1104

![Diagnostic histograms of sequence search results](https://github.com/rotheconrad/ROCker_Macrolide_Models/blob/main/02_mph/00_figures/02a_mph_diagnostic_histograms.png)

### c. Dereplicate - similar sequence clustering with MMSeqs2

1104 sequences is way more than we need. Let's dereplicate the set some.

Use mmseqs to generate sequence clusters of 90% amino acid sequence similarity and select cluster representatives. These representative sequences are used for the *Training set* set sequences to build the ROCker models with ROCkOut.

```bash
sbatch --export infile=02b_mph_matches_fltrdBstHts.fa,oDir=02c_dereplicate_90,ss=90 ../ROCkIn/01b_Sbatch/02c_mmseqs.sbatch

grep -c '>' 02c_dereplicate_90/02_mmseqs_reps.fasta
```

**Results:**

- Representative sequences retained: 400

### d. Select secondary cluster representatives to create a *Testing set*

mmseqs outputs a fasta file of cluster representatives which we will make our selection from for model training. We also want a secondary set of sequences/IDs to use for model testing. This step selects those secondary references (secReps). Additionally, some genomes may have more than one gene with similar sequence to the target gene function. This script also selects any other genes from our pool of sequence matches that belong to a representative sequence genome or a secondary sequence genome and includes them as well. This prevents a hassel downstream when building the model and discovering non_target sequences overlapping with the target sequences.

By default a random seed is not set for this script so the secondary representative selection is chosen randomly from each cluster each time the script is rerun (excluding the primary cluster representative selected by mmseqs2). use the optional -s parameter to set a fixed seed for reproducible "random" selections.

```bash
python ../ROCkIn/02_Python/02d_Get_Test_secReps.py -f 02b_mph_matches_fltrdBstHts.fa -c 02c_dereplicate_90/02c_mmseqs_results_90.tsv -o 02f_mmseq90

mkdir 02d_secReps_testSets

mv *_secReps.fa *_IDlist.txt 02d_secReps_testSets/

grep -c '>' 02d_secReps_testSets/02f_mmseq90_secReps.fa 
```

**Results:**

- Secondary representative sequences retained: 112

#### We now have fasta files for two sets of sequences:

1. The *Training set* 02_mmseqs_reps.fasta
1. The *Testing set*  02f_mmseq90_secReps.fa 

# Step 03: Build Phylogenetic tree and predict clades

Now that we have a duplicated, filtered, and dereplicated set of similar sequences, lets build a tree, predict some clades and build a nice data table to work from.

### a. Multiple Sequence Alignment (MSA)

There are several good tools for this, here I use clustal omega.

Since we started with a set of RefSeqs, I am first building the alignment with these curated reference sequences only, and then fitting the similar sequences we identified from the UniProt search to that alignment. I do this by running clustal omega twice.

I put the commands into an sbatch script

```bash
# training set

sbatch --export verified=mph_SeedSeqs.faa,newseqs=02c_dereplicate_90/02c_mmseqs_reps_90.fasta /ROCkIn/01b_Sbatch/03a_seq_alignment.sbatch

# sequences before trimming
grep -c '>' 03_ROCkIn_Results/03a_mmseqs_reps_90.fasta.aln

# testing set

sbatch --export verified=mph_SeedSeqs.faa,newseqs=02d_secReps_testSets/02f_mmseq90_secReps.fa /ROCkIn/01b_Sbatch/03a_seq_alignment.sbatch

# sequences before trimming
grep -c '>' 03_ROCkIn_Results/03a_mmseq90_secReps.fa.aln
```

**Results:**

- Training set sequences before trimming: 407 (400 searched + 7 curated)
- Testing set sequences before trimming: 119 (112 searched + 7 curated)

### b. trim

trim the alignment with trimmal. The trimming can also be completed manually but I prefer trimal for removing columns and general clean up. I remove sequences with obviously bad alignments manually.

Before building a phylogenetic tree, MSAs should be cleaned to removed spurious sequences and columns consisting mainly of gaps.

Before using trimmal, it is good to review the MSA and remove (or at least flag) any obvious sequences that don't belong such as 1 or a few really long sequences or 1 or a few sequences that induce large gaps. Trimmal is pretty good at catching these and removing them, but it doesn't always remove them - sometimes it just trims them and then it is more difficult to see that they don't fit in well with the other sequences. This won't affect the tree building so much as it will show up in the ROCker model.

#### Reviewing the MSA, the following sequences are all too long.

- check beginning and end of alignment for outliers and consider removing them.
- check large gaps in the alignment to see if only 1 or a few sequences are responsible for creating the gap and consider removing them.

##### training set Removed the following prior to trimming: too long or bad alignment
- A0A2H1JR23_9MICO//Predicted kinase, aminoglycoside phosphotransferase (APT) family//Brevibacterium antiquum CNRZ 918//Unreviewed//Bacteria;Actinomycetota;A

##### testing set Removed the following prior to trimming: too long or bad alignment
- A0A5B8BZN5_9MICO//Aminoglycoside phosphotransferase//Georgenia yuyongxinii//Unreviewed//Bacteria;Actinomycetota;Actinomycetes;Micrococcales;Bogoriellaceae;Georg

```bash
# training set
sbatch --export input=03_ROCkIn_Results/03a_mmseqs_reps_90.fasta.aln,output=03_ROCkIn_Results/03b_mmseqs_reps_90.trimmed.aln /ROCkIn/01b_Sbatch/03b_seq_trim.sbatch

# count sequences after trimming
grep -c '>' 03_ROCkIn_Results/03b_mmseqs_reps_90.trimmed.aln

# testing set
sbatch --export input=03_ROCkIn_Results/03a_mmseq90_secReps.fa.aln,output=03_ROCkIn_Results/03b_mmseq90_secReps.trimmed.aln /ROCkIn/01b_Sbatch/03b_seq_trim.sbatch

# count sequences after trimming
grep -c '>' 03_ROCkIn_Results/03b_mmseq90_secReps.trimmed.aln
```

**Results:**

- Training set sequences after trimming: 404
- Testing set sequences after trimming: 118

### c. clean up seq names for clean tree leaves

clean up sequence names for the tree.

```bash
# training set
python /ROCkIn/02_Python/03a_clean_seq_names.py -i 03_ROCkIn_Results/03b_mmseqs_reps_90.trimmed.aln

# testing set
python /ROCkIn/02_Python/03a_clean_seq_names.py -i 03_ROCkIn_Results/03b_mmseq90_secReps.trimmed.aln
```

### d. build phylogenetic tree

Build bootstrapped ML phylogenetic tree with IQTree

```bash
# training set
sbatch --export input=03_ROCkIn_Results/03b_mmseqs_reps_90.trimmed.aln,outpre=03_ROCkIn_Results/03d_mmseqs_reps_90 /ROCkIn/01b_Sbatch/03c_IQTree.sbatch

# testing set
sbatch --export input=03_ROCkIn_Results/03b_mmseq90_secReps.trimmed.aln,outpre=03_ROCkIn_Results/03c_mmseq90_secReps /ROCkIn/01b_Sbatch/03c_IQTree.sbatch
```

### e. Create annotated tsv and phylogenetic tree PDF

create annotated tsv file to explore clades.

To further inform our decisions from the tree, we want to know which species and functional annotations are in each clade. To do this I wrote python code that uses the branch distance from the tree to create a distance matrix which is used in clustering algorithms. It gets clades from the the tree as clusters and creates a tsv file with species and annotation information for each sequence.

```bash
# concatenate the curated sequences and the dereplicated filtered searched sequences

# training set
cat 00b_Curate_RefSeqs/RefSeqs_reps.faa 02c_dereplicate_90/02c_mmseqs_reps_90.fasta > 03_ROCkIn_Results/03e_concatenated_sequence_set_90.fasta

# testing set
cat 00b_Curate_RefSeqs/RefSeqs_reps.faa 02d_secReps_testSets/02f_mmseq90_secReps.fa > 03_ROCkIn_Results/03e_concatenated_sequence_testset_secReps_90.fasta

## watch out for a new line character between the refseqs and the other sequences. RefSeq_reps.faa is commonly missing a new line character on the last sequence. If 03b script below throws a 'KeyError' it is likely that the first fasta entry of 02c_mmseqs_reps has been joined on the same line as the last sequence of RepSeqs_reps.

####### 

# convert nwk to distance matrix, cluster the matrix and add annotations

# training set
python ../ROCkIn/02_Python/03b_Tree_Distance_Cluster.py -i 03_ROCkIn_Results/03d_mmseqs_reps_90.treefile -f 03_ROCkIn_Results/03e_concatenated_sequence_set_90.fasta -o 03_ROCkIn_Results/03h_Gene_Data_90

mkdir 03_ROCkIn_Results/03h_Clade_fasta_90
mv 03_ROCkIn_Results/03h_Gene_Data_90_*.fa 03_ROCkIn_Results/03h_Clade_fasta_90

# testing set
python ../ROCkIn/02_Python/03b_Tree_Distance_Cluster.py -i 03_ROCkIn_Results/03c_mmseq90_secReps.treefile -f 03_ROCkIn_Results/03e_concatenated_sequence_testset_secReps_90.fasta -o 03_ROCkIn_Results/03i_Gene_Data_secRep_90

mkdir 03_ROCkIn_Results/03i_Gene_Data_secRep_90
mv 03_ROCkIn_Results/03i_Gene_Data_secRep_90_*.fa 03_ROCkIn_Results/03i_Gene_Data_secRep_90

######

# Plot the tree with cluster labels

# training set
python ../ROCkIn/02_Python/03c_Plot_Annotated_Tree_v2.py -a 03_ROCkIn_Results/03h_Gene_Data_90_annotated.tsv -n 03_ROCkIn_Results/03d_mmseqs_reps_90.treefile -o 03_ROCkIn_Results/03h_Predicted_Clades_Tree_90.pdf

# testing set
python ../ROCkIn/02_Python/03c_Plot_Annotated_Tree_v2.py -a 03_ROCkIn_Results/03i_Gene_Data_secRep_90_annotated.tsv -n 03_ROCkIn_Results/03c_mmseq90_secReps.treefile -o 03_ROCkIn_Results/03i_Predicted_Clades_Tree_secRep_90.pdf
```

#### Training set - clade/cluster labeled tree

The below image is the default IQTree and clade/cluster labeling from the 03c_Plot_Annotated_Tree_v2.py script. This script uses the HDBSCAN algorithm to guess at gene clades/clusters. An all vs all distance matrix of tree branch lengths is computed and used as input to HDBSCAN. This labeling is intended as a predicted starting point and is used to label and order the genes in the corresponding output 03h_Gene_Data_90_annotated.tsv data file. The researcher can change the labels in the data file and rerun the 03c plotting script as needed.

![Training set tree](https://github.com/rotheconrad/ROCker_Macrolide_Models/blob/main/02_mph/00_figures/03a_mph_training_set_tree.png)

#### Testing set - clade/cluster labeled tree

The below image is the default IQTree and clade/cluster labeling from the 03c_Plot_Annotated_Tree_v2.py script. This script uses the HDBSCAN algorithm to guess at gene clades/clusters. An all vs all distance matrix of tree branch lengths is computed and used as input to HDBSCAN. This labeling is intended as a predicted starting point and is used to label and order the genes in the corresponding output 03h_Gene_Data_90_annotated.tsv data file. The researcher can change the labels in the data file and rerun the 03c plotting script as needed.

![Testing set tree](https://github.com/rotheconrad/ROCker_Macrolide_Models/blob/main/02_mph/00_figures/03b_mph_testing_set_tree.png)

Review the Data Table and the Phylogenetic tree to make positive and negative UniProt ID selections to give to ROCkOut. Clade/Cluster labels can be changed in the tsv file and the tree can be replotted.

- for the test set be sure to remove all the RepSeqs from the UniProt ID list.
- for the training set remove any Gene Names that aren't a UniProt ID.
- for the training set, look at 02b_filter/02b_mph_matches_pID100.fa and include some sequences that 100% similar to any RepSeqs without a UniProt ID. You only need 1 per RepSeq. This will include some different genomes in the training set that have the identical copies.

We decided to build two models:

1. A model using all the genes in the *Training set* tree as positive references to make a model that captures all mph genes. This model is useful for broad surveillance mph gene functions in microbial communities. We refer to this model as mphAll. Looking at the multiple sequence alignment, phylogenetic tree, gene annotations and taxonomic classifications we did not see any cleary distinguished outgroup form in our results. The mphB reference sequence falls close to one end of our tree and the mphG and mphE reference sequences fall at the other end. To the best of our knowledge, all genes retrieved are bioinformatically likely to perform the mph gene function, and no clear closely related sequences of diverging function were identified in the similar sequence space.
1. A model focusing specifically on the mphA gene clade of the *Training set* using only the sequences from this clade as positive references and using all other sequences as negative references. We refer to this model as mphA. The mphA clade may be of particlar intrest to some ARG researchers and formed a clearly distinguished clade in sequence space from the other mph gene clades.

We make two similar selections for the *Testing set*:

1. We select all genes in the *Testing set* to use as positive references. The genomes containing these genes will be used to create a mock metagenome to challenge the mphAll ROCker model.
1. We select genes in the *Testing set* specifically from the mphA gene clade. The genomes containing these genes will be used to create a mock metagenome to challenge the mphA ROCker model.

**Results**

- mphAll_train_positive.txt
- mphAll_test_positive.txt
- mphA_train_pos.txt, mphA_train_neg.txt
- mphA_test_pos.txt, mphA_test_neg.txt

# Step 04: Build ROCker models

*I don't have the final figures yet. Waiting on Kenji to finalize the ROCkOut code.*

This is a paragraph about any modification to the input UniProt ID lists when building the ROCker model

This is an iterative process of building a model, investigating the results, and tracking down all the peculiarities to arrive at a final postive and/or negative UniProt ID set.

```bash
# mphAll model
mkdir 01a_mphAll_train 01b_mphAll_test 01c_alignment_label_testing

# training set
cd 01a_mphAll_train
sbatch --export odir=model,pos=mphAll_train_positive.txt /ROCkOut/00b_sbatch/ROCkOut_posOnly.sbatch

# testing set
cd 01b_mphAll_test
sbatch --export odir=model,pos=mphAll_test_positive.txt /ROCkOut/00b_sbatch/ROCkOut_posOnly.sbatch

######

# mphA model
mkdir 02a_mphA_train 02b_mphA_test 02c_alignment_label_testing

# training set
cd 02a_mphA_train
sbatch --export odir=model,pos=mphA_train_pos.txt,neg=mphA_train_neg.txt /ROCkOut/00b_sbatch/ROCkOut.sbatch

# testing set
cd 02b_mphA_test
sbatch --export odir=model,pos=mphA_test_pos.txt,neg=mphA_test_neg.txt /ROCkOut/00b_sbatch/ROCkOut.sbatch
```

**Results**

That's it. The models are built. I'll include final figures here once I have them.

ROCkOut as an "-extract" feature that collects the labeled, simulated, short reads into fasta files. We will use these fasta files from the *Testing set* to execute ROCkOut's align, filter, and place functions against the *Training set* models.

# Step 05: Test ROCker models

We build models for the *Testing set* because 1) this is a convenient way to use the ROCkOut figures and interactive plots to fine tune the sequence selection and 2) ROCkOut downloads genomes the genes are part of, simulates short sequence reads for the genomes, and labels the reads and positive, negative, target, or non_target. We can collect these reads into a mock metagenome that is labeled so we can track them and use them to score the results.

We test 4 different read lengths for each model.

Scripts for this section are in the 00_scripts directory of this GitHub Repo. Scripts for all previous sections were from either the ROCkIn or ROCkOut repos.

### a. Run mock metagenomes against models using ROCkOut

```bash
# mphAll model

cd 01a_mphAll_train

sbatch --export infile=../01b_mphAll_test/simulated_reads/simread_100_raw_reads.fasta,model=model,outpre=test_score_100 ../../00b_sbatch/ROCkOut_posOnly_score.sbatch

sbatch --export infile=../01b_mphAll_test/simulated_reads/simread_150_raw_reads.fasta,model=model,outpre=test_score_150 ../../00b_sbatch/ROCkOut_posOnly_score.sbatch

sbatch --export infile=../01b_mphAll_test/simulated_reads/simread_250_raw_reads.fasta,model=model,outpre=test_score_250 ../../00b_sbatch/ROCkOut_posOnly_score.sbatch

sbatch --export infile=../01b_mphAll_test/simulated_reads/simread_300_raw_reads.fasta,model=model,outpre=test_score_300 ../../00b_sbatch/ROCkOut_posOnly_score.sbatch

######

# mphA model

cd 02a_mphA_train

sbatch --export infile=../02b_mphA_test/simulated_reads/simread_100_raw_reads.fasta,model=model,outpre=test_score_100 ../../00b_sbatch/ROCkOut_score.sbatch

sbatch --export infile=../02b_mphA_test/simulated_reads/simread_150_raw_reads.fasta,model=model,outpre=test_score_150 ../../00b_sbatch/ROCkOut_score.sbatch

sbatch --export infile=../02b_mphA_test/simulated_reads/simread_250_raw_reads.fasta,model=model,outpre=test_score_250 ../../00b_sbatch/ROCkOut_score.sbatch

sbatch --export infile=../02b_mphA_test/simulated_reads/simread_300_raw_reads.fasta,model=model,outpre=test_score_300 ../../00b_sbatch/ROCkOut_score.sbatch
```

### b. Create and search custom HMM model

Use [HMMER](http://hmmer.org/) to build a model from the multiple alignment created by ROCkOut. This multiple alignment is the same genes and alignment ROCkOut uses to create it's model, so we use it to create an HMM model for comparison.

We convert the mock metagenome (nucleotide fasta) to all 6 reading frames of amino acid sequence because this is what BLASTx (or diamond BLASTx) does internal to align reads. In this way the results our HMM search will be comparable to BLASTx. We do this for all 4 read lengths.

We use [HMMER](http://hmmer.org/) to search our amino acid sequences against the model and we use a custom script to filter results for best match because there is a low rate of more than one reading frame per read finding a match.

```bash
# mphAll model

# build model (hmmbuild)
hmmbuild mphAll_model-final.hmm model/shared_files/multiple_alignment/complete_multiple_alignment_aa.fasta

# convert mock metagenome nucleotides to 6 aa frames
python 00_scripts/nuc_fasta_6aa_frame_fasta.py -i ../01b_mphAll_test/simulated_reads/simread_100_raw_reads.fasta -o test_score_100/simread_100_raw_reads_6frames.faa

python 00_scripts/nuc_fasta_6aa_frame_fasta.py -i ../01b_mphAll_test/simulated_reads/simread_150_raw_reads.fasta -o test_score_150/simread_150_raw_reads_6frames.faa

python 00_scripts/nuc_fasta_6aa_frame_fasta.py -i ../01b_mphAll_test/simulated_reads/simread_250_raw_reads.fasta -o test_score_250/simread_250_raw_reads_6frames.faa

python 00_scripts/nuc_fasta_6aa_frame_fasta.py -i ../01b_mphAll_test/simulated_reads/simread_300_raw_reads.fasta -o test_score_300/simread_300_raw_reads_6frames.faa

# Map Reads (hmmsearch)
sbatch --export model=mphAll_model-final.hmm,infile=test_score_100/simread_100_raw_reads_6frames.faa,outfile=test_score_100/mphAll-final_hmmSearch_100.tsv 00_scripts/hmmsearch.sbatch

sbatch --export model=mphAll_model-final.hmm,infile=test_score_150/simread_150_raw_reads_6frames.faa,outfile=test_score_150/mphAll-final_hmmSearch_150.tsv 00_scripts/hmmsearch.sbatch

sbatch --export model=mphAll_model-final.hmm,infile=test_score_250/simread_250_raw_reads_6frames.faa,outfile=test_score_250/mphAll-final_hmmSearch_250.tsv 00_scripts/hmmsearch.sbatch

sbatch --export model=mphAll_model-final.hmm,infile=test_score_300/simread_300_raw_reads_6frames.faa,outfile=test_score_300/mphAll-final_hmmSearch_300.tsv 00_scripts/hmmsearch.sbatch


# Filter for best hit and score
python 00_scripts/besthit_filter_hmm.py -i test_score_100/mphAll-final_hmmSearch_100.tsv -o test_score_100/mphAll-final_hmmSearch_100_filtered.tsv

python 00_scripts/besthit_filter_hmm.py -i test_score_150/mphAll-final_hmmSearch_150.tsv -o test_score_150/mphAll-final_hmmSearch_150_filtered.tsv

python 00_scripts/besthit_filter_hmm.py -i test_score_250/mphAll-final_hmmSearch_250.tsv -o test_score_250/mphAll-final_hmmSearch_250_filtered.tsv

python 00_scripts/besthit_filter_hmm.py -i test_score_300/mphAll-final_hmmSearch_300.tsv -o test_score_300/mphAll-final_hmmSearch_300_filtered.tsv

######

# mphA model

# build model (hmmbuild)
hmmbuild mphA_model-final.hmm model/shared_files/multiple_alignment/complete_multiple_alignment_aa.fasta

# convert mock metagenome nucleotides to 6 aa frames
python 00_scripts/nuc_fasta_6aa_frame_fasta.py -i ../02b_mphA_test/simulated_reads/simread_100_raw_reads.fasta -o test_score_100/simread_100_raw_reads_6frames.faa

python 00_scripts/nuc_fasta_6aa_frame_fasta.py -i ../02b_mphA_test/simulated_reads/simread_150_raw_reads.fasta -o test_score_150/simread_150_raw_reads_6frames.faa

python 00_scripts/nuc_fasta_6aa_frame_fasta.py -i ../02b_mphA_test/simulated_reads/simread_250_raw_reads.fasta -o test_score_250/simread_250_raw_reads_6frames.faa

python 00_scripts/nuc_fasta_6aa_frame_fasta.py -i ../02b_mphA_test/simulated_reads/simread_300_raw_reads.fasta -o test_score_300/simread_300_raw_reads_6frames.faa

# Map Reads (hmmsearch)
sbatch --export model=mphA_model-final.hmm,infile=test_score_100/simread_100_raw_reads_6frames.faa,outfile=test_score_100/mphA-final_hmmSearch_100.tsv 00_scripts/hmmsearch.sbatch

sbatch --export model=mphA_model-final.hmm,infile=test_score_150/simread_150_raw_reads_6frames.faa,outfile=test_score_150/mphA-final_hmmSearch_150.tsv 00_scriptsh/hmmsearch.sbatch

sbatch --export model=mphA_model-final.hmm,infile=test_score_250/simread_250_raw_reads_6frames.faa,outfile=test_score_250/mphA-final_hmmSearch_250.tsv 00_scripts/hmmsearch.sbatch

sbatch --export model=mphA_model-final.hmm,infile=test_score_300/simread_300_raw_reads_6frames.faa,outfile=test_score_300/mphA-final_hmmSearch_300.tsv 00_scripts/hmmsearch.sbatch

# Filter for best hit and score
python 00_scripts/besthit_filter_hmm.py -i test_score_100/mphA-final_hmmSearch_100.tsv -o test_score_100/mphA-final_hmmSearch_100_filtered.tsv

python 00_scripts/besthit_filter_hmm.py -i test_score_150/mphA-final_hmmSearch_150.tsv -o test_score_150/mphA-final_hmmSearch_150_filtered.tsv

python 00_scripts/besthit_filter_hmm.py -i test_score_250/mphA-final_hmmSearch_250.tsv -o test_score_250/mphA-final_hmmSearch_250_filtered.tsv

python 00_scripts/besthit_filter_hmm.py -i test_score_300/mphA-final_hmmSearch_300.tsv -o test_score_300/mphA-final_hmmSearch_300_filtered.tsv
```

**Results:**

#### mphAll - 100 bp reads
- Total number of entries in hmm file: 27806
- Number of duplicate hmm matches: 0
- Number of best hit entries written to new file: 27806

#### mphAll - 150 bp reads
- Total number of entries in hmm file: 26264
- Number of duplicate hmm matches: 2
- Number of best hit entries written to new file: 26262 

#### mphAll - 250 bp reads
- Total number of entries in hmm file: 19582
- Number of duplicate hmm matches: 11
- Number of best hit entries written to new file: 19571 

#### mphAll - 300 bp reads
- Total number of entries in hmm file: 17942
- Number of duplicate hmm matches: 13
- Number of best hit entries written to new file: 17929

#### mphA - 100 bp reads
- Total number of entries in hmm file: 16628
- Number of duplicate hmm matches: 0
- Number of best hit entries written to new file: 16628

#### mphA - 150 bp reads
- Total number of entries in hmm file: 17457
- Number of duplicate hmm matches: 2
- Number of best hit entries written to new file: 17455

#### mphA - 250 bp reads
- Total number of entries in hmm file: 13195
- Number of duplicate hmm matches: 0
- Number of best hit entries written to new file: 13195

#### mphA - 300 bp reads
- Total number of entries in hmm file: 11702
- Number of duplicate hmm matches: 4
- Number of best hit entries written to new file: 11698 

### c. Compile, score, and plot results

input rockout score results and hmm results and build bar plot.

For hmm results need to get total positive reads from mock metagenome. Find positive reads not mapped by hmm vs mapped by hmm and find non-positive reads mapped by hmm.
two files 1) mock metagenome 2) hmmsearch result

rocker results same thing but blastx results are split into passing and failing files.
three files 1) mock metagenome 2) passing 3) failing
from: /storage/home/hcoda1/9/rconrad6/scratch/ROCkOut/02_mph/02a_mphA_train

```bash
# mphAll model

python 00_scripts/score_rocker_model.py -mm ../01b_mphAll_test/simulated_reads/simread_100_raw_reads.fasta -bx test_score_100/alignments/simread_100_raw_reads.fasta_ROCkOut_alignments.txt -hm test_score_100/mphAll-final_hmmSearch_100_filtered.tsv -rp test_score_100/ROCkOut_passing_alignments/simread_100_raw_reads.ROCkOut_passing.txt -rf test_score_100/ROCkOut_failing_alignments/simread_100_raw_reads.ROCkOut_failing.txt -o test_score_100/test_scores_100

python 00_scripts/00c_scripts/score_rocker_model.py -mm ../01b_mphAll_test/simulated_reads/simread_150_raw_reads.fasta -bx test_score_150/alignments/simread_150_raw_reads.fasta_ROCkOut_alignments.txt -hm test_score_150/mphAll-final_hmmSearch_150_filtered.tsv -rp test_score_150/ROCkOut_passing_alignments/simread_150_raw_reads.ROCkOut_passing.txt -rf test_score_150/ROCkOut_failing_alignments/simread_150_raw_reads.ROCkOut_failing.txt -o test_score_150/test_scores_150

python 00_scripts/00c_scripts/score_rocker_model.py -mm ../01b_mphAll_test/simulated_reads/simread_250_raw_reads.fasta -bx test_score_250/alignments/simread_250_raw_reads.fasta_ROCkOut_alignments.txt -hm test_score_250/mphAll-final_hmmSearch_250_filtered.tsv -rp test_score_250/ROCkOut_passing_alignments/simread_250_raw_reads.ROCkOut_passing.txt -rf test_score_250/ROCkOut_failing_alignments/simread_250_raw_reads.ROCkOut_failing.txt -o test_score_250/test_scores_250

python 00_scripts/00c_scripts/score_rocker_model.py -mm ../01b_mphAll_test/simulated_reads/simread_300_raw_reads.fasta -bx test_score_300/alignments/simread_300_raw_reads.fasta_ROCkOut_alignments.txt -hm test_score_300/mphAll-final_hmmSearch_300_filtered.tsv -rp test_score_300/ROCkOut_passing_alignments/simread_300_raw_reads.ROCkOut_passing.txt -rf test_score_300/ROCkOut_failing_alignments/simread_300_raw_reads.ROCkOut_failing.txt -o test_score_300/test_scores_300

# mphA model

python 00_scripts/score_rocker_model.py -mm ../02b_mphA_test/simulated_reads/simread_100_raw_reads.fasta -bx test_score_100/alignments/simread_100_raw_reads.fasta_ROCkOut_alignments.txt -hm test_score_100/mphA-final_hmmSearch_100_filtered.tsv -rp test_score_100/ROCkOut_passing_alignments/simread_100_raw_reads.ROCkOut_passing.txt -rf test_score_100/ROCkOut_failing_alignments/simread_100_raw_reads.ROCkOut_failing.txt -o test_score_100/test_scores_100

python 00_scripts/score_rocker_model.py -mm ../02b_mphA_test/simulated_reads/simread_150_raw_reads.fasta -bx test_score_150/alignments/simread_150_raw_reads.fasta_ROCkOut_alignments.txt -hm test_score_150/mphA-final_hmmSearch_150_filtered.tsv -rp test_score_150/ROCkOut_passing_alignments/simread_150_raw_reads.ROCkOut_passing.txt -rf test_score_150/ROCkOut_failing_alignments/simread_150_raw_reads.ROCkOut_failing.txt -o test_score_150/test_scores_150

python 00_scripts/score_rocker_model.py -mm ../02b_mphA_test/simulated_reads/simread_250_raw_reads.fasta -bx test_score_250/alignments/simread_250_raw_reads.fasta_ROCkOut_alignments.txt -hm test_score_250/mphA-final_hmmSearch_250_filtered.tsv -rp test_score_250/ROCkOut_passing_alignments/simread_250_raw_reads.ROCkOut_passing.txt -rf test_score_250/ROCkOut_failing_alignments/simread_250_raw_reads.ROCkOut_failing.txt -o test_score_250/test_scores_250

python 00_scripts/score_rocker_model.py -mm ../02b_mphA_test/simulated_reads/simread_300_raw_reads.fasta -bx test_score_300/alignments/simread_300_raw_reads.fasta_ROCkOut_alignments.txt -hm test_score_300/mphA-final_hmmSearch_300_filtered.tsv -rp test_score_300/ROCkOut_passing_alignments/simread_300_raw_reads.ROCkOut_passing.txt -rf test_score_300/ROCkOut_failing_alignments/simread_300_raw_reads.ROCkOut_failing.txt -o test_score_300/test_scores_300
```

**Results:**

#### mphAll - 100 bp reads
- Positive read count from metagenome: 42985
- Positive read count from Diamond/BlastX alignment: 39714
- Positive reads not aligned by BlastX: 3271
- Positive read count from ROCkOut: 39714
- Positive reads not reported by ROCkOut: 3271
- Positive read count from hmmsearch: 26080
- Positive reads not aligned by hmmsearch: 16905

#### mphAll - 150 bp reads
- Positive read count from metagenome: 29734
- Positive read count from Diamond/BlastX alignment: 27851
- Positive reads not aligned by BlastX: 1883
- Positive read count from ROCkOut: 27851
- Positive reads not reported by ROCkOut: 1883
- Positive read count from hmmsearch: 23504
- Positive reads not aligned by hmmsearch: 6230

#### mphAll - 250 bp reads
- Positive read count from metagenome: 19185
- Positive read count from Diamond/BlastX alignment: 18119
- Positive reads not aligned by BlastX: 1066
- Positive read count from ROCkOut: 18119
- Positive reads not reported by ROCkOut: 1066
- Positive read count from hmmsearch: 16551
- Positive reads not aligned by hmmsearch: 2634

#### mphAll - 300 bp reads
- Positive read count from metagenome: 16227
- Positive read count from Diamond/BlastX alignment: 15486
- Positive reads not aligned by BlastX: 741
- Positive read count from ROCkOut: 15486
- Positive reads not reported by ROCkOut: 741
- Positive read count from hmmsearch: 14317
- Positive reads not aligned by hmmsearch: 1910

#### mphA - 100 bp reads
- Positive read count from metagenome: 15548
- Positive read count from Diamond/BlastX alignment: 15396
- Positive reads not aligned by BlastX: 152
- Positive read count from ROCkOut: 15396
- Positive reads not reported by ROCkOut: 152
- Positive read count from hmmsearch: 8633
- Positive reads not aligned by hmmsearch: 6915

#### mphA - 150 bp reads
- Positive read count from metagenome: 10724
- Positive read count from Diamond/BlastX alignment: 10622
- Positive reads not aligned by BlastX: 102
- Positive read count from ROCkOut: 10622
- Positive reads not reported by ROCkOut: 102
- Positive read count from hmmsearch: 8751
- Positive reads not aligned by hmmsearch: 1973

#### mphA - 250 bp reads
- Positive read count from metagenome: 6698
- Positive read count from Diamond/BlastX alignment: 6599
- Positive reads not aligned by BlastX: 99
- Positive read count from ROCkOut: 6599
- Positive reads not reported by ROCkOut: 99
- Positive read count from hmmsearch: 5835
- Positive reads not aligned by hmmsearch: 863

#### mphA - 300 bp reads
- Positive read count from metagenome: 5660
- Positive read count from Diamond/BlastX alignment: 5593
- Positive reads not aligned by BlastX: 67
- Positive read count from ROCkOut: 5593
- Positive reads not reported by ROCkOut: 67
- Positive read count from hmmsearch: 5089
- Positive reads not aligned by hmmsearch: 571

### d. Positive read count sanity check

This section was performed during the development of ROCkOut but may still serve useful those with extra curiosity.

blast, diamond, and hmm align fewer reads than rocker labels as positive during
the read simulation step.

The purpose of this script is to investigate that.

It plots histograms for the simulated reads labeled as positive split between reads that align and read reads that do not align by blastx or hmm or don't pass the rocker filter.

I just show results for the mphAll model here, but the results were and can be repeated with mphA or any other model.

```bash
python 00_scripts/investigate_pos_read_discrepancy.py -mm ../01b_mphAll_test/simulated_reads/simread_100_raw_reads.fasta -bx test_score_100/alignments/simread_100_raw_reads.fasta_ROCkOut_alignments.txt -hm test_score_100/mphAll-final_hmmSearch_100_filtered.tsv -rp test_score_100/ROCkOut_passing_alignments/simread_100_raw_reads.ROCkOut_passing.txt -rf test_score_100/ROCkOut_failing_alignments/simread_100_raw_reads.ROCkOut_failing.txt -o test_score_100/pos_read_test_100

python 00_scripts/investigate_pos_read_discrepancy.py -mm ../01b_mphAll_test/simulated_reads/simread_150_raw_reads.fasta -bx test_score_150/alignments/simread_150_raw_reads.fasta_ROCkOut_alignments.txt -hm test_score_150/mphAll-final_hmmSearch_150_filtered.tsv -rp test_score_150/ROCkOut_passing_alignments/simread_150_raw_reads.ROCkOut_passing.txt -rf test_score_150/ROCkOut_failing_alignments/simread_150_raw_reads.ROCkOut_failing.txt -o test_score_150/pos_read_test_150

python 00_scripts/investigate_pos_read_discrepancy.py -mm ../01b_mphAll_test/simulated_reads/simread_250_raw_reads.fasta -bx test_score_250/alignments/simread_250_raw_reads.fasta_ROCkOut_alignments.txt -hm test_score_250/mphAll-final_hmmSearch_250_filtered.tsv -rp test_score_250/ROCkOut_passing_alignments/simread_250_raw_reads.ROCkOut_passing.txt -rf test_score_250/ROCkOut_failing_alignments/simread_250_raw_reads.ROCkOut_failing.txt -o test_score_250/pos_read_test_250

python 00_scripts/00c_scripts/investigate_pos_read_discrepancy.py -mm ../01b_mphAll_test/simulated_reads/simread_300_raw_reads.fasta -bx test_score_300/alignments/simread_300_raw_reads.fasta_ROCkOut_alignments.txt -hm test_score_300/mphAll-final_hmmSearch_300_filtered.tsv -rp test_score_300/ROCkOut_passing_alignments/simread_300_raw_reads.ROCkOut_passing.txt -rf test_score_300/ROCkOut_failing_alignments/simread_300_raw_reads.ROCkOut_failing.txt -o test_score_300/pos_read_test_300
```

**Results:**

#### mphAll - 100 bp reads
- Total Positives reads created: 42985
- Total reads not aligned: 3271
- Total reads not in rocker filter output: 3271
- Total reads not in hmm search results: 16905

#### mphAll - 150 bp reads
- Total Positives reads created: 29734
- Total reads not aligned: 1883
- Total reads not in rocker filter output: 1883
- Total reads not in hmm search results: 6230

#### mphAll - 250 bp reads
- Total Positives reads created: 19185
- Total reads not aligned: 1066
- Total reads not in rocker filter output: 1066
- Total reads not in hmm search results: 2634

#### mphAll - 300 bp reads
- Total Positives reads created: 16227
- Total reads not aligned: 741
- Total reads not in rocker filter output: 741
- Total reads not in hmm search results: 1910

There are plots too. I'll put some examples. This section and plots have mainly been to double check ROCkOut and track down discrepencies ect. during development.

### e. Verify TP FP TN FN labels and scoring analysis (additional visuals)

The score_rocker_model.py script writes out a fasta directory with file names labeled for filter method and failed, passed, TP, FP, TN, or FN.

The goal of this section is to test the results of mapping the reads in nucleotide format to the genes in nucleotide with regular blastn vs. mapping the the reads in nucleotides to the genes in amino acids with blastx.

First step is to create gene databases with nucl and prot, then to map all the different read files in the fasta directory created by score_rocker_model.py.

The second step will be to create some figures to visualize the results.

I just show results for the mphAll model here, but the results were and can be repeated with mphA or any other model.

```bash
cd 01c_alignment_label_testing
mkdir 00a_log 00b_alignment_dbs 01a_alignments

# copy files:
cat 01a_mphAll_train/model/positive/*/target_protein/*_nt.fasta > 00b_alignment_dbs/combined_proteins_nt.fasta
cat 01a_mphAll_train/model/positive/*/target_protein/*_AA.fasta > 00b_alignment_dbs/combined_proteins_aa.fasta

# create alignment databases
# using blast-plus v2.10.1
makeblastdb -dbtype prot -in 00b_alignment_dbs/combined_proteins_aa.fasta
makeblastdb -dbtype nucl -in 00b_alignment_dbs/combined_proteins_nt.fasta
# using diamond v2.0.15
diamond makedb --in 00b_alignment_dbs/combined_proteins_aa.fasta -d 00b_alignment_dbs/combined_proteins_aa

# Run blastn and blastx for each mock metagenome
# use sbatch script. script should include best hit filter post alignments.

# blastn 100
sbatch -e 00a_log/mphAll_100_BLASTn.err -o 00a_log/mphAll_100_BLASTn.out --export ref=00b_alignment_dbs/combined_proteins_nt.fasta,qry=../01b_mphAll_test/simulated_reads/simread_100_raw_reads.fasta,odir=01a_alignments 00_scripts/blastn.sbatch

# blastx 100
sbatch -e 00a_log/mphAll_100_BLASTx.err -o 00a_log/mphAll_100_BLASTx.out --export ref=00b_alignment_dbs/combined_proteins_aa.fasta,qry=../01b_mphAll_test/simulated_reads/simread_100_raw_reads.fasta,odir=01a_alignments 00_scripts/blastx.sbatch 

# diamondx 100
sbatch -e 00a_log/mphAll_100_DIAMONDx.err -o 00a_log/mphAll_100_DIAMONDx.out --export ref=00b_alignment_dbs/combined_proteins_aa,qry=../01b_mphAll_test/simulated_reads/simread_100_raw_reads.fasta,odir=01a_alignments 00_scripts/diamondx.sbatch 

# blastn 150
sbatch -e 00a_log/mphAll_150_BLASTn.err -o 00a_log/mphAll_150_BLASTn.out --export ref=00b_alignment_dbs/combined_proteins_nt.fasta,qry=../01b_mphAll_test/simulated_reads/simread_150_raw_reads.fasta,odir=01a_alignments 00_scripts/blastn.sbatch

# blastx 150
sbatch -e 00a_log/mphAll_150_BLASTx.err -o 00a_log/mphAll_150_BLASTx.out --export ref=00b_alignment_dbs/combined_proteins_aa.fasta,qry=../01b_mphAll_test/simulated_reads/simread_150_raw_reads.fasta,odir=01a_alignments 00_scripts/blastx.sbatch 

# diamondx 150
sbatch -e 00a_log/mphAll_150_DIAMONDx.err -o 00a_log/mphAll_150_DIAMONDx.out --export ref=00b_alignment_dbs/combined_proteins_aa,qry=../01b_mphAll_test/simulated_reads/simread_150_raw_reads.fasta,odir=01a_alignments 00_scripts/diamondx.sbatch 

# blastn 250
sbatch -e 00a_log/mphAll_250_BLASTn.err -o 00a_log/mphAll_250_BLASTn.out --export ref=00b_alignment_dbs/combined_proteins_nt.fasta,qry=../01b_mphAll_test/simulated_reads/simread_250_raw_reads.fasta,odir=01a_alignments 00_scripts/blastn.sbatch

# blastx 250
sbatch -e 00a_log/mphAll_250_BLASTx.err -o 00a_log/mphAll_250_BLASTx.out --export ref=00b_alignment_dbs/combined_proteins_aa.fasta,qry=../01b_mphAll_test/simulated_reads/simread_250_raw_reads.fasta,odir=01a_alignments 00_scripts/blastx.sbatch 

# diamondx 250
sbatch -e 00a_log/mphAll_250_DIAMONDx.err -o 00a_log/mphAll_250_DIAMONDx.out --export ref=00b_alignment_dbs/combined_proteins_aa,qry=../01b_mphAll_test/simulated_reads/simread_250_raw_reads.fasta,odir=01a_alignments 00_scripts/diamondx.sbatch 

# blastn 300
sbatch -e 00a_log/mphAll_300_BLASTn.err -o 00a_log/mphAll_300_BLASTn.out --export ref=00b_alignment_dbs/combined_proteins_nt.fasta,qry=../01b_mphAll_test/simulated_reads/simread_300_raw_reads.fasta,odir=01a_alignments 00_scripts/blastn.sbatch

# blastx 300
sbatch -e 00a_log/mphAll_300_BLASTx.err -o 00a_log/mphAll_300_BLASTx.out --export ref=00b_alignment_dbs/combined_proteins_aa.fasta,qry=../01b_mphAll_test/simulated_reads/simread_300_raw_reads.fasta,odir=01a_alignments 00_scripts/blastx.sbatch 

# diamondx 300
sbatch -e 00a_log/mphAll_300_DIAMONDx.err -o 00a_log/mphAll_300_DIAMONDx.out --export ref=00b_alignment_dbs/combined_proteins_aa,qry=../01b_mphAll_test/simulated_reads/simread_300_raw_reads.fasta,odir=01a_alignments 00_scripts/diamondx.sbatch 

# visualize results:

# wrote python script to read blastn and blastx results and to use the fasta files to label the passing, failing, TP, FP, TN, FN reads for each.

python 00_scripts/alignment_label_test.py -bn 01a_alignments/simread_100_raw_reads.blastn -bx 01a_alignments/simread_100_raw_reads.blastx -dx 01a_alignments/simread_100_raw_reads.dmndx -fd ../01a_mphAll_train/test_score_100/fastas -od 02a_plots_100

python 00_scripts/alignment_label_test.py -bn 01a_alignments/simread_150_raw_reads.blastn -bx 01a_alignments/simread_150_raw_reads.blastx -dx 01a_alignments/simread_150_raw_reads.dmndx -fd ../01a_mphAll_train/test_score_150/fastas -od 02b_plots_150

python 00_scripts/alignment_label_test.py -bn 01a_alignments/simread_250_raw_reads.blastn -bx 01a_alignments/simread_250_raw_reads.blastx -dx 01a_alignments/simread_250_raw_reads.dmndx -fd ../01a_mphAll_train/test_score_250/fastas -od 02c_plots_250

python 00_scripts/alignment_label_test.py -bn 01a_alignments/simread_300_raw_reads.blastn -bx 01a_alignments/simread_300_raw_reads.blastx -dx 01a_alignments/simread_300_raw_reads.dmndx -fd ../01a_mphAll_train/test_score_300/fastas -od 02d_plots_300
```

This creates a lot of figures. I'll put some of them here and I'll put the rest in their own directory of the GitHub repo. I've ran this a few times already during the development of ROCkOut but I'll wait to put the final ones here once Kenji has the code all sorted.

### f. Build pie trees with pplacer and itol.

ROCkOut already includes a place function. This step is really just retrieve the right file and using the [iTol](https://itol.embl.de/) website.

I'll put the figure here once I have the final ROCkOut results after Kenji and I get the code sorted.

### g. Make final phylogenetic tree figure

1. We want to create a nice phylogenetic tree for the publication figure
2. ROCkOut isn't able to retrieve complete data for all input UniProt IDs so we want to create the phylogenetic tree from only the sequences ROCkOut is able to use to build the model.
3. We want to include our initial curated sequences in this tree.
4. I'll use clustal omega for the alignment, ali view and trimal to clean it up, and IQTree to compute a maximum likelihood tree with bootstrap support. 

```bash
mkdir final_tree

# copy the sequences from the mphAll model
cp model/shared_files/combined_genomes/combined_proteins_aa.fasta final_tree/

# Clean seq names
python 00_scripts/clean_seq_names.py -i final_tree/combined_proteins_aa.fasta

# Create multiple alignment by first aligning the mph SeedSeqs and then adding the combined proteins to the SeedSeq alignment.
# using clustal omega version 1.2.4
sbatch --export verified=mph_SeedSeqs.faa,newseqs=combined_proteins_aa.fasta,output=final_tree_seqs.aln 00_scripts/seq_alignment.sbatch 

# use Aliview to view and edit the multiple alignment.
# Use Trimal to trim up the alignment.
sbatch --export input=final_tree_seqs.aln,output=final_tree_seqs_trimmed.aln 00_scripts/seq_trim.sbatch

# Build tree using IQ-Tree
sbatch --export input=final_tree_seqs_trimmed.aln,outpre=final_tree 00_scripts//IQTree.sbatch
```

Upload the final_tree.treefile to iTol and customize as necessary.

again I'll put the final figure here once Kenji and I get the ROCkOut code sorted and finalized.

