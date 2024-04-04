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

The initial mphAll training and testing models had a few Non_target genes up in the positive target space. We checked them out and made the following decisions before rebuilding the model:

##### For the training model:

- CP009687.1 - Clostridium aceticum strain DSM 1496 - This genome has a gene "A0A0D8IAI0_9CLOT" annotated as Aminoglycoside phosphotransferase that is highly similar to other macrolide 2'- phosphotransferase and mphB genes. This gene was picked up by ROCkIn and added to the positive references for them mphAll model. This genome also has a second gene "A0A0D8I8G2_9CLOT" annotated as hypothetical protein on NCBI and as Aminoglycoside phosphotransferase domain-containing protein on UniProt. It is highly similar to other Mph(B) family macrolide 2'-phosphotransferase genes. We added this second gene to the positive reference sequences for the mphAll model. This genome also has a third gene "" annotated as

- CP017269.1 - Geosporobacter ferrireducens strain IRF9 - this genome has two "macrolide 2'-phosphotransferase" genes adjacent to each other with 45.36% sequence identity. "A0A1D8GF31_9CLOT" is in the original positive set and we've added "A0A1D8GF44_9CLOT" to the positive set. The second gene also showed high sequence similarity with many other macrolide phosphotransferase and mphB annotated genes.

- CP017603.1 - Clostridium formicaceticum strain ATCC 27076, complete genome - This genome has one Macrolide 2'-phosphotransferase gene picked up by ROCkIn and included in the positive reference set "A0A1D9FGX1_9CLOT" It also has a second degraded partial gene with internal stop codon. This partial gene has no match in the UniProt database and so this non_target reads can't be classified with ROCkOut. I removed this Uniprot ID from the positive reference list.

- CP054705.1 - Salicibibacter cibarius strain NKC5-3 chromosome, complete genome - This genome has two copies of genes labelled as phosphotransferase that are 43.33% sequence identity with eachother but both share high sequence identity with other genes annotated as Macrolide 2'-phosphotransferase genes. The "A0A7T6Z3M7_9BACI" gene was picked up and included in the initial ROCkIn run. We added the second gene "A0A7T6Z7R1_9BACI" to the positive reference set.

- KY399978.1 - Vibrio cholerae O139 strain ICDC-211 plasmid pVC211 - this plasmid has an "mphK" gene originally included in the positive set "A0A1Y0F4Y8_VIBCL" it also has an "mph2" gene that is 100% similar to genes in many other genomes and plasmids typically annotated as an "mphE." We added the new uniprot ID "A0A1Y0F4S0_VIBCL" to positive sequences. This plasmid also has an mphR gene that looks short but if it pops up in the model again check it out.

- MH491967.2 - Proteus mirabilis strain HFK418 plasmid pHFK418-NDM - this plasmid appears to have two full copies and two short fragment copies of genes annotated as "Aldehyde dehydrogenase." The positive reference gene associated with this genome "A0A346FVF1_PROMI" is annotated as an Aldehyde dehydrogenase which is something we saw frequently blended in with mph clades in the ROCkIn results. Blast results on UniProt show these genes with high sequence similarity to mph genes and the MSA shows conserved sites. We added the second copy "A0A2I6SHX8_PROMI" to the positive references. This copy has many 100% identical matches with mphE and mph2 annotated genes in genomes of other species. this plasmid copy gets around. Possible misannotation as Aldehyde dehydrogenase. The A0A2I6SHX8_PROMI ID appears to be linked to a 100% identical gene in a different genome leaving no way to classify these reads with ROCkOut. I've removed these IDs from the positive reference set.

##### For the testing model:

- CP017269.1 - Geosporobacter ferrireducens strain IRF9 - this genome is listed above. Apparently one copy made it in the training set and one in the testing set. Both genes are retained in the training set and removed from the testing set.

- CP054705.1 - Salicibibacter cibarius strain NKC5-3 - this genome is listed above. Apparently one copy made it in the training set and one in the testing set. Both genes are retained in the training set and removed from the testing set.

- KF648874.1 - Exiguobacterium sp. S3-2 plasmid pMC1 - this plasmid has two genes annotated as macrolide 2-phosphotransferase that have 44.37% sequence identity. The "V9Z4V9_9BACL" gene was included initially from ROCkIn and the "V9ZAR1_9BACL" gene was added to the positive reference set in this refinement round.

- MK450539.1 - Pseudomonas aeruginosa DIM-1 - this genome has two mph genes annotated as Mph(F) family macrolide 2'-phosphotransferase and Mph(E) or mph2 family macrolide 2'-phosphotransferase that have 37.63% sequence identity. The mphF "A0A5P4S8P6_PSEAI" gene was included initially from ROCkIn and the "A0A1I9WCL7_PSEAI" gene was added to the positive reference set in this refinement round. This additional gene shares 100% sequence identity with many genes in many species' genomes. The "A0A1I9WCL7_PSEAI" uniprot ID appears to link to these two plasmids KU984333.1, KF512014.1 and not the second copy in the MK450539.1 genome. I'm removing these two IDs from the set.

*note: The second genes are typically picked up by ROCkIn in the initial blast search of the uniprot database but get dereplicated with mmseqs at 90% sequence similarity with other genes. When a genome has two genes but only one gene is selected as a cluster representative, the other gene ID is missed, but shows up in the simulated reads during model building. This is why we investigate them and add them back in when necessary.*

Finally, we chose the FP+ setting with ROCkOut refine.
parameter -c fp+

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

![MCR ROCker Models](https://github.com/rotheconrad/ROCker_Macrolide_Models/blob/main/02_mph/00_figures/04a-mph-ROCker-model-250bp.png)

# Step 05: Test ROCker models

ROCkOut as an "-extract" feature that collects the labeled, simulated, short reads into fasta files. We will use these fasta files from the *Testing set* to execute ROCkOut's align, filter, and place functions against the *Training set* models.

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
- Total number of entries in hmm file: 5914
- Number of duplicate hmm matches: 0
- Number of best hit entries written to new file: 5914

#### mphAll - 150 bp reads
- Total number of entries in hmm file: 4620
- Number of duplicate hmm matches: 0
- Number of best hit entries written to new file: 4620 

#### mphAll - 250 bp reads
- Total number of entries in hmm file: 3395
- Number of duplicate hmm matches: 0
- Number of best hit entries written to new file: 3395 

#### mphAll - 300 bp reads
- Total number of entries in hmm file: 2966
- Number of duplicate hmm matches: 2
- Number of best hit entries written to new file: 2964

#### mphA - 100 bp reads
- Total number of entries in hmm file: 16570
- Number of duplicate hmm matches: 0
- Number of best hit entries written to new file: 16570

#### mphA - 150 bp reads
- Total number of entries in hmm file: 13389
- Number of duplicate hmm matches: 1
- Number of best hit entries written to new file: 13388

#### mphA - 250 bp reads
- Total number of entries in hmm file: 9440
- Number of duplicate hmm matches: 1
- Number of best hit entries written to new file: 9439

#### mphA - 300 bp reads
- Total number of entries in hmm file: 8647
- Number of duplicate hmm matches: 3
- Number of best hit entries written to new file: 8644 

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
- Positive read count from metagenome: 6872
- Positive read count from Diamond/BlastX alignment: 6674
- Positive reads not aligned by BlastX: 198
- Positive read count from ROCkOut: 6674
- Positive reads not reported by ROCkOut: 198
- Positive read count from hmmsearch: 5738
- Positive reads not aligned by hmmsearch: 1134

#### mphAll - 150 bp reads
- Positive read count from metagenome: 4868
- Positive read count from Diamond/BlastX alignment: 4735
- Positive reads not aligned by BlastX: 133
- Positive read count from ROCkOut: 4735
- Positive reads not reported by ROCkOut: 133
- Positive read count from hmmsearch: 4272
- Positive reads not aligned by hmmsearch: 596

#### mphAll - 250 bp reads
- Positive read count from metagenome: 3311
- Positive read count from Diamond/BlastX alignment: 3182
- Positive reads not aligned by BlastX: 129
- Positive read count from ROCkOut: 3182
- Positive reads not reported by ROCkOut: 129
- Positive read count from hmmsearch: 2892
- Positive reads not aligned by hmmsearch: 419

#### mphAll - 300 bp reads
- Positive read count from metagenome: 2603
- Positive read count from Diamond/BlastX alignment: 2516
- Positive reads not aligned by BlastX: 87
- Positive read count from ROCkOut: 2516
- Positive reads not reported by ROCkOut: 87
- Positive read count from hmmsearch: 2356
- Positive reads not aligned by hmmsearch: 247

#### mphA - 100 bp reads
- Positive read count from metagenome: 1331
- Positive read count from Diamond/BlastX alignment: 1331
- Positive reads not aligned by BlastX: 0
- Positive read count from ROCkOut: 1331
- Positive reads not reported by ROCkOut: 0
- Positive read count from hmmsearch: 981
- Positive reads not aligned by hmmsearch: 350

#### mphA - 150 bp reads
- Positive read count from metagenome: 1019
- Positive read count from Diamond/BlastX alignment: 1019
- Positive reads not aligned by BlastX: 0
- Positive read count from ROCkOut: 1019
- Positive reads not reported by ROCkOut: 0
- Positive read count from hmmsearch: 859
- Positive reads not aligned by hmmsearch: 160

#### mphA - 250 bp reads
- Positive read count from metagenome: 684
- Positive read count from Diamond/BlastX alignment: 681
- Positive reads not aligned by BlastX: 3
- Positive read count from ROCkOut: 681
- Positive reads not reported by ROCkOut: 3
- Positive read count from hmmsearch: 586
- Positive reads not aligned by hmmsearch: 98

#### mphA - 300 bp reads
- Positive read count from metagenome: 560
- Positive read count from Diamond/BlastX alignment: 558
- Positive reads not aligned by BlastX: 2
- Positive read count from ROCkOut: 558
- Positive reads not reported by ROCkOut: 2
- Positive read count from hmmsearch: 502
- Positive reads not aligned by hmmsearch: 58

![MCR Test Set Scoring results](https://github.com/rotheconrad/ROCker_Macrolide_Models/blob/main/02_mph/00_figures/04b-mphMockMetaScores.png)

### d. Build pie trees with pplacer and itol.

ROCkOut already includes the pie tree output with its "place" function. This step is really just retrieve the right file and using the [iTol](https://itol.embl.de/) website.

##### Read placement tree for mcr All model.

![Read placement tree for mcr All model](https://github.com/rotheconrad/ROCker_Macrolide_Models/blob/main/02_mph/00_figures/04c-mphall-phylo-placement.png)

##### Read placement tree for mcr 1 model.

![Read placement tree for mcr 1 model](https://github.com/rotheconrad/ROCker_Macrolide_Models/blob/main/02_mph/00_figures/04d-mphA-phylo-placement.png)

### e. Verify TP FP TN FN labels and scoring analysis (additional visuals)

The score_rocker_model.py script writes out a fasta directory with file names labeled for filter method and failed, passed, TP, FP, TN, or FN.

The goal of this section is to test the results of mapping the reads in nucleotide format to the genes in nucleotide with regular blastn vs. mapping the the reads in nucleotides to the genes in amino acids with blastx.

First step is to create gene databases with nucl and prot, then to map all the different read files in the fasta directory created by score_rocker_model.py.

The second step will be to create some figures to visualize the results.

I just show results for the mphAll model here, but the results were and can be repeated with mphA or any other model.

```bash
> python ../../00c_scripts/ROCkOut_filter_visualize.py -p test_score_100/ROCkOut_passing_alignments/simread_100_raw_reads.ROCkOut_passing.txt -f test_score_100/ROCkOut_failing_alignments/simread_100_raw_reads.ROCkOut_failing.txt -o test_score_100/ROCkOut_Filter_Viz_100 -a test_score_100/alignments/simread_100_raw_reads.fasta_ROCkOut_alignments.txt

> python ../../00c_scripts/ROCkOut_filter_visualize.py -p test_score_150/ROCkOut_passing_alignments/simread_150_raw_reads.ROCkOut_passing.txt -f test_score_150/ROCkOut_failing_alignments/simread_150_raw_reads.ROCkOut_failing.txt -o test_score_150/ROCkOut_Filter_Viz_150 -a test_score_150/alignments/simread_150_raw_reads.fasta_ROCkOut_alignments.txt

> python ../../00c_scripts/ROCkOut_filter_visualize.py -p test_score_250/ROCkOut_passing_alignments/simread_250_raw_reads.ROCkOut_passing.txt -f test_score_250/ROCkOut_failing_alignments/simread_250_raw_reads.ROCkOut_failing.txt -o test_score_250/ROCkOut_Filter_Viz_250 -a test_score_250/alignments/simread_250_raw_reads.fasta_ROCkOut_alignments.txt

> python ../../00c_scripts/ROCkOut_filter_visualize.py -p test_score_300/ROCkOut_passing_alignments/simread_300_raw_reads.ROCkOut_passing.txt -f test_score_300/ROCkOut_failing_alignments/simread_300_raw_reads.ROCkOut_failing.txt -o test_score_300/ROCkOut_Filter_Viz_300 -a test_score_300/alignments/simread_300_raw_reads.fasta_ROCkOut_alignments.txt
```

##### mcr All model read mapping and classification for bitscore.

![mcr All model read mapping and classification for bitscore](https://github.com/rotheconrad/ROCker_Macrolide_Models/blob/main/02_mph/00_figures/04e_verify_bitscore.png)

##### mcr All model read mapping and classification for percent sequence identity.

![mcr All model read mapping and classification for percent sequence identity](https://github.com/rotheconrad/ROCker_Macrolide_Models/blob/main/02_mph/00_figures/04f_verify_pid.png)





