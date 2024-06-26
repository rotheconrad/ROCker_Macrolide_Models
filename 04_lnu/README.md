# Working notes for development of ROCker models for lnu genes

These notes contain the step-by-step details using ROCkIn and ROCkOut to curate training and testing sequence, to build the models, and to test the models for lincosamide nucleotidyltransferase genes. ROCkIn and ROCkOut and all their dependencies are required to reproduce the commands.

lnu - lincosamide nucleotidyltransferase
Previously named as linF linA etc.

Resistance to the lincosamide antibiotic by ATP-dependent modification of the 3' and/or 4'-hydroxyl groups of the methylthiolincosamide sugar.

Lincosamides (e.g. lincomycin, clindamycin) are a class of drugs which bind to the 23s portion of the 50S subunit of bacterial ribosomes. This interaction inhibits early elongation of peptide chains by inhibiting the transpeptidase reaction, acting similarly to macrolides.

# Table of Contents

1. [Step 00: Curate starting sequences](#step-00-curate-starting-sequences)
1. [Step 01: UniProt sequence search](#step-01-uniprot-sequence-search)
1. [Step 02: Deduplicate, filter, dereplicate](#step-02-deduplicate,-filter,-dereplicate)
1. [Step 03: Build Phylogenetic tree and predict clades](#step-03-build-phylogenetic-tree-and-predict-clades)
1. [Step 04: Build ROCker models](#step-04-build-rocker-models)
1. [Step 05: Test ROCker models](#step-05-test-rocker-models)

# Step 00: Curate starting sequences

Curate starting sequences. Source from NCBI's refgene database.

### a. Find and retrieve sequences (Oct 9th, 2023)

There are multiple classes of lnu gene: A, B, C, D, F, G, H.

Search for lnu genes at ncbi refgene database returned 23 results
https://www.ncbi.nlm.nih.gov/pathogens/refgene/#lnu

Use the Download button, select "Dataset (.zip)" from the drop down menu, and check the "Reference protein" box to download (reference_protein.faa)

Renamed fasta sequences with Sublime text editor find and replace regular expressions to create shorter format ex: >WP_000700648.1_LnuA_Bacteria

Download fasta sequences and rename files. Collect the set of renamed fasta formatted sequences into a single file for amino acid sequence. These are referred to as the curated or reference sequences (lnu_SeedSeqs.faa)

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

#### Multiple sequence alignment of lnu seed sequences.

![Multiple sequence alignment of lnu seed sequences.](https://github.com/rotheconrad/ROCker_Macrolide_Models/blob/main/04_lnu/00_figures/00a_lnu_SeedSeqs_msa.png)

Multiple sequence alignment image produced with AliView highlighting majority rule concensus characters.

#### Neighbor joining phylogenetic tree of lnu seed sequences.

![Neighbor joining phylogenetic tree of lnu seed sequences.](https://github.com/rotheconrad/ROCker_Macrolide_Models/blob/main/04_lnu/00_figures/00b_lnu_SeedSeqs_tree.png)

Neighbor joining phylogenetic tree image produced with FigTree default settings.

#### Percent sequence identity matrix hierarchical clustered heatmap of lnu seed sequences.

```bash
python /ROCkIn/02_Python/00a_PIM_clustered_heatmap.py -i lnu_SeedSeqs.faa.aln.pim -o lnu_SeedSeqs.faa.aln.pim.pdf
```

![Multiple sequence alignment of lnu seed sequences.](https://github.com/rotheconrad/ROCker_Macrolide_Models/blob/main/04_lnu/00_figures/00c_lnu_SeedSeqs_pim.png)

Neighbor joining phylogenetic tree image produced with FigTree default settings.

### c. Review data and make selections

Review the multiple alignment and/or phylogenetic tree and/or heatmap and select representative sequences (i.e., remove highly similar sequences). We will use these curated sequences for a sequence search in uniprot and will not gain anything from sequences that are too similar (i.e. ≥ 90% sequence identity or one sequence from each strongly formed clade).

The lnu genes have a clear bifurcation into two divergent clades with <30% amino acid sequence identity shared between them. We chose to focus on the clade containing lnuF and so we discarded the lnuA, lnuE, lnuC, lnuD, lnuP, and lnuAN2 sequences removing them from the remainder of the analysis and model building. Our lnu model will focus on the lnu clade with the lnuF, lnuB, lnuG, and lnuH sequences.

Create a new fasta file with selected representative sequences (NCBI_lnu_reps.faa). See file for selection. Using NCBI_lnu_proteins.faa.pim.pdf, I selected 3 lnuF sequences, 1 lnuB sequence, and the lnuG, and lnuH sequences.

# Step 01: UniProt sequence search

Run a search of the curated lnu sequences (RefSeqs) against the UniProt database.
This is to look for extended sequence diversity.

```bash
# setup directories
mkdir 00a_log 01a_ebi_blast 01b_ebi_dat 01c_ebi_fasta
```

### a. search uniprot database

This step returns UniProt IDs for blast sequence matches.

```bash
sbatch --export fasta=lnu_SeedSeqs.faa /ROCkIn/01b_Sbatch/01a_ebi_blast.sbatch

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
cat 01c_ebi_fasta/*.fasta >> 01d_lnu_all_ebi_matches.fa

# count 'em
grep -c '>' 01d_lnu_all_ebi_matches.fa
```

**Results:**

- 5,615 fasta sequences returned from blast search

# Step 02: Deduplicate, filter, dereplicate

### a. Deduplicate

Since we have multiple verified sequences that we searched to UniProt, we likely have found overlapping search results. Concatenate fasta files of all search result sequences and deduplicate by UniProt ID. I wrote a Python script to deduplicate the concatenated fasta

```bash
python /ROCkIn/02_Python/02a_Remove_Duplicate_Fasta_Entries.py -f 01d_lnu_all_ebi_matches.fa -o 02a_lnu_matches_dedup.fa
```

**Results:**

- Total sequences in file: 5615
- Duplicates Removed: 3936
- Unique sequences retained: 1679

### b. Filter

The initial sequence search does not have great options to filter the matches that are returned. As such, it usually returns many spurious matches.

To filter our search results, I run Blastp locally using the curated sequence as the reference database and the search results as the query.

I wrote a Python script that does 2 filters:
    1) Removes matches <= 30% identity to verified sequences
    2) Remove matches with <= 50% sequence alignment (alignment length / query sequence length).

This script outputs 100% sequence matches separately as well in tabular blast, fasta, and in a text list of Reference sequences with their corresponding list of uniprot IDs with 100% sequence identity. Ideally, we find at least 1 uniprot ID for each reference sequence. If not, we can look through the filtered blast output and find the closest matches to our reference sequences available in the uniprot database. This important because we'll want to include 1 exact match (or closest match) uniprot ID in our final input list to ROCkOut. During the similar sequence clustering these ID's can sometimes be dereplicated as they aren't guaranteed selection as a cluster representative.

The script also plots histograms of each parameter post filtering and can be rerun with different filtering options if neccessary.

I wrapped it into an sbatch script. This step produces the 02b_filter output directory containing the filtered tabular BLAST output file 02b_lnu_matches_fltrdBstHts.fa which is used downstream. It also contains several diagnostic histograms which provide a visual representation of the sequence search results and can be used to modify the filtering settings if necessary.

The 02b_filter output directory also contains the file 02b_lnu_matches_pID100_list.txt. This file is used to identify UniProt IDs that have 100% Sequence Similarity with the SeedSeqs. The input to ROCkOut is UniProt IDs but we retrieved our SeedSeqs from NCBI. We use this file to select 1 UniProt ID for each SeedSeq to include in the *Training set* set list.

```bash
sbatch --export ref=lnu_SeedSeqs.faa,qry=02a_lnu_matches_dedup.fa,out=02b_lnu_matches /ROCkIn/01b_Sbatch/02b_Blastp.sbatch

cat 00a_log/02b_BlastP.out
```

**Results:**

- Total number of entries in blast file: 11457
- Number of entries failing the filters: 8781
- Number of entries passing the filters: 2676
- Number of duplicate blast matches passing filter to remove: 2092
- Number of best hit entries written to new file: 584

![Diagnostic histograms of sequence search results](https://github.com/rotheconrad/ROCker_Macrolide_Models/blob/main/04_lnu/00_figures/02a_lnu_diagnostic_histograms.png)

### c. Dereplicate - similar sequence clustering with MMSeqs2

5427 sequences is way more than we need. Let's dereplicate the set some.

Use mmseqs to generate sequence clusters of 90% amino acid sequence similarity and select cluster representatives. These representative sequences are used for the *Training set* set sequences to build the ROCker models with ROCkOut.

```bash
sbatch --export infile=02b_lnu_matches_fltrdBstHts.fa,oDir=02c_dereplicate_90,ss=90 ../ROCkIn/01b_Sbatch/02c_mmseqs.sbatch

grep -c '>' 02c_dereplicate_90/02_mmseqs_reps.fasta
```

**Results:**

- Representative sequences retained: 269

### d. Select secondary cluster representatives to create a *Testing set*

mmseqs outputs a fasta file of cluster representatives which we will make our selection from for model training. We also want a secondary set of sequences/IDs to use for model testing. This step selects those secondary references (secReps). Additionally, some genomes may have more than one gene with similar sequence to the target gene function. This script also selects any other genes from our pool of sequence matches that belong to a representative sequence genome or a secondary sequence genome and includes them as well. This prevents a hassel downstream when building the model and discovering non_target sequences overlapping with the target sequences.

By default a random seed is not set for this script so the secondary representative selection is chosen randomly from each cluster each time the script is rerun (excluding the primary cluster representative selected by mmseqs2). use the optional -s parameter to set a fixed seed for reproducible "random" selections.

```bash
python ../ROCkIn/02_Python/02d_Get_Test_secReps.py -f 02b_lnu_matches_fltrdBstHts.fa -c 02c_dereplicate_90/02c_mmseqs_results_90.tsv -o 02f_mmseq90

mkdir 02d_secReps_testSets

mv *_secReps.fa *_IDlist.txt 02d_secReps_testSets/

grep -c '>' 02d_secReps_testSets/02f_mmseq90_secReps.fa 
```

**Results:**

- Secondary representative sequences retained: 49

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

sbatch --export verified=lnu_SeedSeqs.faa,newseqs=02c_dereplicate_90/02c_mmseqs_reps_90.fasta /ROCkIn/01b_Sbatch/03a_seq_alignment.sbatch

# sequences before trimming
grep -c '>' 03_ROCkIn_Results/03a_mmseqs_reps_90.fasta.aln

# testing set

sbatch --export verified=lnu_SeedSeqs.faa,newseqs=02d_secReps_testSets/02f_mmseq90_secReps.fa /ROCkIn/01b_Sbatch/03a_seq_alignment.sbatch

# sequences before trimming
grep -c '>' 03_ROCkIn_Results/03a_mmseq90_secReps.fa.aln
```

**Results:**

- Training set sequences before trimming: 275 (269 searched + 6 curated)
- Testing set sequences before trimming: 55 (49 searched + 6 curated)

### b. trim

trim the alignment with trimmal. The trimming can also be completed manually but I prefer trimal for removing columns and general clean up. I remove sequences with obviously bad alignments manually.

Before building a phylogenetic tree, MSAs should be cleaned to removed spurious sequences and columns consisting mainly of gaps.

Before using trimmal, it is good to review the MSA and remove (or at least flag) any obvious sequences that don't belong such as 1 or a few really long sequences or 1 or a few sequences that induce large gaps. Trimmal is pretty good at catching these and removing them, but it doesn't always remove them - sometimes it just trims them and then it is more difficult to see that they don't fit in well with the other sequences. This won't affect the tree building so much as it will show up in the ROCker model.

#### Reviewing the MSA, the following sequences are all too long.

- check beginning and end of alignment for outliers and consider removing them.
- check large gaps in the alignment to see if only 1 or a few sequences are responsible for creating the gap and consider removing them.

##### training set Removed the following prior to trimming: too long, too short, or bad alignment
- A0A3A5WST9_9BACE//Nucleotidyltransferase domain-containing protein//Bacteroides sp. AF34-31BH//Unreviewed//Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;Bacteroidaceae;Bacteroides
- A0A5P2ET54_9CLOT//Nucleotidyltransferase domain-containing protein//Clostridium diolis//Unreviewed//Bacteria;Bacillota;Clostridia;Eubacteriales;Clostridiaceae;Clostridium
- A0A1G0WS67_9BACT//Nucleotidyltransferase//Ignavibacteria bacterium RIFOXYC2_FULL_35_21//Unreviewed//Bacteria;Ignavibacteriota;Ignavibacteria
- A0A832DL39_9BACT//Nucleotidyltransferase domain-containing protein//Ignavibacterium album//Unreviewed//Bacteria;Ignavibacteriota;Ignavibacteria;Ignavibacteriales;Ignavibacteriaceae;Ignavibacterium
- A0A349G9H8_9FIRM//Nucleotidyltransferase//Clostridiales bacterium UBA8960//Unreviewed//Bacteria;Bacillota;Clostridia;Eubacteriales;Eubacteriales;Family;XII;Incertae;Sedis
- A0A8T4J2P5_9ACTN//Nucleotidyltransferase domain-containing protein//Streptomyces daliensis//Unreviewed//Bacteria;Actinomycetota;Actinomycetes;Kitasatosporales;Streptomycetaceae;Streptomyces
- A0A0E9LST6_9BACT//Predicted nucleotidyltransferases//Geofilum rubicundum JCM 15548//Unreviewed//Bacteria;Bacteroidota;Bacteroidia;Marinilabiliales;Marinilabiliaceae;Geofilum
- A0A350V2R4_UNCBA//Nucleotidyltransferase//Bacteroidota bacterium//Unreviewed//Bacteria;Bacteroidota
- E1YM00_9BACT//Polymerase beta nucleotidyltransferase domain-containing protein//uncultured Desulfobacterium sp//Unreviewed//Bacteria;Thermodesulfobacteriota;Desulfobacteria;Desulfobacterales;Desulfobacteriaceae;Desulfobacterium;environmental;samples
- A0A970V553_9BACT//Nucleotidyltransferase domain-containing protein//candidate division WS1 bacterium//Unreviewed//Bacteria;candidate;division;WS1
- X1ARS9_9ZZZZ//Polymerase nucleotidyl transferase domain-containing protein//marine sediment metagenome//Unreviewed//unclassified;sequences;metagenomes;ecological;metagenomes
- A0A662W4L5_9EURY//Nucleotidyltransferase//Archaeoglobales archaeon//Unreviewed//Archaea;Euryarchaeota;Archaeoglobi;Archaeoglobales
- A0A1G1HVE5_9BACT//Polymerase nucleotidyl transferase domain-containing protein//Nitrospirae bacterium RBG_16_64_22//Unreviewed//Bacteria;Nitrospirota
- A0A328CAL1_9DELT//Uncharacterized protein//Lujinxingia litoralis//Unreviewed//Bacteria;Deltaproteobacteria;Bradymonadales;Lujinxingiaceae;Lujinxingia
- A0A2H0AI31_9DELT//Nucleotidyltransferase//Deltaproteobacteria bacterium CG23_combo_of_CG06-09_8_20_14_all_51_20//Unreviewed//Bacteria;Deltaproteobacteria
- A0A9N7JZP7_TETHN//Lincosamide nucleotidyltransferase//NBRC 12172) (Pediococcus halophilus)//Unreviewed//Bacteria;Bacillota;Bacilli;Lactobacillales;Enterococcaceae;Tetragenococcus
- A0A7Z7PM48_9BACT//Nudix hydrolase domain-containing protein//Mesotoga infera//Unreviewed//Bacteria;Thermotogota;Thermotogae;Kosmotogales;Kosmotogaceae;Mesotoga
- A0A1C4E2C9_9BACI//Polymerase beta nucleotidyltransferase domain-containing protein//Bacillus wiedmannii//Unreviewed//Bacteria;Bacillota;Bacilli;Bacillales;Bacillaceae;Bacillus;Bacillus;cereus;group
- A0A7L5DI87_9BACT//Nucleotidyltransferase domain-containing protein//Spirosoma rhododendri//Unreviewed//Bacteria;Bacteroidota;Cytophagia;Cytophagales;Cytophagaceae;Spirosoma
- A0A4Y1ZH60_9BACL//Uncharacterized protein//Sporolactobacillus inulinus//Unreviewed//Bacteria;Bacillota;Bacilli;Bacillales;Sporolactobacillaceae;Sporolactobacillus
- A0A3C2D8E3_9BACT//Polymerase nucleotidyl transferase domain-containing protein//Prolixibacteraceae bacterium//Unreviewed//Bacteria;Bacteroidota;Bacteroidia;Marinilabiliales;Prolixibacteraceae
- A0A661AXG8_UNCXX//Nucleotidyltransferase domain-containing protein//bacterium//Unreviewed//Bacteria
- A0A1S6J291_9ACTN//Uncharacterized protein//Streptomyces pactum//Unreviewed//Bacteria;Actinomycetota;Actinomycetes;Kitasatosporales;Streptomycetaceae;Streptomyces
- A0A2S2FU38_9ACTN//Nucleotidyltransferase domain-containing protein//Streptomyces sp. SM17//Unreviewed//Bacteria;Actinomycetota;Actinomycetes;Kitasatosporales;Streptomycetaceae;Streptomyces
- A0A4R0ZE44_9BACL//Nucleotidyltransferase domain-containing protein//Exiguobacterium sp. SH5S13//Unreviewed//Bacteria;Bacillota;Bacilli;Bacillales;Bacillales;Family;XII;Incertae;Sedis;Exiguobacterium
- A0A7Z9WFY3_UNCSY//Nucleotidyltransferase domain-containing protein//Syntrophaceae bacterium//Unreviewed//Bacteria;Thermodesulfobacteriota;Syntrophia;Syntrophales;Syntrophaceae
- A0A5M8QNQ9_9BACT//Nucleotidyltransferase domain-containing protein//Dyadobacter flavalbus//Unreviewed//Bacteria;Bacteroidota;Cytophagia;Cytophagales;Spirosomaceae;Dyadobacter
- A0A4Y1ZGK0_9BACL//Polymerase beta nucleotidyltransferase domain-containing protein//Sporolactobacillus inulinus//Unreviewed//Bacteria;Bacillota;Bacilli;Bacillales;Sporolactobacillaceae;Sporolactobacillus
- A0A2H6CW78_TETHA//Putative lincosamide nucleotidyltransferase//Tetragenococcus halophilus subsp. halophilus//Unreviewed//Bacteria;Bacillota;Bacilli;Lactobacillales;Enterococcaceae;Tetragenococcus
- A0A937WPD0_UNCPO//Nucleotidyltransferase domain-containing protein//Poribacteria bacterium//Unreviewed//Bacteria;Candidatus;Poribacteria
- A0A497NIV7_9ARCH//Nucleotidyltransferase domain-containing protein//Candidatus Bathyarchaeota archaeon//Unreviewed//Archaea;Candidatus;Bathyarchaeota
- A0A958CG54_9CHLR//Nucleotidyltransferase domain-containing protein//Anaerolineae bacterium//Unreviewed//Bacteria;Chloroflexota;Anaerolineae
- A0A971VH29_9BACT//Nucleotidyltransferase domain-containing protein//Bacteroidales bacterium//Unreviewed//Bacteria;Bacteroidota;Bacteroidia;Bacteroidales
- A0A101I8I5_9BACT//LinF//Marinimicrobia bacterium 46_43//Unreviewed//Bacteria;Candidatus;Marinimicrobia
- A0A101JMS4_9ACTN//DNA polymerase subunit beta//Actinoplanes awajinensis subsp. mycoplanecinus//Unreviewed//Bacteria;Actinomycetota;Actinomycetes;Micromonosporales;Micromonosporaceae;Actinoplanes
- A0A972S7A6_9BACT//Polymerase beta nucleotidyltransferase domain-containing protein//Campylobacterota bacterium//Unreviewed//Bacteria;Campylobacterota
- A0A7X3P1G3_UNCPO//Nucleotidyltransferase domain-containing protein//Poribacteria bacterium//Unreviewed//Bacteria;Candidatus;Poribacteria
- A0A5J4L2W2_9ZZZZ//Polymerase beta nucleotidyltransferase domain-containing protein//hot springs metagenome//Unreviewed//unclassified;sequences;metagenomes;ecological;metagenomes
- A0A3A0A7Q1_UNCCH//Nucleotidyltransferase domain-containing protein//Chloroflexota bacterium//Unreviewed//Bacteria;Chloroflexota
- A0A7V2RQ93_UNCBA//DinB family protein//Bacteroidota bacterium//Unreviewed//Bacteria;Bacteroidota
- A0A661UD38_UNCCO//Polymerase beta nucleotidyltransferase domain-containing protein//Coatesbacteria bacterium//Unreviewed//Bacteria;Candidatus;Coatesbacteria
- A0A1W1WTW9_9BACT//Nucleotidyltransferase domain-containing protein//Nitratiruptor tergarcus DSM 16512//Unreviewed//Bacteria;Campylobacterota;Epsilonproteobacteria;Nautiliales;Nitratiruptoraceae;Nitratiruptor
- A0A1M2ZI67_9SPHN//Uncharacterized protein//Sphingomonas sp. 66-10//Unreviewed//Bacteria;Pseudomonadota;Alphaproteobacteria;Sphingomonadales;Sphingomonadaceae;Sphingomonas
- A0A419F9G8_9BACT//Nucleotidyltransferase domain-containing protein//Candidatus Abyssubacteria bacterium SURF_17//Unreviewed//Bacteria;Candidatus;Abyssubacteria
- A0A1F9FX31_9DELT//Polymerase beta nucleotidyltransferase domain-containing protein//Deltaproteobacteria bacterium RIFCSPHIGHO2_02_FULL_44_16//Unreviewed//Bacteria;Deltaproteobacteria
- A0A6N7MMI2_UNCBA//Nucleotidyltransferase domain-containing protein//Bacteroidota bacterium//Unreviewed//Bacteria;Bacteroidota
- A0A1I5RY18_9BACT//Nucleotidyltransferase domain-containing protein//Pseudarcicella hirudinis//Unreviewed//Bacteria;Bacteroidota;Cytophagia;Cytophagales;Spirosomaceae;Pseudarcicella
- A0A660YE10_UNCLA//Nucleotidyltransferase domain-containing protein//Latescibacteria bacterium//Unreviewed//Bacteria;Candidatus;Latescibacteria
- A0A5A7RT87_UNCTM//Nucleotidyltransferase domain-containing protein//Thermoplasmata archaeon//Unreviewed//Archaea;Candidatus;Thermoplasmatota;Thermoplasmata
- A0A926BY26_9BACT//Nucleotidyltransferase domain-containing protein//Fibrella sp//Unreviewed//Bacteria;Bacteroidota;Cytophagia;Cytophagales;Spirosomaceae;Fibrella
- A0A350LYM3_UNCXX//Nucleotidyltransferase//bacterium//Unreviewed//Bacteria
- A0A9D6YFG9_9BACT//Nucleotidyltransferase domain-containing protein//Candidatus Rokubacteria bacterium//Unreviewed//Bacteria;Candidatus;Rokubacteria
- A0A957QLW1_9CHLR//Nucleotidyltransferase domain-containing protein//Anaerolineales bacterium//Unreviewed//Bacteria;Chloroflexota;Anaerolineae;Anaerolineales
- A0A7W1CZ29_9ACTN//Uncharacterized protein//Rubrobacteraceae bacterium//Unreviewed//Bacteria;Actinomycetota;Rubrobacteria;Rubrobacterales;Rubrobacteraceae
- A0A3B9ZRT1_9BACT//Polymerase nucleotidyl transferase domain-containing protein//Prolixibacteraceae bacterium//Unreviewed//Bacteria;Bacteroidota;Bacteroidia;Marinilabiliales;Prolixibacteraceae
- A0A5C9EDU4_9ARCH//Putative Nucleotidyltransferase domain protein//Candidatus Lokiarchaeota archaeon//Unreviewed//Archaea;Asgard;group;Candidatus;Lokiarchaeota
- A0A2N2JSL5_9DELT//Polymerase beta nucleotidyltransferase domain-containing protein//Deltaproteobacteria bacterium HGW-Deltaproteobacteria-11//Unreviewed//Bacteria;Deltaproteobacteria
- A0A522UQF2_9EURY//Nucleotidyltransferase domain-containing protein//Candidatus Methanoperedens sp//Unreviewed//Archaea;Euryarchaeota;Stenosarchaea;group;Methanomicrobia;Methanosarcinales;Candidatus;Methanoperedenaceae;Candidatus;Methanoperedens
- A0A1P9X1P8_9BACT//Polymerase nucleotidyl transferase domain-containing protein//Spirosoma montaniterrae//Unreviewed//Bacteria;Bacteroidota;Cytophagia;Cytophagales;Cytophagaceae;Spirosoma
- A0A925R8W8_9FIRM//Nucleotidyltransferase domain-containing protein//Vallitaleaceae bacterium//Unreviewed//Bacteria;Bacillota;Clostridia;Eubacteriales;Vallitaleaceae
- A0A537J182_9BACT//Nucleotidyltransferase domain-containing protein//Terrabacteria group bacterium ANGP1//Unreviewed//Bacteria
- A0A959S175_9BACT//Nucleotidyltransferase domain-containing protein//Ignavibacteriota bacterium//Unreviewed//Bacteria;Ignavibacteriota
- A0A925NW77_9BACT//Nucleotidyltransferase domain-containing protein//Prolixibacteraceae bacterium//Unreviewed//Bacteria;Bacteroidota;Bacteroidia;Marinilabiliales;Prolixibacteraceae
- A0A923UL60_9BACT//Nucleotidyltransferase domain-containing protein//Arcicella sp//Unreviewed//Bacteria;Bacteroidota;Cytophagia;Cytophagales;Spirosomaceae;Arcicella
- A0A202DFR6_UNCAD//Polymerase beta nucleotidyltransferase domain-containing protein//archaeon D22//Unreviewed//Archaea
- A0A0S8GVD3_UNCZI//DNA polymerase//candidate division Zixibacteria bacterium SM23_73//Unreviewed//Bacteria
- A0A519TH72_9BACT//Nucleotidyltransferase domain-containing protein//Hymenobacter sp//Unreviewed//Bacteria;Bacteroidota;Cytophagia;Cytophagales;Hymenobacteraceae;Hymenobacter
- A0A7C1GPN9_9BACT//NUDIX domain-containing protein//Mesotoga infera//Unreviewed//Bacteria;Thermotogota;Thermotogae;Kosmotogales;Kosmotogaceae;Mesotoga
- A0A357BP05_UNCDE//Nucleotidyltransferase//Deltaproteobacteria bacterium//Unreviewed//Bacteria;Deltaproteobacteria
- A0A0W8FMW7_9ZZZZ//Putative nucleotidyltransferase//hydrocarbon metagenome//Unreviewed//unclassified;sequences;metagenomes;ecological;metagenomes
- A0A2W2CW86_9ACTN//Polymerase nucleotidyl transferase domain-containing protein//Jiangella anatolica//Unreviewed//Bacteria;Actinomycetota;Actinomycetes;Jiangellales;Jiangellaceae;Jiangella
- A0A1H6SQE5_9BACT//Nucleotidyltransferase domain-containing protein//Dyadobacter sp. SG02//Unreviewed//Bacteria;Bacteroidota;Cytophagia;Cytophagales;Spirosomaceae;Dyadobacter
- A0A1G1DXP6_9BACT//Polymerase beta nucleotidyltransferase domain-containing protein//Nitrospinae bacterium RIFCSPLOWO2_12_39_16//Unreviewed//Bacteria;Nitrospinae/Tectomicrobia;group;Nitrospinota
- A0A1I5WL75_9FIRM//Nucleotidyltransferase domain-containing protein//Caldicoprobacter faecalis//Unreviewed//Bacteria;Bacillota;Clostridia;Eubacteriales;Caldicoprobacteraceae;Caldicoprobacter
- A0A521VJ66_9BACT//Nucleotidyltransferase domain-containing protein//Saprospiraceae bacterium//Unreviewed//Bacteria;Bacteroidota;Saprospiria;Saprospirales;Saprospiraceae
- A0A1I2BND6_9BACT//Nucleotidyltransferase domain-containing protein//Thermoflexibacter ruber//Unreviewed//Bacteria;Bacteroidota;Cytophagia;Cytophagales;Thermoflexibacteraceae;Thermoflexibacter
- A0A7V9VD77_9CHLR (BAD ALIGNMENT - LONG BRANCH)
- A0A969FB27_9CHLR (BAD ALIGNMENT - LONG BRANCH)

##### testing set Removed the following prior to trimming: too long too short, or bad alignment
- A0A971TZX4_9THEM//NUDIX domain-containing protein//Thermotogaceae bacterium//Unreviewed//Bacteria;Thermotogota;Thermotogae;Thermotogales;Thermotogaceae
- A0A124FYT1_9BACT//LinF lincosamide nucleotidyltransferase LinF//Mesotoga prima//Unreviewed//Bacteria;Thermotogota;Thermotogae;Kosmotogales;Kosmotogace
- A0A174JHG6_BACUN//Nucleotidyltransferase domain-containing protein//Bacteroides uniformis//Unreviewed//Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;Bacteroidaceae;Bacteroides
- A0A1J5EGN9_9BACT//Nucleotidyltransferase//Desulfobacteraceae bacterium CG2_30_51_40//Unreviewed//Bacteria;Thermodesulfobacteriota;Desulfobacteria;Desulfobacterales;Desulfobacteraceae
- A0A957G6F0_9CHLR//Nucleotidyltransferase domain-containing protein//Anaerolineales bacterium//Unreviewed//Bacteria;Chloroflexota;Anaerolineae;Anaerolineales
- A0A640WB02_UNCTM//Nucleotidyltransferase domain-containing protein//Thermoplasmata archaeon//Unreviewed//Archaea;Candidatus;Thermoplasmatota;Thermoplasmata
- A0A418QE28_9DEIO//Uncharacterized protein//Deinococcus sp. RM//Unreviewed//Bacteria;Deinococcota;Deinococci;Deinococcales;Deinococcaceae;Deinococcus
- A0A7G1GY40_9BACT//Polymerase beta nucleotidyltransferase domain-containing protein//Dissulfurispira thermophila//Unreviewed//Bacteria;Nitrospirota;Thermodesulfovibrionia;Thermodesulfovibrionales;Dissulfurispiraceae;Dissulfurispira
- A0A0K2MI31_CLOBE//Nucleotidyltransferase domain-containing protein//Clostridium beijerinckii NRRL B-598//Unreviewed//Bacteria;Bacillota;Clostridia;Eubacteriales;Clostridiaceae;Clostridium

*Removed the fasta sequences from 02c_mmseqs_reps_90.fasta (training set) and 02f_mmseq90_secReps.fa (testing set) and reran the alignment in step A after removing the above sequences from the fasta files and then trimmed after realignment.*

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

- Training set sequences after trimming: 180
- Testing set sequences after trimming: 50

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

##### Found two long branches in the tree results
- A0A7V9VD77_9CHLR
- A0A969FB27_9CHLR

Removed from the trimmed alignment and rebuilt tree.

#### Training set - clade/cluster labeled tree

The below image is the default IQTree and clade/cluster labeling from the 03c_Plot_Annotated_Tree_v2.py script. This script uses the HDBSCAN algorithm to guess at gene clades/clusters. An all vs all distance matrix of tree branch lengths is computed and used as input to HDBSCAN. This labeling is intended as a predicted starting point and is used to label and order the genes in the corresponding output 03h_Gene_Data_90_annotated.tsv data file. The researcher can change the labels in the data file and rerun the 03c plotting script as needed.

![Training set tree](https://github.com/rotheconrad/ROCker_Macrolide_Models/blob/main/04_lnu/00_figures/03a_lnu_training_set_tree.png)

#### Testing set - clade/cluster labeled tree

The below image is the default IQTree and clade/cluster labeling from the 03c_Plot_Annotated_Tree_v2.py script. This script uses the HDBSCAN algorithm to guess at gene clades/clusters. An all vs all distance matrix of tree branch lengths is computed and used as input to HDBSCAN. This labeling is intended as a predicted starting point and is used to label and order the genes in the corresponding output 03h_Gene_Data_90_annotated.tsv data file. The researcher can change the labels in the data file and rerun the 03c plotting script as needed.

![Testing set tree](https://github.com/rotheconrad/ROCker_Macrolide_Models/blob/main/04_lnu/00_figures/03b_lnu_testing_set_tree.png)

Review the Data Table and the Phylogenetic tree to make positive and negative UniProt ID selections to give to ROCkOut. Clade/Cluster labels can be changed in the tsv file and the tree can be replotted.

- for the test set be sure to remove all the RepSeqs from the UniProt ID list.
- for the training set remove any Gene Names that aren't a UniProt ID.
- for the training set, look at 02b_filter/02b_lnu_matches_pID100.fa and include some sequences that 100% similar to any RepSeqs without a UniProt ID. You only need 1 per RepSeq. This will include some different genomes in the training set that have the identical copies.

**Results**

We decided to build three models:

1. A model capturing all lnu genes (F, B, G, H) used in the *Training set* tree as positive references to make a model that captures all of these lnu genes. This model also uses a negative references from the polymerase clades x-x in the default chart. This model is useful for broad surveillance of lnu F, B, G, H gene functions in microbial communities. We refer to this model as lnuAll.
1. A model focusing specifically on the lnuF gene clade of the *Training set* using only the sequences from this clade as positive references and using all other sequences as negative references. We refer to this model as lnuF.
1. A model focusing specifically on the lnuB/G gene clade of the *Training set* using only the sequences from this clade as positive references and using all other sequences as negative references. We refer to this model as lnuB/G.

We make three similar similar selections for the *Testing set*:

1. We mirrored the gene selections from the lnuAll model with the *Testing set* to use as positive and negative references. The genomes containing these genes will be used to create a mock metagenome to challenge the lnuAll ROCker model.
1. We select genes in the *Testing set* specifically from the lnuF gene clade. The genomes containing these genes will be used to create a mock metagenome to challenge the lnuF ROCker model.
1. We select genes in the *Testing set* specifically from the lnuB/G gene clade. The genomes containing these genes will be used to create a mock metagenome to challenge the lnuF ROCker model.

# Step 04: Build ROCker models

This is an iterative process of building a model, investigating the results, and tracking down all the peculiarities to arrive at a final postive and/or negative UniProt ID set.

The lnu genes share some sequence similarity with dna polymerase sub unit beta.
ROCkIn picks up a large clade of these genes that apppear distinct from the
lnu genes. Using ROCkIn, we decided to use this clade as negative references for the lnu all model.

##### notes on building the lnuAll model:

In the training set model there was one non_target gene that generated reads (orange dots in model figure) right at the bitscore cutoff line.

- CP047640.1 - Proteus sp. ZN5 plasmid pZN5-103kb - this plasmid has two lnu genes annotated as one lnuF and one lnuG that have 34.85% sequence identity. The lnuF gene "A0A6G6TM98_9GAMM" was included in the positive reference set intitially from ROCkIn and the the lnuG "A0A6G6TMC0_9GAMM" has been added to the positive reference set in this refinement round. Both the lnuG and lnuF genes have 100% or very high sequence identity to genes in genomes of many other species. These genes or this plasmid gets around.

In the testing set model there was one non_target gene that generated reads (orange dots in model figure) in sort of a double v pattern crossing the bitscore threshold.

- MW940622.1 - Escherichia sp. strain SH9W plasmid pSH9W-tetX4 - this plasmid seems to contain a degraded psuedo-gene that is split into two fragments by gene prediction. The "A0A8E6P1C7_9ESCH" gene is the first half of an lnuF gene and was intitially picked up by ROCkIn and included in the positive reference set. The "A0A8E6UEX3_9ESCH" gene is the second half of the fragmented lnuF gene. Both gene fragments share 100% sequence identity with numerous genes from genomes of several other species. We added the second gene fragment UniProt ID to the positive reference set in this refinement round, but there is still a third fragment that maps to the end of the gene without a UniProt ID so I removed these IDs.

##### notes on building the lnuF model:

For this model we use only the lnuF clade as positive references and all other
clades as negative references.

The testing set for the lnuF model had the same gene fragment issue as the lnuAll testing model with the lnuF gene on the plasmid:

- MW940622.1 - Escherichia sp. strain SH9W plasmid pSH9W-tetX4 - this plasmid seems to contain a degraded psuedo-gene that is split into two fragments by gene prediction. The "A0A8E6P1C7_9ESCH" gene is the first half of an lnuF gene and was intitially picked up by ROCkIn and included in the positive reference set. The "A0A8E6UEX3_9ESCH" gene is the second half of the fragmented lnuF gene. Both gene fragments share 100% sequence identity with numerous genes from genomes of several other species. We added the second gene fragment UniProt ID to the positive reference set in this refinement round.

```bash
# lnuAll model
mkdir 01a_lnuAll_train 01b_lnuAll_test 01c_alignment_label_testing

# training set
cd 01a_lnuAll_train
sbatch --export odir=model,pos=lnuAll_train_positive.txt,neg=lnuAll_train_negative.txt /ROCkOut/00b_sbatch/ROCkOut.sbatch

# testing set
cd 01b_lnuAll_test
sbatch --export odir=model,pos=lnuAll_test_positive.txt,pos=lnuAll_test_negative.txt /ROCkOut/00b_sbatch/ROCkOut.sbatch

######

# lnuF model

mkdir 02a_lnuF_train 02b_lnuF_test 02c_alignment_label_testing

# training set
cd 02a_lnuF_train
sbatch --export odir=model,pos=lnuF_train_pos.txt,neg=lnuF_train_neg.txt /ROCkOut/00b_sbatch/ROCkOut.sbatch

# testing set
cd 02b_lnuF_test
sbatch --export odir=model,pos=lnuF_test_pos.txt,neg=lnuF_test_neg.txt /ROCkOut/00b_sbatch/ROCkOut.sbatch

######

# lnuB/G model

mkdir 04a_lnuBG_train 04_lnuBG_test

# upload lnuBG_train_pos.txt and lnuBG_train_neg.txt to 04a_lnuBG_train
# upload lnuBG_test_pos.txt and lnuBG_test_neg.txt to 04_lnuBG_test

## training set
cd 04a_lnuBG_train
sbatch --export odir=model,pos=lnuBG_train_pos.txt,neg=lnuBG_train_neg.txt ../../00b_sbatch/ROCkOut.sbatch

## testing set
cd 04b_lnuBG_test
sbatch --export odir=model,pos=lnuBG_test_pos.txt,neg=lnuBG_test_neg.txt ../../00b_sbatch/ROCkOut.sbatch
```

![mph ROCker Models](https://github.com/rotheconrad/ROCker_Macrolide_Models/blob/main/04_lnu/00_figures/04a-lnu-ROCker-model-250bp.png)

Top panel (A) is the ROCker model for the lnu all gene set. Middle panel (B) is the ROCker model for the lnuF gene set. Bottom panel (C) is the ROCker model for the lnuB/G gene set.

# Step 05: Test ROCker models

ROCkOut as an "-extract" feature that collects the labeled, simulated, short reads into fasta files. We will use these fasta files from the *Testing set* to execute ROCkOut's align, filter, and place functions against the *Training set* models.

We build models for the *Testing set* because 1) this is a convenient way to use the ROCkOut figures and interactive plots to fine tune the sequence selection and 2) ROCkOut downloads genomes the genes are part of, simulates short sequence reads for the genomes, and labels the reads and positive, negative, target, or non_target. We can collect these reads into a mock metagenome that is labeled so we can track them and use them to score the results.

We test 4 different read lengths for each model.

Scripts for this section are in the 00_scripts directory of this GitHub Repo. Scripts for all previous sections were from either the ROCkIn or ROCkOut repos.

### a. Run mock metagenomes against models using ROCkOut

```bash
# lnuAll model

cd 01a_lnuAll_train

sbatch --export infile=../01b_lnuAll_test/simulated_reads/simread_100_raw_reads.fasta,model=model,outpre=test_score_100 ../../00b_sbatch/ROCkOut_posOnly_score.sbatch

sbatch --export infile=../01b_lnuAll_test/simulated_reads/simread_150_raw_reads.fasta,model=model,outpre=test_score_150 ../../00b_sbatch/ROCkOut_posOnly_score.sbatch

sbatch --export infile=../01b_lnuAll_test/simulated_reads/simread_250_raw_reads.fasta,model=model,outpre=test_score_250 ../../00b_sbatch/ROCkOut_posOnly_score.sbatch

sbatch --export infile=../01b_lnuAll_test/simulated_reads/simread_300_raw_reads.fasta,model=model,outpre=test_score_300 ../../00b_sbatch/ROCkOut_posOnly_score.sbatch

######

# lnuF model

cd 02a_lnuF_train

sbatch --export infile=../02b_lnuF_test/simulated_reads/simread_100_raw_reads.fasta,model=model,outpre=test_score_100 ../../00b_sbatch/ROCkOut_score.sbatch

sbatch --export infile=../02b_lnuF_test/simulated_reads/simread_150_raw_reads.fasta,model=model,outpre=test_score_150 ../../00b_sbatch/ROCkOut_score.sbatch

sbatch --export infile=../02b_lnuF_test/simulated_reads/simread_250_raw_reads.fasta,model=model,outpre=test_score_250 ../../00b_sbatch/ROCkOut_score.sbatch

sbatch --export infile=../02b_lnuF_test/simulated_reads/simread_300_raw_reads.fasta,model=model,outpre=test_score_300 ../../00b_sbatch/ROCkOut_score.sbatch

######

# lnuB/G model

cd 04a_lnuBG_train

sbatch --export infile=../04b_lnuBG_test/simulated_reads/simread_100_raw_reads.fasta,model=model,outpre=test_score_100 ../../00b_sbatch/ROCkOut_score.sbatch

sbatch --export infile=../04b_lnuBG_test/simulated_reads/simread_150_raw_reads.fasta,model=model,outpre=test_score_150 ../../00b_sbatch/ROCkOut_score.sbatch

sbatch --export infile=../04b_lnuBG_test/simulated_reads/simread_250_raw_reads.fasta,model=model,outpre=test_score_250 ../../00b_sbatch/ROCkOut_score.sbatch

sbatch --export infile=../04b_lnuBG_test/simulated_reads/simread_300_raw_reads.fasta,model=model,outpre=test_score_300 ../../00b_sbatch/ROCkOut_score.sbatch
```

### b. Create and search custom HMM model

Use [HMMER](http://hmmer.org/) to build a model from the multiple alignment created by ROCkOut. This multiple alignment is the same genes and alignment ROCkOut uses to create it's model, so we use it to create an HMM model for comparison.

We convert the mock metagenome (nucleotide fasta) to all 6 reading frames of amino acid sequence because this is what BLASTx (or diamond BLASTx) does internal to align reads. In this way the results our HMM search will be comparable to BLASTx. We do this for all 4 read lengths.

We use [HMMER](http://hmmer.org/) to search our amino acid sequences against the model and we use a custom script to filter results for best match because there is a low rate of more than one reading frame per read finding a match.

```bash
# lnuAll model

# build model (hmmbuild)
hmmbuild lnuAll_model-final.hmm model/shared_files/multiple_alignment/complete_multiple_alignment_aa.fasta

# convert mock metagenome nucleotides to 6 aa frames
python 00_scripts/nuc_fasta_6aa_frame_fasta.py -i ../01b_lnuAll_test/simulated_reads/simread_100_raw_reads.fasta -o test_score_100/simread_100_raw_reads_6frames.faa

python 00_scripts/nuc_fasta_6aa_frame_fasta.py -i ../01b_lnuAll_test/simulated_reads/simread_150_raw_reads.fasta -o test_score_150/simread_150_raw_reads_6frames.faa

python 00_scripts/nuc_fasta_6aa_frame_fasta.py -i ../01b_lnuAll_test/simulated_reads/simread_250_raw_reads.fasta -o test_score_250/simread_250_raw_reads_6frames.faa

python 00_scripts/nuc_fasta_6aa_frame_fasta.py -i ../01b_lnuAll_test/simulated_reads/simread_300_raw_reads.fasta -o test_score_300/simread_300_raw_reads_6frames.faa

# Map Reads (hmmsearch)
sbatch --export model=lnuAll_model-final.hmm,infile=test_score_100/simread_100_raw_reads_6frames.faa,outfile=test_score_100/lnuAll-final_hmmSearch_100.tsv 00_scripts/hmmsearch.sbatch

sbatch --export model=lnuAll_model-final.hmm,infile=test_score_150/simread_150_raw_reads_6frames.faa,outfile=test_score_150/lnuAll-final_hmmSearch_150.tsv 00_scripts/hmmsearch.sbatch

sbatch --export model=lnuAll_model-final.hmm,infile=test_score_250/simread_250_raw_reads_6frames.faa,outfile=test_score_250/lnuAll-final_hmmSearch_250.tsv 00_scripts/hmmsearch.sbatch

sbatch --export model=lnuAll_model-final.hmm,infile=test_score_300/simread_300_raw_reads_6frames.faa,outfile=test_score_300/lnuAll-final_hmmSearch_300.tsv 00_scripts/hmmsearch.sbatch


# Filter for best hit and score
python 00_scripts/besthit_filter_hmm.py -i test_score_100/lnuAll-final_hmmSearch_100.tsv -o test_score_100/lnuAll-final_hmmSearch_100_filtered.tsv

python 00_scripts/besthit_filter_hmm.py -i test_score_150/lnuAll-final_hmmSearch_150.tsv -o test_score_150/lnuAll-final_hmmSearch_150_filtered.tsv

python 00_scripts/besthit_filter_hmm.py -i test_score_250/lnuAll-final_hmmSearch_250.tsv -o test_score_250/lnuAll-final_hmmSearch_250_filtered.tsv

python 00_scripts/besthit_filter_hmm.py -i test_score_300/lnuAll-final_hmmSearch_300.tsv -o test_score_300/lnuAll-final_hmmSearch_300_filtered.tsv

######

# lnuF model

# build model (hmmbuild)
hmmbuild lnuF_model-final.hmm model/shared_files/multiple_alignment/complete_multiple_alignment_aa.fasta

# convert mock metagenome nucleotides to 6 aa frames
python 00_scripts/nuc_fasta_6aa_frame_fasta.py -i ../02b_lnuF_test/simulated_reads/simread_100_raw_reads.fasta -o test_score_100/simread_100_raw_reads_6frames.faa

python 00_scripts/nuc_fasta_6aa_frame_fasta.py -i ../02b_lnuF_test/simulated_reads/simread_150_raw_reads.fasta -o test_score_150/simread_150_raw_reads_6frames.faa

python 00_scripts/nuc_fasta_6aa_frame_fasta.py -i ../02b_lnuF_test/simulated_reads/simread_250_raw_reads.fasta -o test_score_250/simread_250_raw_reads_6frames.faa

python 00_scripts/nuc_fasta_6aa_frame_fasta.py -i ../02b_lnuF_test/simulated_reads/simread_300_raw_reads.fasta -o test_score_300/simread_300_raw_reads_6frames.faa

# Map Reads (hmmsearch)
sbatch --export model=lnuF_model-final.hmm,infile=test_score_100/simread_100_raw_reads_6frames.faa,outfile=test_score_100/lnuF-final_hmmSearch_100.tsv 00_scripts/hmmsearch.sbatch

sbatch --export model=lnuF_model-final.hmm,infile=test_score_150/simread_150_raw_reads_6frames.faa,outfile=test_score_150/lnuF-final_hmmSearch_150.tsv 00_scriptsh/hmmsearch.sbatch

sbatch --export model=lnuF_model-final.hmm,infile=test_score_250/simread_250_raw_reads_6frames.faa,outfile=test_score_250/lnuF-final_hmmSearch_250.tsv 00_scripts/hmmsearch.sbatch

sbatch --export model=lnuF_model-final.hmm,infile=test_score_300/simread_300_raw_reads_6frames.faa,outfile=test_score_300/lnuF-final_hmmSearch_300.tsv 00_scripts/hmmsearch.sbatch

# Filter for best hit and score
python 00_scripts/besthit_filter_hmm.py -i test_score_100/lnuF-final_hmmSearch_100.tsv -o test_score_100/lnuF-final_hmmSearch_100_filtered.tsv

python 00_scripts/besthit_filter_hmm.py -i test_score_150/lnuF-final_hmmSearch_150.tsv -o test_score_150/lnuF-final_hmmSearch_150_filtered.tsv

python 00_scripts/besthit_filter_hmm.py -i test_score_250/lnuF-final_hmmSearch_250.tsv -o test_score_250/lnuF-final_hmmSearch_250_filtered.tsv

python 00_scripts/besthit_filter_hmm.py -i test_score_300/lnuF-final_hmmSearch_300.tsv -o test_score_300/lnuF-final_hmmSearch_300_filtered.tsv

######

# lnuB/G model

# build model (hmmbuild)
hmmbuild lnuBG_model-final.hmm model/shared_files/multiple_alignment/complete_multiple_alignment_aa.fasta

# convert mock metagenome nucleotides to 6 aa frames
python ../../00c_scripts/nuc_fasta_6aa_frame_fasta.py -i ../04b_lnuBG_test/simulated_reads/simread_100_raw_reads.fasta -o test_score_100/simread_100_raw_reads_6frames.faa

python ../../00c_scripts/nuc_fasta_6aa_frame_fasta.py -i ../04b_lnuBG_test/simulated_reads/simread_150_raw_reads.fasta -o test_score_150/simread_150_raw_reads_6frames.faa

python ../../00c_scripts/nuc_fasta_6aa_frame_fasta.py -i ../04b_lnuBG_test/simulated_reads/simread_250_raw_reads.fasta -o test_score_250/simread_250_raw_reads_6frames.faa

python ../../00c_scripts/nuc_fasta_6aa_frame_fasta.py -i ../04b_lnuBG_test/simulated_reads/simread_300_raw_reads.fasta -o test_score_300/simread_300_raw_reads_6frames.faa


# Map Reads (hmmsearch)
sbatch --export model=lnuBG_model-final.hmm,infile=test_score_100/simread_100_raw_reads_6frames.faa,outfile=test_score_100/lnuBG-final_hmmSearch_100.tsv ../../00b_sbatch/hmmsearch.sbatch

sbatch --export model=lnuBG_model-final.hmm,infile=test_score_150/simread_150_raw_reads_6frames.faa,outfile=test_score_150/lnuBG-final_hmmSearch_150.tsv ../../00b_sbatch/hmmsearch.sbatch

sbatch --export model=lnuBG_model-final.hmm,infile=test_score_250/simread_250_raw_reads_6frames.faa,outfile=test_score_250/lnuBG-final_hmmSearch_250.tsv ../../00b_sbatch/hmmsearch.sbatch

sbatch --export model=lnuBG_model-final.hmm,infile=test_score_300/simread_300_raw_reads_6frames.faa,outfile=test_score_300/lnuBG-final_hmmSearch_300.tsv ../../00b_sbatch/hmmsearch.sbatch

# Filter for best hit and score
python ../../00c_scripts/besthit_filter_hmm.py -i test_score_100/lnuBG-final_hmmSearch_100.tsv -o test_score_100/lnuBG-final_hmmSearch_100_filtered.tsv

python ../../00c_scripts/besthit_filter_hmm.py -i test_score_150/lnuBG-final_hmmSearch_150.tsv -o test_score_150/lnuBG-final_hmmSearch_150_filtered.tsv

python ../../00c_scripts/besthit_filter_hmm.py -i test_score_250/lnuBG-final_hmmSearch_250.tsv -o test_score_250/lnuBG-final_hmmSearch_250_filtered.tsv

python ../../00c_scripts/besthit_filter_hmm.py -i test_score_300/lnuBG-final_hmmSearch_300.tsv -o test_score_300/lnuBG-final_hmmSearch_300_filtered.tsv
```

**Results:**

#### lnuAll - 100 bp reads
- Total number of entries in hmm file: 3077
- Number of duplicate hmm matches: 0
- Number of best hit entries written to new file: 3077

#### lnuAll - 150 bp reads
- Total number of entries in hmm file: 2348
- Number of duplicate hmm matches: 1
- Number of best hit entries written to new file: 2347 

#### lnuAll - 250 bp reads
- Total number of entries in hmm file: 1532
- Number of duplicate hmm matches: 0
- Number of best hit entries written to new file: 1532 

#### lnuAll - 300 bp reads
- Total number of entries in hmm file: 1340
- Number of duplicate hmm matches: 0
- Number of best hit entries written to new file: 1340

#### lnuF - 100 bp reads
- Total number of entries in hmm file: 5551
- Number of duplicate hmm matches: 0
- Number of best hit entries written to new file: 5551

#### lnuF - 150 bp reads
- Total number of entries in hmm file: 4429
- Number of duplicate hmm matches: 1
- Number of best hit entries written to new file: 4428

#### lnuF - 250 bp reads
- Total number of entries in hmm file: 2944
- Number of duplicate hmm matches: 0
- Number of best hit entries written to new file: 2944

#### lnuF - 300 bp reads
- Total number of entries in hmm file: 2519
- Number of duplicate hmm matches: 1
- Number of best hit entries written to new file: 2519 

#### lnuB/G - 100 bp reads
- Total number of entries in hmm file: 2668
- Number of duplicate hmm matches: 0
- Number of best hit entries written to new file: 2668

#### lnuB/G - 150 bp reads
- Total number of entries in hmm file: 2070
- Number of duplicate hmm matches: 1
- Number of best hit entries written to new file: 2069

#### lnuB/G - 250 bp reads
- Total number of entries in hmm file: 1336
- Number of duplicate hmm matches: 0
- Number of best hit entries written to new file: 1336

#### lnuB/G - 300 bp reads
- Total number of entries in hmm file: 1173
- Number of duplicate hmm matches: 0
- Number of best hit entries written to new file: 1173 

### c. Compile, score, and plot results

input rockout score results and hmm results and build bar plot.

For hmm results need to get total positive reads from mock metagenome. Find positive reads not mapped by hmm vs mapped by hmm and find non-positive reads mapped by hmm.
two files 1) mock metagenome 2) hmmsearch result

rocker results same thing but blastx results are split into passing and failing files.
three files 1) mock metagenome 2) passing 3) failing
from: /storage/home/hcoda1/9/rconrad6/scratch/ROCkOut/02_lnu/02a_lnuF_train

```bash
# lnuAll model

python 00_scripts/score_rocker_model.py -mm ../01b_lnuAll_test/simulated_reads/simread_100_raw_reads.fasta -bx test_score_100/alignments/simread_100_raw_reads.fasta_ROCkOut_alignments.txt -hm test_score_100/lnuAll-final_hmmSearch_100_filtered.tsv -rp test_score_100/ROCkOut_passing_alignments/simread_100_raw_reads.ROCkOut_passing.txt -rf test_score_100/ROCkOut_failing_alignments/simread_100_raw_reads.ROCkOut_failing.txt -o test_score_100/test_scores_100

python 00_scripts/00c_scripts/score_rocker_model.py -mm ../01b_lnuAll_test/simulated_reads/simread_150_raw_reads.fasta -bx test_score_150/alignments/simread_150_raw_reads.fasta_ROCkOut_alignments.txt -hm test_score_150/lnuAll-final_hmmSearch_150_filtered.tsv -rp test_score_150/ROCkOut_passing_alignments/simread_150_raw_reads.ROCkOut_passing.txt -rf test_score_150/ROCkOut_failing_alignments/simread_150_raw_reads.ROCkOut_failing.txt -o test_score_150/test_scores_150

python 00_scripts/00c_scripts/score_rocker_model.py -mm ../01b_lnuAll_test/simulated_reads/simread_250_raw_reads.fasta -bx test_score_250/alignments/simread_250_raw_reads.fasta_ROCkOut_alignments.txt -hm test_score_250/lnuAll-final_hmmSearch_250_filtered.tsv -rp test_score_250/ROCkOut_passing_alignments/simread_250_raw_reads.ROCkOut_passing.txt -rf test_score_250/ROCkOut_failing_alignments/simread_250_raw_reads.ROCkOut_failing.txt -o test_score_250/test_scores_250

python 00_scripts/00c_scripts/score_rocker_model.py -mm ../01b_lnuAll_test/simulated_reads/simread_300_raw_reads.fasta -bx test_score_300/alignments/simread_300_raw_reads.fasta_ROCkOut_alignments.txt -hm test_score_300/lnuAll-final_hmmSearch_300_filtered.tsv -rp test_score_300/ROCkOut_passing_alignments/simread_300_raw_reads.ROCkOut_passing.txt -rf test_score_300/ROCkOut_failing_alignments/simread_300_raw_reads.ROCkOut_failing.txt -o test_score_300/test_scores_300

# lnuF model

python 00_scripts/score_rocker_model.py -mm ../02b_lnuF_test/simulated_reads/simread_100_raw_reads.fasta -bx test_score_100/alignments/simread_100_raw_reads.fasta_ROCkOut_alignments.txt -hm test_score_100/lnuF-final_hmmSearch_100_filtered.tsv -rp test_score_100/ROCkOut_passing_alignments/simread_100_raw_reads.ROCkOut_passing.txt -rf test_score_100/ROCkOut_failing_alignments/simread_100_raw_reads.ROCkOut_failing.txt -o test_score_100/test_scores_100

python 00_scripts/score_rocker_model.py -mm ../02b_lnuF_test/simulated_reads/simread_150_raw_reads.fasta -bx test_score_150/alignments/simread_150_raw_reads.fasta_ROCkOut_alignments.txt -hm test_score_150/lnuF-final_hmmSearch_150_filtered.tsv -rp test_score_150/ROCkOut_passing_alignments/simread_150_raw_reads.ROCkOut_passing.txt -rf test_score_150/ROCkOut_failing_alignments/simread_150_raw_reads.ROCkOut_failing.txt -o test_score_150/test_scores_150

python 00_scripts/score_rocker_model.py -mm ../02b_lnuF_test/simulated_reads/simread_250_raw_reads.fasta -bx test_score_250/alignments/simread_250_raw_reads.fasta_ROCkOut_alignments.txt -hm test_score_250/lnuF-final_hmmSearch_250_filtered.tsv -rp test_score_250/ROCkOut_passing_alignments/simread_250_raw_reads.ROCkOut_passing.txt -rf test_score_250/ROCkOut_failing_alignments/simread_250_raw_reads.ROCkOut_failing.txt -o test_score_250/test_scores_250

python 00_scripts/score_rocker_model.py -mm ../02b_lnuF_test/simulated_reads/simread_300_raw_reads.fasta -bx test_score_300/alignments/simread_300_raw_reads.fasta_ROCkOut_alignments.txt -hm test_score_300/lnuF-final_hmmSearch_300_filtered.tsv -rp test_score_300/ROCkOut_passing_alignments/simread_300_raw_reads.ROCkOut_passing.txt -rf test_score_300/ROCkOut_failing_alignments/simread_300_raw_reads.ROCkOut_failing.txt -o test_score_300/test_scores_300

# lnuB/G model

python ../../00c_scripts/score_rocker_model.py -mm ../04b_lnuBG_test/simulated_reads/simread_100_raw_reads.fasta -bx test_score_100/alignments/simread_100_raw_reads.fasta_ROCkOut_alignments.txt -hm test_score_100/lnuBG-final_hmmSearch_100_filtered.tsv -rp test_score_100/ROCkOut_passing_alignments/simread_100_raw_reads.ROCkOut_passing.txt -rf test_score_100/ROCkOut_failing_alignments/simread_100_raw_reads.ROCkOut_failing.txt -o test_score_100/test_scores_100

python ../../00c_scripts/score_rocker_model.py -mm ../04b_lnuBG_test/simulated_reads/simread_150_raw_reads.fasta -bx test_score_150/alignments/simread_150_raw_reads.fasta_ROCkOut_alignments.txt -hm test_score_150/lnuBG-final_hmmSearch_150_filtered.tsv -rp test_score_150/ROCkOut_passing_alignments/simread_150_raw_reads.ROCkOut_passing.txt -rf test_score_150/ROCkOut_failing_alignments/simread_150_raw_reads.ROCkOut_failing.txt -o test_score_150/test_scores_150

python ../../00c_scripts/score_rocker_model.py -mm ../04b_lnuBG_test/simulated_reads/simread_250_raw_reads.fasta -bx test_score_250/alignments/simread_250_raw_reads.fasta_ROCkOut_alignments.txt -hm test_score_250/lnuBG-final_hmmSearch_250_filtered.tsv -rp test_score_250/ROCkOut_passing_alignments/simread_250_raw_reads.ROCkOut_passing.txt -rf test_score_250/ROCkOut_failing_alignments/simread_250_raw_reads.ROCkOut_failing.txt -o test_score_250/test_scores_250

python ../../00c_scripts/score_rocker_model.py -mm ../04b_lnuBG_test/simulated_reads/simread_300_raw_reads.fasta -bx test_score_300/alignments/simread_300_raw_reads.fasta_ROCkOut_alignments.txt -hm test_score_300/lnuBG-final_hmmSearch_300_filtered.tsv -rp test_score_300/ROCkOut_passing_alignments/simread_300_raw_reads.ROCkOut_passing.txt -rf test_score_300/ROCkOut_failing_alignments/simread_300_raw_reads.ROCkOut_failing.txt -o test_score_300/test_scores_300
```

**Results:**

#### lnuAll - 100 bp reads
- Positive read count from metagenome: 2462
- Positive read count from Diamond/BlastX alignment: 2432
- Positive reads not aligned by BlastX: 30
- Positive read count from ROCkOut: 2432
- Positive reads not reported by ROCkOut: 30
- Positive read count from hmmsearch: 2057
- Positive reads not aligned by hmmsearch: 405

#### lnuAll - 150 bp reads
- Positive read count from metagenome: 1760
- Positive read count from Diamond/BlastX alignment: 1717
- Positive reads not aligned by BlastX: 43
- Positive read count from ROCkOut: 1717
- Positive reads not reported by ROCkOut: 43
- Positive read count from hmmsearch: 1552
- Positive reads not aligned by hmmsearch: 208

#### lnuAll - 250 bp reads
- Positive read count from metagenome: 1201
- Positive read count from Diamond/BlastX alignment: 1141
- Positive reads not aligned by BlastX: 60
- Positive read count from ROCkOut: 1141
- Positive reads not reported by ROCkOut: 60
- Positive read count from hmmsearch: 1027
- Positive reads not aligned by hmmsearch: 174

#### lnuAll - 300 bp reads
- Positive read count from metagenome: 965
- Positive read count from Diamond/BlastX alignment: 926
- Positive reads not aligned by BlastX: 39
- Positive read count from ROCkOut: 926
- Positive reads not reported by ROCkOut: 39
- Positive read count from hmmsearch: 873
- Positive reads not aligned by hmmsearch: 92

#### lnuF - 100 bp reads
- Positive read count from metagenome: 4126
- Positive read count from Diamond/BlastX alignment: 4060
- Positive reads not aligned by BlastX: 66
- Positive read count from ROCkOut: 4060
- Positive reads not reported by ROCkOut: 66
- Positive read count from hmmsearch: 3232
- Positive reads not aligned by hmmsearch: 894

#### lnuF - 150 bp reads
- Positive read count from metagenome: 2998
- Positive read count from Diamond/BlastX alignment: 2943
- Positive reads not aligned by BlastX: 55
- Positive read count from ROCkOut: 2943
- Positive reads not reported by ROCkOut: 55
- Positive read count from hmmsearch: 2611
- Positive reads not aligned by hmmsearch: 387

#### lnuF - 250 bp reads
- Positive read count from metagenome: 2004
- Positive read count from Diamond/BlastX alignment: 1972
- Positive reads not aligned by BlastX: 32
- Positive read count from ROCkOut: 1972
- Positive reads not reported by ROCkOut: 32
- Positive read count from hmmsearch: 1775
- Positive reads not aligned by hmmsearch: 229

#### lnuF - 300 bp reads
- Positive read count from metagenome: 1627
- Positive read count from Diamond/BlastX alignment: 1604
- Positive reads not aligned by BlastX: 23
- Positive read count from ROCkOut: 1604
- Positive reads not reported by ROCkOut: 23
- Positive read count from hmmsearch: 1479
- Positive reads not aligned by hmmsearch: 148

#### lnuB/G - 100 bp reads
- Positive read count from metagenome: 1643
- Positive read count from Diamond/BlastX alignment: 1605
- Positive reads not aligned by BlastX: 38
- Positive read count from ROCkOut: 1605
- Positive reads not reported by ROCkOut: 38
- Positive read count from hmmsearch: 1335
- Positive reads not aligned by hmmsearch: 308

#### lnuB/G - 150 bp reads
- Positive read count from metagenome: 1187
- Positive read count from Diamond/BlastX alignment: 1146
- Positive reads not aligned by BlastX: 41
- Positive read count from ROCkOut: 1146
- Positive reads not reported by ROCkOut: 41
- Positive read count from hmmsearch: 1028
- Positive reads not aligned by hmmsearch: 159

#### lnuB/G - 250 bp reads
- Positive read count from metagenome: 793
- Positive read count from Diamond/BlastX alignment: 757
- Positive reads not aligned by BlastX: 36
- Positive read count from ROCkOut: 757
- Positive reads not reported by ROCkOut: 36
- Positive read count from hmmsearch: 666
- Positive reads not aligned by hmmsearch: 127

#### lnuB/G - 300 bp reads
- Positive read count from metagenome: 626
- Positive read count from Diamond/BlastX alignment: 606
- Positive reads not aligned by BlastX: 20
- Positive read count from ROCkOut: 606
- Positive reads not reported by ROCkOut: 20
- Positive read count from hmmsearch: 574
- Positive reads not aligned by hmmsearch: 52

![lnu Test Set Scoring results](https://github.com/rotheconrad/ROCker_Macrolide_Models/blob/main/04_lnu/00_figures/04b-lnuMockMetaScores.png)

Top panel is for the lnu all model, middle panel is for the lnuF model, and bottom panel is for the lnuB/G model.

### d. Build pie trees with pplacer and itol.

ROCkOut already includes the pie tree output with its "place" function. This step is really just retrieve the right file and using the [iTol](https://itol.embl.de/) website.

##### Read placement tree for lnu All model.

![Read placement tree for lnu All model](https://github.com/rotheconrad/ROCker_Macrolide_Models/blob/main/04_lnu/00_figures/04c-lnuall-phylo-placement.png)

##### Read placement tree for lnuF model.

![Read placement tree for lnuF model](https://github.com/rotheconrad/ROCker_Macrolide_Models/blob/main/04_lnu/00_figures/04d-lnuF-phylo-placement.png)

##### Read placement tree for lnuB/G model.

![Read placement tree for lnuB/G model](https://github.com/rotheconrad/ROCker_Macrolide_Models/blob/main/04_lnu/00_figures/04e-lnuBG-phylo-placement.png)

### f. Verify TP FP TN FN labels and scoring analysis (additional visuals)

The score_rocker_model.py script writes out a fasta directory with file names labeled for filter method and failed, passed, TP, FP, TN, or FN.

The goal of this section is to test the results of mapping the reads in nucleotide format to the genes in nucleotide with regular blastn vs. mapping the the reads in nucleotides to the genes in amino acids with blastx.

First step is to create gene databases with nucl and prot, then to map all the different read files in the fasta directory created by score_rocker_model.py.

The second step will be to create some figures to visualize the results.

I just show results for the lnuAll model here, but the results were and can be repeated with lnuF or any other model.

```bash
> python ../../00c_scripts/ROCkOut_filter_visualize.py -p test_score_100/ROCkOut_passing_alignments/simread_100_raw_reads.ROCkOut_passing.txt -f test_score_100/ROCkOut_failing_alignments/simread_100_raw_reads.ROCkOut_failing.txt -o test_score_100/ROCkOut_Filter_Viz_100 -a test_score_100/alignments/simread_100_raw_reads.fasta_ROCkOut_alignments.txt

> python ../../00c_scripts/ROCkOut_filter_visualize.py -p test_score_150/ROCkOut_passing_alignments/simread_150_raw_reads.ROCkOut_passing.txt -f test_score_150/ROCkOut_failing_alignments/simread_150_raw_reads.ROCkOut_failing.txt -o test_score_150/ROCkOut_Filter_Viz_150 -a test_score_150/alignments/simread_150_raw_reads.fasta_ROCkOut_alignments.txt

> python ../../00c_scripts/ROCkOut_filter_visualize.py -p test_score_250/ROCkOut_passing_alignments/simread_250_raw_reads.ROCkOut_passing.txt -f test_score_250/ROCkOut_failing_alignments/simread_250_raw_reads.ROCkOut_failing.txt -o test_score_250/ROCkOut_Filter_Viz_250 -a test_score_250/alignments/simread_250_raw_reads.fasta_ROCkOut_alignments.txt

> python ../../00c_scripts/ROCkOut_filter_visualize.py -p test_score_300/ROCkOut_passing_alignments/simread_300_raw_reads.ROCkOut_passing.txt -f test_score_300/ROCkOut_failing_alignments/simread_300_raw_reads.ROCkOut_failing.txt -o test_score_300/ROCkOut_Filter_Viz_300 -a test_score_300/alignments/simread_300_raw_reads.fasta_ROCkOut_alignments.txt
```

##### lnu All model read mapping and classification for bitscore.

![lnu All model read mapping and classification for bitscore](https://github.com/rotheconrad/ROCker_Macrolide_Models/blob/main/04_lnu/00_figures/04f_verify_bitscore.png)

##### lnu All model read mapping and classification for percent sequence identity.

![lnu All model read mapping and classification for percent sequence identity](https://github.com/rotheconrad/ROCker_Macrolide_Models/blob/main/04_lnu/00_figures/04g_verify_pid.png)







