# Construction an HMM for the BPTI Kunitz-Protease Inhibitor Domain
This repository implements a pipeline to detect the BPTI/Kunitz protease inhibitor domain (PF00014) using profile HMMs. It integrates structural bioinformatics and sequence modeling to build a structure-guided HMM and compare its performance to a standard MSA-based HMM, testing whether structure improves detection. in real proteins across datasets

# Project Overview
This repository contains an original implementation of a pipeline for detecting the BPTI/Kunitz-type protease inhibitor domain (Pfam: PF00014) using profile Hidden
Markov Models (HMMs). this project combines structural bioinformatics with sequencebased modeling to differentiate proteins containing the Kunitz domain from those that
do not. The primary objective is to design and validate a structure-guided Hidden Markov Model (HMM), and to evaluate its predictive performance against a conventional sequence-based HMM constructed from multiple sequence alignment (MSA). This comparative approach aims to assess whether incorporating structural
information enhances the specificity and sensitivity of Kunitz domain protein detection


# HMMs for Protein Domains
Protein domains like Kunitz often preserve their core structure and functional residues
even across distantly related sequences. Profile HMMs model this variability by
capturing:
• Conserved residues (emissions)

• Insertions/deletions (gaps)

• Probabilities across positions

By training an HMM on multiple structurally aligned domain instances, we can more
accurately detect related proteins—even when sequence identity is low.


# Required Tools
Install these via conda:

conda install -c bioconda cd-hit hmmer blast muscle

conda install -c conda-forge biopython

<img width="760" height="280" alt="image" src="https://github.com/user-attachments/assets/04c74ece-ea24-46a4-a6c8-ccb955864e31" />


# Pipeline Steps
1. Extract Representative Kunitz Sequences from PDB
   
• Use advanced search: Data Collection Resolution <= 3.5 AND ( Identifier = "PF00014"
AND Annotation Type = "Pfam" ) AND Polymer Entity Sequence Length <= 80 AND
Polymer Entity Sequence Length >= 45

• Press the custom report with the following flags:
• Entry ID

• PDB ID

• Entity ID

• Auth Asym ID

• Sequence

• Annotation Identifier

• Data collection resolution

This will output a .CSV file. Then execute the following bash script:
Code: bash script_recover_representative_kunitz.sh

This will do the following operations:

• Extract sequences of PF00014 domains from the PDB custom report

• Cluster the sequences using CD-HIT at 90% identity threshold

• Extract the most representative ID from each cluster

• Retrieve the sequences of the representative IDs and store them in a new FASTA file.

• Generate the PDBefold file (convert IDs to the expected format:PDB:CHAIN): tmp_pdb_efold_ids.txt

• BEFORE GOING TO THE NEXT STEP: check in the tmp_pdb_efold_ids.txt file.
sequences that are too long and remove them.


# 2. Download the full Swiss-Prot protein dataset in FASTA format, used to extract negative sets and full Kunitz reference
• Go to https://www.uniprot.org/ and save the FASTA file as uniprot_sprot.fasta

. Download all Kunitz proteins in FASTA format

• Go to https://www.uniprot.org/ and save the FASTA file as all_kunitz.fasta

# 3. Structural Alignment

Use PDBeFold :

• Input: tmp_pdb_efold_ids.txt (PDB IDs)

• Output: Aligned FASTA file pdb_kunitz_rp.ali

# 4. Build HMM from Structural Alignment

bash create_hmm_str.sh

bash create_testing_sets.sh

• Build a structural HMM from the PDBeFold structural alignment

• Remove training sequences from the full dataset and generate random subsets of
positive and negative sequences, which will be used to create test sets

• Automatically identify the optimal E-value thresholds via MCC evaluation (2-fold CV)

• Perform evaluation on: Set 1, using Set 2’s threshold; Set 2, using Set 1’s threshold; Combined Set 1 + Set 2, using both thresholds for overall assessment

• Report MCC, precision, recall, and identify false positives and false negatives Write detailed results to hmm_results_strali.txt

# 5. Evaluate Performance

Use the provided script performance.py to compute:

• MCC (Matthews Correlation Coefficient)

• Accuracy (Q2)

• True Positive Rate (TPR)

• Precision (PPV)

# 6. Plot the confusion matrices for each run resulting in the hmm_results_strali.txt file
Execute the python script confusion_matrix.py , put in the values for FN (false negatives),FP (false positives),TN (true negatives) and TP (true positives) variables. The
script will output a .png image with the plotted confusion matrix. The script has to be done for each set of values.


# Output Files Summary

• hmm_results_strali.txt and hmm_results_seqali.txt contain:

• The best E-value thresholds selected by maximizing the Matthews Correlation Coefficient
(MCC) -Performance metrics for each test set and overall, calculated using the E-value that yielded the highest MCC, based on either full sequence or best single domain
evaluations -Lists of false positives and false negatives

• neg_1.fasta and neg_2.fasta: FASTA files of non-Kunitz sequences used respectively as
the negative set 1 and set 2.

• pos_1.fasta and pos_2.fasta: FASTA files of Kunitz (positive) sequences used in set 1
and set 2.

• pdb_kunitz_rp_clean.fasta: Cleaned representative Kunitz sequences used for both
structural and sequence alignments (after filtering for length).

• pdb_kunitz_rp_seqali.fasta and pdb_kunitz_rp_seqali.hmm: The multiple sequence
alignment and resulting HMM model produced using MUSCLE.

• pdb_kunitz_rp_strali.fasta and pdb_kunitz_rp_strali.hmm: The multiple structure
alignment and resulting HMM model built from PDBeFold.
• set_1_strali.class and set_2_strali.class: Classification results from hmmsearch for
set 1 and set 2 (structure-based model).

• temp_overall.class: Classification results for the combined dataset (positive +
negative), used to estimate overall performance for each threshold

