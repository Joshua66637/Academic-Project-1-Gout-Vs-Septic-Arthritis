# Project Title: Differential Gene Expression Analysis
# Course: Statistics and Data Analysis for Bioinformatics (MSc Bioinformatics)

## Description
Arthritis refers to any disorder causing inflammation in the joints of patients. Gout is a common type of arthritis that occurs due to the precipitation of monosodium urate in the joints as a result of hyperuricemia, leading to an autoimmune response and inflammation. The main genes involved in gout are related to the innate inflammatory pathways and neutrophil extracellular trapping. 
Septic arthritis occurs when the joint inflammation is due to bacterial invasion of the synovium and joint space. Neutrophils are  crucial for engulfing the pathogens but can also cause tissue destruction. Synovial fluid analysis and Gram staining are common diagnostic methods.

**Research Question:** Gout and Septic Arthritis exhibit many clinical symptoms in common, such as inflammation, fever, redness, and reduced movement in the affected joints. Even Gram staining and the presence of monosodium urate crystals cannot reliably exclude infection due to their high false positive rate. Given the rapid spread and potential lethality of septic arthritis, it is advisable to consider alternative diagnostic methods, such as the analysis of bulk RNAseq data. This simple project performs a differential gene expression analysis to differentiate Gout and Septic Arthritis (SA).

## Data Files
The following data files are required for this analysis:
- `Annotations.csv`: Contains gene annotations.
- `DE_GOUT_vs_HC.csv`: Differential expression results for Gout vs Healthy Control.
- `DE_SA_vs_HC.csv`: Differential expression results for SA vs Healthy Control.
- `Sample_Information.csv`: Information about the samples used.
- `Expression_Table.csv`: Gene expression matrix.

## Results

- Sample summary: 9 samples per group; sex ratio ≈ 1:1. Mean neutrophils — **Sepsis: 12.78**, **Gout: 7.48**, **Healthy: 5.04**.

- Significant differentially expressed genes: **15** (Gout), **296** (Septic Arthritis).  
- Shared significant genes (8): **MYO3B, SPP1, GATD3A, PCP2, KLHDC7A, EGFL6*** and ***RPL35P5**.  
  
- General Linear models showed that there was no significant difference in neutrophil count and patient sex across disease groups.  
- A hypergeometric test (p= **5.2e-13**) showed significant overepresentation of shared genes between the two diseases.  
 
- Expression patterns: Two putative biomarkers: **SNORA73B** (Up in SA, Down in Gout) and **MYO3B** (Down in SA, Up in Gout). Both genes have an intermediate expression value in healthy controls, thus helping to diffrentiate all three groups. SNORA73(small nucleolar RNA73) aids resistance of cells to oxidative stress, inflammation and cell death. MYO3B (class III myosin) is an actin-based motor protein which binds to cell-surface proteins. Differential expression of myosins has been linked to autoimmune diseases.
