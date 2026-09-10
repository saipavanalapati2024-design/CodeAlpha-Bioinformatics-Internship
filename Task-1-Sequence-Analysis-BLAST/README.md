# Task 1 - DNA/Protein Sequence Analysis and BLAST

## Objective

Retrieve a biological sequence and use BLAST to identify homologous sequences.

## Sequence

- Protein: Cellular tumor antigen p53
- Gene: TP53
- Organism: Homo sapiens
- UniProt accession: P04637
- Sequence length: 393 amino acids
- Source: UniProt
- FASTA file: [sequences/TP53_p04637.fasta](sequences/TP53_p04637.fasta)

## BLAST Analysis

- Program: BLASTP
- Database: NCBI non-redundant protein database (nr)
- Top hit: Cellular tumor antigen p53 isoform X1, *Gorilla gorilla gorilla*
- Query coverage: 100%
- Identity: 392/393 (99.75%)
- Similarity/positives: 392/393 (99.75%)
- Alignment score: 814 bits
- Gaps: 0/393 (0%)
- E-value: 0.0

The results show that TP53 is highly conserved between human and gorilla. The
complete query coverage, very high identity, strong alignment score, and zero
gaps support a close homologous relationship.

## Evidence

- [UniProt sequence screenshot](Screenshots/01_Uniport_TP53.png)
- [Downloaded FASTA screenshot](Screenshots/02_TP53_FASTA.png)
- [BLAST input screenshot](Screenshots/03_BLAST_Input.png)
- [BLAST results summary](Screenshots/04_BLAST_Results_summary.png)
- [Top BLAST hit summary](Screenshots/05_Top_BLAST_summary.png)
- [Pairwise alignment](Screenshots/06_BLAST_Alignment.png)

## Report

The detailed 2-3 page analysis report is available here:

[Task_1_BLAST_Analysis_Report.pdf](report/Task_1_BLAST_Analysis_Report.pdf)