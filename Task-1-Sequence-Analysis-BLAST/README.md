 Task 1 – DNA/Protein Sequence Analysis and BLAST Homology Search

Objective

The objective of this task is to retrieve a protein sequence from UniProt and perform a BLASTP search to identify homologous protein sequences and analyze their sequence similarity.

 Protein Information

| Parameter | Details |
|---|---|
| Protein | Cellular tumor antigen p53 |
| Gene | TP53 |
| Organism | Homo sapiens |
| UniProt Accession | P04637 |
| Molecule Type | Protein |
| Sequence Length | 393 amino acids |
| BLAST Program | BLASTP |
| Database | NCBI nr |

Methodology

1. The human TP53 protein sequence was retrieved from UniProt using accession P04637.
2. The protein sequence was downloaded in FASTA format.
3. BLASTP was performed using the NCBI BLAST web server.
4. The NCBI non-redundant (nr) protein database was used.
5. The top homologous sequences were analyzed based on query coverage, percentage identity and E-value.
6. Pairwise sequence alignments were examined to identify conserved regions.

 Top BLAST Hits

| Hit Protein | Organism | Query Cover | E-value | Identity | Accession |
|---|---|---:|---:|---:|---|
| Cellular tumor antigen p53 isoform X1 | Gorilla gorilla gorilla | 100% | 0.0 | 99.75% | XP_063556659.1 |
| Chain K, Cellular tumor antigen p53 | Homo sapiens | 100% | 0.0 | 100.00% | 7XZZ_K |
| Chain C, Cellular tumor antigen p53 | Homo sapiens | 100% | 0.0 | 100.00% | 8R1F_C |
| Cellular tumor antigen p53 isoform a | Homo sapiens | 100% | 0.0 | 100.00% | NP_000537.3 |

Results

The BLASTP analysis identified highly similar homologous sequences for the human TP53 protein. The top hit from Gorilla gorilla gorilla showed 100% query coverage, 99.75% sequence identity and an E-value of 0.0.

Conclusion

The BLASTP analysis demonstrated that TP53 is a highly conserved protein. The high sequence identity and complete query coverage observed among the top hits support the evolutionary conservation of TP53.

Files Included

- `sequence/TP53_P04637.fasta` – Protein FASTA sequence
- `screenshots/` – BLAST and sequence analysis screenshots
- `report/Task_1_BLAST_Analysis_Report.pdf` – Detailed analysis report
