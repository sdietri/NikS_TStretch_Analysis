# NikS_TStretch_Analysis

Script for the NikS promotor T-stretch length analysis as used by Eisenbart et al 2020.

Requirements:
- NCBI Blast suite 2.8+
- Python 3.5+, seaborn 0.10.1, pandas 1.0.3, BioPython 1.76


Input required:
- Your Email address to use for Entrez access
- SRA sequences from SRP162088 "Helicobacter pylori genome evolution within the human stomach", converted to FASTA. File names needs to follow the scheme "SRR-id.fas" to allow for Entrez information lookup
- Start sequence of NikS gene for identification of upstream promotor region

Command:

python3 NiksSearch.py "email" "Folder to SRR fastas" "Start sequence NikS fasta"
