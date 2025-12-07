# CRISPR Guide RNA Designer

This small Python toolkit finds PAM sites (default SpCas9 NGG), generates candidate gRNAs, computes simple on-target heuristics (GC content) and a basic off-target penalty against a provided background sequence (e.g., a chromosome or whole-genome FASTA). It outputs an optimality score combining on-target and off-target metrics.

Usage examples:

```bash
# Simple usage with a single sequence:
python run_designer.py --sequence "ATGCGT..." --pam NGG

# Provide a background FASTA to evaluate off-targets (will scan windows for potential off-target matches):
python run_designer.py --sequence-file gene.fasta --background background.fasta
```

Requirements: see requirements.txt
