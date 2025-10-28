# Taxonomic overmatching analysis 
This repository documents the taxonomic overmatching analysis for the preprint and submission of "Environmental DNA Reveals Reykjavík’s Human and Ecological History"; doi: https://doi.org/10.1101/2025.10.08.681091. The analysis includes remapping and taxonomic identification of eleven different genera within the two plant families Rosaceae and Poaceae. The goal was to explain the presence in our data of Iceland's non-native detected taxonomic groups and highlight the database limitations we encountered. The following analysis can be split into three parts. 

## 1. Database construction 
We downloaded all available references for the two families, Poaceae and Rosaceae, with a custom Python script ```scripts/download_all_for_fam.py```. The list of NCBI accession numbers: ```Poaceae_db/poa_accession.txt``` and Rosaceae_db/rosa_accession.txt```. 

## 2. Read extraction and remapping
We extracted all reads that were initially classified to the 11 different genera and remapped them to their corresponding family database (```scripts/get_reads.sh```). The mapping scripts are found in ```scripts/map_bowtie.sh``` and ```scripts/map_botwie2P.sh```. Merged data was analysed with metaDMG LCA using ```scripts/metaDMG_run.sh```. 

## 3. Data visualisation 
All data was visualised using custom R scripts found in ```scripts/plotting.R```. The plots used in the publication can be found in the main repository. 
