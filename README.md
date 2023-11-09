# arfAdapt
The data and code made available here is for the analysis of all omics data acquired and presented within the manuscript.

## Code 

- The following two code files were used to perform analysis on RNA and ATAC sequencing data according to standard published pipelines. Here in this repository, they are presented only as code files, and are not directly functional within this directory.
    - DESeq2_RG_SDELCISAR_13_rnaseq.Enrichment.Rmd
    - ATACseq_RG_SDELCISAR_11_Enrichment.Rmd

- The code file adhoc_omics_integration.R is an ad hoc integration approach developed for this project in the lab, and the code is presented within the repository along with all its data dependencies and outputs.

## Data

All data except the following were generated as part of the project and raw files have been made available in the GEO repository associated with the manuscript.
- UCSC data set called Data_UCSC/hgTables.txt

While using data from the GEO repository, or from supplementary data files published with the manuscript, please observe the following nomencleture: All ATAC data is labelled with the letter A and all RNA data is labelled with the letter R. A or R is followed by either N, 1 or 6 which correspond to shNT treatement, shASAP1 treatment or shARF6 treatment respectively. Any other treatments or conditions should be ignored. This is followed by 03, 07, 14 or 21 which represent the time points and finally 01, 02, or 03, which represents biological replicates. An example of RNA-seq data for shARF6 treatment for day 3 replicate 1 should read R60301.


**Note**: ChEA3 data available within the RNA folder was generated manually using their [online tool](https://maayanlab.cloud/chea3/) and downloaded.
