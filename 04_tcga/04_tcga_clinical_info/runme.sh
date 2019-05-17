#!/bin/bash -l

Rscript ./01_retrieve_clinical_info.r
Rscript ./03_add_sample_key.r
Rscript ./04_annotate_clinical_to_trm.r
Rscript ./05_beta_pau_vs_stypes_boxplot_query_genes.r
