# hnRNPL project code
This repository contains code used in the analyses for the article "A conserved role of hnRNPL in regulating alternative splicing of transcriptional regulators necessary for B cell activation".

The RNA-seq datasets were processed using [this RNA-seq pipeline](https://github.com/tellyalogicalguy/RNAseq_pipeline).

The output files produced by the pipeline were compared using R scripts within the `Compilation` folder.

### Gene expression
`10_compilation_expr_FDR_pt1_L2FC_inc.R` contains code to compare gene expression changes between datasets.  
`2_gseaForAllDatasets.R` is to do GSEA analysis using the `fgsea` R package for differentially expressed genes for all datasets.  
### Splicing
`1_compilation_splicSE_FDR_lt_pt2.R` contains code for identifying conserved skipped exon alternative splicing events between mouse and human datasets.  
`8_covFilt_datatable.R` is for adding coverage and inclusion level filters for the comparative splicing analysis.  
`9_gprofilerForSplicingChanges_covFilt_datatable.R` is for doing Gene Ontology enrichment analyses for splicing changes using `gprofiler2` R package.  
