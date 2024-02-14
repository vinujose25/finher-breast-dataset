R script structure
==================

When running scripts, keep the order of the scripts in mind, as any artifacts generated in one script are used in other scripts following it.

1. initialize.R
1.1. The main script from which other scripts are managed

2. functions.R
2.1. Supporting functions for other scripts.

3. functions_cpp.R
3.1. R interface for detection calling functions implemented in C++

4. sample_characteristics_curation.R
4.1. Metadata curation of the expression profiling assay
4.2. Clinical data curation

5. expression_data_processing_and_qc.R
5.1. Expression data processing, as described in Figure 2
5.2. Comparison of MAS and RMA normalization methods.
5.3. Sample and probeset QC before and after filtering.

6. ffpe_norm_comparison_supplementary.R
6.1. Comparison of MAS and RMA normalization methods using independent datasets.

7. ffpe_norm_comparison_finher_cocorrelation.R
7.1. ER/HER gene co-correlation compared using FinHER data

8. technical_validation.
8.1. IHC vs PAM50 subtyping comparison
8.2. ER/HER2 IHC status prediction from respective mRNAs
8.3. Technical replicate agreement.

8. geo_GSE47994_updation.R
8.1. Update existing GSE47994 with two technical replicates missing from the original submission.
8.2. Use RMA normalized expression data for backward compatibility, as the original submission used RMA normalized data.
8.3. Update GSE47994 with experimental metadata missing from the original submission

9. geo_reprocessed_submission.R
9.1. Prepare gene expression matrices and sample QC data to include in GSE47998 as supplemental data.
