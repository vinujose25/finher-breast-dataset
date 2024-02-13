
R script structure
==================

When rerunning scripts, keep the order of the scripts in mind, as any artifacts generated in one script may be used in other scripts following it.

1. initialize.R
1.1. The main script from which other scripts are managed

2. functions.R
2.1. Supporting functions for  other scripts.

3. functions_cpp.R
3.1. R interface for detection calling functions implemented in C++

4. sample_characteristics_curation.R
4.1. Metadata curation of the expression profiling assay
4.2. Clinical data curation

5. expression_data_processing_and_qc.R
5.1. Expression data processing, as described in Figure 2
5.2. Comparison of MAS and RMA normalization methods.
5.3. Sample and probeset QC.before and after filtering.

6. supplementary_ffpe_norm_comparison.R
6.1. Comparison of MAS and RMA normalization methods using independent datasets.
6.2. ER/HER gene co-correlation structure comparison using FinHER data

7. technical_validation.
7.1. Technical replicate agreement.

8. figures_tables_data.R
8.1. Generate journal-specific figures and tables.

9. geo_GSE47994_updation.R
9.1. Update existing GSE47994 with two technical replicates missing in original submission.
9.2. Use RMA normalize expression data for backward compatibility, as original submission used RMA normalized data.
9.3. Update GSE47994 with experimental metadata missing in original submission

10. geo_reprocessed_submission.R
10.1. Prepares the following gene expression matrices and sample qc data to include in GSE47998 as supplementary.
10.1.1. RMA and MAS normalized data - all samples + all probesets
10.1.2. RMA and MAS normalized data - all samples + probeset filtering
10.1.3. RMA and MAS normalized data - sample filtering + probeset filtering



