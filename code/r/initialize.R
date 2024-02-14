# initialize.R



# Script structure
# >>>>>>>>>>>>>>>>
# 1. Directory settings: Directory structure for result storage.
# 2. Library: Loading of depended R packages and local functions
# 3. Summary of the purpose of each scripts and their ordinal positions




# 1. Directory settings: Directory structure for result storage.
# ==============================================================================

# Output directory
if(!dir.exists("results/figures")){
  dir.create(path = "results/figures", recursive = TRUE, mode = "0700")
}
if(!dir.exists("results/data")){
  dir.create(path = "results/data", recursive = TRUE, mode = "0700")
}


out_figures <- "results/figures/"
out_data <- "results/data/"


#
# ==============================================================================




# 2. Library: Loading of depended R packages and local functions
# ==============================================================================

# Public
library(affy)
library(tidyverse) # tidyverse_packages(), List all packages in the tidyverse
library(sva) # ComBat
library(pROC) # roc, auc
library(irr) # kappa2: Cohen's kappa
# library(hablar) # retype: Transforms all elements into simple classes
# library(DescTools) # RobScale: Robust Scaling With Median and Mad
library(rmarkdown) # render: To render the Rmarkdown file into specific output format
library(readxl) # read_excel()
library(writexl) # write_xlsx()
library(ggpubr) # ggarrange()
library(GGally) # ggpairs(), ggcatmat()
library(Biobase) # rowMedians()


# Private
source("code/r/functions_cpp.R") # R interface for detection calling functions implemented in cpp
source("code/r/functions.R")


#
# ==============================================================================




# 3. Summary of the purpose of each scripts and their ordinal positions
# ==============================================================================


# 1. initialize.R:
#       a. Initialize folder structure to store artefacts fo the analysis
#       b. Load necessary R packages
# 2. functions.R:
#       a. Supporting functions for other scripts sourced under "Library" section.
# 3. functions_cpp.R
#       a. R interface for detection calling functions implemented in C++
# 4. experiment_metadata_curation.R
#       a. Curation of experimental data associated with each samples.
# 5. expression_data_processing_and_qc.R:
#       a. Normalization and filtering of samples with percent-present values
#          - less than 10%.
#       b. The feasibility RMA and MAS normalization of FFPE derived
#          - expression profiles were compared.
#       c. QC stats and plots generated.
# 6. ffpe_norm_comparison_supplementary.R
#       a. MAS and RMA normalization comparison on independent datasets.
# 7. ffpe_norm_comparison_finher_cocorrelation.R
#       a. MAS-RMA normalization comparison w.r.t ER pathway gene co-correlation structure.
# 7. technical_validation.R
#       a. IHC vs PAM50 subtyping comparison
#       b. ER/HER2 IHC status prediction from respective mRNAs
#       c. Technical validation using technical replicates.
# 9. geo_GSE47994_updation.R
#       a. Update existing GSE47994 with two technical replicates missing in
#          - original submission.
#       b. Use RMA normalize expression data for backward compatibility,
#          as original submission used RMA normalized data
#       b. Update GSE47994 with experimental metadata missing in original submission
# 10. geo_reprocessed_submission.R
#       a. Prepares the following gene expression matrices and sample qc data
#          - to include in GSE47998 as supplementary.
#           a.1 RMA and MAS normalized data - all samples
#           a.2 RMA and MAS normalized data - sample filtered
#
# ==============================================================================



