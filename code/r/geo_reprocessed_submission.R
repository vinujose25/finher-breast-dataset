# geo_reprocessed_submission.R


# What this script does?
# >>>>>>>>>>>>>>>>>>>>>>
#
# This script prepares the following gene expression matrices and sample qc data
# to include in GSE47998 as supplementary
# 1) RMA and MAS normalized data - all samples + all probesets + combat
# 2) RMA and MAS normalized data - sample filtered + probeset filtered + combat

# ! Note MAS normalized data with "sample filtered + probeset filtered + combat" is
# recommended for downstream analysis.




# Script structure
# >>>>>>>>>>>>>>>>
#
# 1. Load data
# 2. Update inhouse to GSE47994 map with replicates geo sample accession
# 3. Format inhouse data for GEO, make it congruent to GSE47994, and write out
#   3.1. RMA and MAS normalized data - all samples
#   3.2. RMA and MAS normalized data - sample filtered




# 1. Load data
# ==============================================================================

# FinHER inhouse matrix - All (TN/HER2) and filtered
load("results/data/finher_expr_norm1.RData")
load("results/data/finher_expr_norm2.RData")

# FinHER inhouse sample qc - All (TN/HER2) and filtered
load("results/data/sample_qc_norm1.RData")
load("results/data/sample_qc_norm2.RData")

# FinHER inhouse probeset qc - All (TN/HER2) and filtered
load("results/data/pset_qc_norm1.RData")
load("results/data/pset_qc_norm2.RData")

# FinHER inhouse to GSE47994 mapping
celmap <- read_tsv(file = "results/data/GSE47994_to_inhouse_sample_map.tsv")
dim(celmap) # 337 3


names(finher_expr_norm2) # rma mas
names(finher_expr_norm2$mas) # "sum"    "pvac"   "maxvar" "combat"

names(sample_qc_norm2) # rma mas
names(sample_qc_norm2$mas) # dataframe of 15 columns
# [1] "CEL_filename" "PP"           "Neg_Avg"      "Pos_Avg"      "GE_Avg"       "PC1"         
# [7] "PC2"          "PC3"          "PC4"          "PC1_Var_Prop" "PC2_Var_Prop" "PC3_Var_Prop"
# [13] "PC4_Var_Prop" "Subject_Id"   "Batch" 

names(pset_qc_norm2) # rma mas
names(pset_qc_norm2$mas) # dataframe of 8 columns
# [1] "Probeset_Id"        "Pset_GE_Var"        "Pset_GE_Avg"        "Pset_PP_Cor"
# [5] "Pset_PP_Cor_Pval"   "Probe_PC1_Var_Prop" "PVAC_Status"        "Probeset_Type"

names(celmap)
# [1] "inhouse"  "GSE47994" "inhouse2"

# ==============================================================================


# 2. Update inhouse to GSE47994 map with replicate's geo accession
# ==============================================================================

# Replicates
celmap %>% dplyr::filter(str_detect(inhouse, "_"))

#   inhouse            GSE47994       inhouse2
# 1 A1764-009_a.CEL.gz GSM1164303     A1764-009_a
# 2 A1764-040_a.CEL.gz GSM1164334     A1764-040_a
# 3 A1764-009_b.CEL.gz rep_GSM1164303 A1764-009_b
# 4 A1764-040_b.CEL.gz rep_GSM1164334 A1764-040_b

# New geo accessions of replicates added

# Ref: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE47994
# GSM8011739	Breast Cancer 26 b
# GSM8011740	Breast Cancer 113 b

# Get sample id for correct mapping
sample_qc_norm1$mas %>% dplyr::filter(str_detect(CEL_filename, "_")) %>%
  dplyr::select(CEL_filename, Subject_Id, Batch)
#   CEL_filename Subject_Id Batch
# 1 A1764-009_a  26         Batch1
# 2 A1764-040_a  113        Batch1
# 3 A1764-009_b  26         Batch1
# 4 A1764-040_b  113        Batch1

# Correct mapping
# GSM8011739	Breast Cancer 26 b  A1764-009_b
# GSM8011740	Breast Cancer 113 b A1764-040_b

# Update celmap
celmap <- celmap %>%
  dplyr::mutate(GSE47994 = case_when(
    inhouse2 == "A1764-009_b" ~ "GSM8011739",
    inhouse2 == "A1764-040_b" ~ "GSM8011740",
    TRUE ~ GSE47994))

celmap %>% dplyr::filter(str_detect(inhouse, "_"))
#   inhouse            GSE47994   inhouse2
# 1 A1764-009_a.CEL.gz GSM1164303 A1764-009_a
# 2 A1764-040_a.CEL.gz GSM1164334 A1764-040_a
# 3 A1764-009_b.CEL.gz GSM8011739 A1764-009_b
# 4 A1764-040_b.CEL.gz GSM8011740 A1764-040_b

# ==============================================================================



# 3. Format inhouse data for GEO, make it congruent to GSE47994, and write out
# ==============================================================================

# Rename matrix samples with GEO accession
finher_expr_norm1$rma$sum <- format_finher_matrix_for_geo(x = finher_expr_norm1$rma$sum, celmap)
finher_expr_norm1$rma$combat <- format_finher_matrix_for_geo(x = finher_expr_norm1$rma$combat, celmap)

finher_expr_norm1$mas$sum <- format_finher_matrix_for_geo(x = finher_expr_norm1$mas$sum, celmap)
finher_expr_norm1$mas$combat <- format_finher_matrix_for_geo(x = finher_expr_norm1$mas$combat, celmap)



finher_expr_norm2$rma$sum <- format_finher_matrix_for_geo(x = finher_expr_norm2$rma$sum, celmap)
finher_expr_norm2$rma$combat <- format_finher_matrix_for_geo(x = finher_expr_norm2$rma$combat, celmap)

finher_expr_norm2$mas$sum <- format_finher_matrix_for_geo(x = finher_expr_norm2$mas$sum, celmap)
finher_expr_norm2$mas$combat <- format_finher_matrix_for_geo(x = finher_expr_norm2$mas$combat, celmap)


# Subset sample qc metrics and include sample GEO accession
nme = c("CEL_filename", "PP", "Neg_Avg","Pos_Avg", "GE_Avg",
        "PC1", "PC2", "PC3", "PC4",
        "PC1_Var_Prop", "PC2_Var_Prop", "PC3_Var_Prop", "PC4_Var_Prop",
        "Subject_Id", "Batch")

sample_qc_norm1$rma <- format_finher_sample_qc_for_geo(x = sample_qc_norm1$rma, nme = nme, celmap = celmap)
sample_qc_norm1$mas <- format_finher_sample_qc_for_geo(x = sample_qc_norm1$mas, nme = nme, celmap = celmap)

sample_qc_norm2$rma <- format_finher_sample_qc_for_geo(x = sample_qc_norm2$rma, nme = nme, celmap = celmap)
sample_qc_norm2$mas <- format_finher_sample_qc_for_geo(x = sample_qc_norm2$mas, nme = nme, celmap = celmap)


# Agreement between probeset data between PVAC, Sample QC, and Probeset QC dataframes
# Norm1 - RMA
x_pvac = finher_expr_norm1$rma$pvac
x_sampleqc = sample_qc_norm1$rma
x_psetqc = pset_qc_norm1$rma

x_pvac %>% glimpse()
x_sampleqc %>% glimpse()
x_psetqc %>% glimpse()

identical(x_pvac$Probeset_Id, x_psetqc$Probeset_Id) # TRUE




nme = c("Probeset_Id", "PVAC_Status", "Probeset_Type")

# ! Note: P value <0.000001 is set to zero
pset_qc_norm1$rma <- format_pset_qc_for_geo(x = pset_qc_norm1$rma, nme)
pset_qc_norm1$mas <- format_pset_qc_for_geo(x = pset_qc_norm1$mas, nme)

pset_qc_norm2$rma <- format_pset_qc_for_geo(x = pset_qc_norm2$rma, nme)
pset_qc_norm2$mas <- format_pset_qc_for_geo(x = pset_qc_norm2$mas, nme)



# Write out GEO formatted data to include as supplementary to GSE47994
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


# # Norm1 file names to write out
# 
# expr_all_sample_all_pset
# sample_qc
# pset_pvac_filtering
# pset_qc
# annot_maxvar_after_pset_filtering
# expr_maxvar_combat_after_pset_filtering

 
# # Norm2 file names to write out
# 
# expr_sample_filtered_all_pset
# sample_qc
# pset_pvac_filtering
# pset_qc
# annot_maxvar_after_sample_and_pset_filtering
# expr_maxvar_combat_after_sample_and_pset_filtering



# "GSE47994_supp_all_samples/mas/"
nme = str_c(out_data,"GSE47994_supp_all_samples/mas/")
dir.create(path = nme, recursive = TRUE, mode = "0700")

write_tsv(x= finher_expr_norm1$mas$sum,
          file = str_c(nme, "expr_all_sample_all_pset.tsv"))
write_tsv(x= sample_qc_norm1$mas,
          file = str_c(nme, "sample_qc.tsv"))
write_tsv(x= finher_expr_norm1$mas$pvac,
          file = str_c(nme, "pset_pvac_filtering.tsv"))
write_tsv(x= pset_qc_norm1$mas,
          file = str_c(nme, "pset_qc.tsv"))
write_tsv(x= finher_expr_norm1$mas$maxvar,
          file = str_c(nme, "annot_maxvar_after_pset_filtering.tsv"))
write_tsv(x= finher_expr_norm1$mas$combat,
          file = str_c(nme, "expr_maxvar_combat_after_pset_filtering.tsv"))



# "GSE47994_supp_all_samples/rma/"
nme = str_c(out_data,"GSE47994_supp_all_samples/rma/")
dir.create(path = nme, recursive = TRUE, mode = "0700")

write_tsv(x= finher_expr_norm1$rma$sum,
          file = str_c(nme, "expr_all_sample_all_pset.tsv"))
write_tsv(x= sample_qc_norm1$rma,
          file = str_c(nme, "sample_qc.tsv"))
write_tsv(x= finher_expr_norm1$rma$pvac,
          file = str_c(nme, "pset_pvac_filtering.tsv"))
write_tsv(x= pset_qc_norm1$rma,
          file = str_c(nme, "pset_qc.tsv"))
write_tsv(x= finher_expr_norm1$rma$maxvar,
          file = str_c(nme, "annot_maxvar_after_pset_filtering.tsv"))
write_tsv(x= finher_expr_norm1$rma$combat,
          file = str_c(nme, "expr_maxvar_combat_after_pset_filtering.tsv"))



# "GSE47994_supp_sample_filtered/mas/"
nme = str_c(out_data,"GSE47994_supp_sample_filtered/mas/")
dir.create(path = nme, recursive = TRUE, mode = "0700")

write_tsv(x= finher_expr_norm2$mas$sum,
          file = str_c(nme, "expr_sample_filtered_all_pset.tsv"))
write_tsv(x= sample_qc_norm2$mas,
          file = str_c(nme, "sample_qc.tsv"))
write_tsv(x= finher_expr_norm2$mas$pvac,
          file = str_c(nme, "pset_pvac_filtering.tsv"))
write_tsv(x= pset_qc_norm2$mas,
          file = str_c(nme, "pset_qc.tsv"))
write_tsv(x= finher_expr_norm2$mas$maxvar,
          file = str_c(nme, "annot_maxvar_after_sample_and_pset_filtering.tsv"))
write_tsv(x= finher_expr_norm2$mas$combat,
          file = str_c(nme, "expr_maxvar_combat_after_sample_and_pset_filtering.tsv"))



# "GSE47994_supp_sample_filtered/rma/"
nme = str_c(out_data,"GSE47994_supp_sample_filtered/rma/")
dir.create(path = nme, recursive = TRUE, mode = "0700")

write_tsv(x= finher_expr_norm2$rma$sum,
          file = str_c(nme, "expr_sample_filtered_all_pset.tsv"))
write_tsv(x= sample_qc_norm2$rma,
          file = str_c(nme, "sample_qc.tsv"))
write_tsv(x= finher_expr_norm2$rma$pvac,
          file = str_c(nme, "pset_pvac_filtering.tsv"))
write_tsv(x= pset_qc_norm2$rma,
          file = str_c(nme, "pset_qc.tsv"))
write_tsv(x= finher_expr_norm2$rma$maxvar,
          file = str_c(nme, "annot_maxvar_after_sample_and_pset_filtering.tsv"))
write_tsv(x= finher_expr_norm2$rma$combat,
          file = str_c(nme, "expr_maxvar_combat_after_sample_and_pset_filtering.tsv"))


# ==============================================================================


