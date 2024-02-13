# technical_validation.R

# Objective
# >>>>>>>>>
#
# Biological and technical validation of FinHER dataset


# Script structure
# >>>>>>>>>>>>>>>>
#
# Load data
# IHC vs PAM50 subtyping
# ER/HER2 IHC status prediction (AUC)
# Technical replicate comparison



# ==============================================================================
# Load data
# ==============================================================================

load("results/data/finher_expr_norm2.RData")
load("results/data/finher_assay.RData")
# load("data/finher_clin.RData") # Legally protected data !!!

# FinHER cleaned expression
finher_expr <- list()

finher_expr[["rma"]] <- finher_expr_norm2$rma$maxvar %>%
  dplyr::left_join(finher_expr_norm2$rma$combat, by = "Probeset_Id") %>%
  dplyr::select(-Probeset_Id, - Gene_Symbol) %>%
  dplyr::rename(Ncbi_gene_id = "Entrez_Gene") %>%
  dplyr::mutate(Ncbi_gene_id = str_c("ncbi_", Ncbi_gene_id))

finher_expr[["mas"]] <- finher_expr_norm2$mas$maxvar %>%
  dplyr::left_join(finher_expr_norm2$mas$combat, by = "Probeset_Id") %>%
  dplyr::select(-Probeset_Id, - Gene_Symbol) %>%
  dplyr::rename(Ncbi_gene_id = "Entrez_Gene") %>%
  dplyr::mutate(Ncbi_gene_id = str_c("ncbi_", Ncbi_gene_id))


# Integrated clinical data + assay data
finher_clin <- tibble(CEL_filename = names(finher_expr$rma)[-1]) %>% # For sample ordering
  dplyr::left_join(finher_assay %>%
                     dplyr::mutate(CEL_filename = str_c(CEL_filename, ".CEL.gz")),
                   by = "CEL_filename") %>%
  dplyr::left_join(finher_clin, by = "Subject_Id")

identical(finher_clin$CEL_filename, names(finher_expr$rma)[-1]) # TRUE


# PAM50 centroids
pam50 <- list()
pam50[["centroid"]] <- read_tsv(file = "data/pam50/pam50_centroids.txt")
pam50[["annot"]] <- read_tsv(file = "data/pam50/pam50_annotation.txt")

pam50 <- pam50$centroid %>%
  dplyr::rename(pcrID = "...1") %>%
  dplyr::left_join(pam50$annot, by = "pcrID") %>%
  dplyr::mutate(EntrezGene = str_c("ncbi_", EntrezGene)) %>%
  dplyr::rename(Ncbi_gene_id = "EntrezGene")
nrow(pam50) # 50


# ER/HER2 genes
"ncbi_2099" %in% finher_expr$rma$Ncbi_gene_id # ER: TRUE
"ncbi_2064" %in% finher_expr$rma$Ncbi_gene_id # HER2: TRUE

"ncbi_2099" %in% finher_expr$mas$Ncbi_gene_id # ER: TRUE
"ncbi_2064" %in% finher_expr$mas$Ncbi_gene_id # HER2: TRUE

module_list <- list(
  Er_gene = tibble(
    Ncbi_gene_id = "ncbi_2099",
    Hugo_gene_symbol = "ESR1",
    Direction = 1
  ),
  Her2_gene = tibble(
    Ncbi_gene_id = "ncbi_2064",
    Hugo_gene_symbol = "ERBB2",
    Direction = 1
  )
)



#
# ==============================================================================




# ==============================================================================
# IHC vs PAM50 subtyping
# ==============================================================================

finher_clin %>% glimpse()


# IHC subtyping
# >>>>>>>>>>>>>

# Note !!!!
# IHC subtyping did not depend on molecular data


finher_clin %>%
  dplyr::group_by(ER_IHC, PR_IHC, HR_IHC, HER2_IHC_CISH, Subtype_IHC) %>%
  dplyr::summarise(N = n())
#   ER_IHC   PR_IHC   HR_IHC   HER2_IHC_CISH Subtype_IHC     N
# 1 Negative Negative Negative Negative      TN            120
# 2 Negative Negative Negative Positive      HR-HER2+       89
# 3 Negative Positive Positive Positive      HR+HER2+        4
# 4 Positive Negative Positive Positive      HR+HER2+       27
# 5 Positive Positive Positive Positive      HR+HER2+       60


finher_clin <- finher_clin %>%
  dplyr::mutate(Subtype_IHC_2 = if_else(str_detect(Subtype_IHC, "HER2+"),
                                        "HER2", Subtype_IHC) %>%
                  str_replace(pattern = "HR+HER2-", replacement = "HR")
  )

finher_clin %>%
  dplyr::group_by(Subtype_IHC, Subtype_IHC_2) %>%
  dplyr::summarise(N = n())
#   Subtype_IHC Subtype_IHC_2     N
# 1 HR-HER2+    HER2             89
# 2 HR+HER2+    HER2             91
# 3 TN          TN              120


# PAM50 subtyping
# >>>>>>>>>>>>>>>

# RMA
xx <- get_pam50_subtype(x = finher_expr$rma, pam50 = pam50) %>%
  dplyr::mutate(
    Pam50_subtype_2 = purrr::map_chr(
      Pam50_subtype,
      function(x){
        case_when(
          x == "Basal" ~ "TN",
          x == "Her2" ~ "HER2",
          x == "LumB" ~ "HR",
          x == "LumA" ~ "HR",
          TRUE ~ NA_character_ # Normals will be set as NA
        )
      }
    ),
    Pam50_subtype_3 = purrr::map_chr(
      Pam50_subtype,
      function(x){
        case_when(
          x == "Basal" ~ "TN",
          x == "Her2" ~ "HR-HER2+",
          x == "LumB" ~ "HR+HER2+",
          x == "LumA" ~ "HR+HER2-",
          TRUE ~ NA_character_ # Normals will be set as NA
        )
      }
    )
  ) %>%
  dplyr::rename_if(
    startsWith(names(.), "Pam50"),
    ~(str_replace(string = .x, pattern = "Pam50", replacement = "Rma_pam50"))
  )
# The below code will also work
# dplyr::rename_all(
#   ~(str_replace(string = .x, pattern = "Pam50", replacement = "Rma_pam50"))
# )

# "25 pam50 genes were missing in the dataset"


# Merge with finher_clin
finher_clin <- finher_clin %>%
  dplyr::left_join(xx, by = "CEL_filename")


# MAS
xx <- get_pam50_subtype(x = finher_expr$mas, pam50 = pam50) %>%
  dplyr::mutate(
    Pam50_subtype_2 = purrr::map_chr(
      Pam50_subtype,
      function(x){
        case_when(
          x == "Basal" ~ "TN",
          x == "Her2" ~ "HER2",
          x == "LumB" ~ "HR",
          x == "LumA" ~ "HR",
          TRUE ~ NA_character_ # Normals will be set as NA
        )
      }
    ),
    Pam50_subtype_3 = purrr::map_chr(
      Pam50_subtype,
      function(x){
        case_when(
          x == "Basal" ~ "TN",
          x == "Her2" ~ "HR-HER2+",
          x == "LumB" ~ "HR+HER2+",
          x == "LumA" ~ "HR+HER2-",
          TRUE ~ NA_character_ # Normals will be set as NA
        )
      }
    )
  ) %>%
  dplyr::rename_if(
    startsWith(names(.), "Pam50"),
    ~(str_replace(string = .x, pattern = "Pam50", replacement = "Mas_pam50"))
  )
# The below code will also work
# dplyr::rename_all(
#   ~(str_replace(string = .x, pattern = "Pam50", replacement = "Rma_pam50"))
# )
 
#  "25 pam50 genes were missing in the dataset"

# Merge with finher_clin
finher_clin <- finher_clin %>%
  dplyr::left_join(xx, by = "CEL_filename")


# Version 1
# >>>>>>>>>

finher_clin %>%
  dplyr::group_by(Subtype_IHC_2, Rma_pam50_subtype_2, Mas_pam50_subtype_2) %>%
  dplyr::summarise(N = n())
#   Subtype_IHC_2 Rma_pam50_subtype_2 Mas_pam50_subtype_2     N
# 1 HER2          HER2                HER2                   57
# 2 HER2          HER2                HR                      1
# 3 HER2          HER2                TN                      1
# 4 HER2          HER2                NA                      1
# 5 HER2          HR                  HER2                    7
# 6 HER2          HR                  HR                     75
# 7 HER2          HR                  NA                      1
# 8 HER2          TN                  TN                      5
# 9 HER2          NA                  HR                      3
# 10 HER2          NA                  TN                      1
# 11 HER2          NA                  NA                     28
# 12 TN            HER2                HER2                    6
# 13 TN            HR                  HR                     21
# 14 TN            HR                  NA                      1
# 15 TN            TN                  TN                     63
# 16 TN            TN                  NA                      7
# 17 TN            NA                  TN                      1
# 18 TN            NA                  NA                     21

kappa2(ratings = finher_clin[,c("Subtype_IHC_2", "Rma_pam50_subtype_2")])
# Cohen's Kappa for 2 Raters (Weights: unweighted)
#  Subjects = 246
#    Raters = 2
# Kappa = 0.342
# z = 10.4
# p-value = 0

kappa2(ratings = finher_clin[,c("Subtype_IHC_2", "Mas_pam50_subtype_2")])
#  Cohen's Kappa for 2 Raters (Weights: unweighted)
# Subjects = 246
# Raters = 2
# Kappa = 0.338
# z = 10
# p-value = 0


# Version 2 !!!!!
# >>>>>>>>>>>>>>>

finher_clin %>%
  dplyr::group_by(Subtype_IHC, Rma_pam50_subtype_3, Mas_pam50_subtype_3) %>%
  dplyr::summarise(N = n())
# Subtype_IHC Rma_pam50_subtype_3 Mas_pam50_subtype_3     N
# 1 HR-HER2+    HR-HER2+            HR-HER2+               49
# 2 HR-HER2+    HR-HER2+            TN                      1
# 3 HR-HER2+    HR+HER2-            HR-HER2+                1
# 4 HR-HER2+    HR+HER2-            HR+HER2-               10
# 5 HR-HER2+    HR+HER2+            HR-HER2+                2
# 6 HR-HER2+    HR+HER2+            HR+HER2+                1
# 7 HR-HER2+    TN                  TN                      5
# 8 HR-HER2+    NA                  HR+HER2-                1
# 9 HR-HER2+    NA                  TN                      1
# 10 HR-HER2+    NA                  NA                     18
# # … with 20 more rows

kappa2(ratings = finher_clin[,c("Subtype_IHC", "Rma_pam50_subtype_3")])
# Cohen's Kappa for 2 Raters (Weights: unweighted)
# Subjects = 246
# Raters = 2
# Kappa = 0.456
# z = 14.1
# p-value = 0

kappa2(ratings = finher_clin[,c("Subtype_IHC", "Mas_pam50_subtype_3")])
#  Cohen's Kappa for 2 Raters (Weights: unweighted)
# Subjects = 240
# Raters = 2
# Kappa = 0.445
# z = 13.6
# p-value = 0

#
# ==============================================================================



# ==============================================================================
# ER/HER2 IHC status prediction (AUC) from ER/HER2 gene
# ==============================================================================


# RMA
# >>>

# converting x to genes on column and samples on rows
x = t_tibble(finher_expr$rma , names_x_desc = "CEL_filename")
score <- get_module_score(x = x, module_list = module_list, by = "Ncbi_gene_id") %>%
  dplyr::rename_if(
    str_detect(names(.), "gene"),
    ~(str_c("Rma_", .x %>% str_to_lower()))
  )
# x: A tibble with genes on column and samples on rows.
#     1st column is considered as sample names.

finher_clin <- finher_clin %>%
  left_join(score, by = "CEL_filename")

pROC::roc(formula = ER_IHC  ~ Rma_er_gene, data = finher_clin)$auc
# Area under the curve: 0.896
pROC::roc(formula = HER2_IHC_CISH  ~ Rma_her2_gene, data = finher_clin)$auc
# Area under the curve: 0.9279



# MAS
# >>>

# converting x to genes on column and samples on rows
x = t_tibble(finher_expr$mas , names_x_desc = "CEL_filename")
score <- get_module_score(x = x, module_list = module_list, by = "Ncbi_gene_id") %>%
  dplyr::rename_if(
    str_detect(names(.), "gene"),
    ~(str_c("Mas_", .x %>% str_to_lower()))
  )
# x: A tibble with genes on column and samples on rows.
#     1st column is considered as sample names.

finher_clin <- finher_clin %>%
  left_join(score, by = "CEL_filename")


pROC::roc(formula = ER_IHC  ~ Mas_er_gene, data = finher_clin)$auc
# Area under the curve: 0.8664
pROC::roc(formula = HER2_IHC_CISH  ~ Mas_her2_gene, data = finher_clin)$auc
# Area under the curve: 0.9264


#
# ==============================================================================



# ==============================================================================
# Technical replicate comparison
# ==============================================================================

load("results/data/sample_qc_norm1.RData")
load("results/data/finher_expr_norm1.RData")

replicates <- c("A1764-009_a","A1764-009_b",
                "A1764-040_a","A1764-040_b")

# Sample QC metrics
sample_qc_norm1$rma %>%
  dplyr::filter(CEL_filename %in% replicates) %>%
  dplyr::arrange(CEL_filename) %>%
  dplyr::select(CEL_filename, PP, Neg_Avg, Pos_Avg, GE_Avg)
#   CEL_filename    PP Neg_Avg Pos_Avg GE_Avg
# 1 A1764-009_a   53.1    4.54    7.98   4.32
# 2 A1764-009_b   55.3    4.47    8.05   4.33
# 3 A1764-040_a   37.2    4.84    7.34   4.30
# 4 A1764-040_b   45.0    4.73    7.40   4.30
sample_qc_norm1$mas %>%
  dplyr::filter(CEL_filename %in% replicates) %>%
  dplyr::arrange(CEL_filename) %>%
  dplyr::select(CEL_filename, PP, Neg_Avg, Pos_Avg, GE_Avg)
#   CEL_filename    PP Neg_Avg Pos_Avg GE_Avg
# 1 A1764-009_a   53.1    5.69    8.83   5.99
# 2 A1764-009_b   55.3    5.50    8.82   5.90
# 3 A1764-040_a   37.2    6.00    8.12   5.81
# 4 A1764-040_b   45.0    5.72    8.06   5.71


# Expression correlation
# >>>>>>>>>>>>>>>>>>>>>>

# RMA
x <- finher_expr_norm1$rma$combat %>%
  dplyr::select(Probeset_Id, str_c(replicates,".CEL.gz"))


cor(x[,-1])
#                    A1764-009_a.CEL.gz A1764-009_b.CEL.gz A1764-040_a.CEL.gz A1764-040_b.CEL.gz
# A1764-009_a.CEL.gz          1.0000000          0.9938489          0.7133188          0.7093989
# A1764-009_b.CEL.gz          0.9938489          1.0000000          0.7096301          0.7064587
# A1764-040_a.CEL.gz          0.7133188          0.7096301          1.0000000          0.9824050
# A1764-040_b.CEL.gz          0.7093989          0.7064587          0.9824050          1.0000000


# MAS
x <- finher_expr_norm1$mas$combat %>%
  dplyr::select(Probeset_Id, str_c(replicates,".CEL.gz"))


cor(x[,-1])
#                    A1764-009_a.CEL.gz A1764-009_b.CEL.gz A1764-040_a.CEL.gz A1764-040_b.CEL.gz
# A1764-009_a.CEL.gz          1.0000000          0.9907463          0.6784602          0.6757887
# A1764-009_b.CEL.gz          0.9907463          1.0000000          0.6738596          0.6710225
# A1764-040_a.CEL.gz          0.6784602          0.6738596          1.0000000          0.9728422
# A1764-040_b.CEL.gz          0.6757887          0.6710225          0.9728422          1.0000000

#
# ==============================================================================
