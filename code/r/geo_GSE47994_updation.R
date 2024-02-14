# geo_GSE47994_updation.R

# Why this script?
# >>>>>>>>>>>>>>>>
#
# The original FinHER data available from GSE47994 
# 1) did not include the two technical replicates
# 2) did not include experimental metadata
# 3) did not considered degraded FFPE RNA specific data processing protocols.
# The original FinHER data available from GSE47994 is RMA normalized dataset
# without sample and probeset filtering.
#
#
# This script does the following
# 1) identify the two samples from GSE47994 for which technical replicates
# are missing.
# 2) prepare the RMA reprocessed expression matrix with replicates to update GSE47994
#    (Note that the manuscript recommends MAS normalization for expression profiles
#    from degraded RNA. Updation of GSE47994 with RMA normalization is for backward
#    compatibility with studies already used the previous RMA normalized expression
#    matrix without replicates)
# 3) prepare metadata matrix to update GSE47994



# Script structure
# ================
#
# 1. Load data
# 2. Map inhouse sample id to existing FinHER data in GSE47994
# 3. Identify replicates from GSE47994
# 4. Write out tsv files of RMA normalized data matrix with replicates and
#    - experiment metadata


# 1. Load data
# ==============================================================================

# FinHER inhouse - All TN and HER2
load("results/data/finher_expr_norm1.RData")
inhouse <- finher_expr_norm1$rma$sum
dim(inhouse) # 49386   338


# FinHER GSE47994 - All TN and HER2 (without replicates!)
GSE47994_expr <- read_tsv("data/geo_existing_finher_data/GSE47994_expr.csv")
GSE47994_expr <- GSE47994_expr %>% dplyr::rename(Probeset_Id = ID_REF)
dim(GSE47994_expr) # 49386   336


# FinHER inhouse - All metadata
load("results/data/finher_assay.RData")
inhouse_metadata <- finher_assay
dim(inhouse_metadata) # 337  19

# ==============================================================================



# 2. Map inhouse sample id to existing FinHER data in GSE47994
# ==============================================================================


# inhouse to GSE47994 co-correlation
cor_inhouse_GSE47994 = cor(
  inhouse %>% dplyr::select(-Probeset_Id),
  tibble(Probeset_Id = inhouse$Probeset_Id) %>%
    dplyr::left_join(GSE47994_expr, by = "Probeset_Id") %>%
    dplyr::select(-Probeset_Id)
)
dim(cor_inhouse_GSE47994) # 337 335; inhouse on rows & gse47994 on columns

# mapping
celmap <- tibble(inhouse = names(inhouse %>% dplyr::select(-Probeset_Id)))

# GSE47994 to inhouse map
idx = purrr::map_int(
      as_tibble(cor_inhouse_GSE47994),
      ~(which(.x == max(.x))) # index of inhouse sample highly correlated GSE47994
    )
xmap = tibble(GSE47994 = colnames(cor_inhouse_GSE47994),
              inhouse = rownames(cor_inhouse_GSE47994)[idx])
celmap <- celmap %>% dplyr::left_join(xmap, by = "inhouse")

# verify mapping
x <- na.omit(celmap)
cor_x = cor(
  inhouse %>%
    dplyr::select(x$inhouse),
  tibble(Probeset_Id = inhouse$Probeset_Id) %>%
    dplyr::left_join(GSE47994_expr, by = "Probeset_Id") %>%
    dplyr::select(x$GSE47994)
)

diag(cor_x) %>% summary()
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 1       1       1       1       1       1

# Note that the original GSE47994 contains RMA normalized dataset before
# - sample/probeset filtering and without experimental batch correction.


# Double check whether GSE47994 contains data before batch correction
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 
# Combat correction
# ! In finher_expr_norm2$mas$combat, ComBat is done on maxvar after pset filtering 
xx <-  inhouse[,-1] %>%
  as.matrix()
rownames(xx) <- inhouse$Probeset_Id

batch <- finher_assay$Batch
names(batch) <- str_c(finher_assay$CEL_filename, ".CEL.gz")
batch <- batch[colnames(xx)]

identical(colnames(xx), names(batch)) # TRUE

inhouse2 <- sva::ComBat(
  dat = xx,
  batch = batch,
  mean.only = FALSE # default is FALSE
) %>%
  as_tibble(rownames = "Probeset_Id")


identical(names(inhouse), names(inhouse2)) # TRUE

# verify mapping with combat corrected inhouse
x <- na.omit(celmap)
cor_x2 = cor(
  inhouse2 %>%
    dplyr::select(x$inhouse),
  tibble(Probeset_Id = inhouse2$Probeset_Id) %>%
    dplyr::left_join(GSE47994_expr, by = "Probeset_Id") %>%
    dplyr::select(x$GSE47994)
)

diag(cor_x2) %>% summary()
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.9825  0.9981  0.9993  0.9969  0.9995  0.9996 

# Reduced correlation of GSE47994 with combat corrected inhouse data 
# -comapred to perfect correlation in non-corrected inhouse data.

rm(inhouse2)


# ==============================================================================




# 3. Identify samples having replicate from GSE47994
# ==============================================================================

# Replicates
celmap %>% dplyr::filter(str_detect(inhouse, "_"))

#   inhouse            GSE47994
# 1 A1764-009_a.CEL.gz GSM1164303
# 2 A1764-040_a.CEL.gz GSM1164334
# 3 A1764-009_b.CEL.gz NA
# 4 A1764-040_b.CEL.gz NA

# Note that there is no replicates in GSE47994


# Give new name for replicates
celmap <- celmap %>%
  dplyr::mutate(GSE47994 = purrr::map2_chr(
    inhouse, GSE47994,
    ~switch(.x,
            "A1764-009_b.CEL.gz" = "rep_GSM1164303",
            "A1764-040_b.CEL.gz" = "rep_GSM1164334",
            .y)
  ))

celmap <- celmap %>%
  dplyr::mutate(inhouse2 = str_replace(string = inhouse, pattern = ".CEL.gz", replacement = ""))


# Replicates
celmap %>% dplyr::filter(str_detect(inhouse, "_"))

#   inhouse            GSE47994       inhouse2
# 1 A1764-009_a.CEL.gz GSM1164303     A1764-009_a
# 2 A1764-040_a.CEL.gz GSM1164334     A1764-040_a
# 3 A1764-009_b.CEL.gz rep_GSM1164303 A1764-009_b
# 4 A1764-040_b.CEL.gz rep_GSM1164334 A1764-040_b

write_tsv(celmap, file = str_c(out_data,"GSE47994_to_inhouse_sample_map.tsv"))


# ==============================================================================



# 4. Write out tsv files of RMA renormalized data matrix and experiment metadata
# ==============================================================================

# Keeping the RMA normalized data in GSE47994 to make it backward compatible
# or in other words to make third party reanalysis using GSE47994 intact.

inhouse %>% dim() # 49386   338; RMA normalized
inhouse_metadata %>% dim() # 337  19


# Format RMA re-normalized data matrix for GSE47994

identical(names(inhouse)[-1], celmap$inhouse) # TRUE

inhouse_x <- inhouse %>%
  dplyr::rename_with(
    .fn = function(x){
      purrr::map_chr(x, function(x, celmap){
        celmap$GSE47994[which(x==celmap$inhouse)]
      },
      celmap = celmap)
    },
  .cols = celmap$inhouse)

identical(names(inhouse_x)[-1], celmap$GSE47994) # TRUE # verification


# Format experiment metadata for GSE47994

identical(inhouse_metadata$CEL_filename, celmap$inhouse2) # TRUE

inhouse_metadata_x <- celmap %>%
  dplyr::rename(Sample_id = "GSE47994",
                Filename = "inhouse2") %>%
  dplyr::select(-inhouse) %>%
  dplyr::left_join(
    inhouse_metadata %>%
      dplyr::rename(Filename = "CEL_filename"),
    by = "Filename"
  )

# to make it congruent to existing GSE47994
inhouse_metadata_x <- inhouse_metadata_x %>%
  dplyr::mutate(Sample_title = str_c("Breast Cancer ", Subject_Id),
                CEL_file = str_c(Sample_id, "_", Subject_Id, ".CEL.gz"))


identical(inhouse_metadata_x$Sample_id, names(inhouse_x)[-1]) # TRUE


# Write out in original order as in GSE47994

inhouse_x <- inhouse_x %>%
  dplyr::rename(ID_REF = "Probeset_Id") %>%
  dplyr::select(ID_REF,
                names(GSE47994_expr)[-1],
                "rep_GSM1164303", "rep_GSM1164334")

inhouse_metadata_x <- tibble(Sample_id = names(inhouse_x)[-1]) %>%
  dplyr::left_join(inhouse_metadata_x, by = "Sample_id") %>%
  dplyr::rename_with(.fn = function(x){str_c("characteristics: ", x)})

identical(names(inhouse_x)[-1], inhouse_metadata_x$`characteristics: Sample_id`) # TRUE

# inhouse_metadata_x[,1:2] %>% as.data.frame()

# write_tsv(x=inhouse_x %>% dplyr::mutate(across(-"ID_REF",~round(.x, digits = 5))),
#           file = str_c(out_data,"GSE47994_rma_renormalized_with_replicates_rounded.tsv"))
# Samples were missing due to manual conversion of tsv to excel format.
# Hence writing directly in excel format.
# Another replicate only matrix is also generated as requested by GEO team during the update process.

write_xlsx(
  x = inhouse_x %>%
    dplyr::mutate(across(-"ID_REF",~round(.x, digits = 5))),
  path = str_c(out_data,"GSE47994_rma_renormalized_with_replicates.xlsx"),
  format_headers = FALSE,
  use_zip64 = FALSE
)

write_xlsx(
  x = inhouse_x %>%
    dplyr::select("ID_REF", "rep_GSM1164303", "rep_GSM1164334") %>%
    dplyr::mutate(across(-"ID_REF",~round(.x, digits = 5))),
  path = str_c(out_data,"GSE47994_rma_renormalized_replicates_only.xlsx"),
  format_headers = FALSE,
  use_zip64 = FALSE
) # GEO team asked for expression data only for missing replicates !!!!!!!!!!!!

write_tsv(x=inhouse_metadata_x,
          file = str_c(out_data,"GSE47994_metadata_updated.tsv"))


# ==============================================================================




