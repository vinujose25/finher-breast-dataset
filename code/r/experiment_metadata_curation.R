# sample_characteristics_curation.R


# What the script does?
# >>>>>>>>>>>>>>>>>>>>>
#
# This script curates experimental data associated with samples.


# Script structure
# >>>>>>>>>>>>>>>>
#
# 1. Experimental data curation
# 2. Clinical data curation



# ==============================================================================
# 1. Experimental data curation
# ==============================================================================

path_files  <- paste0("data/cel/loi2013/batch", c(1,2,3,6))

# Note on FFPE expression profiling assay by Genentech:
#
# Genentech profiles FFPE samples from FinHER, Pregnancy and Pilot study together
# in 6 batches. Batches 4 and 5 contains only Pregnancy project samples.
# Batches 2 and 3 contains only FinHER samples. Barthes 1 and 6 contains samples
# from multiple projects. Only "CEL" files of FinHER samples were considered
# in this project.
#
# Details on batches containing FinHER samples.
# 1. In addition to 87 FinHER samples, Batch1 worksheet contains details for
#   8 samples from Pilot project whose "CEL" files were removed from the
#   respective folder. Rows of these 8 samples must be removed.
# 2. Batch2 contains only 95 FinHER samples.
# 3. Batch3 contains only 95 FinHER samples.
# 4. In addition to 60 FinHER samples, Batch6 worksheet contains details for
#   8 samples from Pilot project and 4 samples from Pregnancy project
#   whose "CEL" files are not present in the respective folder.
#   Rows for these 12 samples must be removed.
# 5. All the above batches contain details of Control samples, whose "CEL" are
#   not present in the respective folder.
#
name_files <- list.files(path = path_files,
                         pattern = ".xls",
                         full.names = TRUE)
finher_assay <- purrr::map(name_files,
                             readxl::read_excel,
                             sheet = "cDNA Labeling",
                             skip = 6,
                             n_max = 97)
# Appending batch to tibble list
finher_assay <- purrr::map2(finher_assay,
                              name_files,
                              function(x, y) {
                                x$Batch <- stringr::str_split_fixed(string = y,
                                                                    pattern = "/",
                                                                    n = 5)[,4]
                                x
                              }
)


finher_assay <- dplyr::bind_rows(finher_assay)

# Column "...14" contains specific comments on few cel files, hennce removed.
finher_assay <- dplyr::select(finher_assay, -"...14")

finher_assay <- dplyr::rename(finher_assay,
                                "CEL_filename" = "A-nr:",
                                "Plate_Pos" = "Plate pos.",
                                "Subject_Id" = "Sample ID:",
                                "CRNA_Conc_ng_Per_ul" = "cRNA Conc. ng/µl",
                                "Ratio_260_280" = "Ratio 260/ 280",
                                "Ratio_260_230" = "Ratio 260/ 230",
                                "Conc_ug_Per_ul" = "Conc. ug/ul",
                                "Final_Vol_ul_1" = "Final Vol. ul...8",
                                "CRNA_Yield_ug" = "cRNA yield ug",
                                "CRNA_Used_ug" = "cRNA used µg",
                                "Vol_CRNA_Used_ul" ="Vol. of cRNA used ul",
                                "Final_Vol_ul_2" = "Final Vol. ul...12",
                                "RNase_Free_Water_ul" = "RNase-free Water ul",
                                "CDNA_Conc_ng_Per_ul" = "cDNA Conc. ng/µl",
                                "CDNA_Yield_ug" = "cDNA yield ug",
                                "CDNA_Used_ug" = "cDNA used µg",
                                "Vol_CDNA_Used_ul" = "Vol. of cDNA used ul")

# NOTE:
# The following colnames only present in metadata from batch1
# c("CRNA_Conc_ng_Per_ul", "CRNA_Yield_ug", "CRNA_Used_ug", "Vol_CRNA_Used_ul")
# All the remaining batches (i.e batch2:batch6) have the follwing column names
# c("CDNA_Conc_ng_Per_ul", "CDNA_Yield_ug", "CDNA_Used_ug", "Vol_CDNA_Used_ul")
#
# Seems like a typo from batch1, but kept as is since this is not verified with the experimenter


# "Control" samples present in each batch with same name.
# ReadAffy() will throw error if identical cel names are present.
finher_assay <- finher_assay %>% dplyr::filter(CEL_filename != "Control")
# Removing gap filling data
finher_assay <- finher_assay %>% dplyr::filter(CEL_filename != "0")

# To distinguish samples associated with different projects
finher_assay$Dataset <- "FinHER"
finher_assay$Dataset[grepl(finher_assay$Subject_Id, pattern = "FFPE" )] <- "Pilot"
finher_assay$Dataset[grepl(finher_assay$Subject_Id, pattern = "ARN GBC" )] <- "Pregnancy"

# Fixing non congruent cel file name format
finher_assay$CEL_filename[finher_assay$CEL_filename == "A1764_040_a"] <- "A1764-040_a"


# Discarding non-FinHER samples !!!!!!!!!
finher_assay$Dataset %>% table()
# FinHER     Pilot Pregnancy
# 337        16         4
finher_assay <- finher_assay %>%
  dplyr::filter(Dataset == "FinHER")




# Check variable data-types and recode variables if needed
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

finher_assay[,1:7] %>% glimpse()
finher_assay[,8:14] %>% glimpse() # Recode Batch
finher_assay[,15:ncol(finher_assay)] %>% glimpse()

# Batch
finher_assay <- finher_assay %>%
  mutate_at("Batch",~str_to_title(Batch))

dim(finher_assay)
# [1] 337  19



# Save RData !!
#
# save(finher_assay, file=paste0(out_data,"finher_assay.RData"))

#
# ==============================================================================


