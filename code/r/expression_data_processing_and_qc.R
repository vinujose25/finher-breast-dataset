# expression_data_processing.R


# What the script does?
# >>>>>>>>>>>>>>>>>>>>>
#
# This script loads raw expression profiles of HER2 and TN samples from the
# FinHER trial into an AffyBatch object. FinHER samples' expression profiles
# were already published (GSE47994 - TN+HER2; GSE65095 - HER2 only - redundant!)
# but without stringent sample and probeset QC.
# The present script uses the inhouse version (identical to published version)
# of raw CEL files, instead of the published versions.
# Mapping of inhouse-id to geo-sample-accession from GSE47994 has been done
# in geo_GSE47994_updation.R.
#
# This script normalizes the raw expression data using
# RMA and MAS algorithms. This is to check the validity of both methods to
# normalize degraded FFPE expression profiles. Plots are generated.
#
# Normalization is done before sample filtering and after filtering of samples
# with percent-present cutoff value of 10%.
# Samples with less than 10% percent-present values were filtered out.
#


# Script structure
# >>>>>>>>>>>>>>>>
#
# 1. Load non-expression data.
# 2. Load raw expression data and compute percent-present calls from raw data.
# 3. Normalization level-1 (before sample filtering).
# 4. QC level-1 (before sample filtering).
# 5. Sample QC metric distribution, PCA plot, Probeset QC plot
# 6. Sample filtering
# 7. Normalization level-2 (after sample filtering).
# 8. QC level-2 (after sample filtering).
# 9. Normalization QC plot
# 10. RLE plot



# 1. Load non-expression data.
# ==============================================================================

load("data/hgu219probe.gcIndexList.RData") # !!!! include function to prepare this data frame
load("data/hgu219_annot36.RData")
load("data/hgu133plus2.hgu219_genomicPosControl.RData")

hgu219_annot36 <- hgu219_annot36 %>%
  tibble::as_tibble() %>%
  dplyr::select(Probe.Set.ID, Gene.Symbol, Entrez.Gene) %>%
  dplyr::rename(Probeset_Id = "Probe.Set.ID",
         Gene_Symbol = "Gene.Symbol",
         Entrez_Gene = "Entrez.Gene")


load("results/data/finher_assay.RData") # n = 337, contains 2 technical replicates

#
# ==============================================================================




# 2. Load raw expression data and compute percent-present calls from raw data.
# ==============================================================================


# Hydra module used: "module load R-bundle-Bioconductor/3.10-foss-2019b"

path_files  <- paste0("data/cel/loi2013/", unique(finher_assay$Batch))
# "finher_assay$Batch" contains tilte case values (Batch1, Bacth2 etc)
# While path file needs lower case
path_files <- str_to_lower(path_files)

name_files <- lapply(unique(finher_assay$Batch),
                     function(b, finher_assay){
                       finher_assay$CEL_filename[finher_assay$Batch == b]  %>%
                         paste0(".CEL.gz")
                     },
                     finher_assay
)


xx <- map2(path_files,
           name_files,
           function(x, y){
             ReadAffy(celfile.path = x,
                      filenames = y,
                      compress = T)
           }
)


finher_expr_raw <- list() # to keep intermediate files

finher_expr_raw[["raw"]] <- merge.AffyBatch(xx[[1]], xx[[2]])
finher_expr_raw[["raw"]] <- merge.AffyBatch(finher_expr_raw[["raw"]], xx[[3]])
finher_expr_raw[["raw"]] <- merge.AffyBatch(finher_expr_raw[["raw"]], xx[[4]])
rm(xx)
gc()


# Percent-present (pp) values
finher_expr_raw[["pp"]] <- getu219dcal(aBatch = finher_expr_raw[["raw"]],
                                       gcIndexList = hgu219probe.gcIndexList,
                                       type = "medGC")
# Converting pp to tibble
xx <- as_tibble(finher_expr_raw$pp)
xx <- bind_cols(tibble(Probeset_Id = rownames(finher_expr_raw$pp)), xx)
finher_expr_raw$pp <- xx
# Note that pp values for 23 antigenomic probesets were missing,
#  as they were used to estimate background. Hence, pp values for only
#  49386-23 = 49363 probesets were present.


identical(sampleNames(finher_expr_raw$raw), names(finher_expr_raw$pp)[-1])# TRUE

# # Save RData
# save(finher_expr_raw, file=paste0(out_data,"finher_expr_raw.RData"))


#
# ==============================================================================




# 3. Normalization level-1 (before sample filtering).
# ==============================================================================

finher_expr_norm1 <- vector(mode = "list", length = 2)
names(finher_expr_norm1) <- c("rma", "mas")


# RMA
# >>>

finher_expr_norm1$rma[["norm"]] <- finher_expr_raw[["raw"]] %>%
  bg.correct(method = "rma") %>%
  normalize.AffyBatch.quantiles(type = "pmonly")


finher_expr_norm1$rma[["sum"]] <- rma(finher_expr_raw[["raw"]])
xx <- as_tibble(exprs(finher_expr_norm1$rma[["sum"]]))
xx <- bind_cols(tibble(Probeset_Id = rownames(exprs(finher_expr_norm1$rma[["sum"]]))), xx)
finher_expr_norm1$rma[["sum"]] <- xx



# PVAC probeset filtering
# (Hydra module used: "module load R-bundle-Bioconductor/3.10-foss-2019b")
probe_index_list <- lapply(hgu219probe.gcIndexList, function(x){ x$indexIntensity })
probe_index_list <- probe_index_list[hgu219_annot36$Probeset_Id]
# Explicitly use negative anti-genomic probesets.
# Anti-genomic probeset name has "GC" in it.
neg_probesets=hgu219_annot36$Probeset_Id[grepl(hgu219_annot36$Probeset_Id,pattern="GC",fixed=T)]
# Use of negative probesets is simpler compared to use of all absent probesets.
# All absent probesets method requires detection call to be generated.
# With PM only array detection calling is not straight forward,
# hence better to use affymetrix antigenomic probsets as negative probesets.
finher_expr_norm1$rma[["pvac"]] <- pvacFilter(
  abatch = finher_expr_norm1$rma[["norm"]], # probe level data
  probeIndexList=probe_index_list,
  dCal = NULL,
  neg.probeSets = neg_probesets
)
xx <- as_tibble(finher_expr_norm1$rma[["pvac"]]$pvac.df)
xx$pc1.varProp_cutoff <- finher_expr_norm1$rma[["pvac"]]$pvac.cutoff
xx <- xx %>% rename(Probeset_Id = "probeSet.id")
finher_expr_norm1$rma[["pvac"]]  <- xx


# Maxvar probeset collapsing
identical(finher_expr_norm1$rma$sum$Probeset_Id,
          finher_expr_norm1$rma$pvac$Probeset_Id) # TRUE
finher_expr_norm1$rma[["maxvar"]] <-  get_max_var_annot(
  ge = finher_expr_norm1$rma$sum[finher_expr_norm1$rma$pvac$Is_Selected_Probeset, ],
  xannot = hgu219_annot36,
  pset_colname = "Probeset_Id",
  gene_colname = "Entrez_Gene"
)


# Combat
xx <-  finher_expr_norm1$rma[["maxvar"]] %>%
  left_join(finher_expr_norm1$rma[["sum"]], by = "Probeset_Id") %>%
  dplyr::select(-c(1:3)) %>%
  as.matrix()
rownames(xx) <- finher_expr_norm1$rma[["maxvar"]]$Probeset_Id

batch <- finher_assay$Batch
names(batch) <- str_c(finher_assay$CEL_filename, ".CEL.gz")
batch <- batch[colnames(xx)]

identical(colnames(xx), names(batch)) # TRUE

finher_expr_norm1$rma[["combat"]] <- sva::ComBat(
  dat = xx,
  batch = batch,
  mean.only = FALSE # default is FALSE
) %>%
  as_tibble(rownames = "Probeset_Id")




# MAS
# >>>

finher_expr_norm1$mas[["norm"]] <- finher_expr_raw[["raw"]] %>%
  bg.correct(method = "mas") %>%
  normalize.AffyBatch.constant()


finher_expr_norm1$mas[["sum"]] <- expresso(
  afbatch = finher_expr_raw[["raw"]],
  bgcorrect.method = "mas",
  normalize.method = "constant",
  pmcorrect.method = "pmonly",
  summary.method = "mas"
)
xx <- as_tibble(log2(exprs(finher_expr_norm1$mas[["sum"]]) + 1))
xx <- bind_cols(tibble(Probeset_Id = rownames(exprs(finher_expr_norm1$mas[["sum"]]))), xx)
finher_expr_norm1$mas[["sum"]] <- xx


# PVAC probeset filtering
# (Hydra module used: "module load R-bundle-Bioconductor/3.10-foss-2019b")
probe_index_list <- lapply(hgu219probe.gcIndexList, function(x){ x$indexIntensity })
probe_index_list <- probe_index_list[hgu219_annot36$Probeset_Id]
# Explicitly use negative anti-genomic probesets.
# Anti-genomic probeset name has "GC" in it.
neg_probesets=hgu219_annot36$Probeset_Id[grepl(hgu219_annot36$Probeset_Id,pattern="GC",fixed=T)]
# Use of negative probesets is simpler compared to use of all absent probesets.
# All absent probesets method requires detection call to be generated.
# With PM only array detection calling is not straight forward,
# hence better to use affymetrix antigenomic probesets as negative probesets.
finher_expr_norm1$mas[["pvac"]] <- pvacFilter(
  abatch = finher_expr_norm1$mas[["norm"]], # probe level data
  probeIndexList=probe_index_list,
  dCal = NULL,
  neg.probeSets = neg_probesets
)
xx <- as_tibble(finher_expr_norm1$mas[["pvac"]]$pvac.df)
xx$pc1.varProp_cutoff <- finher_expr_norm1$mas[["pvac"]]$pvac.cutoff
xx <- xx %>% rename(Probeset_Id = "probeSet.id")
finher_expr_norm1$mas[["pvac"]]  <- xx


# Maxvar probeset collapsing
identical(finher_expr_norm1$mas$sum$Probeset_Id,
          finher_expr_norm1$mas$pvac$Probeset_Id) # TRUE
finher_expr_norm1$mas[["maxvar"]] <-  get_max_var_annot(
  ge = finher_expr_norm1$mas$sum[finher_expr_norm1$mas$pvac$Is_Selected_Probeset, ],
  xannot = hgu219_annot36,
  pset_colname = "Probeset_Id",
  gene_colname = "Entrez_Gene"
)


# Combat
xx <-  finher_expr_norm1$mas[["maxvar"]] %>%
  left_join(finher_expr_norm1$mas[["sum"]], by = "Probeset_Id") %>%
  dplyr::select(-c(1:3)) %>%
  as.matrix()
rownames(xx) <- finher_expr_norm1$mas[["maxvar"]]$Probeset_Id

batch <- finher_assay$Batch
names(batch) <- str_c(finher_assay$CEL_filename, ".CEL.gz")
batch <- batch[colnames(xx)]

identical(colnames(xx), names(batch)) # TRUE

finher_expr_norm1$mas[["combat"]] <- sva::ComBat(
  dat = xx,
  batch = batch,
  mean.only = FALSE # default is FALSE
) %>%
  as_tibble(rownames = "Probeset_Id")




# Discarding AffyBatch object of normalized expression data at probe level
finher_expr_norm1$rma <- finher_expr_norm1$rma[c("sum", "pvac", "maxvar", "combat")]
finher_expr_norm1$mas <- finher_expr_norm1$mas[c("sum", "pvac", "maxvar", "combat")]

# Impose variable naming convention (except for sample names); Xyz_Abc_Pqr
names(finher_expr_raw$pp) # Following variable naming convention
names(finher_expr_norm1$rma$sum) # Following variable naming convention
names(finher_expr_norm1$rma$pvac) # Not following variable naming convention

finher_expr_norm1$rma$pvac <- finher_expr_norm1$rma$pvac %>%
  rename(
    "Is_Neg_Probeset" = "isNeg.probeSet",
    "Is_Selected_Probeset" = "isSelected.probeSet",
    "PC1_Var" = "pc1.var",
    "PC2_Var" = "pc2.var",
    "Total_Var" = "total.var",
    "PC1_Var_Prop" = "pc1.varProp",
    "PC2_Var_Prop" = "pc2.varProp",
    "PC1_Var_Prop_Cutoff" = "pc1.varProp_cutoff"
  )
finher_expr_norm1$mas$pvac <- finher_expr_norm1$mas$pvac %>%
  rename(
    "Is_Neg_Probeset" = "isNeg.probeSet",
    "Is_Selected_Probeset" = "isSelected.probeSet",
    "PC1_Var" = "pc1.var",
    "PC2_Var" = "pc2.var",
    "Total_Var" = "total.var",
    "PC1_Var_Prop" = "pc1.varProp",
    "PC2_Var_Prop" = "pc2.varProp",
    "PC1_Var_Prop_Cutoff" = "pc1.varProp_cutoff"
  )



# # Save RData @@@@@@@@
# save(finher_expr_norm1, file = paste0(out_data, "finher_expr_norm1.RData"))

#
# ==============================================================================




# 4. QC level-1 (before sample filtering).
# ==============================================================================

# QC of finher_expr_norm1
sample_qc_norm1 <- vector(mode = "list", length = 2)
pset_qc_norm1 <- vector(mode = "list", length = 2)
names(sample_qc_norm1) <- c("rma", "mas")
names(pset_qc_norm1) <- c("rma", "mas")


# Sample QC plots
# >>>>>>>>>>>>>>>

# load(paste0(out_data,"finher_expr_raw.RData"))
# load(paste0(out_data,"finher_expr_norm1.RData"))


# RMA
sample_qc_norm1[["rma"]] <- get_sample_qc_data(
  ge = finher_expr_norm1$rma$sum,
  pp = finher_expr_raw$pp,
  neg_pset = grep(finher_expr_norm1$rma$sum$Probeset_Id,
                  pattern = "GC",
                  fixed = T,
                  value = T),
  pos_pset =  hgu133plus2.hgu219_genomicPosControl$hgu219
)


# MAS
sample_qc_norm1[["mas"]] <- get_sample_qc_data(
  ge = finher_expr_norm1$mas$sum,
  pp = finher_expr_raw$pp,
  neg_pset = grep(finher_expr_norm1$mas$sum$Probeset_Id,
                  pattern = "GC",
                  fixed = T,
                  value = T),
  pos_pset =  hgu133plus2.hgu219_genomicPosControl$hgu219
)


# Integrate sample qc data with assay and clin
sample_qc_norm1[["rma"]] <- sample_qc_norm1[["rma"]] %>%
  dplyr::mutate(CEL_filename = str_replace(CEL_filename, "\\.CEL\\.gz", "")) %>%
  left_join(finher_assay %>% select(CEL_filename, Subject_Id, Batch),
            by = "CEL_filename") 
# %>%
#   left_join(finher_clin, by = "Subject_Id")

sample_qc_norm1[["mas"]] <- sample_qc_norm1[["mas"]] %>%
  dplyr::mutate(CEL_filename = str_replace(CEL_filename, "\\.CEL\\.gz", "")) %>%
  left_join(finher_assay %>% select(CEL_filename, Subject_Id, Batch),
            by = "CEL_filename") 
# %>%
#   left_join(finher_clin, by = "Subject_Id")



# # Save RData
# save(sample_qc_norm1, file = paste0(out_data,"sample_qc_norm1.RData"))


# Plotting

load(paste0(out_data,"sample_qc_norm1.RData"))

get_sample_qc_plots(
  qc_data = sample_qc_norm1$rma,
  outdir = out_figures,
  filename = "Level1_RMA"
)

get_sample_qc_plots(
  qc_data = sample_qc_norm1$mas,
  outdir = out_figures,
  filename = "Level1_MAS"
)




# Probeset QC plots
# >>>>>>>>>>>>>>>>>

# RMA
pset_qc_norm1[["rma"]] <- get_pset_qc_data(
  ge = finher_expr_norm1$rma$sum,
  pp = finher_expr_raw$pp,
  pvac = finher_expr_norm1$rma$pvac,
  pos_pset =  hgu133plus2.hgu219_genomicPosControl$hgu219
)

# MAS
pset_qc_norm1[["mas"]] <- get_pset_qc_data(
  ge = finher_expr_norm1$mas$sum,
  pp = finher_expr_raw$pp,
  pvac = finher_expr_norm1$mas$pvac,
  pos_pset =  hgu133plus2.hgu219_genomicPosControl$hgu219
)


# # Save RData
# save(pset_qc_norm1, file = paste0(out_data,"pset_qc_norm1.RData"))


# plotting
load(paste0(out_data,"pset_qc_norm1.RData"))

get_pset_qc_plots(
  qc_data = pset_qc_norm1$rma,
  outdir = out_figures,
  filename = "Level1_RMA"
)

get_pset_qc_plots(
  qc_data = pset_qc_norm1$mas,
  outdir = out_figures,
  filename = "Level1_MAS"
)

#
# ==============================================================================




# 5. Sample QC metric distribution, PCA plot, Probeset QC plot
# ==============================================================================


# MAS
# >>>>>>>>>>

# Pos, Neg and Genomic average distribution

p1 <- sample_qc_norm1$mas %>%
  dplyr::mutate(Neg = Neg_Avg,
                Pos = Pos_Avg) %>%
  tidyr::gather("Key", "Value", "Neg_Avg", "GE_Avg", "Pos_Avg") %>%
  dplyr::mutate(
    Key = purrr::map(
      Key,
      ~switch (.x,
               Neg_Avg = "Nagative\ncontrol",
               GE_Avg = "Genomic\nprobeset",
               Pos_Avg = "Positive\ncontrol")
    ),
    Key = factor(Key, levels = c("Nagative\ncontrol",
                                 "Positive\ncontrol",
                                 "Genomic\nprobeset"))
  ) %>%
  ggplot(aes(x= Key, y=Value)) +
  geom_line(aes(color = PP > 10, linetype = Neg < Pos, group = CEL_filename))+
  geom_boxplot(fill = NA) +
  scale_linetype_manual(values = c("TRUE" = "solid", "FALSE" = "dashed")) +
  scale_color_manual(values = c("TRUE" = "gray50", "FALSE" = "goldenrod")) +
  theme_bw() +
  labs(y = "Mean expression") +
  theme(axis.title.x = element_blank()) +
  theme(panel.grid = element_blank())



# Percent-present distribution

p2 <- sample_qc_norm1$mas %>%
  ggplot(aes(x= "FinHER", y = PP)) +
  geom_boxplot(fill = NA) +
  geom_jitter(aes(color = PP > 10), shape = 1, width = .25, show.legend = FALSE) +
  scale_color_manual(values = c("TRUE" = "gray50", "FALSE" = "goldenrod")) +
  theme_bw() +
  labs(y = "Percent-Present (PP)") +
  theme(axis.title.x = element_blank()) +
  theme(panel.grid = element_blank())



# PCA plot

p3 <- sample_qc_norm1$mas %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(color = PP > 10), shape = 1) +
  geom_hline(yintercept = 0, color = "gray20") +
  geom_vline(xintercept = 0, color = "gray20") +
  scale_color_manual(values = c("TRUE" = "gray50", "FALSE" = "goldenrod")) +
  theme_bw() +
  labs(x = str_c("PC1 (",unique(sample_qc_norm1$mas$PC1_Var_Prop) %>% round(digits = 1),"%)"),
       y = str_c("PC2 (",unique(sample_qc_norm1$mas$PC2_Var_Prop) %>% round(digits = 1),"%)")) +
  theme(panel.grid = element_blank())



# Probeset QC plot

# https://www.researchgate.net/figure/Mean-variance-relationship-Here-we-show-the-sample-variance-across-lanes-in-the-liver_fig2_24282377
p4 <- pset_qc_norm1$mas %>%
  ggplot(aes(x = Pset_GE_Avg, y = Pset_GE_Var)) +
  geom_point(aes(color = PVAC_Status), shape = 1) +
  geom_smooth(aes(linetype = PVAC_Status), color = "gray20", method = "lm", se = FALSE) +
  theme_bw() +
  labs(x = "Probeset mean", y = "Probeset varaiance") +
  scale_linetype_manual(values = c("Selected" = "solid", "Discarded" = "dashed")) +
  scale_color_manual(values = c("Selected" = "gray50", "Discarded" = "goldenrod")) +
  guides(linetype = "none", color = guide_legend(title = "Probeset\nPVAC status")) +
  theme(panel.grid = element_blank())


# Integrated plot
 
pdf(file = paste0(out_figures, "sample_pset_qc_mas.pdf"),
    width = 7.5,
    height = 5,
    title = "Sample and Probeset QC")

ggpubr::ggarrange(
  
  ggpubr::ggarrange( #ggdplyr::arrange(
    p1, p3 + guides(color = FALSE) ,
    labels = c("a","b"),
    ncol = 2, nrow = 1, widths = c(.6,.4)#, vjust = 1.1
  ),
  
  ggpubr::ggarrange( #ggdplyr::arrange(
    p2, p4,
    labels = c("c","d"),
    ncol = 2, nrow = 1, widths = c(.3,.7)#, vjust = 1.1
  ),
  
  ncol = 1, nrow = 2
)

dev.off()



# RMA
# >>>>>>>>>>>>>>>


# Pos, Neg and Genomic average distribution

p1 <- sample_qc_norm1$rma %>%
  dplyr::mutate(Neg = Neg_Avg,
                Pos = Pos_Avg) %>%
  tidyr::gather("Key", "Value", "Neg_Avg", "GE_Avg", "Pos_Avg") %>%
  dplyr::mutate(
    Key = purrr::map(
      Key,
      ~switch (.x,
               Neg_Avg = "Nagative\ncontrol",
               GE_Avg = "Genomic\nprobeset",
               Pos_Avg = "Positive\ncontrol")
    ),
    Key = factor(Key, levels = c("Nagative\ncontrol",
                                 "Positive\ncontrol",
                                 "Genomic\nprobeset"))
  ) %>%
  ggplot(aes(x= Key, y=Value)) +
  geom_line(aes(color = PP > 10, linetype = Neg < Pos, group = CEL_filename))+
  geom_boxplot(fill = NA) +
  scale_linetype_manual(values = c("TRUE" = "solid", "FALSE" = "dashed")) +
  scale_color_manual(values = c("TRUE" = "gray50", "FALSE" = "goldenrod")) +
  theme_bw() +
  labs(y = "Mean expression") +
  theme(axis.title.x = element_blank()) +
  theme(panel.grid = element_blank())



# Percent-present distribution

p2 <- sample_qc_norm1$rma %>%
  ggplot(aes(x= "FinHER", y = PP)) +
  geom_boxplot(fill = NA) +
  geom_jitter(aes(color = PP > 10), shape = 1, width = .25, show.legend = FALSE) +
  scale_color_manual(values = c("TRUE" = "gray50", "FALSE" = "goldenrod")) +
  theme_bw() +
  labs(y = "Percent-Present (PP)") +
  theme(axis.title.x = element_blank()) +
  theme(panel.grid = element_blank())



# PCA plot

p3 <- sample_qc_norm1$rma %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(color = PP > 10), shape = 1) +
  geom_hline(yintercept = 0, color = "gray20") +
  geom_vline(xintercept = 0, color = "gray20") +
  scale_color_manual(values = c("TRUE" = "gray50", "FALSE" = "goldenrod")) +
  theme_bw() +
  labs(x = str_c("PC1 (",unique(sample_qc_norm1$mas$PC1_Var_Prop) %>% round(digits = 1),"%)"),
       y = str_c("PC2 (",unique(sample_qc_norm1$mas$PC2_Var_Prop) %>% round(digits = 1),"%)")) +
  theme(panel.grid = element_blank())



# Probeset QC plot

# https://www.researchgate.net/figure/Mean-variance-relationship-Here-we-show-the-sample-variance-across-lanes-in-the-liver_fig2_24282377
p4 <- pset_qc_norm1$rma %>%
  ggplot(aes(x = Pset_GE_Avg, y = Pset_GE_Var)) +
  geom_point(aes(color = PVAC_Status), shape = 1) +
  geom_smooth(aes(linetype = PVAC_Status), color = "gray20", method = "lm", se = FALSE) +
  theme_bw() +
  labs(x = "Probeset mean", y = "Probeset varaiance") +
  scale_linetype_manual(values = c("Selected" = "solid", "Discarded" = "dashed")) +
  scale_color_manual(values = c("Selected" = "gray50", "Discarded" = "goldenrod")) +
  guides(linetype = "none", color = guide_legend(title = "Probeset\nPVAC status")) +
  theme(panel.grid = element_blank())


# Integrated plot

pdf(file = paste0(out_figures, "sample_pset_qc_rma.pdf"),
    width = 7.5,
    height = 5,
    title = "Sample and Probeset QC")

ggpubr::ggarrange(
  
  ggpubr::ggarrange( #ggdplyr::arrange(
    p1, p3 + guides(color = FALSE) ,
    labels = c("a","b"),
    ncol = 2, nrow = 1, widths = c(.6,.4)#, vjust = 1.1
  ),
  
  ggpubr::ggarrange( #ggdplyr::arrange(
    p2, p4,
    labels = c("c","d"),
    ncol = 2, nrow = 1, widths = c(.3,.7)#, vjust = 1.1
  ),
  
  ncol = 1, nrow = 2
)

dev.off()


#
# ==============================================================================




# 6. Sample filtering
# ==============================================================================

# Sample filtering (PP > 10), replicates removal ("A1764-009_b", "A1764-040_b")
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

sample_filtered <- sample_qc_norm1$rma$CEL_filename[sample_qc_norm1$rma$PP <= 10]
sample_filtered <- c(sample_filtered, "A1764-009_b", "A1764-040_b") %>%
  str_c(".CEL.gz")

# [1] "A1764-010.CEL.gz"   "A1764-011.CEL.gz"   "A1764-014.CEL.gz"
# [4] "A1764-035.CEL.gz"   "A1764-039.CEL.gz"   "A1764-052.CEL.gz"
# [7] "A1764-055.CEL.gz"   "A1764-064.CEL.gz"   "A1764-083.CEL.gz"
# [10] "A1764-086.CEL.gz"   "A1764-088.CEL.gz"   "A1764-101.CEL.gz"
# [13] "A1764-145.CEL.gz"   "A1764-151.CEL.gz"   "A1764-158.CEL.gz"
# [16] "A1764-176.CEL.gz"   "A1764-195.CEL.gz"   "A1764-201.CEL.gz"
# [19] "A1764-204.CEL.gz"   "A1764-219.CEL.gz"   "A1764-467.CEL.gz"
# [22] "A1764-470.CEL.gz"   "A1764-473.CEL.gz"   "A1764-475.CEL.gz"
# [25] "A1764-485.CEL.gz"   "A1764-488.CEL.gz"   "A1764-489.CEL.gz"
# [28] "A1764-490.CEL.gz"   "A1764-491.CEL.gz"   "A1764-502.CEL.gz"
# [31] "A1764-512.CEL.gz"   "A1764-526.CEL.gz"   "A1764-528.CEL.gz"
# [34] "A1764-529.CEL.gz"   "A1764-531.CEL.gz"   "A1764-009_b.CEL.gz"
# [37] "A1764-040_b.CEL.gz"

#
# ==============================================================================




# 7. Normalization level-2 (after sample filtering).
# ==============================================================================

finher_expr_norm2 <- vector(mode = "list", length = 2)
names(finher_expr_norm2) <- c("rma", "mas")

# Sample filtering index
index <- !(sampleNames(finher_expr_raw$raw) %in% sample_filtered)
table(index)
# FALSE  TRUE
# 37   300


# RMA
# >>>

finher_expr_norm2$rma[["norm"]] <- finher_expr_raw[["raw"]][ , index] %>%
  bg.correct(method = "rma") %>%
  normalize.AffyBatch.quantiles(type = "pmonly")


finher_expr_norm2$rma[["sum"]] <- rma(finher_expr_raw[["raw"]][ , index])
xx <- as_tibble(exprs(finher_expr_norm2$rma[["sum"]]))
xx <- bind_cols(tibble(Probeset_Id = rownames(exprs(finher_expr_norm2$rma[["sum"]]))), xx)
finher_expr_norm2$rma[["sum"]] <- xx



# PVAC probeset filtering
# (Hydra module used: "module load R-bundle-Bioconductor/3.10-foss-2019b")
probe_index_list <- lapply(hgu219probe.gcIndexList, function(x){ x$indexIntensity })
probe_index_list <- probe_index_list[hgu219_annot36$Probeset_Id]
# Explicitly use negative anti-genomic probesets.
# Anti-genomic probeset name has "GC" in it.
neg_probesets=hgu219_annot36$Probeset_Id[grepl(hgu219_annot36$Probeset_Id,pattern="GC",fixed=T)]
# Use of negative probesets is simpler compared to use of all absent probsets.
# All absent probsets method requires detection call to be generated.
# With PM only array detection calling is not straight forward,
# hence better to use affymetrix antigenomic probsets as negative probesets.
finher_expr_norm2$rma[["pvac"]] <- pvacFilter(
  abatch = finher_expr_norm2$rma[["norm"]], # probe level data
  probeIndexList=probe_index_list,
  dCal = NULL,
  neg.probeSets = neg_probesets
)
xx <- as_tibble(finher_expr_norm2$rma[["pvac"]]$pvac.df)
xx$pc1.varProp_cutoff <- finher_expr_norm2$rma[["pvac"]]$pvac.cutoff
xx <- xx %>% rename(Probeset_Id = "probeSet.id")
finher_expr_norm2$rma[["pvac"]]  <- xx



# Maxvar probeset collapsing
identical(finher_expr_norm2$rma$sum$Probeset_Id,
          finher_expr_norm2$rma$pvac$Probeset_Id) # TRUE
finher_expr_norm2$rma[["maxvar"]] <-  get_max_var_annot(
  ge = finher_expr_norm2$rma$sum[finher_expr_norm2$rma$pvac$Is_Selected_Probeset, ],
  xannot = hgu219_annot36,
  pset_colname = "Probeset_Id",
  gene_colname = "Entrez_Gene"
)


# Combat
xx <-  finher_expr_norm2$rma[["maxvar"]] %>%
  left_join(finher_expr_norm2$rma[["sum"]], by = "Probeset_Id") %>%
  dplyr::select(-c(1:3)) %>%
  as.matrix()
rownames(xx) <- finher_expr_norm2$rma[["maxvar"]]$Probeset_Id

batch <- finher_assay$Batch
names(batch) <- str_c(finher_assay$CEL_filename, ".CEL.gz")
batch <- batch[colnames(xx)]

identical(colnames(xx), names(batch)) # TRUE

finher_expr_norm2$rma[["combat"]] <- sva::ComBat(
  dat = xx,
  batch = batch,
  mean.only = FALSE # default is FALSE
) %>%
  as_tibble(rownames = "Probeset_Id")



# MAS
# >>>

finher_expr_norm2$mas[["norm"]] <- finher_expr_raw[["raw"]][ , index] %>%
  bg.correct(method = "mas") %>%
  normalize.AffyBatch.constant()


finher_expr_norm2$mas[["sum"]] <- expresso(
  afbatch = finher_expr_raw[["raw"]][ , index],
  bgcorrect.method = "mas",
  normalize.method = "constant",
  pmcorrect.method = "pmonly",
  summary.method = "mas"
)
xx <- as_tibble(log2(exprs(finher_expr_norm2$mas[["sum"]]) + 1))
xx <- bind_cols(tibble(Probeset_Id = rownames(exprs(finher_expr_norm2$mas[["sum"]]))), xx)
finher_expr_norm2$mas[["sum"]] <- xx


# PVAC probeset filtering
# (Hydra module used: "module load R-bundle-Bioconductor/3.10-foss-2019b")
probe_index_list <- lapply(hgu219probe.gcIndexList, function(x){ x$indexIntensity })
probe_index_list <- probe_index_list[hgu219_annot36$Probeset_Id]
# Explicitly use negative anti-genomic probesets.
# Anti-genomic probeset name gas "GC" in it.
neg_probesets=hgu219_annot36$Probeset_Id[grepl(hgu219_annot36$Probeset_Id,pattern="GC",fixed=T)]
# Use of negative probesets is simpler compared to use of all absent probsets.
# All absent probsets method requires detection call to be generated.
# With PM only array detection calling is not straight forward,
# hence better to use affymetrix antigenomic probsets as negative probesets.
finher_expr_norm2$mas[["pvac"]] <- pvacFilter(
  abatch = finher_expr_norm2$mas[["norm"]], # probe level data
  probeIndexList=probe_index_list,
  dCal = NULL,
  neg.probeSets = neg_probesets
)
xx <- as_tibble(finher_expr_norm2$mas[["pvac"]]$pvac.df)
xx$pc1.varProp_cutoff <- finher_expr_norm2$mas[["pvac"]]$pvac.cutoff
xx <- xx %>% rename(Probeset_Id = "probeSet.id")
finher_expr_norm2$mas[["pvac"]]  <- xx



# Maxvar probeset collapsing
identical(finher_expr_norm2$mas$sum$Probeset_Id,
          finher_expr_norm2$mas$pvac$Probeset_Id) # TRUE
finher_expr_norm2$mas[["maxvar"]] <-  get_max_var_annot(
  ge = finher_expr_norm2$mas$sum[finher_expr_norm2$mas$pvac$Is_Selected_Probeset, ],
  xannot = hgu219_annot36,
  pset_colname = "Probeset_Id",
  gene_colname = "Entrez_Gene"
)


# Combat
xx <-  finher_expr_norm2$mas[["maxvar"]] %>%
  left_join(finher_expr_norm2$mas[["sum"]], by = "Probeset_Id") %>%
  dplyr::select(-c(1:3)) %>%
  as.matrix()
rownames(xx) <- finher_expr_norm2$mas[["maxvar"]]$Probeset_Id

batch <- finher_assay$Batch
names(batch) <- str_c(finher_assay$CEL_filename, ".CEL.gz")
batch <- batch[colnames(xx)]

identical(colnames(xx), names(batch)) # TRUE

finher_expr_norm2$mas[["combat"]] <- sva::ComBat(
  dat = xx,
  batch = batch,
  mean.only = FALSE # default is FALSE
) %>%
  as_tibble(rownames = "Probeset_Id")



# Discarding AffyBatch object of normalized expression data at probe level
finher_expr_norm2$rma <- finher_expr_norm2$rma[c("sum", "pvac", "maxvar", "combat")]
finher_expr_norm2$mas <- finher_expr_norm2$mas[c("sum", "pvac", "maxvar", "combat")]

# Impose variable naming convension (except for sample names); Xyz_Abc_Pqr
names(finher_expr_raw$pp) # Following variable naming convention
names(finher_expr_norm2$rma$sum) # Following variable naming convention
names(finher_expr_norm2$rma$pvac) # Not following variable naming convention

finher_expr_norm2$rma$pvac <- finher_expr_norm2$rma$pvac %>%
  rename(
    "Is_Neg_Probeset" = "isNeg.probeSet",
    "Is_Selected_Probeset" = "isSelected.probeSet",
    "PC1_Var" = "pc1.var",
    "PC2_Var" = "pc2.var",
    "Total_Var" = "total.var",
    "PC1_Var_Prop" = "pc1.varProp",
    "PC2_Var_Prop" = "pc2.varProp",
    "PC1_Var_Prop_Cutoff" = "pc1.varProp_cutoff"
  )
finher_expr_norm2$mas$pvac <- finher_expr_norm2$mas$pvac %>%
  rename(
    "Is_Neg_Probeset" = "isNeg.probeSet",
    "Is_Selected_Probeset" = "isSelected.probeSet",
    "PC1_Var" = "pc1.var",
    "PC2_Var" = "pc2.var",
    "Total_Var" = "total.var",
    "PC1_Var_Prop" = "pc1.varProp",
    "PC2_Var_Prop" = "pc2.varProp",
    "PC1_Var_Prop_Cutoff" = "pc1.varProp_cutoff"
  )


# # Save RData @@@@@@@@
# save(finher_expr_norm2, file = paste0(out_data, "finher_expr_norm2.RData"))

#
# ==============================================================================




# 8. QC level-2 (after sample filtering).
# ==============================================================================

# QC of finher_expr_norm2
sample_qc_norm2 <- vector(mode = "list", length = 2)
pset_qc_norm2 <- vector(mode = "list", length = 2)
names(sample_qc_norm2) <- c("rma", "mas")
names(pset_qc_norm2) <- c("rma", "mas")


# Sample QC plots
# >>>>>>>>>>>>>>>

load(paste0(out_data,"finher_expr_raw.RData"))
load(paste0(out_data,"finher_expr_norm2.RData"))


# RMA
sample_qc_norm2[["rma"]] <- get_sample_qc_data(
  ge = finher_expr_norm2$rma$sum,
  pp = finher_expr_raw$pp[, names(finher_expr_norm2$rma$sum)], # due to sample filtering
  neg_pset = grep(finher_expr_norm2$rma$sum$Probeset_Id,
                  pattern = "GC",
                  fixed = T,
                  value = T),
  pos_pset =  hgu133plus2.hgu219_genomicPosControl$hgu219
)


# MAS
sample_qc_norm2[["mas"]] <- get_sample_qc_data(
  ge = finher_expr_norm2$mas$sum,
  pp = finher_expr_raw$pp[, names(finher_expr_norm2$mas$sum)], # due to sample filtering
  neg_pset = grep(finher_expr_norm2$mas$sum$Probeset_Id,
                  pattern = "GC",
                  fixed = T,
                  value = T),
  pos_pset =  hgu133plus2.hgu219_genomicPosControl$hgu219
)


# Integrate sample qc data with assay and clin
sample_qc_norm2[["rma"]] <- sample_qc_norm2[["rma"]] %>%
  dplyr::mutate(CEL_filename = str_replace(CEL_filename, "\\.CEL\\.gz", "")) %>%
  left_join(finher_assay %>% select(CEL_filename, Subject_Id, Batch),
            by = "CEL_filename") 
# %>%
#   left_join(finher_clin, by = "Subject_Id")

sample_qc_norm2[["mas"]] <- sample_qc_norm2[["mas"]] %>%
  dplyr::mutate(CEL_filename = str_replace(CEL_filename, "\\.CEL\\.gz", "")) %>%
  left_join(finher_assay %>% select(CEL_filename, Subject_Id, Batch),
            by = "CEL_filename") 
# %>%
#   left_join(finher_clin, by = "Subject_Id")



# # Save RData
# save(sample_qc_norm2, file = paste0(out_data,"sample_qc_norm2.RData"))


# Plotting

load(paste0(out_data,"sample_qc_norm2.RData"))

get_sample_qc_plots(
  qc_data = sample_qc_norm2$rma,
  outdir = out_figures,
  filename = "Level2_RMA"
)

get_sample_qc_plots(
  qc_data = sample_qc_norm2$mas,
  outdir = out_figures,
  filename = "Level2_MAS"
)




# Probeset QC plots
# >>>>>>>>>>>>>>>>>

# RMA
pset_qc_norm2[["rma"]] <- get_pset_qc_data(
  ge = finher_expr_norm2$rma$sum,
  pp = finher_expr_raw$pp[, names(finher_expr_norm2$rma$sum)], # due to sample filtering
  pvac = finher_expr_norm2$rma$pvac,
  pos_pset =  hgu133plus2.hgu219_genomicPosControl$hgu219
)

# MAS
pset_qc_norm2[["mas"]] <- get_pset_qc_data(
  ge = finher_expr_norm2$mas$sum,
  pp = finher_expr_raw$pp[, names(finher_expr_norm2$mas$sum)], # due to sample filtering
  pvac = finher_expr_norm2$mas$pvac,
  pos_pset =  hgu133plus2.hgu219_genomicPosControl$hgu219
)


# # Save RData
# save(pset_qc_norm2, file = paste0(out_data,"pset_qc_norm2.RData"))


# plotting
load(paste0(out_data,"pset_qc_norm2.RData"))

get_pset_qc_plots(
  qc_data = pset_qc_norm2$rma,
  outdir = out_figures,
  filename = "Level2_RMA"
)

get_pset_qc_plots(
  qc_data = pset_qc_norm2$mas,
  outdir = out_figures,
  filename = "Level2_MAS"
)

#
# ==============================================================================




# 9. Normalization QC plot
# ==============================================================================

# MAS
# >>>

# Correaltion
p1_annot <- cor(sample_qc_norm2$mas[,"PP"], sample_qc_norm2$mas[,c("Neg_Avg", "Pos_Avg", "GE_Avg")]) %>%
  as_tibble() %>%
  tidyr::gather("Key", "Value") %>%
  dplyr::mutate(
    Key = purrr::map_chr(
      Key,
      ~(switch (.x,
                Neg_Avg = "Negative control",
                Pos_Avg = "Positive control",
                GE_Avg = "Genomic probeset"))
    ) %>%
      factor(levels = c("Negative control",  "Positive control", "Genomic probeset")),
    
    Value = str_c("Pearson: ", round(Value, digits = 2))
  )


p1 <- sample_qc_norm2$mas %>%
  dplyr::mutate(Neg = Neg_Avg,
                Pos = Pos_Avg) %>%
  tidyr::gather("Key", "Value", "Neg_Avg", "GE_Avg", "Pos_Avg") %>%
  dplyr::mutate(
    Key = purrr::map(
      Key,
      ~switch (.x,
               Neg_Avg = "Negative control",
               GE_Avg = "Genomic probeset",
               Pos_Avg = "Positive control")
    ),
    Key = factor(Key, levels = c("Negative control",
                                 "Positive control",
                                 "Genomic probeset"))
  ) %>%
  ggplot(aes(x = Value, y = PP)) +
  geom_point(shape = 1, color = "gray50") + # "gray20"
  geom_smooth(method = "lm", se = FALSE, color = "goldenrod") + # "gray20"
  geom_text(data = p1_annot,
            aes(x = Inf, y = -Inf, label = Value),
            vjust = -1.5, hjust = 1.1, color = "gray20",
            fontface = "bold.italic", size = 3) +
  theme_bw() +
  labs(x = "Mean expression", y = "Percentage") +
  # theme(axis.title = element_blank()) +
  facet_grid("Percent-present" ~ Key, scales = "free", switch = "y")


# RMA
# >>>

# Correaltion
p2_annot <- cor(sample_qc_norm2$rma[,"PP"], sample_qc_norm2$rma[,c("Neg_Avg", "Pos_Avg", "GE_Avg")]) %>%
  as_tibble() %>%
  tidyr::gather("Key", "Value") %>%
  dplyr::mutate(
    Key = purrr::map_chr(
      Key,
      ~(switch (.x,
                Neg_Avg = "Negative control",
                Pos_Avg = "Positive control",
                GE_Avg = "Genomic probeset"))
    ) %>%
      factor(levels = c("Negative control",  "Positive control", "Genomic probeset")),
    
    Value = str_c("Pearson: ", round(Value, digits = 2))
  )


p2 <- sample_qc_norm2$rma %>%
  dplyr::mutate(Neg = Neg_Avg,
                Pos = Pos_Avg) %>%
  tidyr::gather("Key", "Value", "Neg_Avg", "GE_Avg", "Pos_Avg") %>%
  dplyr::mutate(
    Key = purrr::map(
      Key,
      ~switch (.x,
               Neg_Avg = "Negative control",
               GE_Avg = "Genomic probeset",
               Pos_Avg = "Positive control")
    ),
    Key = factor(Key, levels = c("Negative control",
                                 "Positive control",
                                 "Genomic probeset"))
  ) %>%
  ggplot(aes(x = Value, y = PP)) +
  geom_point(shape = 1, color = "gray50") + # "gray20"
  geom_smooth(method = "lm", se = FALSE, color = "goldenrod") + # "gray20"
  geom_text(data = p2_annot,
            aes(x = Inf, y = -Inf, label = Value),
            vjust = -1.5, hjust = 1.1, color = "gray20",
            fontface = "bold.italic", size = 3) +
  theme_bw() +
  labs(x = "Mean expression", y = "Percentage") +
  # theme(axis.title = element_blank()) +
  facet_grid("Percent-present" ~ Key, scales = "free", switch = "y")



# Sample and probeset QC plot
pdf(file = paste0(out_figures, "norm_comparison_after_samp_filtering.pdf"),
    width = 7.5,
    height = 5,
    title = "Normalization comparison")

ggpubr::ggarrange( #ggdplyr::arrange(
  p1, p2,
  labels = c("a","b"),
  ncol = 1, nrow = 2#, widths = c(.6,.4)#, vjust = 1.1
)


dev.off()


#
# ==============================================================================



# 10. RLE plot
# ==============================================================================



# MAS full pset rle
# >>>>>>>>>>>>>>>>>

x1 <- finher_expr_norm2$mas$sum
# x2 as combat corrected x1

# Combat correction
# ! In finher_expr_norm2$mas$combat, ComBat is done on maxvar after pset filtering 
xx <-  finher_expr_norm2$mas[["sum"]][,-1] %>%
  as.matrix()
rownames(xx) <- finher_expr_norm2$mas[["sum"]]$Probeset_Id

batch <- finher_assay$Batch
names(batch) <- str_c(finher_assay$CEL_filename, ".CEL.gz")
batch <- batch[colnames(xx)]

identical(colnames(xx), names(batch)) # TRUE

x2 <- sva::ComBat(
  dat = xx,
  batch = batch,
  mean.only = FALSE # default is FALSE
) %>%
  as_tibble(rownames = "Probeset_Id")


rle1 <- get_rle(ge = x1)
rle2 <- get_rle(ge = x2)


df <-  bind_rows(
  rle1 %>%
    dplyr::select(-1) %>%
    purrr::map_dfr(~(summary(.x))) %>%
    dplyr::rename(Min = "Min.",
                  First_Quartile = "1st Qu.",
                  Third_Quartile = "3rd Qu.",
                  Max = "Max.") %>%
    dplyr::mutate(
      CEL_filename = names(rle1)[-1] %>%
        str_replace("\\.CEL\\.gz",""),
      Type = "Before ComBat"
    ) %>%
    dplyr::left_join(finher_assay, by = "CEL_filename"),
  
  rle2 %>%
    dplyr::select(-1) %>%
    purrr::map_dfr(~(summary(.x))) %>%
    dplyr::rename(Min = "Min.",
                  First_Quartile = "1st Qu.",
                  Third_Quartile = "3rd Qu.",
                  Max = "Max.") %>%
    dplyr::mutate(
      CEL_filename = names(rle2)[-1] %>%
        str_replace("\\.CEL\\.gz",""),
      Type = "After ComBat"
    ) %>%
    dplyr::left_join(finher_assay, by = "CEL_filename")
)


p1 <- df %>%
  tidyr::gather("Key", "Value", "First_Quartile", "Median", "Third_Quartile") %>%
  dplyr::mutate(Type = factor(Type, levels = c("Before ComBat", "After ComBat")),
                Batch = str_replace(Batch, "Batch6", "Batch4")) %>%
  # Note the orginal profiling expresimnet contains samples from multiple datasets
  # That is the reasnon for Batch6.
  ggplot(aes(x= CEL_filename, y = Value, group = Key)) +
  geom_line(aes(color = Batch)) +
  facet_wrap(facets = ~ Type, nrow = 2, ncol = 1, scales = "free") +
  labs(x = "Samples", y = "Relative log expression (RLE)") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom")


# Sample and probeset QC plot
pdf(file = paste0(out_figures, "rle_mas_all_pset.pdf"),
    width = 7.5,
    height = 5,
    title = "RLE plot")

print(p1)

dev.off()



# MAS PVAC pset rle
# >>>>>>>>>>>>>>>>>

x1 <- finher_expr_norm2$mas$maxvar %>%
  dplyr::select("Probeset_Id") %>%
  dplyr::left_join(finher_expr_norm2$mas$sum, by = "Probeset_Id")
x2 <- finher_expr_norm2$mas$combat

rle1 <- get_rle(ge = x1)
rle2 <- get_rle(ge = x2)


df <-  bind_rows(
  rle1 %>%
    dplyr::select(-1) %>%
    purrr::map_dfr(~(summary(.x))) %>%
    dplyr::rename(Min = "Min.",
                  First_Quartile = "1st Qu.",
                  Third_Quartile = "3rd Qu.",
                  Max = "Max.") %>%
    dplyr::mutate(
      CEL_filename = names(rle1)[-1] %>%
        str_replace("\\.CEL\\.gz",""),
      Type = "Before ComBat"
    ) %>%
    dplyr::left_join(finher_assay, by = "CEL_filename"),
  
  rle2 %>%
    dplyr::select(-1) %>%
    purrr::map_dfr(~(summary(.x))) %>%
    dplyr::rename(Min = "Min.",
                  First_Quartile = "1st Qu.",
                  Third_Quartile = "3rd Qu.",
                  Max = "Max.") %>%
    dplyr::mutate(
      CEL_filename = names(rle2)[-1] %>%
        str_replace("\\.CEL\\.gz",""),
      Type = "After ComBat"
    ) %>%
    dplyr::left_join(finher_assay, by = "CEL_filename")
)


p1 <- df %>%
  tidyr::gather("Key", "Value", "First_Quartile", "Median", "Third_Quartile") %>%
  dplyr::mutate(Type = factor(Type, levels = c("Before ComBat", "After ComBat")),
                Batch = str_replace(Batch, "Batch6", "Batch4")) %>%
  # Note the orginal profiling expresimnet contains samples from multiple datasets
  # That is the reasnon for Batch6.
  ggplot(aes(x= CEL_filename, y = Value, group = Key)) +
  geom_line(aes(color = Batch)) +
  facet_wrap(facets = ~ Type, nrow = 2, ncol = 1, scales = "free") +
  labs(x = "Samples", y = "Relative log expression (RLE)") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom")


# Sample and probeset QC plot
pdf(file = paste0(out_figures, "rle_mas_pvac_pset.pdf"),
    width = 7.5,
    height = 5,
    title = "RLE plot")

print(p1)

dev.off()




# RMA full pset rle
# >>>>>>>>>>>>>>>>>

x1 <- finher_expr_norm2$rma$sum
# x2 as combat corrected x1

# Combat correction
# ! In finher_expr_norm2$rma$combat, ComBat is done on maxvar after pset filtering 
xx <-  finher_expr_norm2$rma[["sum"]][,-1] %>%
  as.matrix()
rownames(xx) <- finher_expr_norm2$rma[["sum"]]$Probeset_Id

batch <- finher_assay$Batch
names(batch) <- str_c(finher_assay$CEL_filename, ".CEL.gz")
batch <- batch[colnames(xx)]

identical(colnames(xx), names(batch)) # TRUE

x2 <- sva::ComBat(
  dat = xx,
  batch = batch,
  mean.only = FALSE # default is FALSE
) %>%
  as_tibble(rownames = "Probeset_Id")


rle1 <- get_rle(ge = x1)
rle2 <- get_rle(ge = x2)


df <-  bind_rows(
  rle1 %>%
    dplyr::select(-1) %>%
    purrr::map_dfr(~(summary(.x))) %>%
    dplyr::rename(Min = "Min.",
                  First_Quartile = "1st Qu.",
                  Third_Quartile = "3rd Qu.",
                  Max = "Max.") %>%
    dplyr::mutate(
      CEL_filename = names(rle1)[-1] %>%
        str_replace("\\.CEL\\.gz",""),
      Type = "Before ComBat"
    ) %>%
    dplyr::left_join(finher_assay, by = "CEL_filename"),
  
  rle2 %>%
    dplyr::select(-1) %>%
    purrr::map_dfr(~(summary(.x))) %>%
    dplyr::rename(Min = "Min.",
                  First_Quartile = "1st Qu.",
                  Third_Quartile = "3rd Qu.",
                  Max = "Max.") %>%
    dplyr::mutate(
      CEL_filename = names(rle2)[-1] %>%
        str_replace("\\.CEL\\.gz",""),
      Type = "After ComBat"
    ) %>%
    dplyr::left_join(finher_assay, by = "CEL_filename")
)


p1 <- df %>%
  tidyr::gather("Key", "Value", "First_Quartile", "Median", "Third_Quartile") %>%
  dplyr::mutate(Type = factor(Type, levels = c("Before ComBat", "After ComBat")),
                Batch = str_replace(Batch, "Batch6", "Batch4")) %>%
  # Note the orginal profiling expresimnet contains samples from multiple datasets
  # That is the reasnon for Batch6.
  ggplot(aes(x= CEL_filename, y = Value, group = Key)) +
  geom_line(aes(color = Batch)) +
  facet_wrap(facets = ~ Type, nrow = 2, ncol = 1, scales = "free") +
  labs(x = "Samples", y = "Relative log expression (RLE)") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom")


# Sample and probeset QC plot
pdf(file = paste0(out_figures, "rle_rma_all_pset.pdf"),
    width = 7.5,
    height = 5,
    title = "RLE plot")

print(p1)

dev.off()



# RMA PVAC pset rle
# >>>>>>>>>>>>>>>>>

x1 <- finher_expr_norm2$rma$maxvar %>%
  dplyr::select("Probeset_Id") %>%
  dplyr::left_join(finher_expr_norm2$rma$sum, by = "Probeset_Id")
x2 <- finher_expr_norm2$rma$combat

rle1 <- get_rle(ge = x1)
rle2 <- get_rle(ge = x2)


df <-  bind_rows(
  rle1 %>%
    dplyr::select(-1) %>%
    purrr::map_dfr(~(summary(.x))) %>%
    dplyr::rename(Min = "Min.",
                  First_Quartile = "1st Qu.",
                  Third_Quartile = "3rd Qu.",
                  Max = "Max.") %>%
    dplyr::mutate(
      CEL_filename = names(rle1)[-1] %>%
        str_replace("\\.CEL\\.gz",""),
      Type = "Before ComBat"
    ) %>%
    dplyr::left_join(finher_assay, by = "CEL_filename"),
  
  rle2 %>%
    dplyr::select(-1) %>%
    purrr::map_dfr(~(summary(.x))) %>%
    dplyr::rename(Min = "Min.",
                  First_Quartile = "1st Qu.",
                  Third_Quartile = "3rd Qu.",
                  Max = "Max.") %>%
    dplyr::mutate(
      CEL_filename = names(rle2)[-1] %>%
        str_replace("\\.CEL\\.gz",""),
      Type = "After ComBat"
    ) %>%
    dplyr::left_join(finher_assay, by = "CEL_filename")
)


p1 <- df %>%
  tidyr::gather("Key", "Value", "First_Quartile", "Median", "Third_Quartile") %>%
  dplyr::mutate(Type = factor(Type, levels = c("Before ComBat", "After ComBat")),
                Batch = str_replace(Batch, "Batch6", "Batch4")) %>%
  # Note the orginal profiling expresimnet contains samples from multiple datasets
  # That is the reasnon for Batch6.
  ggplot(aes(x= CEL_filename, y = Value, group = Key)) +
  geom_line(aes(color = Batch)) +
  facet_wrap(facets = ~ Type, nrow = 2, ncol = 1, scales = "free") +
  labs(x = "Samples", y = "Relative log expression (RLE)") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom")


# Sample and probeset QC plot
pdf(file = paste0(out_figures, "rle_rma_pvac_pset.pdf"),
    width = 7.5,
    height = 5,
    title = "RLE plot")

print(p1)

dev.off()


#
# ==============================================================================

