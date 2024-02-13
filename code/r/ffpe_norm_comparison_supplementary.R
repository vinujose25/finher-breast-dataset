# ffpe_norm_comparison_supplementary.R


# What this script does?
# >>>>>>>>>>>>>>>>>>>>>>
#
# To compare the feasibility of scale (MAS) and quantile (RMA) normalization
# for FFPE derived degraded RNA based expression profiles.
#
# In comparable expression profiles (i.e. similar quality / percent-present values),
# both methods seems feasible. However for non-comparable expression profiles
# (i.e. variable quality / wide range percent-present values),
# the feasibility of both methods are unknown.
#
# This script compare MAS and RMA algorithm in two datasets derived from
# frozen / ffpe samples profiled on u219(PM-only) arrays
# with respect to correlation between percent-present values and
# average negative/positive/genomic expression values.
# High pp value per samples implies low negative, high positive and high genomic 
# expression in the respective sample.
# Hence for ideal normalization method there should be positive correlation
# between pp values and positive/genomic average expression and negative correlation
# between pp values and negative average expression.

# The two datasets used in this comparison are
# colon: u219 + frozen
# pancreas: u219 + ffpe


# Script structure
# >>>>>>>>>>>>>>>>
# 1. Load supporting data
# 2. Prepare clinical data
# 3. Load raw expression data (cel files) and compute percent-present (pp) values
# 4. Normalization
# 5. QC metric computation
# 6. QC plots




# ==============================================================================
# 1. Load supporting data
# ==============================================================================

load("data/hgu219probe.gcIndexList.RData") # !!!! include function to prepare this data frame
# load("data/hgu219_annot36.RData")
load("data/hgu133plus2.hgu219_genomicPosControl.RData")

#
# ==============================================================================




# ==============================================================================
# 2. Prepare clinical data
# ==============================================================================

clin <- vector(mode = "list", length = 2)
names(clin) <- c("colon", "pancreas")


# colon frozen u219, n = 98 (98 Mucosa samples from 246 paired samples)
clin$colon <- read_tsv(file = "data/GSE44076_metadata.txt") %>%
  dplyr::mutate(
    Cel_filename = str_split_fixed(
      string =  Sample_supplementary_file,
      pattern = "suppl\\/",
      n = 2)[, 2]
  ) %>%
  dplyr::filter(Sample_type2 == "Tumor")  # n = 98


# pancreas ffpe u219, n = 80 (all samples)
clin$pancreas <- read_tsv(file = "data/GSE85916_metadata.txt") %>%
  dplyr::mutate(
    Cel_filename = str_split_fixed(
      string =  Sample_supplementary_file,
      pattern = "suppl\\/",
      n = 2)[, 2]
  )


# Save RData !!!
#
# supp_clin <- clin
# save(supp_clin, file = str_c(out_data,"supp_clin.RData"))


#
# ==============================================================================



# ==============================================================================
# 3. Load raw expression data (cel files) and compute percent-present (pp) values
# ==============================================================================

expr_raw <- vector(mode = "list", length = 2)
names(expr_raw) <- c("colon", "pancreas")
expr_raw <- purrr::map(expr_raw,~(list(raw = NULL, pp = NULL)))


# Load cel files
# >>>>>>>>>>>>>>

# colon frozen u219, n = 98 (98 Mucosa samples from 246 paired samples)
expr_raw$colon$raw <- ReadAffy(
  celfile.path = "data/GSE44076_RAW/",
  filenames = supp_clin$colon$Cel_filename,
  compress = T
)

# pancreas ffpe u219, n = 80 (all samples)
expr_raw$pancreas$raw <-  ReadAffy(
  celfile.path = "data/GSE85916_RAW/",
  filenames = supp_clin$pancreas$Cel_filename,
  compress = T
)



# Detection calling
# >>>>>>>>>>>>>>>>>


# colon frozen u219, n = 98 (98 Mucosa samples from 246 paired samples)
expr_raw$colon$pp <- getu219dcal(
  aBatch = expr_raw$colon$raw,
  gcIndexList = hgu219probe.gcIndexList,
  type = "medGC"
)
# Converting pp to tibble
expr_raw$colon$pp <- as_tibble(expr_raw$colon$pp, rownames = "Probeset_Id")
# Note that pp values for 23 antigenomic probesets were missing,
#  as they were used to estimate background. Hence, pp values for only
#  49386-23 = 49363 probesets were present.
identical(sampleNames(expr_raw$colon$raw),
          names(expr_raw$colon$pp)[-1])# TRUE



# pancreas ffpe u219, n = 80 (all samples)
expr_raw$pancreas$pp <-  getu219dcal(
  aBatch = expr_raw$pancreas$raw,
  gcIndexList = hgu219probe.gcIndexList,
  type = "medGC"
)
# Converting pp to tibble
expr_raw$pancreas$pp <- as_tibble(expr_raw$pancreas$pp, rownames = "Probeset_Id")
# Note that pp values for 23 antigenomic probesets were missing,
#  as they were used to estimate background. Hence, pp values for only
#  49386-23 = 49363 probesets were present.
identical(sampleNames(expr_raw$pancreas$raw),
          names(expr_raw$pancreas$pp)[-1])# TRUE



# Save RData !!!
#
# supp_expr_raw <- expr_raw
# save(supp_expr_raw, file = str_c(out_data,"supp_expr_raw.RData"))

#
# ==============================================================================



# ==============================================================================
# 4. Normalization
# ==============================================================================

expr_norm <- vector(mode = "list", length = 2)
names(expr_norm) <- c("colon", "pancreas")
expr_norm <- purrr::map(expr_norm,~(list(mas = NULL, rma = NULL)))

# MAS
# >>>


# colon frozen u219, n = 98 (98 Mucosa samples from 246 paired samples)
expr_norm$colon$mas <- expresso(
  afbatch = supp_expr_raw$colon$raw,
  bgcorrect.method = "mas",
  normalize.method = "constant",
  pmcorrect.method = "pmonly",
  summary.method = "mas"
)
expr_norm$colon$mas <- as_tibble(log2(exprs(expr_norm$colon$mas) + 1),
                                 rownames = "Probeset_Id")


# pancreas ffpe u219, n = 80 (all samples)
expr_norm$pancreas$mas <-  expresso(
  afbatch = supp_expr_raw$pancreas$raw,
  bgcorrect.method = "mas",
  normalize.method = "constant",
  pmcorrect.method = "pmonly",
  summary.method = "mas"
)
expr_norm$pancreas$mas <- as_tibble(log2(exprs(expr_norm$pancreas$mas) + 1),
                                    rownames = "Probeset_Id")



# RMA
# >>>


# colon frozen u219, n = 98 (98 Mucosa samples from 246 paired samples)
expr_norm$colon$rma <- rma(supp_expr_raw$colon$raw)
expr_norm$colon$rma <- as_tibble(exprs(expr_norm$colon$rma),
                                 rownames = "Probeset_Id")


# pancreas ffpe u219, n = 80 (all samples)
expr_norm$pancreas$rma <-  rma(supp_expr_raw$pancreas$raw)
expr_norm$pancreas$rma <- as_tibble(exprs(expr_norm$pancreas$rma),
                                    rownames = "Probeset_Id")


# Save RData !!!
#
# supp_expr_norm <- expr_norm
# save(supp_expr_norm, file = str_c(out_data,"supp_expr_norm.RData"))



# # Id congruence checking: clin vs raw vs norm
# # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# purrr::map2(supp_expr_raw, supp_expr_norm,
#             function(x, y){
#               c(identical(names(x$pp), names(y$mas)),
#                 identical(names(x$pp), names(y$rma)) )
#             })
# # All true
# purrr::map2(supp_clin, supp_expr_norm,
#             function(x, y){
#               c(identical(x$Cel_filename, names(y$mas)[-1]),
#                 identical(x$Cel_filename, names(y$rma)[-1]) )
#             })
# # All true

#
# ==============================================================================



# ==============================================================================
# 5. QC metric computation
# ==============================================================================

norm_qc <- vector(mode = "list", length = 2)
names(norm_qc) <- c("colon", "pancreas")

norm_qc <- purrr::map(
  c("colon", "pancreas"),
  function(nme, expr_norm, expr_raw, neg_pset, pos_pset){

    print(nme)

     # colon and pancreas datasets are u219
      list(
        mas = get_norm_qc_data(ge = expr_norm[[nme]]$mas,
                               pp = expr_raw[[nme]]$pp,
                               neg_pset = neg_pset$u219,
                               pos_pset = pos_pset$hgu219),

        rma = get_norm_qc_data(ge = expr_norm[[nme]]$rma,
                               pp = expr_raw[[nme]]$pp,
                               neg_pset = neg_pset$u219,
                               pos_pset = pos_pset$hgu219)

      )
    

  },
  expr_norm = supp_expr_norm,
  expr_raw = supp_expr_raw,
  neg_pset =  list(u219 = grep(supp_expr_norm$colon$mas$Probeset_Id,
                               pattern = "GC",
                               fixed = T,
                               value = T)),
  pos_pset = hgu133plus2.hgu219_genomicPosControl)

names(norm_qc) <- c("colon", "pancreas")

# Save RData !!!
#
# supp_norm_qc <- norm_qc
# save(supp_norm_qc, file = str_c(out_data,"supp_norm_qc.RData"))

#
# ==============================================================================



# ==============================================================================
# 6. QC plots
# ==============================================================================

plot_list <- purrr::map(
  supp_norm_qc,
  function(df_list){

    purrr::map(
      df_list,
      function(df){

        # Correaltion
        plot_annot <- cor(df[,"PP"], df[,c("Neg_Avg", "Pos_Avg", "GE_Avg")]) %>%
          as_tibble() %>%
          tidyr::gather("Key", "Value") %>%
          dplyr::mutate(
            Key = purrr::map_chr(
              Key,
              ~(switch (.x,
                        Neg_Avg = "Negative control",
                        Pos_Avg = "Positive control",
                        GE_Avg = "Genomic avg."))
            ) %>%
              factor(levels = c("Negative control",  "Positive control", "Genomic avg.")),

            Value = str_c("Pearson: ", round(Value, digits = 2))
          )

        # Plot
        df %>%
          tidyr::gather(
            key = "Key", value = "Value", "Neg_Avg", "Pos_Avg", "GE_Avg") %>%
          dplyr::mutate(
            Key = purrr::map_chr(
              Key,
              ~(switch (.x,
                        Neg_Avg = "Negative control",
                        Pos_Avg = "Positive control",
                        GE_Avg = "Genomic avg."))
            ) %>%
              factor(levels = c("Negative control",  "Positive control", "Genomic avg."))
          ) %>%
          ggplot(aes(x = Value, y = PP)) +
          geom_point(shape = 1, color = "gray50") + # "gray20"
          geom_smooth(method = "lm", se = FALSE, color = "goldenrod") +
          geom_text(data = plot_annot,
                    aes(x = Inf, y = -Inf, label = Value),
                    vjust = -1.5, hjust = 1.1,
                    fontface = "bold.italic", size = 3) +
          theme_bw() +
          theme(axis.title = element_blank()) +
          facet_grid("Percent-present" ~ Key, scales = "free", switch = "y")

      }
    )
  }
)


# Convert multi-level list into single level list to use with ggarrange()
p <- purrr::map2(
  names(plot_list),
  plot_list,
  function(nme, x){
    names(x) <- str_c(nme, "_", names(x))
    x
  }
)
p <- c(p[[1]], #colon 
       p[[2]]) #pancreas



# Normalization comparison plot
pdf(file = paste0(out_figures, "supp_ffpe_norm_comparison.pdf"),
    width = 7.5,
    height = 8,
    title = "Normalization comparison")

ggpubr::ggarrange( #ggdplyr::arrange(
  plotlist = p,
  labels = c("a","b","c","d"),
  ncol = 2, nrow = 4#, widths = c(.4,.6), vjust = 1.1
)

dev.off()

#
# ==============================================================================




