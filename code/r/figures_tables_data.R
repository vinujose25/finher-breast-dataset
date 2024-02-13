# s5_figures_tables_data.R


# Objective:
# >>>>>>>>>>
#
# Final figures, tables, and data for publication


# Script structure
# >>>>>>>>>>>>>>>>
#
# Load data
# Fig1: Sample QC metric distribution, PCA plot, and Probeset QC plot
# Fig2: Normalization QC plot
# Fig3: RLE plot
# Supp.Fig1: Normalization comparison plot
# Curated clinical and expression data in table format



# ==============================================================================
# Load data
# ==============================================================================

load("results/data/sample_qc_norm1.RData")
load("results/data/sample_qc_norm2.RData")
load("results/data/pset_qc_norm1.RData")
load("results/data/pset_qc_norm2.RData")

load("results/data/finher_expr_norm1.RData")
load("results/data/finher_expr_norm2.RData")
load("results/data/finher_assay.RData")

load("results/data/supp_norm_qc.RData")

#
# ==============================================================================



# ==============================================================================
# Fig1: Sample QC metric distribution, PCA plot, Probeset QC plot
# ==============================================================================


# Pos, Neg and Genomic average distribution
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

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
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>

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
# >>>>>>>>

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
# >>>>>>>>>>>>>>>>

# https://www.researchgate.net/figure/Mean-variance-relationship-Here-we-show-the-sample-variance-across-lanes-in-the-liver_fig2_24282377
p4 <- pset_qc_norm2$mas %>%
  ggplot(aes(x = Pset_GE_Avg, y = Pset_GE_Var)) +
  geom_point(aes(color = PVAC_Status), shape = 1) +
  geom_smooth(aes(linetype = PVAC_Status), color = "gray20", method = "lm", se = FALSE) +
  theme_bw() +
  labs(x = "Probeset mean", y = "Probeset varaiance") +
  scale_linetype_manual(values = c("Selected" = "solid", "Discarded" = "dashed")) +
  scale_color_manual(values = c("Selected" = "gray50", "Discarded" = "goldenrod")) +
  guides(linetype = "none", color = guide_legend(title = "Probeset\nPVAC status")) +
  theme(panel.grid = element_blank())


# Sample and probeset QC plot
# pdf(file = paste0(out_figures, "fig1_sample_pset_qc.pdf"),
pdf(file = paste0(out_figures, "fig1_sample_pset_qc_v2.pdf"),
    width = 7.5,
    height = 5,
    title = "Figure 1: Sample QC metric distribution")

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


# png(file = paste0(out_figures, "fig1_sample_pset_qc.png"),
#     width = 7.5,
#     height = 5,
#     units = "in",
#     res = 150,
#     title = "Figure1: Sample QC metric distribution")
#
# ggpubr::ggarrange(
#
#   ggpubr::ggarrange( #ggdplyr::arrange(
#     p1, p3 + guides(color = FALSE) ,
#     labels = c("A","B"),
#     ncol = 2, nrow = 1, widths = c(.6,.4)#, vjust = 1.1
#   ),
#
#   ggpubr::ggarrange( #ggdplyr::arrange(
#     p2, p4,
#     labels = c("C","D"),
#     ncol = 2, nrow = 1, widths = c(.3,.7)#, vjust = 1.1
#   ),
#
#   ncol = 1, nrow = 2
# )
#
# dev.off()

#
# ==============================================================================



# ==============================================================================
# Fig2: Normalization QC plot
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
pdf(file = paste0(out_figures, "fig2_norm_comparison.pdf"),
    width = 7.5,
    height = 5,
    title = "Figure2: Normalization QC metric comparison")

ggpubr::ggarrange( #ggdplyr::arrange(
  p1, p2,
  labels = c("a","b"),
  ncol = 1, nrow = 2#, widths = c(.6,.4)#, vjust = 1.1
)


dev.off()


# png(file = paste0(out_figures, "fig2_norm_comparison.png"),
#      width = 7.5,
#      height = 5,
#      units = "in",
#      res = 150,
#      title = "Figure2: Normalization QC metric comparison")
#
#   ggpubr::ggarrange( #ggdplyr::arrange(
#     p1, p2,
#     labels = c("A","B"),
#     ncol = 1, nrow = 2#, widths = c(.6,.4)#, vjust = 1.1
#   )
#
# dev.off()


#
# ==============================================================================



# ==============================================================================
# Fig3: RLE plot
# ==============================================================================



# MAS full pset rle
# >>>>>>>>>>>>>>>>>

x1 <- finher_expr_norm2$mas$sum
# x2 as combat corrected x1

# Combat
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

# To test RLE plot with only pvac probesets
# identical(finher_expr_norm2$mas$sum$Probeset_Id,
#           finher_expr_norm2$mas$pvac$Probeset_Id) # TRUE
# identical(finher_expr_norm2$mas$combat$Probeset_Id,
#           finher_expr_norm2$mas$pvac$Probeset_Id) # TRUE
# index <- finher_expr_norm2$mas$pvac$Is_Selected_Probeset
# rle1 <- get_rle(ge = finher_expr_norm2$mas$sum[index, ])
# rle2 <- get_rle(ge = finher_expr_norm2$mas$combat[index, ])


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
pdf(file = paste0(out_figures, "fig3_rle.pdf"),
    width = 7.5,
    height = 5,
    title = "Figure3: RLE plot")

print(p1)

dev.off()


png(file = paste0(out_figures, "fig3_rle.png"),
    width = 7.5,
    height = 5,
    units = "in",
    res = 150,
    title = "Figure3: RLE plot")

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

# To test RLE plot with only pvac probesets
# identical(finher_expr_norm2$mas$sum$Probeset_Id,
#           finher_expr_norm2$mas$pvac$Probeset_Id) # TRUE
# identical(finher_expr_norm2$mas$combat$Probeset_Id,
#           finher_expr_norm2$mas$pvac$Probeset_Id) # TRUE
# index <- finher_expr_norm2$mas$pvac$Is_Selected_Probeset
# rle1 <- get_rle(ge = finher_expr_norm2$mas$sum[index, ])
# rle2 <- get_rle(ge = finher_expr_norm2$mas$combat[index, ])


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
pdf(file = paste0(out_figures, "fig3_rle_pvac.pdf"),
    width = 7.5,
    height = 5,
    title = "Figure3: RLE plot")

print(p1)

dev.off()


png(file = paste0(out_figures, "fig3_rle_pvac.png"),
    width = 7.5,
    height = 5,
    units = "in",
    res = 150,
    title = "Figure3: RLE plot")

print(p1)

dev.off()


#
# ==============================================================================



# ==============================================================================
# Supp.Fig1: Normalization comparison plot
# ==============================================================================


plot_list <- purrr::map(
  supp_norm_qc[c("colon","pancreas")],
  # Note:
  # Breast 1&2 u133plus2 datasets are only used for inhouse exploratory analysis
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
p <- c(p[[1]], p[[2]]) # "colon"    "pancreas"



# Normalization comparison plot
pdf(file = paste0(out_figures, "sfig1_normalization_comparison.pdf"),
    width = 7.5,
    height = 4.5,
    title = "Normalization comparison")

ggpubr::ggarrange( #ggdplyr::arrange(
  plotlist = p,
  labels = c("a","b","c","d"),
  #  "colon_mas"    "colon_rma"    "pancreas_mas" "pancreas_rma"
  #  colon: frozen, pancreas: ffpe
  ncol = 2, nrow = 2#, widths = c(.4,.6), vjust = 1.1
)

dev.off()


# png(file = paste0(out_figures, "sfig1_normalization_comparison.png"),
#     width = 7.5,
#     height = 4.5,
#     units = "in",
#     res = 150,
#     title = "Normalization comparison")
#
# ggpubr::ggarrange( #ggdplyr::arrange(
#   plotlist = p,
#   labels = c("A","B","C","D"),
#   #  "colon_mas"    "colon_rma"    "pancreas_mas" "pancreas_rma"
#   #  colon: frozen, pancreas: ffpe
#   ncol = 2, nrow = 2#, widths = c(.4,.6), vjust = 1.1
# )
#
# dev.off()


#
# ==============================================================================



# ==============================================================================
# Curated clinical and expression data in table format
# ==============================================================================

# Expression and clinical
# >>>>>>>>>>>>>>>>>>>>>>>

identical(finher_clin$CEL_filename, names(finher_expr_norm2$mas$combat)[-1])
identical(finher_clin$CEL_filename, names(finher_expr_norm2$rma$combat)[-1])
# TRUE

write_delim(x = finher_clin, path = str_c(out_tables,"finher_clin.tsv"), delim = "\t")
write_delim(x = finher_expr_norm2$mas$combat, path = str_c(out_tables,"finher_expr_mas.tsv"), delim = "\t")
write_delim(x = finher_expr_norm2$rma$combat, path = str_c(out_tables,"finher_expr_rma.tsv"), delim = "\t")


# Annot
# >>>>>

identical(finher_expr_mas$Probeset_Id, finher_expr_norm2$mas$maxvar$Probeset_Id)
identical(finher_expr_rma$Probeset_Id, finher_expr_norm2$rma$maxvar$Probeset_Id)
# TRUE

write_delim(x = finher_expr_norm2$mas$maxvar, path = str_c(out_tables,"finher_expr_mas_annot.tsv"), delim = "\t")
write_delim(x = finher_expr_norm2$rma$maxvar, path = str_c(out_tables,"finher_expr_rma_annot.tsv"), delim = "\t")

#
# ==============================================================================
