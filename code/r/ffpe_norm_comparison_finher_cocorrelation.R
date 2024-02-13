# ffpe_norm_comparison_finher_cocorrelation.R

# What this script does?
# >>>>>>>>>>>>>>>>>>>>>>
#
# This script compare MAS and RMA algorithm in terms of co correlation between 
# genes from ER pathway in MAS and RMA normalized FinHER dataset.
# A better co correlation structure means better normalization

# Script structure
# >>>>>>>>>>>>>>>>
# 
# 1. Load (maxvar probeset) expression data
# 2. Load ER signature
# 3. Co-correlation within up-regulated and down-regulated genes
# 4. Heat map


# 1. Load (maxvar probeset) expression data
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# source("src/r/functions.R")
# load("results/data/finher_expr_norm2.RData")

dmas <- left_join(finher_expr_norm2$mas$maxvar, finher_expr_norm2$mas$combat, by = "Probeset_Id")
drma <- left_join(finher_expr_norm2$rma$maxvar, finher_expr_norm2$rma$combat, by = "Probeset_Id")

dmas <- dmas %>%
  dplyr::select(-c("Probeset_Id", "Gene_Symbol")) %>%
  dplyr::rename(Ncbi_gene_id = "Entrez_Gene") %>%
  dplyr::mutate(Ncbi_gene_id = str_c("ncbi_", Ncbi_gene_id))


drma <- drma %>%
  dplyr::select(-c("Probeset_Id", "Gene_Symbol")) %>%
  dplyr::rename(Ncbi_gene_id = "Entrez_Gene") %>%
  dplyr::mutate(Ncbi_gene_id = str_c("ncbi_", Ncbi_gene_id))

dmas <- dmas %>% t_tibble(names_x_desc = "Sample_name")
drma <- drma %>% t_tibble(names_x_desc = "Sample_name")




# Load ER signature
# >>>>>>>>>>>>>>>>>

# ER signature from Desmedt et al. 2008 CCR doi: 10.1158/1078-0432.CCR-07-4756

er_desmedt2008 <- readxl::read_excel(path = "data/gene_modules_erher2/Desmedt2008_Immune-Proliferation_Table_S1_clean_er.xlsx")

# module cleaning
er_desmedt2008 <- er_desmedt2008 %>%
  dplyr::mutate(
    Ncbi_gene_id = str_c("ncbi_", EntrezGene.ID),
    Direction = if_else(coefficient < 0 ,-1, 1)
  ) %>%
  dplyr::select(Ncbi_gene_id, Direction) %>%
  na.omit()


# gene filtering according to mas/rma finher data
er_desmedt2008 <- er_desmedt2008 %>%
  dplyr::filter(Ncbi_gene_id %in% names(drma)) %>%
  dplyr::filter(Ncbi_gene_id %in% names(dmas)) %>%
  dplyr::distinct(Ncbi_gene_id, .keep_all = TRUE)

er_up = er_desmedt2008 %>%
  dplyr::filter(Direction == 1) %>%
  dplyr::select(Ncbi_gene_id) %>%
  tibble::deframe()
er_down = er_desmedt2008 %>%
  dplyr::filter(Direction == -1) %>%
  dplyr::select(Ncbi_gene_id) %>%
  tibble::deframe()


# 3. Co-correlation within up-regulated and down-regulated genes
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# All-er
drma %>% dplyr::select(all_of(er_up),all_of(er_down)) %>% cor() %>% mean()
# [1] 0.1120924
dmas %>% dplyr::select(all_of(er_up),all_of(er_down)) %>% cor() %>% mean()
# [1] 0.2009971

# Up-er
drma %>% dplyr::select(er_up) %>% cor() %>% mean()
# [1] 0.3284034
dmas %>% dplyr::select(er_up) %>% cor() %>% mean()
# [1] 0.3957734

# Down-er
drma %>% dplyr::select(er_down) %>% cor() %>% mean()
# [1] 0.258074
dmas %>% dplyr::select(er_down) %>% cor() %>% mean()
# [1] 0.3217032



# 4. Heat map
# >>>>>>>>>>>>

# Co-correlation heat map of unregulated genes

x <- drma %>%
  dplyr::select(er_up) %>%
  cor()
x <- as.dendrogram(hclust(d = dist(x)))
rma_cluster_order <- order.dendrogram(x)


# ER genes in RMA
ggdf <- drma %>%
  dplyr::select(all_of(er_up)) %>%
  dplyr::select(all_of(rma_cluster_order)) %>%
  cor() %>% # default pearson
  as_tibble(rownames = "Row") %>%
  tidyr::gather(key = "Column", value = "Pearson", -1) %>%
  dplyr::mutate(Pearson = round(Pearson, digits = 1),
                Row = factor(Row, levels = unique(Row)),
                Column = factor(Column, levels = levels(Row)))


p1 <- ggplot(data = ggdf, aes(x = Column, y = Row)) +
  geom_raster(aes(fill = Pearson)) +
  scale_fill_gradient2(low = "#053061",
                       mid = "#f7f7f7",
                       high = "#67001f",
                       midpoint = 0,
                       limits = c(-1,1))+
  # geom_text(aes(label = Pearson), size = 5) +
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0))+
  theme(legend.position = "bottom",
        axis.text = element_blank(),
        axis.ticks = element_blank()) + #axis.title = element_blank()
  labs(x = "Genes upregulated in ER signaling", y = "Genes upregulated in ER signaling",
       title = "RMA normalized dataset - average co-correlation: 0.33 ")

# ER genes in MAS
ggdf <- dmas %>%
  dplyr::select(all_of(er_up)) %>%
  dplyr::select(all_of(rma_cluster_order)) %>%
  cor() %>% # default pearson
  as_tibble(rownames = "Row") %>%
  tidyr::gather(key = "Column", value = "Pearson", -1) %>%
  dplyr::mutate(Pearson = round(Pearson, digits = 1),
                Row = factor(Row, levels = unique(Row)),
                Column = factor(Column, levels = levels(Row)))


p2 <- ggplot(data = ggdf, aes(x = Column, y = Row)) +
  geom_raster(aes(fill = Pearson)) +
  scale_fill_gradient2(low = "#053061",
                       mid = "#f7f7f7",
                       high = "#67001f",
                       midpoint = 0,
                       limits = c(-1,1))+
  # geom_text(aes(label = Pearson), size = 5) +
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0))+
  theme(legend.position = "bottom",
        axis.text = element_blank(),
        axis.ticks = element_blank()) + #axis.title = element_blank()
  labs(x = "Genes upregulated in ER signaling", y = "Genes upregulated in ER signaling",
       title = "MAS normalized dataset - average co-correlation: 0.40 ")



# plot printin
# library(ggpubr) # ggarrange: Arrange nultiple ggplots
# # annotate_figure: Annotate figures including: i) ggplots, ii) arranged ggplots from ggarrange(), grid.arrange() and plot_grid().

# png(filename = str_c(out_figures, "ER_upgene_co_correlation.png"),
#     width = 6, height = 7.5, units= "in", res =300)
pdf(file = str_c(out_figures, "ER_upgene_co_correlation.pdf"),
    width = 6, height = 7.5)
print(
  ggarrange( p1, p2,
             ncol = 1,
             nrow = 2,
             labels = c("a", "b"),
             align = "hv",
             widths = 1,
             # heights = c(.61,.39),
             legend = "right",
             common.legend = T
  ),
  newpage = F
)
dev.off()


# Different strategy to categorizes er genes based on high/mid/low expression
# irrespective of their up/down regulation status.
#
# Common low/mid/high
x <- drma %>% dplyr::select(er_desmedt2008$Ncbi_gene_id) %>%
  purrr::map_dbl(mean) %>% summary()
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 2.524   4.755   5.574   5.547   6.354   9.433
dmas %>% dplyr::select(er_desmedt2008$Ncbi_gene_id) %>%
  purrr::map_dbl(mean) %>% summary()
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 4.592   5.934   6.537   6.568   7.133   9.898

er_gene_class <- er_desmedt2008 %>%
  dplyr::mutate(
    avg_expr_rma = drma %>%
      dplyr::select(er_desmedt2008$Ncbi_gene_id) %>%
      purrr::map_dbl(mean),
    avg_expr_mas = dmas %>%
      dplyr::select(er_desmedt2008$Ncbi_gene_id) %>%
      purrr::map_dbl(mean),
    avg_expr_rma_class = purrr::map_chr(
      avg_expr_rma,
      function(x,expr_sum){
        case_when(
          x <= expr_sum["1st Qu."] ~ "low",
          x >= expr_sum["3rd Qu."] ~ "high",
          TRUE ~ "mid"
        )
      },
      expr_sum = summary(avg_expr_rma)
    ),
    avg_expr_mas_class = purrr::map_chr(
      avg_expr_mas,
      function(x,expr_sum){
        case_when(
          x <= expr_sum["1st Qu."] ~ "low",
          x >= expr_sum["3rd Qu."] ~ "high",
          TRUE ~ "mid"
        )
      },
      expr_sum = summary(avg_expr_mas)
    )
  )

er_high <- er_gene_class %>%
  dplyr::filter(avg_expr_rma_class == avg_expr_mas_class) %>%
  dplyr::filter(avg_expr_rma_class == "high") %>%
  dplyr::select(Ncbi_gene_id) %>%
  tibble::deframe()

er_low <- er_gene_class %>%
  dplyr::filter(avg_expr_rma_class == avg_expr_mas_class) %>%
  dplyr::filter(avg_expr_rma_class == "low") %>%
  dplyr::select(Ncbi_gene_id) %>%
  tibble::deframe()

er_mid <- er_gene_class %>%
  dplyr::filter(avg_expr_rma_class == avg_expr_mas_class) %>%
  dplyr::filter(avg_expr_rma_class == "mid") %>%
  dplyr::select(Ncbi_gene_id) %>%
  tibble::deframe()


# Intersection between ER-up regulated and ER genes with high/mid/low expression in dataset

# er-high
drma %>% dplyr::select(intersect(er_high,er_up)) %>% cor() %>% mean()
# [1] 0.380493
dmas %>% dplyr::select(intersect(er_high,er_up)) %>% cor() %>% mean()
# [1] 0.4399851

# er-mid
drma %>% dplyr::select(intersect(er_mid,er_up)) %>% cor() %>% mean()
# [1] 0.3161805
dmas %>% dplyr::select(intersect(er_mid,er_up)) %>% cor() %>% mean()
# [1] 0.386422

# er-low
drma %>% dplyr::select(intersect(er_low,er_up)) %>% cor() %>% mean()
# [1] 0.4012327
dmas %>% dplyr::select(intersect(er_low,er_up)) %>% cor() %>% mean()
# [1] 0.462462


# Intersection between ER-down regulated and ER genes with high/mid/low expression in dataset

# er-high
drma %>% dplyr::select(intersect(er_high,er_down)) %>% cor() %>% mean()
# [1] 0.2428905
dmas %>% dplyr::select(intersect(er_high,er_down)) %>% cor() %>% mean()
# [1] 0.2888708

# er-mid
drma %>% dplyr::select(intersect(er_mid,er_down)) %>% cor() %>% mean()
# [1] 0.3057669
dmas %>% dplyr::select(intersect(er_mid,er_down)) %>% cor() %>% mean()
# [1] 0.3763115

# er-low
drma %>% dplyr::select(intersect(er_low,er_down)) %>% cor() %>% mean()
# [1] 0.2825506
dmas %>% dplyr::select(intersect(er_low,er_down)) %>% cor() %>% mean()
# [1] 0.3408725

