# functions.R


# PVAC probeset filtering
pvacFilter <- function(abatch, probeIndexList, dCal = NULL, neg.probeSets=NULL){
  #abatch: bg corrected normalized abatch
  #probeIndexList: intensity index list
  #dCal: detection call logical matrix (relevant if neg.probeSet is NULL)
  #neg.probeSets: negative probesets


  if(is.null(neg.probeSets)){
    neg.probeSets=apply(dCal,1,function(x){all(!x)})
    neg.probeSets=names(neg.probeSets)[neg.probeSets]
  }

  dat=t(log2(intensity(abatch)))
  dat=scale(x=dat,center=T,scale=T) # base::scale
  pcVarList=lapply(names(probeIndexList),
                   function(x,probeIndexList,dat){
                     index=probeIndexList[[x]]
                     pcVar=(prcomp(x=dat[,index],retx=F,center=F,scale.=F)$sdev)^2
                     #return(pcVar[1]/sum(pcVar))
                     return(list(pc1.var=pcVar[1],pc2.var=pcVar[2],total.var=sum(pcVar)))
                   },probeIndexList,dat)

  pcVarDf=data.frame(probeSet.id=names(probeIndexList),pc1.var=NA,
                     pc2.var=NA,total.var=NA,stringsAsFactors=F)
  rownames(pcVarDf)=pcVarDf$probeSet.id
  pcVarDf$pc1.var=unlist(lapply(pcVarList,function(x){x$pc1.var}))
  pcVarDf$pc2.var=unlist(lapply(pcVarList,function(x){x$pc2.var}))
  pcVarDf$total.var=unlist(lapply(pcVarList,function(x){x$total.var}))
  pcVarDf$pc1.varProp=pcVarDf$pc1.var/pcVarDf$total.var
  pcVarDf$pc2.varProp=pcVarDf$pc2.var/pcVarDf$total.var
  #isNeg.probeSet
  pcVarDf$isNeg.probeSet=FALSE
  index=sapply(neg.probeSets,function(x,id){which(x==id)},id=pcVarDf$probeSet.id)
  pcVarDf$isNeg.probeSet[index]=TRUE
  #isSelected.probeSet
  pvac.neg=pcVarDf$pc1.varProp[pcVarDf$isNeg.probeSet]
  cutoff=min(c(quantile(x=pvac.neg,probs=.99),.5))# setting ctoff at max 0.5
  pcVarDf$isSelected.probeSet=FALSE
  pcVarDf$isSelected.probeSet[pcVarDf$pc1.varProp>cutoff]=TRUE
  pcVarDf=pcVarDf[,c("probeSet.id","isNeg.probeSet","isSelected.probeSet",
                     "pc1.var","pc2.var","total.var","pc1.varProp","pc2.varProp")]

  outList=list(pvac.df=pcVarDf,pvac.cutoff=cutoff)

  return(outList)
}


# Extract/compute normalization (MAS/RMA) qc data
get_norm_qc_data <- function(ge, pp, neg_pset, pos_pset){

  # ge: tibble, gene expression
  # pp: tibble, percent-present
  # neg_pset: char vector, negative probesets
  # pos_pset: char vector, positive probesets

  cat("Preparing qc data...\n")

  qc_data <- tibble(
    CEL_filename = names(ge)[-1],
    PP =  pp %>%
      dplyr::select(-Probeset_Id) %>%
      transmute_all(function(x){x < .01}) %>%
      purrr::map_dbl(function(x){ (sum(x) / length(x)) * 100}),
    Neg_Avg = ge %>%
      filter(Probeset_Id %in% neg_pset) %>%
      dplyr::select(-Probeset_Id) %>%
      purrr::map_dbl(function(x){mean(x)}),
    Pos_Avg = ge %>%
      filter(Probeset_Id %in% pos_pset) %>%
      dplyr::select(-Probeset_Id) %>%
      purrr::map_dbl(function(x){mean(x)}),
    GE_Avg = ge %>%
      dplyr::select(-Probeset_Id) %>%
      purrr::map_dbl(function(x){mean(x)})
  )

  return(qc_data)
}


# Extract/compute sample qc data
get_sample_qc_data <- function(ge, pp, neg_pset, pos_pset){
  # ge: tibble, gene expression
  # pp: tibble, percent-present
  # neg_pset: char vector, negative probesets
  # pos_pset: char vector, positive probesets

  # ge <- data_array$sum
  # pp <- data_array$pp
  # neg_pset <- grep(data_array$sum$Probeset_Id, pattern = "GC", fixed = T, value = T)
  # pos_pset <- hgu133plus2.hgu219_genomicPosControl$hgu219

  cat("Running PCA...\n")
  pca <- ge %>%
    dplyr::select(-Probeset_Id) %>%
    t() %>%
    prcomp(retx = TRUE, center = TRUE, scale. = TRUE)
  pca_var_prop <- (pca$sdev^2 / sum(pca$sdev^2))*100


  cat("Preparing qc data...\n")
  qc_data <- tibble(CEL_filename = names(ge)[-1],
                    PP =  pp %>%
                      dplyr::select(-Probeset_Id) %>%
                      transmute_all(function(x){x < .01}) %>%
                      purrr::map_dbl(function(x){ (sum(x) / length(x)) * 100}),
                    Neg_Avg = ge %>%
                      filter(Probeset_Id %in% neg_pset) %>%
                      dplyr::select(-Probeset_Id) %>%
                      purrr::map_dbl(function(x){mean(x)}),
                    Pos_Avg = ge %>%
                      filter(Probeset_Id %in% pos_pset) %>%
                      dplyr::select(-Probeset_Id) %>%
                      purrr::map_dbl(function(x){mean(x)}),
                    GE_Avg = ge %>%
                      dplyr::select(-Probeset_Id) %>%
                      purrr::map_dbl(function(x){mean(x)}),
                    PC1 = pca$x[ , "PC1"],
                    PC2 = pca$x[ , "PC2"],
                    PC3 = pca$x[ , "PC3"],
                    PC4 = pca$x[ , "PC4"],
                    PC1_Var_Prop = pca_var_prop[1],
                    PC2_Var_Prop = pca_var_prop[2],
                    PC3_Var_Prop = pca_var_prop[3],
                    PC4_Var_Prop = pca_var_prop[4],
  )

  return(qc_data)
}


# Generate sample qc plots
get_sample_qc_plots <- function(qc_data, outdir, filename){
  # qc_data qc_data_tibble
  # outdir output dir
  # filename file name pattern

  # qc_data <-  sample_qc_data
  # outdir <- outdir
  # filename <- "RMA"
  #

  # Main QC plot
  # >>>>>>>>>>>>

  # Pos -Neg plot
  plot_pos_neg <- qc_data %>%
    mutate(Pos_Neg_Filter = map2_lgl(Neg_Avg, Pos_Avg, ~.x < .y)) %>%
    gather(key="key", value="value", Neg_Avg, Pos_Avg) %>%
    mutate_at("key",
              function(x){
                case_when(
                  x == "Neg_Avg" ~ "Negative control\n(background)",
                  x == "Pos_Avg" ~ "Positive control\n(housekeeping)"
                )
              }
    ) %>%
    ggplot(aes(x=key,y=value))+
    theme_bw() +
    geom_line(aes(color= PP>10 , linetype = Pos_Neg_Filter,group = CEL_filename))+
    geom_boxplot(fill=NA) +
    scale_linetype_manual(values = c("TRUE" = "solid", "FALSE" = "dashed")) +
    scale_color_manual(values = c("TRUE" = "gray20", "FALSE" = "goldenrod")) + ##e08214
    guides(linetype = guide_legend(title = "Neg < Pos")) +
    labs(y = "Gene expression") +
    theme(axis.title.x = element_blank(),
          legend.text = element_text(size = 8),
          legend.title =  element_text(size = 8),
          legend.key.size = unit( 1, "lines"),
          legend.background = element_rect(fill = NA),
          legend.margin = margin(r=-5),
          legend.box.margin = margin(),
          legend.box = "horizontal",
          legend.position = c(0.01,0.99),
          legend.justification = c(0,1))


  # Batch plot
  plot_batch <- qc_data %>%
    mutate_at("Batch",~if_else(Batch == "Batch6", "Batch4", Batch)) %>%
    gather(key="key", value="value", PP, Neg_Avg, Pos_Avg, GE_Avg) %>%
    mutate_at(
      "key",
      function(x){
        case_when(
          x == "PP" ~ "Percent-Present",
          x == "Neg_Avg" ~ "Negative control",
          x == "Pos_Avg" ~ "Positive control",
          x == "GE_Avg" ~ "Genomic avarage"
        )
      }
    ) %>%
    mutate_at(
      "key",
      ~factor(.x, levels = c("Percent-Present",
                             "Negative control",
                             "Positive control",
                             "Genomic avarage"))
    ) %>%
    ggplot(aes(x=Batch,y=value))+
    theme_bw() +
    geom_boxplot(aes(fill = Batch)) +
    scale_fill_manual(values = c("#a6cee3","#1f78b4","#b2df8a","#33a02c")) +
    guides(fill = "none") +
    theme(axis.title.y = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    facet_wrap(~key, nrow = 1, scales = "free", strip.position = "left")



  plot_distribution <-  qc_data %>%
    gather(key="key", value="value", Neg_Avg, Pos_Avg, GE_Avg) %>%
    mutate_at(
      "key",
      function(x){
        case_when(
          x == "Neg_Avg" ~ "Negative control",
          x == "Pos_Avg" ~ "Positive control",
          x == "GE_Avg" ~ "Genomic avarage"
        )
      }
    ) %>%
    mutate_at(
      "key",
      ~factor(.x, levels = c("Percent-Present",
                             "Negative control",
                             "Positive control",
                             "Genomic avarage"))
    ) %>%
    ggplot(aes(x=value,y=PP))+
    theme_bw() +
    geom_point( shape = 1, color = "gray20") + # aes(color = Batch),
    # scale_color_manual(values = c("#a6cee3","#1f78b4","#b2df8a","#33a02c")) +
    geom_smooth(method = "loess", se = FALSE, color = "goldenrod") +
    guides(color = "none") +
    labs(x="Gene expression") +
    theme(axis.title.y = element_blank()) +
    facet_grid("Percent-Present"~key, scales = "free", switch = "y")

  # To annotate plot_distribution
  corr_df <- tibble(
    "Precent-Present" = "Percent-Present",
    key = c("Negative control",
            "Positive control",
            "Genomic avarage"),
    Pearson = c(cor(qc_data$PP, qc_data$Neg_Avg),
                cor(qc_data$PP, qc_data$Pos_Avg),
                cor(qc_data$PP, qc_data$GE_Avg))
  ) %>%
    mutate_at("Pearson", ~str_c("Pearson:", round(Pearson, digits = 2))) %>%
    mutate_at("key", ~factor(key, levels = c("Negative control",
                                             "Positive control",
                                             "Genomic avarage")))

  plot_distribution <- plot_distribution +
    geom_text(data = corr_df, aes(x = -Inf, y= +Inf, label = Pearson, hjust = -.1, vjust = 1.5))



  # Main qc plot
  pdf(file = paste0(outdir,filename,"_main_qc.pdf"),
      width = 7.5,
      height = 4,
      title = str_c(filename,"main QC", sep = " "))

  xx <- ggpubr::ggarrange( #ggdplyr::arrange(
    ggpubr::ggarrange( #ggdplyr::arrange(
      plot_pos_neg ,
      plot_distribution,
      labels = c("A","B"), ncol = 2, widths = c(.4,.6), vjust = 1.1
    ),
    plot_batch,
    labels = c("","C"), nrow = 2, ncol = 1 , vjust = 1.1)

  print(xx, newpage = F)

  dev.off()



  # PCA QC plot
  # >>>>>>>>>>>

  pca_qc_plots <- get_pca_qc_plots(
    qc_data = qc_data,
    aes_color = c("PP", "Neg_Avg", "Pos_Avg", "GE_Avg", "Batch",
                  "ER_IHC", "HER2_IHC_CISH", "Subtype_IHC")
  )


  pdf(file = paste0(outdir,filename,'_pca_qc.pdf'),
      width = 7.5,
      height = 8.75,
      title = str_c(filename, "PCA QC", sep = " "))

  # Level 1 grouping to merge legend
  xx <-  map2(pca_qc_plots$PCA_1.2,pca_qc_plots$PCA_3.4,
              # ~ggdplyr::arrange(.x,.y,nrow=1,legend="top", common.legend = T))
              ~ggpubr::ggarrange(.x,.y,nrow=1,legend="top", common.legend = T))
  # Level to grouping to label plots
  # yy <- ggdplyr::arrange(plotlist = xx, nrow = 4, ncol = 2,
  #                        labels = c("A","B","C","D","E","F","G","H"), vjust = 1.1)
  yy <- ggpubr::ggarrange(plotlist = xx, nrow = 4, ncol = 2,
                         labels = c("A","B","C","D","E","F","G","H"), vjust = 1.1)

  print(yy, newpage = F)

  dev.off()



  # # Paired qc plot: data_eblm_array$sum_batch
  # pdf(file = paste0(outdir,'/sample_',filename,'_qc_pair_plot.pdf'),
  #     onefile=TRUE,
  #     width = 7.5,
  #     height = 8.75,
  #     title = filename)
  #
  # xx <- GGally::ggpairs(
  #   data = qc_data,
  #   columns = c("PP","Neg_Avg","Pos_Avg","GE_Avg","PC1","PC2","Batch"),
  #   mapping = aes(color = Batch,
  #                 fill = Batch),
  #   diag = list(continuous = wrap("densityDiag",
  #                                 fill = NA)),
  #   # see ?ggally_densityDiag; prepend "ggally_" to plot constructs, i.e "densityDiag"
  #   upper = list(continuous = wrap("cor",
  #                                  size = 3,
  #                                  alignPercent = 1),
  #                combo =wrap("box_no_facet",
  #                            color = "gray20")),
  #   lower = list(continuous = wrap("points",
  #                                  shape = 1))
  #   # see ?ggally_cor; prepend "ggally_" to plot constructs, i.e "cor"
  #   # see ?ggally_box_no_facet; prepend "ggally_" to plot constructs, i.e "box_no_facet"
  # ) +
  #   theme(axis.text.x = element_text(angle = 45,
  #                                    hjust = 1,
  #                                    vjust = 1))
  #
  # print(xx)
  #
  # dev.off()

}


# Extract/compute probeset qc data
get_pset_qc_data <- function(ge, pp, pvac, pos_pset){ #neg_pset,

  # !!!! neg_pset is not used in this function

  # ge: gene expression
  # pp: percent-present
  # pvac: pvac df
  # neg_pset: negative probesets
  # pose_pset: positive probesets

  # ge <- data_array$sum
  # pp <- data_array$pp
  # pvac <- data_array$pvac
  # neg_pset <- grep(data_array$sum$Probeset_Id, pattern = "GC", fixed = T, value = T)
  # pos_pset =  hgu133plus2.hgu219_genomicPosControl$hgu219


  pp =  pp %>%
    dplyr::select(-Probeset_Id) %>%
    transmute_all(function(x){x < .01}) %>%
    purrr::map_dbl(function(x){ (sum(x) / length(x)) * 100})


  pvac <-   pvac %>%
    mutate(
      Probe_PC1_Var_Prop = PC1_Var_Prop,
      PVAC_Status = if_else(Is_Selected_Probeset, "Selected", "Discarded"),
      Probeset_Type = if_else(
        Is_Neg_Probeset,
        "Anitgenomic",
        if_else(Probeset_Id %in% pos_pset,
                "Housekeeping","Genomic")
      )
    )


  qc_data <- tibble(

    Probeset_Id = ge %>%
      # filter(!Probeset_Id %in% neg_pset) %>%
      dplyr::select(Probeset_Id) %>%
      tibble::deframe(),

    Values = ge %>%
      # filter(!Probeset_Id %in% neg_pset) %>%
      dplyr::select(-Probeset_Id) %>%
      t() %>%
      as.data.frame() %>%
      map_chr(
        function(x,pp){
          str_c(
            c(
              var(x),
              mean(x),
              unlist(cor.test(x, pp)[c("estimate", "p.value")])
            ),
            collapse="///"
          )
        },
        pp
      )
  )

  qc_data <- qc_data %>%
    separate(
      col = "Values",
      into = c("Pset_GE_Var","Pset_GE_Avg", "Pset_PP_Cor","Pset_PP_Cor_Pval"),
      sep= "///",
      convert = TRUE
    ) %>%
    left_join(
      pvac %>%
        dplyr::select(Probeset_Id, Probe_PC1_Var_Prop, PVAC_Status, Probeset_Type),
      by = "Probeset_Id"
    )

  return(qc_data)
}


# Generate probeset qc plots
get_pset_qc_plots <- function(qc_data, outdir, filename){
  # qc_data qc_data_tibble
  # outdir output dir
  # filename file name pattern

  # qc_data <-  pset_qc$level1$rma
  # outdir <- out_figures
  # filename <- "Level1_RMA2"

  # ggpairs plot
  # >>>>>>>>>>>>

  # Integrate PVAC_Status and Probeset_Type for coloring
  qc_data <- qc_data %>%
    mutate(
      Probeset_Status = if_else(
        Probeset_Type == "Genomic",
        str_c("PVAC", str_to_lower(PVAC_Status), sep = " "), Probeset_Type
      )
    ) %>%
    mutate_at(
      "Probeset_Status",
      ~case_when(
        Probeset_Status == "Anitgenomic" ~ "Negative control",
        Probeset_Status == "Housekeeping" ~ "Positive control",
        TRUE ~ Probeset_Status
      )
    ) %>%
    mutate_at(
      "Probeset_Status",
      ~factor(
        .x,
        levels = c("PVAC selected", "PVAC discarded",
                   "Positive control", "Negative control")
      )
    )


  yy <- GGally::ggpairs(
    data = qc_data,
    columns = c("Pset_GE_Var", "Pset_GE_Avg", "Pset_PP_Cor"), # "Probe_PC1_Var_Prop",
    mapping = aes(color = Probeset_Status, fill = Probeset_Status),
    diag = list(continuous = wrap("densityDiag", fill = NA)),
    # see ?ggally_densityDiag; prepend "ggally_" to plot constructs, i.e "densityDiag"
    upper = list(continuous = wrap("cor", alignPercent = 1, displayGrid = FALSE, size = 3)),
    lower = list(continuous = wrap("points", shape = 1))
  ) + theme_bw()

  # Ref: https://stackoverflow.com/questions/34740210/how-to-change-the-color-palette-for-ggallyggpairs
  for(i in 1:yy$nrow) {
    for(j in 1:yy$ncol){
      yy[i,j] <- yy[i,j] +
        scale_fill_manual(values=c("gray20", "gray60", "#33a02c", "#1f78b4")) +
        scale_color_manual(values=c("gray20", "gray60", "#33a02c", "#1f78b4"))
    }
  }


  # Paired qc plot: pset_dc_data
  tiff(filename = paste0(outdir,filename,"_pset_qc.tiff"),
       units = "in",
       width = 5,
       height = 5,
       res = 250)
  print(yy)
  dev.off()

}


# Compute relative log expression
get_rle <- function(ge){
  # ge: tibble gene expression data
  # first column as Probeset_Id

  row_median <- ge %>%
    dplyr::select(-Probeset_Id) %>%
    as.matrix() %>%
    Biobase::rowMedians()

  rle <- ge %>%
    dplyr::select(-Probeset_Id) %>%
    map_df(function(x, row_median){x - row_median}, row_median)

  return(bind_cols(dplyr::select(ge, Probeset_Id), rle))
}


# Estimate PAM50 subtypes
get_pam50_subtype <- function(x, pam50){

  # x :  a gene expression tibble with genes on rows and samples on columns
  #       1st column "Ncbi_gene_id"
  # pam50: a tibble with pam50 centroids and gene annotations

  # x <- finher_expr$rma
  # pam50 <- pam50

  pam50 <- pam50 %>%
    dplyr::filter(Ncbi_gene_id %in% x$Ncbi_gene_id)

  print( str_c(50 - nrow(pam50), " pam50 genes were missing in the dataset"))

  # subsetting
  x <- pam50 %>%
    dplyr::select(Ncbi_gene_id) %>%
    dplyr::left_join(x, by = "Ncbi_gene_id")

  # Median gene centering recommended
  x[, -1] <- x[, -1] %>%
    apply(1,function(xx){xx-median(xx)}) %>%
    # gene centering, the result is transposed
    t() %>%
    as.data.frame()

  # Pearson correlation
  # (Median centering is valid only with pearson correlation.
  # Foe Spearman correlation, median centering has no effect).
  x_corr <- cor(pam50[, c("Basal", "Her2", "LumA", "LumB", "Normal")],
                x[,-1],
                method = "pearson")


  # Max correlated centroid as subtype, but if max corr <.1 the tumor is unclassified
  subtype <- purrr::map_dfr(
    colnames(x_corr),
    function(sample, x_corr, subtype){

      corr <- x_corr[, sample]

      list(
        CEL_filename = sample,
        Pam50_centroid_corr = corr[which(corr == max(corr))],
        Pam50_subtype = if_else(max(corr) < .1,
                                "Unclassified",
                                subtype[which(corr == max(corr))])
      )

    },
    x_corr,
    subtype = rownames(x_corr)
  )

  subtype
}


# Updated t_tibble()
# Transpose gene expression matrix with first column as gene/sample ids
t_tibble <- function(x, names_x_desc = NULL){

  # What this function does?
  # >>>>>>>>>>>>>>>>>>>>>>>>
  # Transposes a gene expression tibble.
  # By definition, row names are not recommended for tibble, instead any
  # rowname must be included as a column in the tibble. For instance the
  # as_tibble() function has "rownames" option to give names for existing rownames
  # so as to include the rownames as the first column of a tibble. However,
  # rownames are conventional and usually included in a tibble as the first column.
  # In this case, the default transpose function may not work properly as
  # transposing a tibble may cause data type conflicts and may convert the
  # entire column of a transposed tibble as the with the type of first column of
  # non-transposed tibble. This will be case with gene expression tibbles with
  # first column as gene symbols, and column names are sample names.
  # Note that this function is defined to use with gene expression tibbles.
  # Gene expression tibbles have the identical data type for all columns except
  # the first column representing rownames.

  # input
  # >>>>>
  # x: A tibble to transpose. The 1st column is considered as rownames.
  # names_x_desc: A character string that describes what names(x) represents.
  #               This character string will be used as the name of the first
  #               column of the trassposed tibble.

  # output
  # >>>>>>
  # A transposed tibble


  # Note: t(x) will set all data to single type

  # If names_x_desc is NULL, it is set as the name of the 1st column of x
  if (is.null(names_x_desc)) {
    names_x_desc <- names(x)[1]
  }

  # Extract rownames of x. This will later used to set column names of t(x)
  rownames_x <- x %>% dplyr::select(1) %>% tibble::deframe()


  tx <- t(x %>% dplyr::select(-1))
  colnames(tx) <- rownames_x

  return(as_tibble(tx, rownames = names_x_desc))
}


# Updated get_max_var_annot()
# Identifying maximum variant probeset per gene
get_max_var_annot <- function(ge, xannot, pset_colname, gene_colname){

  # What this function does?
  # >>>>>>>>>>>>>>>>>>>>>>>>
  # Identifies the maximum variant probe/probeset representing a single gene and
  # its annotation.

  # input
  # >>>>>
  # ge: Gene expression tibble with probe/probeset id as the first column and
  #     column names as sample names (names(ge)).The name of the probe/probeset id column
  #     must be the value of "pset_colname".
  # xannot: Annotation tibble with probe/probeset id and gene id columns named
  #           exactly as the value of "pset_colname" and "gene_colname".
  # pset_colname: Name of column containing probe/probeset ids,
  #                 must present in both ge and xannot.
  # gene_colname: Name of column containing gene ids,
  #                 must present in both ge and xannot.

  # output
  # >>>>>>
  # A tibble representing annotation of maximum variant probeset.


  # Converting all relevant ids to character
  ge <- ge %>%
    dplyr::mutate_at(pset_colname, ~as.character(.x))
  xannot <- xannot %>%
    dplyr::mutate_at(c(pset_colname, gene_colname), ~as.character(.x))


  # print(table((xannot %>% dplyr::select(all_of(pset_colname)) %>% tibble::deframe()) %in% (ge %>% dplyr::select(all_of(pset_colname)) %>% tibble::deframe())))


  pset <- intersect(ge %>% dplyr::select(all_of(pset_colname)) %>% tibble::deframe(),
                    xannot %>% dplyr::select(all_of(pset_colname)) %>% tibble::deframe())

  # # ***debugging
  # print(c(length(pset),nrow(xannot)))
  # print(length(pset)/nrow(xannot))


  # Filter and sort @ by
  # Ref: https://dplyr.tidyverse.org/articles/programming.html
  # Ref: https://stackoverflow.com/questions/26497751/pass-a-vector-of-variable-names-to-dplyr::arrange-in-dplyr
  # Both base::as.name(by) and rlang::sym(by) will take the strings as input and
  # return the value encoded by the string as symbols which can be used as
  # column names
  # is.symbol(by) # FALSE
  # is.symbol(as.name(by)) # TRUE
  # is.symbol(sym(by)) # TRUE
  ge <- ge %>%
    # dplyr::mutate_at(pset_colname, ~as.character(.x)) %>%
    dplyr::filter_at(pset_colname, dplyr::all_vars((. %in% pset))) %>%
    dplyr::arrange(!! as.name(pset_colname))

  xannot <- xannot %>%
    # dplyr::mutate_at(pset_colname, ~as.character(.x)) %>%
    dplyr::filter_at(pset_colname, dplyr::all_vars((. %in% pset))) %>%
    dplyr::arrange(!! as.name(pset_colname))

  # Ref !!: https://community.rstudio.com/t/eval-vs-bang-bang-in-functions-using-dplyr/3977/2



  # Pset variance vector
  pset_var <- ge %>%
    t_tibble() %>% # Transposing gene expression tibble (genes on column)
    dplyr::select(-1) %>% # Removing rownames (samplenames)
    purrr::map_dbl(var)

  # Unique gene ids
  ugid <- xannot %>%
    dplyr::select(all_of(gene_colname)) %>%
    dplyr::distinct_at(gene_colname, .keep_all = FALSE) %>%
    tibble::deframe()

  xannot_gene = xannot %>% dplyr::select(gene_colname) %>% tibble::deframe()

  # # ***debugging
  # print(
  #   c(identical(ge %>% dplyr::select(all_of(pset_colname)) %>% tibble::deframe(),
  #               xannot %>% dplyr::select(all_of(pset_colname)) %>% tibble::deframe()),
  #     identical(names(pset_var),
  #               xannot %>% dplyr::select(all_of(pset_colname)) %>% tibble::deframe())) %>%
  #     all()
  # )
  # # Note: names(pset_var) and xannot_gene are in order !!!!!!!!!!!!!!!!!!!!!!!
  # identical(names(pset_var), xannot[, pset_colname] %>% tibble::deframe())
  # TRUE


  # Max var pset per unique gene id
  max_var_pset <-  purrr::map_chr(
    ugid,
    function(gene, pset_var, xannot_gene ){

      # Per gene annot data frame
      # Note: names(pset_var) and xannot_gene are in order !!!!!!!!!!!!!!!!!!!!
      gene_pset_idx  <- which(xannot_gene == gene)
      gene_pset_var <- pset_var[gene_pset_idx]

      # Get max var pset
      names(gene_pset_var)[gene_pset_var == max(gene_pset_var)][1] #if identical maximums select 1st one
    },
    pset_var,
    xannot_gene
  )

  max_var_annot <- tibble(pset = max_var_pset)
  names(max_var_annot) = c(get("pset_colname"))

  max_var_annot <- left_join(max_var_annot,
                             xannot,
                             by = pset_colname)
  return(max_var_annot)

}


# Module score computation
get_module_score <- function(x, module_list, by = "Entrez_Gene"){

  # What this function does?
  # >>>>>>>>>>>>>>>>>>>>>>>>
  # Compute gene module score as a weighted average

  # input
  # >>>>>
  # x: A tibble with genes on column and samples on rows.
  #     1st column is considered as sample names.
  # module_list: List of gene modules. Each gene module must contain at least
  #               two columns; a gene id column with column name exactly as the
  #               value of the "by" input and a direction column with
  #               column name as "Direction". Further, the gene id
  #               in the gene module must match the gene ids of x, ie. column
  #               names of x. The Direction in the gene module represents the
  #               gene's association with the phenotype that the
  #               gene module represents. For instance, for a gene module
  #               representing TP53 mutation status a direction of +1 suggests
  #               up regulated in mutant phenotype and -1 represents
  #               down regulated in mutant phenotype.
  # by: The column name of gene module which contain gene id. The type gene id
  #     in the gene module must match that of column names of x.

  # output
  # >>>>>>
  # A tibble of module score. Column names represents module names and the 1st
  # column contains sample ids.


  sig_score <- purrr::map_dfr( # map_dfr
    module_list,
    function(sig, xdata, by){

      sig <- sig %>%
        dplyr::filter_at(
          by,
          dplyr::all_vars(. %in% names(xdata))
        )


      xdata <- xdata %>%
        dplyr::select(
          c(names(xdata)[1],
            sig %>% dplyr::select(by) %>% tibble::deframe())
        )


      score <- (xdata %>% dplyr::select(-1) %>% as.matrix()) %*% (sig %>% dplyr::select(Direction) %>% as.matrix())

      # score/nrow(sig) # this old code generates matrix
      as.numeric(score/nrow(sig))
    },
    xdata = x,
    by = by
  )

  sig_score <- dplyr::bind_cols(x %>% dplyr::select(1),
                                sig_score)

}



# To prepare data matrix for GEO submission
format_finher_matrix_for_geo <- function(x, celmap){

  # Rename finher samples to respective sample geo accessions

  celmap2 <- tibble(inhouse = names(x)[-1]) %>%
    dplyr::left_join(celmap, by="inhouse")

  print(identical(names(x)[-1], celmap2$inhouse)) # TRUE

  xx <- x %>%
    dplyr::rename_with(
      .fn = function(x){
        purrr::map_chr(x, function(x, celmap){

          idx <- which(x==celmap$inhouse)

          if(length(idx) == 1){
            return(celmap$GSE47994[idx])
          } else {
            return(x)
          }

        },
        celmap = celmap)
      } #, .cols = names(x)[-1]
    )

  print(identical(names(xx)[-1], celmap2$GSE47994)) # verification, expected TRUE


  # Rename Probeset_Id to ID_REF
  xx <- xx %>%   dplyr::rename(ID_REF = "Probeset_Id")

  # Round numerical values
  xx <- xx %>%
    dplyr::mutate(across(-"ID_REF",~round(.x, digits = 3)))

  return(xx)
}


# To prepare QC data for GEO submission
format_finher_sample_qc_for_geo <- function(x, nme, celmap){

  tibble(Filename = x$CEL_filename) %>%
    dplyr::left_join(
      celmap %>%
        dplyr::select(-inhouse) %>%
        dplyr::rename(Sample_geo_accession = "GSE47994",
                      Filename = "inhouse2"),
      by = "Filename"
    ) %>%
    dplyr::left_join(
      x %>%
        dplyr::select(nme) %>%
        dplyr::rename(Filename = "CEL_filename"),
      by = "Filename"
    )

}


# Format (round) probeset qc metrics
format_pset_qc_for_geo <- function(x, nme){
  
  x %>%
    dplyr::mutate(across(-nme,~round(.x, digits = 5)))
  # P value <0.000001 is set to zero
}

