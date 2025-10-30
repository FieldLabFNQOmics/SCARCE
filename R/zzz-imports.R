#' Internal: centralised imports for SCARCE
#' @name scarce-imports
#' @keywords internal
#' @noRd

#' @importFrom methods as
#' @importFrom utils URLdecode capture.output head tail modifyList
#' @importFrom stats dist hclust fisher.test na.omit p.adjust pnorm quantile setNames
#' @importFrom grDevices colorRampPalette dev.off png

#' @importFrom magrittr %>%
#' @importFrom dplyr filter mutate group_by ungroup summarise across everything arrange desc row_number relocate left_join pull bind_cols bind_rows case_when reframe rowwise ungroup
#' @importFrom tidyr pivot_longer pivot_wider ends_with
#' @importFrom tibble rownames_to_column column_to_rownames enframe deframe
#' @importFrom stringr str_detect str_extract str_split
#' @importFrom rlang .data sym
#' @importFrom foreach %do% foreach

#' @importFrom ggplot2 ggplot aes geom_point geom_bar geom_text position_stack coord_equal geom_histogram
#' @importFrom ggplot2 guides guide_legend theme theme_bw theme_classic element_text labs
#' @importFrom ggplot2 scale_color_manual scale_colour_gradient2 scale_color_discrete scale_fill_manual scale_y_continuous ggtitle ggsave
#' @importFrom patchwork wrap_plots plot_annotation
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation rowAnnotation Legend draw anno_simple
#' @importFrom circlize colorRamp2
#' @importFrom grid unit
#' @importFrom ComplexUpset upset
#' @importFrom networkD3 sankeyNetwork
#' @importFrom RColorBrewer brewer.pal
#' @importFrom SCINA SCINA
#' @importFrom Seurat CreateSeuratObject DefaultAssay<- DimPlot DoHeatmap ElbowPlot Embeddings FeaturePlot FindClusters FindNeighbors GetAssayData Idents Idents<- NormalizeData RunPCA RunUMAP ScaleData
#' @importFrom SeuratObject Features

#' @importFrom data.table fread := setnames rbindlist
#' @importFrom readr write_tsv
#' @importFrom readxl read_xlsx

#' @importFrom rhdf5 h5read
#' @importFrom HDF5Array HDF5Array
#' @importFrom Matrix Matrix

#' @importFrom VariantAnnotation scanVcfHeader samples ScanVcfParam readVcf expand geno ref alt
#' @importFrom Rsamtools TabixFile
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom GenomicRanges start
#' @importFrom GenomeInfoDb seqnames
#' @importFrom S4Vectors elementNROWS

#' @importFrom uwot umap
#' @importFrom dbscan dbscan
#' @importFrom ggnewscale new_scale_color new_scale_fill
#' @importFrom htmlwidgets saveWidget
#' @importFrom scales label_number
#' @importFrom withr with_options

NULL

utils::globalVariables(c(
  "REF","ALT","ref_len","alt_len","complex_snv","complex_snv_pos",
  "simple_ins","complex_ins","complex_del","simple_snv","mnv","skip",
  "var_type","CHROM","START","END","new_REF","new_ALT","var_key","filter",
  "tissue_type","cell_type","cell_name","marker","Symbol","probability",
  "cell_barcode","SCINA_label","UMAP_1","UMAP_2","cluster","deframe",
  "barcode","NGT","reference_allele","alternate_allele","chromosome",
  "start_position","end_position","allele","strand","variant","PICK","ID",
  "SYMBOL","Gene","Feature","HGVSp","MAX_AF","Consequence",
  "Existing_variation","BIOTYPE","SIFT","PolyPhen","plot_ID","priority_flag",
  "alt_cnt_total","alt_proportion_total","data_cnt_total","data_proportion_total",
  "mean_GQ_total","mean_AF_total","mean_AF_mutants","rownumber","var_ngt",
  "mean_GQ","mean_AF","padj","Z_score","OR","sig","feature","count",
  "alt_frequency","data_count","alt_count","variant_order","prep_vcf_dir",
  "alt_cnt_ct","alt_proportion_ct","data_cnt_ct","alt_proportion_in_ct"
))
