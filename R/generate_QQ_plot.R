#' generate_QQ_plot.R
#'
#' Generates a QQ plot and histograms of p values from presumed null variants based
#' on VEP variant annotations.
#'
#' @param results_dir The directory where cell type variant enrichment results
#' are stored, namely `{prefix}`_`{cell_type}`.tsv files
#' @param prefix Prefix used in `find_somatic_variants()` to generate results
#' @param plot_dir Directory to output QQ plot
#' @param results_files Optionally a vector of specific results files to use
#' @return List of plots
#' @export

generate_qq_plot <- function(
    results_dir=NULL,
    prefix=NULL,
    plot_dir,
    results_files=NULL
){
  if(!is.null(results_files)){
    cell_type_files=results_files
  }else{
    all_files <- list.files(results_dir, pattern = prefix, full.names = T)
    
    cell_type_files <- all_files[!grepl("summary|vep|genotypes|variant_counts_per_celltype",all_files)]
    cell_type_files <- cell_type_files[grep(".tsv$",cell_type_files)]
    cell_type_files <- cell_type_files[!grepl("pass_only|priority", cell_type_files)]
  }
  
  message("Pulling p values from ", length(cell_type_files), " results files including:")
  for(i in cell_type_files) message(i)
  
  null_p_values <- foreach(file=cell_type_files, .combine = "bind_rows") %do% {
    fread(file) %>%
      filter(priority_flag==0,
             str_detect(Consequence, "synonymous_variant|intron_variant|intergenic_variant|upstream_gene_variant|downstream_gene_variant|3_prime_UTR_variant|5_prime_UTR_variant") | Consequence=="",
             !str_detect(Consequence, "missense_variant|frameshift_variant|stop_gained|stop_lost|start_lost|splice"),
             #!str_detect(filter,"min_varcount_total|min_varportion_total|max_varportion_total|min_datacount_total|min_dataportion_total|min_quality_mean_total|min_depth_mean_total|min_allelefreq_mean_total|min_allelefreq_mean_mutants|varcount_celltype"
             data_proportion_total >= 0.5,
             alt_proportion_total < 0.01,
             data_cnt_ct >= 10,
             data_cnt_other >= 10,
             ref_cnt_ct >= 5,
             ref_cnt_other >= 5,
      )
  }
  
  n_p <- length(null_p_values$p)
  
  if(n_p==0){
    stop("No p values found after filtering, check input files")
  }else{
    message("Found ", n_p, " presumed null p values")
  }
  
  
  hist(null_p_values$p, breaks = 20)
  
  p_histogram <- recordPlot()

  hist(null_p_values$alt_cnt_ct, breaks = 20)
  
  alt_cnt_ct_histogram <- recordPlot()
  
  p <- null_p_values$p
  p <- p[is.finite(p) & !is.na(p) & p > 0 & p <= 1]
  
  obs <- sort(p)
  m <- length(obs)
  exp <- (1:m) / (m + 1)
  
  plot_file = paste0(plot_dir,"/",prefix,"_QQ_plot.png")
  
  plot(
    -log10(exp),
    -log10(obs),
    pch = 16, cex = 0.5,
    xlab = "Expected -log10(p)",
    ylab = "Observed -log10(p)"
  )
  abline(0, 1, col = "red", lwd = 2)
  
  qq_plot <- recordPlot()
  
  png(plot_file, res = 300, width = 1000, height = 1000)
  print(qq_plot)
  dev.off()
  
  message("Saved QQ plot to ", plot_file)
  
  return(list(qq_plot=qq_plot,p_histogram=p_histogram,alt_cnt_ct_histogram=alt_cnt_ct_histogram))
}


















