#' find_somatic_variants.R
#'
#' Find variants enriched in one or more annotated cell types from a Tapestri
#' H5 or a multi-sample VCF. Exactly one of `h5_in` or `vcf_in` must be supplied.
#' Writes summary tables, per-cell genotypes, and (optionally) per-cell-type
#' enrichment stats to `out_dir` using `file_prefix`.
#'
#' @param h5_in Path to input MissionBio Tapestri `.h5` file. Provide this OR `vcf_in`.
#' @param vcf_in Path to input multi-sample VCF (`.vcf.gz`) with `.tbi` index. If a plain `.vcf`
#'   is provided and `prep_vcf=TRUE`, it will be bgzipped and indexed into `prep_vcf_dir`.
#' @param prep_vcf Logical; if `TRUE` and `vcf_in` is an uncompressed `.vcf`, create a
#'   bgzipped and tabix-indexed copy before processing.
#' @param prep_vcf_dir Directory to write the prepared `*.vcf.gz` and `*.tbi` when
#'   `prep_vcf=TRUE`. Defaults to the directory of `vcf_in` if `NULL`.
#' @param cell_type Character vector of cell-type labels to test for variant enrichment. Use `"all"`
#'   to test every cell type present in `cell_annotations`. Only types present are used.
#' @param cell_annotations Named vector (length = number of cells) mapping barcode
#'   (names) → cell-type (values). Barcodes are matched to VCF/H5 samples; common
#'   prefix/suffix is auto-trimmed if needed.
#' @param out_dir Output directory for all result files.
#' @param file_prefix Prefix used to name output files (e.g., `"<prefix>_summary.tsv"`).
#' @param cell_filter A vector of barcodes/cell_IDs to include
#' @param skip_vep Logical; if `TRUE`, skip VEP annotation and priority variant
#'   detection. Disables `priority_only` filtering.
#' @param priority_only Logical; if `TRUE`, restrict cell-type enrichment testing and
#'   outputs to passing variants flagged as priority by VEP (rare protein-impacting).
#'   VEP is only run on variants passing the variant-level filters. Ignored when `skip_vep=TRUE`.
#' @param pass_only_variants Logical; if `TRUE`, restricts downstream outputs to variants
#'   which pass the variant-level filters to speed up processing time.
#' @param block Integer; number of variants per processing chunk for memory efficiency.
#'
#' @param min_DP_per_call Minimum depth (DP) for a per-cell genotype to be considered;
#'   calls failing are masked (`AF=NA`, `NGT=3`).
#' @param min_GQ_per_call Minimum genotype quality (GQ) for a per-cell genotype to be
#'   considered; calls failing are masked (`AF=NA`, `NGT=3`).
#' @param min_AF_per_call Minimum AF for a mutant call (NGT 1/2); failing calls are
#'   masked as missing (`AF=NA`, `NGT=3`) rather than set to reference. (Default = 0)
#' @param min_altreads_per_call Minimum alt reads (AF * DP) for a mutant call; failing
#'   calls are masked as missing. (Default = 0)
#'
#' @param min_varcount_total Minimum number of mutant cells (NGT %in% 1,2) across all cells.
#' @param min_varportion_total Minimum fraction of mutant among informative cells
#'   (alt_cnt_total / data_cnt_total) required globally.
#' @param max_varportion_total Maximum allowed global mutant fraction (to exclude
#'   likely germline/high-frequency sites).
#' @param min_datacount_total Minimum number of cells with informative data globally.
#' @param min_dataportion_total Minimum fraction of cells with informative data globally.
#' @param min_quality_mean_total Minimum mean GQ across cells with AF>0.
#' @param min_allelefreq_mean_total Minimum mean AF across cells with AF>0.
#' @param min_allelefreq_mean_mutants Minimum mean AF among mutant cells only
#'   (NGT %in% 1,2).
#' @param rare_cutoff MAX_AF threshold used with VEP (e.g., 0.01–0.05). Variants with
#'   MAX_AF < cutoff or with missing MAX_AF are treated as “rare/unreported” for
#'   priority calls.
#'
#' @param min_varcount_celltype Minimum number of mutant cells within the focal cell type.
#' @param min_varportion_celltype Minimum mutant fraction within the focal cell type
#'   (alt_cnt_ct / data_cnt_ct).
#' @param min_datacount_celltype Minimum number of informative cells within the focal type.
#' @param min_varportion_in_celltype Minimum fraction of all mutant cells that fall
#'   inside the focal type (alt_cnt_ct / alt_cnt_total), i.e., enrichment within type.
#'
#' @param max_p_value Maximum adjusted p-value (BH) to keep in cell-type enrichment results.
#' @param odds_ratio_cutoff Minimum odds ratio threshold for enrichment (set negative
#'   to disable).
#' @param zscore_cutoff Minimum Wald Z-score threshold for enrichment (set ≤0 to disable).
#'
#' @param celltype_sort_value Sorting key for per-cell-type result tables; `"p"`
#'   sorts by adjusted p-value ascending (default). `"z"` sorts by descending z-score.
#'   "cooc" sorts by co-occurence score
#'
#' @param overwrite_variant_summary Logical; overwrite existing global summary files
#'   (`*_summary.tsv`, `*_summary_pass_only.tsv`, etc.) if present.
#' @param overwrite_celltype_enrichment Logical; overwrite existing per-cell-type files
#'   if present.
#' @param run_cell_type_enrichment Logical; if `TRUE`, compute per-cell-type counts,
#'   odds ratios, p-values, adjust p-values, and write `"<prefix>_<celltype>.tsv"`
#'   and filtered variants.
#'
#' @param genome_version Genome identifier passed to VariantAnnotation when reading VCF
#'   (e.g., `"hg19"`, `"GRCh38"`). Must match contigs in the input VCF.
#' @param threads Number of threads to use where parallelism is available. (Most
#'   operations here are chunked; heavy parallelism occurs upstream/downstream.)
#' @param use_vep_file Path to existing VEP annotation file to use
#' @param AF_in_percent Force convert AF from 0-100 scale to 0-1 scale
#' @param p_adjust_method P adjustment method passed to stats::p.adjust. 

#' @param run_cooccurrence Logical; score passing variants by how many other variants
#'   share their carrier cells, against a shuffled null. Adds cooc_* columns to the
#'   summary and per-cell-type files. (Default = FALSE)
#' @param cooc_min_af,cooc_min_altreads Thresholds for a confident alt call in
#'   co-occurrence scoring. (Defaults = 0.3, 2)
#' @param cooc_min_dp_ref Minimum depth for a confident reference call in co-occurrence
#'   scoring. (Default = 5)
#' @param cooc_min_carriers Minimum confident carriers for a variant to be scored. (Default = 3)
#' @param cooc_jaccard,cooc_overlap Two variants are partners if their carrier sets have
#'   Jaccard >= cooc_jaccard or overlap coefficient >= cooc_overlap (allowing nested
#'   subclones), computed over jointly covered cells. (Defaults = 0.6, 0.9)
#' @param cooc_min_distance Variants closer than this (bp, same chromosome) are not counted
#'   as partners, avoiding shared-read and same-amplicon artefacts. (Default = 1000)
#' @param cooc_max_reference Maximum variants used as the comparison set; larger sets are
#'   randomly sampled. (Default = 50000)
#' @param cooc_seed Seed for reference sampling and the shuffled null. (Default = 1)
#' @param max_cooc_p Optional maximum cooc_p for per-cell-type results. Variants that
#'   could not be scored are not filtered. (Default = NULL)
#' @export

find_somatic_variants <- function(h5_in=NULL,
                                  vcf_in=NULL,
                                  prep_vcf=F,
                                  prep_vcf_dir=NULL,
                                  cell_type="all",
                                  cell_annotations=NULL,
                                  out_dir,
                                  file_prefix,
                                  skip_vep = F,
                                  priority_only=F,
                                  pass_only_variants=FALSE,
                                  block = 1000,
                                  # Call-Level thresholds
                                  min_DP_per_call= 10, 
                                  min_GQ_per_call= 30,
                                  min_AF_per_call = 0.2,
                                  min_altreads_per_call = 0,
                                  # Variant-level thresholds
                                  min_varcount_total        = 30,
                                  min_varportion_total      = 0,
                                  max_varportion_total      = 0.7,
                                  min_datacount_total       = 100,
                                  min_dataportion_total     = 0.25,
                                  min_quality_mean_total    = 0.3,
                                  min_depth_mean_total      = 5,
                                  min_allelefreq_mean_total = 0,
                                  min_allelefreq_mean_mutants = 0.1,
                                  rare_cutoff = 0.02,
                                  # celltype-specific thresholds
                                  min_varcount_celltype     = 1,
                                  min_varportion_celltype   = 0,
                                  min_datacount_celltype    = 1,
                                  min_varportion_in_celltype = 0,
                                  # Somatic / statistical thresholds
                                  max_p_value = 1,
                                  odds_ratio_cutoff = -1000,
                                  zscore_cutoff     = 0,
                                  # Sorting options
                                  celltype_sort_value = "p",
                                  # other
                                  overwrite_variant_summary=TRUE,
                                  overwrite_celltype_enrichment=TRUE,
                                  run_cell_type_enrichment=TRUE,
                                  genome_version="hg19",
                                  threads=4,
                                  use_vep_file=NULL,
                                  AF_in_percent=NULL,
                                  cell_filter=NULL,
                                  p_adjust_method="BY",
                                  # Co-occurrence scoring
                                  run_cooccurrence = FALSE,
                                  cooc_min_af = 0.3,
                                  cooc_min_altreads = 2,
                                  cooc_min_dp_ref = 5,
                                  cooc_min_carriers = 3,
                                  cooc_jaccard = 0.6,
                                  cooc_overlap = 0.9,
                                  cooc_min_distance = 1000,
                                  cooc_max_reference = 50000,
                                  cooc_seed = 1,
                                  max_cooc_p = NULL
                          ) {

  if(is.null(cell_annotations) & run_cell_type_enrichment){
    stop("ERROR: Cell type enrichment requires cell annotations!")
  }
  
  #########################################################
  # ------------------- Sanity Checks ------------------- #
  #########################################################
  
  if (isTRUE(priority_only) && isTRUE(skip_vep)) {
    warning("priority_only ignored because skip_vep=TRUE")
    priority_only <- FALSE
  }

  if (!is.null(max_cooc_p) && !run_cooccurrence) stop("max_cooc_p requires run_cooccurrence = TRUE")
  if (celltype_sort_value == "cooc" && !run_cooccurrence) stop("celltype_sort_value = 'cooc' requires run_cooccurrence = TRUE")
  
  ############################################################
  # ------------------- Check Cell Types ------------------- #
  ############################################################
  
  if(!is.null(cell_annotations)){
    all_cell_types <- sort(unique(cell_annotations))
    
    if (!identical(cell_type, "all")) {
      analysis_cell_types <- intersect(all_cell_types, cell_type)
    }else{
      analysis_cell_types <- all_cell_types
    }
    
    if (length(analysis_cell_types)==0) stop("Selected cell type doesn't match available data")
  }
  
  #######################################################
  # ------------------- Check Paths ------------------- #
  #######################################################
  
  # stores summary + genotype info for all cell barcodes
  genotype_file <- file.path(out_dir,paste0(file_prefix,"_genotypes.tsv"))
  # stores only the summary data for all variants
  summary_file <- file.path(out_dir,paste0(file_prefix,"_summary.tsv"))
  # stores only the summary data for passing variants
  summary_pass_file <- file.path(out_dir,paste0(file_prefix,"_summary_pass_only.tsv"))
  # stores only the summary data for passing priority variants
  summary_pass_priority_file <- file.path(out_dir,paste0(file_prefix,"_summary_pass_only_priority.tsv"))
  
  counts_per_cell_type_file <- file.path(out_dir,paste0(file_prefix,"_variant_counts_per_celltype.tsv"))
  
  for(path in c(genotype_file,summary_file,summary_pass_file,summary_pass_priority_file,counts_per_cell_type_file)){
    if(file.exists(path)){
      if(overwrite_variant_summary){
        message(paste0("Removing ", path, "..."))
        file.remove(path)
      }else{
        stop(paste0(path," already exists!"))
      }
    }
  }
  
  if(!is.null(cell_annotations)){
    for (ct in analysis_cell_types) {
      # stores cell type-specific stats for all variants
      ct_file <- file.path(out_dir,paste0(file_prefix,"_",ct,".tsv"))
      # stores cell type-specific stats for passing variants
      ct_pass_file <- file.path(out_dir,paste0(file_prefix,"_",ct,"_pass_only.tsv"))
      # stores cell type-specific stats for passing priority variants
      ct_pass_priority_file <- file.path(out_dir,paste0(file_prefix,"_",ct,"_pass_only_priority.tsv"))
      
      for(path in c(ct_file,ct_pass_file,ct_pass_priority_file)){
        if(file.exists(path)){
          if(overwrite_celltype_enrichment){
            message(paste0("Removing ", path, "..."))
            file.remove(path)
          }else{
            stop(paste0(path," already exists!"))
          }
        }
      }
    }
  }
  vep_file <- file.path(out_dir, paste0(file_prefix, "_vep_annotated.vcf"))
  vep_in_file <- file.path(out_dir, paste0(file_prefix, "_vep_in.tsv"))
  
  #############################################################
  # ------------------- Get Cell Barcodes ------------------- #
  #############################################################
  
  if(!is.null(vcf_in)){
    hdr <- scanVcfHeader(vcf_in)
    barcodes_full <- samples(hdr)
    barcodes <- barcodes_full
  }else if(!is.null(h5_in)){
    barcodes_full <- as.vector(rhdf5::h5read(h5_in, "/assays/dna_variants/ra/barcode"))
    barcodes <- barcodes_full
  }
    
  if(!is.null(cell_annotations)){
    if (length(barcodes) != length(cell_annotations)) stop("cell annotations are not the same length as vcf samples")
    
    # attempt to trim barcodes if needed
    if (!identical(sort(names(cell_annotations)),sort(barcodes))) {
      message("variant barcodes don't match cell annotations")
      message("attemtping to find match in substring...")
      barcode1 <- barcodes[1]
      cell1 <- names(cell_annotations)[1]
      match <- F
      for(i in 1:nchar(barcode1)){
        start_sub <- substr(barcode1,1,i)
        end_sub <- substr(barcode1,nchar(barcode1)-i+1,nchar(barcode1))
        if(start_sub==cell1){
          match <- T
          barcode_suffix <- substr(barcode1, i+1, nchar(barcode1))
          barcodes <- gsub(barcode_suffix,"",barcodes)
          if (identical(sort(names(cell_annotations)),sort(barcodes))){
            message(paste0("substring match found, trimming ", barcode_suffix," suffix from barcodes"))
          }
        }else if (end_sub==cell1){
          match <- T
          barcode_prefix <- substr(barcode1, 1, nchar(barcode1)-i)
          barcodes <- gsub(barcode_prefix,"",barcodes,fixed = TRUE)
          if (identical(sort(names(cell_annotations)),sort(barcodes))){
            message(paste0("substring match found, trimming ", barcode_prefix," prefix from barcodes"))
          }
        }
        if(match){
          break
        }
      }
      if(!match){
        stop("no substring match found")
      }
    }
    # put annotations in the same order as the VCF/h5 samples (cells are matched by position below)
    if (!setequal(names(cell_annotations), barcodes)) stop("cell_annotations names don't match the VCF/h5 sample names")
    cell_annotations <- cell_annotations[barcodes]

  }
  
  if(!is.null(cell_filter)){
    if(!all(cell_filter %in% names(cell_annotations))){
      stop("cell_filter vector must match cell names/barcodes in cell_annotations")
    }
    keep_cells <- which(names(cell_annotations) %in% cell_filter)
  } else {
    keep_cells <- seq_along(barcodes)
  }
  
  barcodes <- barcodes[keep_cells]
  cell_annotations <- cell_annotations[keep_cells]
  ct_idx <- split(seq_along(cell_annotations), cell_annotations)
  
  
  ##############################################################
  # ------------------- Normalise Variants ------------------- #
  ##############################################################
  
  message("Processing per-variant info...")
  
  if(!is.null(vcf_in)){
    
    if(prep_vcf){
      message("Generating compressed vcf and indexing...")
      out_vcf <- file.path(prep_vcf_dir,gsub(".vcf",".vcf.gz",basename(vcf_in)))
      #module load htslib/1.16
      system(sprintf("bgzip -c %s > %s",vcf_in,out_vcf))
      system(sprintf("tabix -p vcf %s",out_vcf))
      vcf_in <- out_vcf
    }else if(!str_detect(vcf_in,".vcf.gz")){
      stop("vcf.gz required, use prep_vcf=TRUE")
    }else if(!file.exists(gsub(".vcf.gz|.vcf",".vcf.gz.tbi",vcf_in))){
      stop("vcf index required, use prep_vcf=TRUE")
    }
    
    tabix <- TabixFile(vcf_in, yieldSize = block*100)
    
    param <- ScanVcfParam(
      fixed = c("ALT"),
      info = NA,     
      geno = NA
    )
    
    if(isOpen(tabix)){
      close(tabix)
    }
    
    open(tabix)
    variant_info_chunks <- list()
    repeat {
      v <- readVcf(tabix, genome = genome_version, param = param)
      if (nrow(v) == 0) break
      vcf <- expand(v)
      rr   <- rowRanges(vcf)
      
      # ---- base columns ----
      CHROM  <- as.character(seqnames(rr))
      POS    <- start(rr)
      REF    <- as.character(ref(vcf))
      ALT    <- as.character(alt(vcf))
      #ALT <- vapply(
      #  as.list(alt(vcf)),
      #  function(x) paste(as.character(x), collapse = ","),
      #  character(1)
      #)
      
      ALT[is.na(ALT) | ALT == ""] <- "*"

      variant_info_chunks[[length(variant_info_chunks) + 1]] <- data.frame(CHROM,POS,REF,ALT)
      
    }
    close(tabix)
    
    variant_info <- bind_rows(variant_info_chunks)
    rm(variant_info_chunks)
    
    variants <- .normalise_vcf(variant_info)
    rm(variant_info)

    
    
  }else if(!is.null(h5_in)){
    variant_info <- h5read(h5_in, "/assays/dna_variants/ca")[c("CHROM","POS","REF","ALT")] %>% as.data.frame()
    variants <- .normalise_vcf(variant_info)
    rm(variant_info)
  }
  
  #######################################################################
  # ------------------- Filter Pass 1  (Formatting) ------------------- #
  #######################################################################
  
  # store which variants are being kept
  is_keep <- !variants$filter
  keep_idx <- which(is_keep)
  
  if(any(!is_keep)){
    message(paste0("WARNING: ", length(which(!is_keep)),"/",nrow(variants), " variants are being filtered out due to spurious formatting"))
  }
  
  # filter out unwanted variants
  variants <- variants %>%
    filter(!filter) %>%
    dplyr::select(-filter)

  # Add a plotting ID. If using VEP, this will get overwritten later if variant is protein-coding
  variants <- variants %>%
    mutate(plot_ID=variant)
  
  if(!skip_vep){
    # passing variants are appended to this file inside the loop
    file.create(vep_in_file)
  }
  
  #################################################################
  # ------------------- Process Genotype Info ------------------- #
  #################################################################
  
  message("Processing per-cell info...")
  
  # write header for cell-specific file
  writeLines("#FORMAT=AF;DP;GQ;NGT", genotype_file)
  
  cooc_calls <- list(); cooc_info <- list(); cooc_counts <- list()

  if(!is.null(vcf_in)){
    
    nr <- nrow(variants)
    pass_total <- 0L

    hdr <- scanVcfHeader(vcf_in)
    ad_is_scalar <- geno(hdr)["AD","Number"] == "1"
    has_GQ <- "GQ" %in% rownames(geno(hdr))
    geno_fields <- c("AD","DP","GT", if (has_GQ) "GQ")
    
    tabix <- TabixFile(vcf_in, yieldSize = block)
    
    if(isOpen(tabix)){
      close(tabix)
    }
    
    param <- ScanVcfParam(
      fixed = c("ALT"),
      info = NA,     
      geno = geno_fields
    )
    
    open(tabix)
    
    raw_row_counter <- 0
    variants_row_counter <- 0
    first_block <- TRUE
    last_keep_row <- max(keep_idx)
    repeat {
      if (raw_row_counter >= last_keep_row && variants_row_counter >= nr) break
      # read in vcf block
      v <- readVcf(tabix, genome = genome_version, param = param)
      if (nrow(v) == 0) break
      # expand multi-allelic vcf
      vcf <- expand(v)
      
      # find rows to match variants of interest
      current_rows <- (raw_row_counter+1):(raw_row_counter+nrow(vcf))
      raw_row_counter <- raw_row_counter + nrow(vcf)
      keep_rows <- which(is_keep[current_rows])
      
      if (length(keep_rows) == 0L) {
        if (variants_row_counter >= nr) break
        message(sprintf(
          "Processing variants %d-%d / %d (0 kept in this chunk)",
          min(variants_row_counter + 1L, nr), variants_row_counter, nr
        ))
        next
      }
      
      # find which rows to pull from variants df
      variant_rows <- (variants_row_counter+1):(variants_row_counter+length(keep_rows))
      
      message(paste0("Processing variants ", variants_row_counter+1, "-", min(variants_row_counter + length(keep_rows), nr), " / ", nr))
      
      variants_row_counter <- variants_row_counter + length(keep_rows)
      
      # get variant data for the block
      block_variants <- variants[variant_rows, , drop = FALSE]
      # empty VEP columns, filled in after the loop for passing variants
      if (!skip_vep) block_variants <- .add_empty_vep_cols(block_variants)

      
      if (ad_is_scalar) {
        ADb_alt <- geno(vcf)$AD[keep_rows,,drop=FALSE]
        trueDPb <- geno(vcf)$DP[keep_rows,,drop=FALSE]
        AFb <- ADb_alt/trueDPb
      } else {
        # read in AD array and split to ref/alt matrices
        ADb <- geno(vcf)$AD[keep_rows,,1:2,drop=FALSE]
        ADb_ref <- ADb[,,1, drop=FALSE]
        dim(ADb_ref) <- dim(ADb)[1:2]
        ADb_alt <- ADb[,,2, drop=FALSE]
        dim(ADb_alt) <- dim(ADb)[1:2]
        # calculate DP from AD ref+alt
        trueDPb <- ADb_ref + ADb_alt
        # calculate AF
        AFb <- ADb_alt/trueDPb
      }


      # read in GQ
      if (has_GQ) {
        GQb <- geno(vcf)$GQ[keep_rows,,drop=FALSE]
      } else {
        GQb <- matrix(NA_real_, nrow = length(keep_rows), ncol = length(barcodes_full))
      }

      # read in GT
      GTb <- geno(vcf)$GT[keep_rows,,drop=FALSE]
      
      NGTb <- .pull_expanded_ngt_from_vcf(v)[keep_rows,,drop=FALSE]
      
      dimnames(AFb) <- dimnames(trueDPb) <- dimnames(GQb) <- dimnames(NGTb) <- list(NULL, barcodes_full)
      
      
      AFb     <- AFb[, keep_cells, drop = FALSE]
      trueDPb <- trueDPb[, keep_cells, drop = FALSE]
      GQb     <- GQb[, keep_cells, drop = FALSE]
      NGTb    <- NGTb[, keep_cells, drop = FALSE]
      
      dimnames(AFb) <- dimnames(trueDPb) <- dimnames(GQb) <- dimnames(NGTb) <- list(NULL, barcodes)
      
      # Apply per-call filters
      if (!has_GQ){
          message("Warning: min_GQ_per_call and min_quality_mean_total are being ignored because vcf_in does not contain GQ field")
          min_quality_mean_total <- -Inf
      }
      pass_call <- (trueDPb >= min_DP_per_call) & (if (has_GQ) GQb >= min_GQ_per_call else TRUE)
      low_af_alt <- (NGTb == 1L | NGTb == 2L) & is.finite(AFb) &
                    (AFb < min_AF_per_call | round(AFb * trueDPb) < min_altreads_per_call)
      low_af_alt[is.na(low_af_alt)] <- FALSE
      pass_call <- pass_call & !low_af_alt
      AFb[!pass_call] <- NA_real_
      GQb[!pass_call] <- NA_real_
      NGTb[!pass_call] <- 3L
      
      # calculate variant-level stats
      var_stats <- .calculate_variant_level_stats(NGT = NGTb,AF = AFb, GQ = GQb, DP = trueDPb)

      filter_counts <- c(
        min_varcount_total           = sum(var_stats$alt_cnt_total        < min_varcount_total, na.rm = TRUE),
        min_varportion_total         = sum(var_stats$alt_proportion_total < min_varportion_total, na.rm = TRUE),
        max_varportion_total         = sum(var_stats$alt_proportion_total > max_varportion_total, na.rm = TRUE),
        min_datacount_total          = sum(var_stats$data_cnt_total       < min_datacount_total, na.rm = TRUE),
        min_dataportion_total        = sum(var_stats$data_proportion_total< min_dataportion_total, na.rm = TRUE),
        min_quality_mean_total       = sum(var_stats$mean_GQ_total        < min_quality_mean_total, na.rm = TRUE),
        min_depth_mean_total         = sum(var_stats$mean_DP_total        < min_depth_mean_total, na.rm = TRUE),
        min_allelefreq_mean_total    = sum(var_stats$mean_AF_total        < min_allelefreq_mean_total, na.rm = TRUE),
        min_allelefreq_mean_mutants  = sum(var_stats$mean_AF_mutants      < min_allelefreq_mean_mutants, na.rm = TRUE)
      )

      message("Block filter fails (n=", nrow(var_stats), "): ",
              paste(names(filter_counts)[filter_counts > 0],
                    filter_counts[filter_counts > 0], sep = "=", collapse = ", "))
      
      # assemble lines

      block_summary <- bind_cols(block_variants, var_stats)
      block_summary$filter <- .make_filter(cbind(
        min_varcount_total          = var_stats$alt_cnt_total        < min_varcount_total,
        min_varportion_total        = var_stats$alt_proportion_total < min_varportion_total,
        max_varportion_total        = var_stats$alt_proportion_total > max_varportion_total,
        min_datacount_total         = var_stats$data_cnt_total       < min_datacount_total,
        min_dataportion_total       = var_stats$data_proportion_total< min_dataportion_total,
        min_quality_mean_total      = var_stats$mean_GQ_total        < min_quality_mean_total,
        min_depth_mean_total        = var_stats$mean_DP_total        < min_depth_mean_total,
        min_allelefreq_mean_total   = var_stats$mean_AF_total        < min_allelefreq_mean_total,
        min_allelefreq_mean_mutants = var_stats$mean_AF_mutants      < min_allelefreq_mean_mutants
      ))

      if (run_cooccurrence) block_summary <- .add_empty_cooc_cols(block_summary)

      
      if (isTRUE(pass_only_variants)) {
        pass_idx       <- which(block_summary$filter == ".")
        n_block        <- nrow(block_summary)
        message(length(pass_idx), "/", n_block, " passed variant-level filters in this block")
        pass_total <- pass_total + length(pass_idx)
        message(pass_total, " variants passed so far...")
        block_summary  <- block_summary[pass_idx, , drop = FALSE]
        NGTb           <- NGTb[pass_idx, , drop = FALSE]
        var_stats      <- var_stats[pass_idx, , drop = FALSE]
        block_variants <- block_variants[pass_idx, , drop = FALSE]
        AFb            <- AFb[pass_idx, , drop = FALSE]
        trueDPb        <- trueDPb[pass_idx, , drop = FALSE]
        GQb            <- GQb[pass_idx, , drop = FALSE]
        if (nrow(block_summary) == 0L) next
      }

      if (run_cooccurrence) {
        pr <- which(block_summary$filter == ".")
        if (length(pr)) {
          calls <- .cooc_calls(AFb[pr, , drop = FALSE], trueDPb[pr, , drop = FALSE],
                               cooc_min_af, cooc_min_altreads, cooc_min_dp_ref)
          n_alt <- rowSums(calls == as.raw(2L))
          cooc_counts[[length(cooc_counts) + 1]] <- data.frame(variant = block_summary$variant[pr],
                                                                cooc_n_carriers = n_alt)
          ok <- n_alt >= cooc_min_carriers
          if (any(ok)) {
            cooc_calls[[length(cooc_calls) + 1]] <- calls[ok, , drop = FALSE]
            cooc_info[[length(cooc_info) + 1]] <- block_summary[pr[ok], c("variant", "chromosome", "start_position")]
          }
        }
      }
      
      # write variant info to summary file

      write_tsv(block_summary, summary_file, append = TRUE, col_names = first_block)
      
      write_tsv(block_summary %>%
                  filter(filter==".") %>%
                  dplyr::select(-filter), summary_pass_file, append = TRUE, col_names = first_block)
      
      # VEP input for passing variants (VEP runs once after the loop)
      if (!skip_vep) {
        write_tsv(block_summary %>%
                    filter(filter==".") %>%
                    mutate(allele=paste0(reference_allele,"/",alternate_allele),
                           strand="+") %>%
                    dplyr::select(chromosome, start_position, end_position, allele, strand, variant),
                  vep_in_file, append = TRUE, col_names = FALSE)
      }
      
      # write counts per cell type
      if(!is.null(cell_annotations)){
        alt <- (NGTb == 1L | NGTb == 2L)
        dat <- alt | (NGTb == 0L)
        counts_per_ct <- do.call(cbind, lapply(names(ct_idx), function(ct) {
          j <- ct_idx[[ct]]
          m <- cbind(rowSums(alt[, j, drop = FALSE]), rowSums(dat[, j, drop = FALSE]))
          colnames(m) <- paste0(ct, c("_alt_count", "_data_count"))
          m
        })) %>% as.data.frame()

        counts_per_ct <- bind_cols(block_variants, counts_per_ct)
        write_tsv(counts_per_ct, counts_per_cell_type_file, append = TRUE, col_names = first_block)
      }
            
      
      # get per cell info
      
      block_per_cell_info <- array(
        sprintf("%s;%s;%s;%s", round(as.vector(AFb),3), as.vector(trueDPb), as.vector(GQb), as.vector(NGTb)),
        dim = dim(AFb), dimnames = list(NULL, barcodes)
      )
      block_per_cell_info <- as.data.frame(block_per_cell_info, check.names = FALSE)
      
      block_per_cell_info <- block_summary %>%
        dplyr::select(-dplyr::starts_with("cooc_")) %>%
        bind_cols(block_per_cell_info)
      
      data.table::fwrite(block_per_cell_info, genotype_file, sep = "\t", append = TRUE, col.names = first_block, quote = FALSE)
      
      if(run_cell_type_enrichment){
        
      
       # per cell type stats
        for (ct in analysis_cell_types) {
          
          # stores cell type-specific stats for all variants
          ct_file <- file.path(out_dir,paste0(file_prefix,"_",ct,".tsv"))
          
          ct_stats <- .calculate_cell_type_stats(NGTb,ct,cell_annotations,var_stats)
          
          stats_block_ct <- bind_cols(block_summary, ct_stats)

          add <- .make_filter(cbind(
            varcount_celltype      = ct_stats$alt_cnt_ct           < min_varcount_celltype,
            varportion_celltype    = ct_stats$alt_proportion_ct    < min_varportion_celltype,
            datacount_celltype     = ct_stats$data_cnt_ct          < min_datacount_celltype,
            varportion_in_celltype = ct_stats$alt_proportion_in_ct < min_varportion_in_celltype
          ))

          f <- stats_block_ct$filter
          stats_block_ct$filter <- ifelse(add == ".", f,
                                    ifelse(f == ".", add, paste0(f, ";", add)))
          
          write_tsv(stats_block_ct, ct_file, append = TRUE, col_names = first_block)
          
        }
      }
        
      first_block <- FALSE
    }
    close(tabix)
    
  }else if(!is.null(h5_in)){
    
    pass_total <- 0L
    # working with HDF5Array for memory efficiency 
    
    ds_AF  <- "/assays/dna_variants/layers/AF"
    ds_DP  <- "/assays/dna_variants/layers/DP"
    ds_GQ  <- "/assays/dna_variants/layers/GQ"
    ds_NGT <- "/assays/dna_variants/layers/NGT"
    
    AF_raw <- HDF5Array(h5_in, ds_AF)[keep_idx, keep_cells, drop = FALSE]
    DP     <- HDF5Array(h5_in, ds_DP)[keep_idx, keep_cells, drop = FALSE]
    GQ     <- HDF5Array(h5_in, ds_GQ)[keep_idx, keep_cells, drop = FALSE]
    NGT    <- HDF5Array(h5_in, ds_NGT)[keep_idx, keep_cells, drop = FALSE]
    
    nr <- nrow(variants)
    nc <- ncol(AF_raw)
    
    dimnames(AF_raw) <- dimnames(DP) <- dimnames(GQ) <- dimnames(NGT) <- list(NULL, barcodes)
    
    
    if(is.null(AF_in_percent)){
      # Probe AF to find scale
      probe_rows <- seq_len(min(500, nr))
      probe_cols <- seq_len(min(500, nc))
      af_probe   <- as.matrix(AF_raw[probe_rows, probe_cols, drop = FALSE])
      af_in_percent <- any(is.finite(af_probe) & af_probe > 1)
    }else{
      af_in_percent=AF_in_percent
    }
    
    
    
    # work in blocks for memory efficiency
    message(paste0("Found ", nr, " variants"))
    message(paste0("Working in blocks of ", block, " variants"))
    first_block <- TRUE
    # ---------- Stream row blocks spanning ALL columns ----------
    for (i0 in seq(1L, nr, by = block)) {
      
      message(paste0("Processing variants ", i0, "-", min(i0 + block - 1L, nr), " / ", nr))
      
      rows <- i0:min(i0 + block - 1L, nr)
      
      # get variant data for the block
      block_variants <- variants[rows, , drop = FALSE]
      # empty VEP columns, filled in after the loop for passing variants
      if (!skip_vep) block_variants <- .add_empty_vep_cols(block_variants)

      
      # materialize only this block into RAM
      AFb_raw <- as.matrix(AF_raw[rows, , drop = FALSE])# used for per-cell "AF="
      AFb     <- if (af_in_percent) AFb_raw / 100 else AFb_raw  # 0–1 for stats
      DPb     <- as.matrix(DP[rows, , drop = FALSE])
      GQb     <- as.matrix(GQ[rows, , drop = FALSE])
      NGTb    <- as.matrix(NGT[rows, , drop = FALSE])
      
      
      # double check AF in percent
      
      if(any(is.finite(AFb) & AFb > 1)){
        stop("ERROR: Detected AF > 1 that was missed in initial probe. Please set AF_in_percent=TRUE")
      }
      
      # Apply per-call filters

      pass_call <- (DPb >= min_DP_per_call) & (GQb >= min_GQ_per_call)
      low_af_alt <- (NGTb == 1L | NGTb == 2L) & is.finite(AFb) &
                    (AFb < min_AF_per_call | round(AFb * DPb) < min_altreads_per_call)
      low_af_alt[is.na(low_af_alt)] <- FALSE
      pass_call <- pass_call & !low_af_alt
      AFb[!pass_call] <- NA_real_
      GQb[!pass_call] <- NA_real_
      NGTb[!pass_call] <- 3L
      
      var_stats <- .calculate_variant_level_stats(NGT = NGTb,AF = AFb, GQ = GQb, DP = DPb)
      
      filter_counts <- c(
        min_varcount_total           = sum(var_stats$alt_cnt_total        < min_varcount_total, na.rm = TRUE),
        min_varportion_total         = sum(var_stats$alt_proportion_total < min_varportion_total, na.rm = TRUE),
        max_varportion_total         = sum(var_stats$alt_proportion_total > max_varportion_total, na.rm = TRUE),
        min_datacount_total          = sum(var_stats$data_cnt_total       < min_datacount_total, na.rm = TRUE),
        min_dataportion_total        = sum(var_stats$data_proportion_total< min_dataportion_total, na.rm = TRUE),
        min_quality_mean_total       = sum(var_stats$mean_GQ_total        < min_quality_mean_total, na.rm = TRUE),
        min_depth_mean_total         = sum(var_stats$mean_DP_total        < min_depth_mean_total, na.rm = TRUE),
        min_allelefreq_mean_total    = sum(var_stats$mean_AF_total        < min_allelefreq_mean_total, na.rm = TRUE),
        min_allelefreq_mean_mutants  = sum(var_stats$mean_AF_mutants      < min_allelefreq_mean_mutants, na.rm = TRUE)
      )

      message("Block filter fails (n=", nrow(var_stats), "): ",
              paste(names(filter_counts)[filter_counts > 0],
                    filter_counts[filter_counts > 0], sep = "=", collapse = ", "))

      # assemble lines
      block_summary <- bind_cols(block_variants, var_stats)
      block_summary$filter <- .make_filter(cbind(
        min_varcount_total          = var_stats$alt_cnt_total        < min_varcount_total,
        min_varportion_total        = var_stats$alt_proportion_total < min_varportion_total,
        max_varportion_total        = var_stats$alt_proportion_total > max_varportion_total,
        min_datacount_total         = var_stats$data_cnt_total       < min_datacount_total,
        min_dataportion_total       = var_stats$data_proportion_total< min_dataportion_total,
        min_quality_mean_total      = var_stats$mean_GQ_total        < min_quality_mean_total,
        min_depth_mean_total        = var_stats$mean_DP_total        < min_depth_mean_total,
        min_allelefreq_mean_total   = var_stats$mean_AF_total        < min_allelefreq_mean_total,
        min_allelefreq_mean_mutants = var_stats$mean_AF_mutants      < min_allelefreq_mean_mutants
      ))

      if (run_cooccurrence) block_summary <- .add_empty_cooc_cols(block_summary)
      
      if (isTRUE(pass_only_variants)) {
        pass_idx       <- which(block_summary$filter == ".")
        n_block        <- nrow(block_summary)
        message(length(pass_idx), "/", n_block, " passed variant-level filters in this block")
        pass_total <- pass_total + length(pass_idx)
        message(pass_total, " variants passed so far...")
        block_summary  <- block_summary[pass_idx, , drop = FALSE]
        NGTb           <- NGTb[pass_idx, , drop = FALSE]
        var_stats      <- var_stats[pass_idx, , drop = FALSE]
        block_variants <- block_variants[pass_idx, , drop = FALSE]
        AFb            <- AFb[pass_idx, , drop = FALSE]
        DPb            <- DPb[pass_idx, , drop = FALSE]
        GQb            <- GQb[pass_idx, , drop = FALSE]
        if (nrow(block_summary) == 0L) next
      }
      
      if (run_cooccurrence) {
        pr <- which(block_summary$filter == ".")
        if (length(pr)) {
          calls <- .cooc_calls(AFb[pr, , drop = FALSE], DPb[pr, , drop = FALSE],
                               cooc_min_af, cooc_min_altreads, cooc_min_dp_ref)
          n_alt <- rowSums(calls == as.raw(2L))
          cooc_counts[[length(cooc_counts) + 1]] <- data.frame(variant = block_summary$variant[pr],
                                                                cooc_n_carriers = n_alt)
          ok <- n_alt >= cooc_min_carriers
          if (any(ok)) {
            cooc_calls[[length(cooc_calls) + 1]] <- calls[ok, , drop = FALSE]
            cooc_info[[length(cooc_info) + 1]] <- block_summary[pr[ok], c("variant", "chromosome", "start_position")]
          }
        }
      }


      write_tsv(block_summary, summary_file, append = TRUE, col_names = first_block)
      
      write_tsv(block_summary %>%
                  filter(filter==".") %>%
                  dplyr::select(-filter), summary_pass_file, append = TRUE, col_names = first_block)
      
      # VEP input for passing variants (VEP runs once after the loop)
      if (!skip_vep) {
        write_tsv(block_summary %>%
                    filter(filter==".") %>%
                    mutate(allele=paste0(reference_allele,"/",alternate_allele),
                           strand="+") %>%
                    dplyr::select(chromosome, start_position, end_position, allele, strand, variant),
                  vep_in_file, append = TRUE, col_names = FALSE)
      }
      
      # write counts per cell type
      if(!is.null(cell_annotations)){
        alt <- (NGTb == 1L | NGTb == 2L)
        dat <- alt | (NGTb == 0L)
        counts_per_ct <- do.call(cbind, lapply(names(ct_idx), function(ct) {
          j <- ct_idx[[ct]]
          m <- cbind(rowSums(alt[, j, drop = FALSE]), rowSums(dat[, j, drop = FALSE]))
          colnames(m) <- paste0(ct, c("_alt_count", "_data_count"))
          m
        })) %>% as.data.frame()

        counts_per_ct <- bind_cols(block_variants, counts_per_ct)
        write_tsv(counts_per_ct, counts_per_cell_type_file, append = TRUE, col_names = first_block)
      }
      # get per cell info
      
      block_per_cell_info <- array(
        sprintf("%s;%s;%s;%s", round(as.vector(AFb),3), as.vector(DPb), as.vector(GQb), as.vector(NGTb)),
        dim = dim(AFb), dimnames = list(NULL, barcodes)
      )
      block_per_cell_info <- as.data.frame(block_per_cell_info, check.names = FALSE)
      
      block_per_cell_info <- block_summary %>%
        dplyr::select(-dplyr::starts_with("cooc_")) %>%
        bind_cols(block_per_cell_info)
      
      data.table::fwrite(block_per_cell_info, genotype_file, sep = "\t", append = TRUE, col.names = first_block, quote = FALSE)
      
      if(run_cell_type_enrichment){
        
        
        # per cell type stats
      
        for (ct in analysis_cell_types) {
          ct_file <- file.path(out_dir,paste0(file_prefix,"_",ct,".tsv"))
          
          ct_stats <- .calculate_cell_type_stats(NGTb,ct,cell_annotations,var_stats)
          
          stats_block_ct <- bind_cols(block_summary, ct_stats)

          add <- .make_filter(cbind(
            varcount_celltype      = ct_stats$alt_cnt_ct           < min_varcount_celltype,
            varportion_celltype    = ct_stats$alt_proportion_ct    < min_varportion_celltype,
            datacount_celltype     = ct_stats$data_cnt_ct          < min_datacount_celltype,
            varportion_in_celltype = ct_stats$alt_proportion_in_ct < min_varportion_in_celltype
          ))

          f <- stats_block_ct$filter
          stats_block_ct$filter <- ifelse(add == ".", f,
                                    ifelse(f == ".", add, paste0(f, ";", add)))
          
          write_tsv(stats_block_ct, ct_file, append = TRUE, col_names = first_block)
          
        }
        
      }
      first_block <- FALSE
    }
  }


  ###############################################################
  # ----------------- Co-occurrence Scoring ------------------- #
  ###############################################################

  if (run_cooccurrence) {
    message("Scoring variant co-occurrence...")
    counts_df <- bind_rows(cooc_counts)
    if (length(cooc_calls)) {
      scores <- .cooccurrence_scores(do.call(rbind, cooc_calls), bind_rows(cooc_info),
                                     jaccard_thresh = cooc_jaccard, overlap_thresh = cooc_overlap,
                                     min_distance = cooc_min_distance, max_reference = cooc_max_reference,
                                     seed = cooc_seed)
    } else {
      scores <- data.frame(variant = character(), cooc_partners = integer(), cooc_p = numeric())
    }
    cooc_annot <- counts_df %>%
      left_join(scores, by = "variant") %>%
      mutate(cooc_testable = as.integer(cooc_n_carriers >= cooc_min_carriers))
    message(sum(cooc_annot$cooc_testable), "/", nrow(cooc_annot), " passing variants had enough carriers to score")

    .fill_vep_cols(summary_file, cooc_annot, label = "co-occurrence")
    .fill_vep_cols(summary_pass_file, cooc_annot, label = "co-occurrence")
    if (run_cell_type_enrichment) {
      for (ct in analysis_cell_types) {
        .fill_vep_cols(file.path(out_dir, paste0(file_prefix, "_", ct, ".tsv")), cooc_annot, label = "co-occurrence")
      }
    }
    rm(cooc_calls, cooc_info, cooc_counts, counts_df, scores, cooc_annot)
    invisible(gc())
  }


  ###################################################################
  # ------------------- VEP on Passing Variants ------------------- #
  ###################################################################

  if(!skip_vep && file.size(vep_in_file) == 0){
    message("No variants passed the variant-level filters, skipping VEP")
    skip_vep <- TRUE
  }

  if(!skip_vep){
    if(!is.null(use_vep_file)){
      if(file.exists(use_vep_file)){
        vep_run <- .parse_vep_vcf(use_vep_file)
      }
    }else if(file.exists(vep_file) & !overwrite_variant_summary){
      message(paste0("Reading in existing VEP results from ", vep_file))
      vep_run <- .parse_vep_vcf(vep_file)
    }else{
      # Run VEP
      message("Annotating passing variants using ENSEMBL-vep...")

      vep_run <- run_vep(
        vep_in_file,
        output = vep_file,
        return = "data.frame",
        args = c(
          "--force_overwrite",
          "--biotype","--domains",
          "--max_af","--check_existing","--numbers","--regulatory",
          "--sift","b","--polyphen","b", "--hgvs", "--pick_allele",
          "--no_stats"
        )
      )
    }
    # Wrangle VEP Output and find priority variants

    message("Retrieving VEP annotations and finding priority variants...")

    is_rare_or_unreported <- function(s, cutoff) {
      if (is.na(s) || s == "" || grepl("^;*$", s)) return(TRUE)
      vals <- suppressWarnings(as.numeric(unlist(strsplit(s, ";", fixed = TRUE))))
      if (length(vals) == 0 || all(is.na(vals))) return(TRUE)
      any(vals < cutoff, na.rm = TRUE)
    }

    collapse_terms <- function(x) paste(x, collapse = ";")

    # --pick_allele output has one block per variant and no PICK column;
    # VEP files made with --flag_pick_allele still need the PICK filter
    if ("PICK" %in% names(vep_run)) vep_run <- vep_run %>% filter(PICK==1)

    vep_info_cols <- vep_run %>%
      dplyr::select(ID, SYMBOL, Gene, Feature, HGVSp, MAX_AF, Consequence, Existing_variation, BIOTYPE, SIFT, PolyPhen) %>%
      # generate a user-friendly plot label
      mutate(plot_ID=ifelse(HGVSp!="" & SYMBOL!="", paste0(SYMBOL,gsub(".*:p.","",HGVSp)), ID) ) %>%
      mutate(plot_ID=ifelse(plot_ID %in% plot_ID[which(duplicated(plot_ID))], paste0(plot_ID,"(",str_extract(ID,"([^0-9]+$)", group = 1),")"), plot_ID)) %>%
      mutate(plot_ID=ifelse(nchar(plot_ID)>30,
                            ifelse(str_detect(plot_ID, "\\)$"), paste0(substr(plot_ID,1,30),"...)"), paste0(substr(plot_ID,1,30),"...")),
                            plot_ID)) %>%
      mutate(plot_ID=make.unique(plot_ID)) %>%
      group_by(ID,plot_ID) %>%
      summarise(across(everything(), collapse_terms)) %>%
      # decide on priority variants (rare nonsense/missense/frameshift)
      rowwise() %>%
      mutate(priority_flag = ifelse(
        str_detect(Consequence, "missense_variant|stop_gained|stop_lost|start_lost|frameshift_variant|inframe_(deletion|insertion)|protein_altering_variant")  &
          is_rare_or_unreported(MAX_AF,rare_cutoff),
        1,0
      ))

    if(any(duplicated(vep_info_cols$ID))){
      stop("Problem with VEP annotations: Non-unique rows in output")
    }

    passing_ids <- data.table::fread(vep_in_file, header = FALSE, sep = "\t", colClasses = "character")[[6]]
    message(sum(passing_ids %in% vep_info_cols$ID), "/", length(passing_ids), " passing variants were annotated by VEP")
    message(paste0("Identified ", sum(vep_info_cols$priority_flag == 1), " priority variants"))

    # write the annotations into the empty VEP columns of every file the loop wrote
    vep_annot <- vep_info_cols %>% ungroup() %>% dplyr::rename(variant = ID)
    .fill_vep_cols(summary_file, vep_annot)
    .fill_vep_cols(summary_pass_file, vep_annot)
    .fill_vep_cols(counts_per_cell_type_file, vep_annot)
    .fill_vep_cols(genotype_file, vep_annot, na_string = "", n_skip = 1)
    if(run_cell_type_enrichment){
      for (ct in analysis_cell_types) {
        .fill_vep_cols(file.path(out_dir,paste0(file_prefix,"_",ct,".tsv")), vep_annot)
      }
    }

    # passing priority variants (previously written inside the loop)
    pass_rows <- data.table::fread(summary_pass_file, sep = "\t", header = TRUE, colClasses = "character",
                                   na.strings = NULL, quote = "", strip.white = FALSE, data.table = FALSE)
    pass_rows <- pass_rows[pass_rows$priority_flag == "1", names(pass_rows) != "priority_flag", drop = FALSE]
    data.table::fwrite(pass_rows, summary_pass_priority_file, sep = "\t", quote = FALSE, na = "NA")

    rm(vep_run, vep_info_cols, vep_annot, pass_rows)
    invisible(gc())
  }
  
  if(run_cell_type_enrichment){
    # sorting
    
    ## read in just sorting column
    ## sort
    ## get indexes
    ## rewrite file block-wise using index
    
    # padj
    ## read in p-value col
    ## adjust pvalues
    ## re-write block-wise using index
    
    message("Sorting and Adjusting P Values...")
    
    for (ct in analysis_cell_types) {
      # stores cell type-specific stats for all variants
      ct_file <- file.path(out_dir,paste0(file_prefix,"_",ct,".tsv"))
      # stores cell type-specific stats for passing variants
      ct_pass_file <- file.path(out_dir,paste0(file_prefix,"_",ct,"_pass_only.tsv"))
      # stores cell type-specific stats for passing priority variants
      ct_pass_priority_file <- file.path(out_dir,paste0(file_prefix,"_",ct,"_pass_only_priority.tsv"))
      
      
      full_results <- fread(ct_file)

      # priority_only: test and report passing priority variants only
      if(priority_only){
        full_results <- full_results[priority_flag %in% 1]
      }

      
      tested_idx <- with(full_results, filter == "." & is.finite(p))
      
      # adjust p-values only on passing variants
      full_results[, padj := NA_real_]
      if (any(tested_idx)) {
        full_results$padj[tested_idx] <- p.adjust(full_results$p[tested_idx], method = p_adjust_method)
      }
      
      
      add <- .make_filter(cbind(
        odds_ratio_cutoff = full_results$OR      < odds_ratio_cutoff,
        zscore_cutoff     = full_results$Z_score < zscore_cutoff,
        p_value           = full_results$padj    > max_p_value,
        cooc_support      = if (!is.null(max_cooc_p)) !is.na(full_results$cooc_p) & full_results$cooc_p > max_cooc_p
                            else rep(FALSE, nrow(full_results))
      ))
      f <- full_results$filter
      full_results$filter <- ifelse(add == ".", f,
                              ifelse(f == ".", add, paste0(f, ";", add)))
      
      
      if(celltype_sort_value=="p"){
        full_results <- full_results %>%
          arrange(padj,p)
      }else if(celltype_sort_value=="z"){
        full_results <- full_results %>%
          arrange((desc(Z_score)))
      }else if(celltype_sort_value=="cooc"){
        full_results <- full_results %>%
          arrange(cooc_p, desc(OR))
      }
      
      if(!priority_only){
        write_tsv(full_results, ct_file, append = FALSE, col_names = TRUE)
        
        write_tsv(full_results %>%
                    filter(filter==".") %>%
                    dplyr::select(-filter), ct_pass_file, append = FALSE, col_names = TRUE)
        
        if (isFALSE(skip_vep)) {
          write_tsv(full_results %>%
                      filter(filter==".",priority_flag==1) %>%
                      dplyr::select(-c(filter,priority_flag)), ct_pass_priority_file, append = FALSE, col_names = TRUE)
        }
        
        
        
        
        
      }else{
        write_tsv(full_results, ct_file, append = FALSE, col_names = TRUE)
        write_tsv(full_results %>%
                    filter(filter==".") %>%
                    dplyr::select(-filter), ct_pass_priority_file, append = FALSE, col_names = TRUE)
      }
    }
  }
}

.normalise_vcf <- function(vcf_df){
  out <- vcf_df %>%
    mutate(
      # calculate lengths
      ref_len = nchar(REF),
      alt_len = nchar(ALT),
      # detect INS/DEL/SNV based on ref/alt lengths
      simple_ins = (ref_len == 1) & (alt_len > 1) & (substr(REF,1,1)==substr(ALT,1,1)),
      complex_ins = (ref_len > 1) & (alt_len > ref_len) & (substr(REF,1,1)==substr(ALT,1,1)),
      #ins_len0 = (ref_len == 0) & (alt_len > ref_len),
      simple_del = ((alt_len == 0) & (ref_len > alt_len)) | (ALT=="-" & (ref_len > alt_len)),
      complex_del = (alt_len > 0) & (ref_len > alt_len) & (substr(REF,1,1)==substr(ALT,1,1)),
      simple_snv = (ref_len == 1) & (alt_len == 1) & (REF != ALT),
      complex_snv = (ref_len == alt_len) & (ref_len > 1) & (REF != ALT)) %>%
    mutate(
      complex_snv_pos = .complex_snv_pos(REF, ALT, complex_snv)
    ) %>%
    dplyr::as_tibble() %>%
    mutate(mnv = complex_snv & is.na(complex_snv_pos),
           complex_snv = complex_snv & !is.na(complex_snv_pos)
    ) %>%
    mutate(
      # skip ambiguous calls
      skip = str_detect(ALT, stringr::fixed("*")) | str_detect(REF, stringr::fixed("*")) | (!simple_ins&!simple_del&!complex_ins&!complex_del&!simple_snv&!complex_snv&!mnv),
      # handling positions
      START = case_when(
        simple_ins ~ POS + 1,
        complex_ins ~ POS + 1,
        #ins_len0 ~ POS,
        simple_del ~ POS,
        complex_del ~ POS + 1,
        simple_snv ~ POS,
        complex_snv ~ POS + complex_snv_pos - 1,
        mnv ~ POS
      ),
      END = case_when(
        simple_ins ~ POS,
        complex_ins ~ POS,
        #ins_len0 ~ POS-1,
        simple_del ~ POS + ref_len - 1,
        complex_del ~ POS + ref_len - alt_len,
        simple_snv ~ POS,
        complex_snv ~ POS + complex_snv_pos - 1,
        mnv ~ POS + ref_len - 1
      ),
      var_type = case_when(
        simple_ins ~ "INS", 
        complex_ins ~ "INS", 
        #ins_len0 ~ "INS", 
        simple_del ~ "DEL", 
        complex_del ~ "DEL", 
        simple_snv ~ "SNV", 
        complex_snv ~ "SNV",
        mnv ~ "MNV",
        TRUE ~ NA_character_),
      new_REF = case_when(
        simple_ins ~ "-", 
        complex_ins ~ "-", 
        #ins_len0 ~ "-", 
        simple_del ~ REF, 
        complex_del ~ substr(REF,2,1+ref_len-alt_len),
        simple_snv ~ REF, 
        complex_snv ~ substr(REF,complex_snv_pos,complex_snv_pos),
        mnv ~ REF
      ),
      new_ALT = case_when(
        simple_ins ~ substr(ALT,2,alt_len), 
        complex_ins ~ substr(ALT,2,alt_len-ref_len+1), 
        #ins_len0 ~ ALT, 
        simple_del ~ "-", 
        complex_del ~ "-", # ATGAG ATG
        simple_snv ~ ALT, 
        complex_snv ~ substr(ALT,complex_snv_pos,complex_snv_pos),
        mnv ~ ALT
      ),
      var_key = case_when(
        simple_ins ~ paste0(CHROM,":",START,"+",new_ALT),
        complex_ins ~ paste0(CHROM,":",START,"+",new_ALT), 
        #ins_len0 ~ paste0(CHROM,":",START,"+",new_ALT), 
        simple_del ~ paste0(CHROM,":",START,"-",END,"-",new_REF), 
        complex_del ~ paste0(CHROM,":",START,"-",END,"-",new_REF), 
        simple_snv ~ paste0(CHROM,":",START,new_REF,">",new_ALT), 
        complex_snv ~ paste0(CHROM,":",START,new_REF,">",new_ALT),
        mnv ~ paste0(CHROM,":",START,new_REF,">",new_ALT)
      )
    ) %>%
    mutate(filter=skip | is.na(var_type)) %>%
    dplyr::select(chromosome=CHROM, start_position = START, end_position=END, reference_allele=new_REF, alternate_allele=new_ALT, variant_type=var_type, variant=var_key,filter)
}


.pull_expanded_ngt_from_vcf <- function(v) {
  GT <- geno(v)$GT # get unexpanded genotypes (includes multiallelic sites)
  nvar <- nrow(GT) # number of variants
  nsamp <- ncol(GT) # number of samples
  
  # count number of alt alleles per row/variant
  nALT <- elementNROWS(alt(v))
  
  # get per-sample alleles A and B
  pat <- "^([0-9]+|\\.)[\\|/]([0-9]+|\\.)$"
  # sample allele A
  A <- suppressWarnings(as.integer(sub(pat, "\\1", GT, perl = TRUE))) %>%
    matrix(nrow=nvar, ncol = nsamp)
  # sample allele B
  B <- suppressWarnings(as.integer(sub(pat, "\\2", GT, perl = TRUE))) %>%
    matrix(nrow=nvar, ncol = nsamp)
  A[GT == "./."] <- NA_integer_
  B[GT == "./."] <- NA_integer_
  
  # Build expanded row mapping (same order as expand(v))
  exp_rows <- sum(nALT) # nrows after expansion
  row_map  <- rep.int(seq_len(nvar), nALT)  # which original row for each expanded row
  alt_idx  <- unlist(lapply(nALT, seq_len), use.names = FALSE)
  
  # expand A and B rows
  Aexp <- A[row_map, , drop = FALSE]
  Bexp <- B[row_map, , drop = FALSE]
  # altid matrix map
  J <- matrix(alt_idx, nrow = exp_rows, ncol = nsamp)
  
  # Build NGT (0,1,2,3) per expanded row x sample
  NGT <- matrix(3L, nrow = exp_rows, ncol = nsamp,
                dimnames = list(NULL, colnames(GT)))
  
  # missing stays 3
  not_na <- !(is.na(Aexp) | is.na(Bexp))
  
  # 0/0 → 0
  NGT[not_na & (Aexp == 0L) & (Bexp == 0L)] <- 0L
  # j/j → 2
  NGT[not_na & (Aexp == J)  & (Bexp == J)]  <- 2L
  # 0/j or j/0 → 1
  het <- not_na & (
    ((Aexp == 0L) & (Bexp == J)) |
      ((Aexp == J)  & (Bexp == 0L))
  )
  NGT[het] <- 1L
  
  NGT
}

.calculate_variant_level_stats <- function(NGT,AF,GQ,DP) {
  nc=ncol(NGT)
  # Totals across all cells (for "other" computations)
  alt_cnt_total  <- rowSums(NGT == 1L | NGT == 2L, na.rm = TRUE) # how many cells carry alt allele
  ref_cnt_total  <- rowSums(NGT == 0L, na.rm = TRUE) # how many cells carry ref allele
  data_cnt_total <- alt_cnt_total + ref_cnt_total # how many cells have non-missing data
  alt_proportion_total  <- ifelse(data_cnt_total > 0, alt_cnt_total / data_cnt_total, 0)
  data_proportion_total <- ifelse(data_cnt_total > 0, data_cnt_total/ nc, 0)
  
  ## Calculate means/medians only when AF>0 (mask)
  mask <- is.finite(AF) & (AF > 0)
  denom <- pmax(1L, rowSums(mask))
  mean_AF_total <- rowSums(AF * mask, na.rm = TRUE) / denom
  mean_AF_total[rowSums(mask) == 0] <- NA_real_
  
  AF_na <- AF; AF_na[!mask] <- NA_real_
  median_AF_total <- matrixStats::rowMedians(AF_na, na.rm = TRUE)
  
  ## Same as above but only for cells called as mutatns in NGT
  mask <- NGT %in% c(1L,2L)
  dim(mask) <- dim(NGT)
  denom <- pmax(1L, rowSums(mask))
  mean_AF_mutants <- rowSums(AF * mask, na.rm = TRUE) / denom
  mean_AF_mutants[rowSums(mask) == 0] <- NA_real_
  AF_na <- AF; AF_na[!mask] <- NA_real_
  median_AF_mutants <- matrixStats::rowMedians(AF_na, na.rm = TRUE)
  
  
  # GQ mean where AF>0
  GQuse <- GQ
  GQuse[!mask] <- NA_real_
  mean_GQ_total <- rowMeans(GQuse, na.rm = TRUE)
  
  # DP mean 
  DP[!mask] <- NA_real_
  mean_DP_total <- rowMeans(DP, na.rm = TRUE)
  
  out <- data.frame(
    mean_DP_total,
    mean_GQ_total,
    mean_AF_total,
    median_AF_total,
    alt_cnt_total,
    ref_cnt_total,
    data_cnt_total,
    alt_proportion_total,
    data_proportion_total,
    mean_AF_mutants,
    median_AF_mutants
  )
  return(out)
}

.calculate_cell_type_stats <- function(NGT,cell_type,cell_annotations,variant_stats){
  
  barcode_idx <- which(cell_annotations == cell_type)
  
  alt_cnt_total <- variant_stats$alt_cnt_total
  ref_cnt_total <- variant_stats$ref_cnt_total
  data_cnt_total <- variant_stats$data_cnt_total
  
  # Counts within this celltype
  NGT_ct <- NGT[, barcode_idx, drop = FALSE] # get NGT for celltype
  alt_cnt_ct  <- rowSums(NGT_ct == 1L | NGT_ct == 2L, na.rm = TRUE)
  het_cnt_ct  <- rowSums(NGT_ct == 1L, na.rm = TRUE)
  hom_cnt_ct <- rowSums(NGT_ct == 2L, na.rm = TRUE)
  ref_cnt_ct  <- rowSums(NGT_ct == 0L, na.rm = TRUE)
  no_data_cnt_ct  <- rowSums(NGT_ct == 3L, na.rm = TRUE)
  data_cnt_ct <- alt_cnt_ct + ref_cnt_ct
  alt_proportion_ct  <- ifelse(data_cnt_ct > 0, alt_cnt_ct / data_cnt_ct, 0)
  alt_proportion_in_ct <- ifelse(alt_cnt_total > 0, alt_cnt_ct / alt_cnt_total, 0)
  data_proportion_ct <- ifelse(data_cnt_ct > 0, data_cnt_ct/ length(barcode_idx), 0)
  
  # Counts for "all other" cell types
  alt_cnt_other  <- alt_cnt_total - alt_cnt_ct
  ref_cnt_other <- ref_cnt_total - ref_cnt_ct
  data_cnt_other <- data_cnt_total - data_cnt_ct
  alt_proportion_other <- ifelse(data_cnt_other > 0, alt_cnt_other / data_cnt_other, 0)
  
  a <- alt_cnt_ct
  b <- ref_cnt_ct
  c <- alt_cnt_other
  d <- ref_cnt_other
  
  # Haldane-Anscombe continuity correction only if any cell is zero
  needs_cc <- (a == 0L | b == 0L | c == 0L | d == 0L)
  
  a2 <- ifelse(needs_cc, a + 0.5, a)
  b2 <- ifelse(needs_cc, b + 0.5, b)
  c2 <- ifelse(needs_cc, c + 0.5, c)
  d2 <- ifelse(needs_cc, d + 0.5, d)
  
  # Wald (using corrected counts when needed)
  OR    <- (a2 / b2) / (c2 / d2)
  logOR <- log(OR)
  SE    <- sqrt(1/a2 + 1/b2 + 1/c2 + 1/d2)
  Z_score     <- logOR / SE
  p_wald <- 2 * pnorm(-abs(Z_score))
  ci_low  <- exp(logOR - 1.96 * SE)
  ci_high <- exp(logOR + 1.96 * SE)
  
  # Fisher decision
  has_data    <- (a + b + c + d) > 0L
  need_fisher <- has_data & (pmin(a, b, c, d) < 5L)
  
  p_fisher <- rep(NA_real_, length(a))
  if (any(need_fisher)) {
    idx <- which(need_fisher)
    p_fisher[idx] <- vapply(idx, function(i) {
      m <- matrix(c(a[i], b[i], c[i], d[i]), nrow = 2, byrow = TRUE)
      fisher.test(m)$p.value
    }, numeric(1))
  }
  
  p_final <- ifelse(!is.na(p_fisher), p_fisher, p_wald)
  
  out <- data.frame(
    alt_cnt_ct,
    het_cnt_ct,
    hom_cnt_ct,
    ref_cnt_ct,
    data_cnt_ct,
    alt_proportion_ct,
    alt_proportion_in_ct,
    alt_cnt_other,
    ref_cnt_other,
    data_cnt_other,
    alt_proportion_other,
    OR,
    ci_low = ci_low,
    ci_high = ci_high,
    Z_score,
    p_wald = p_wald,
    p_fisher = p_fisher,
    p = p_final
    
  )
}

.make_filter <- function(mat) {
  out <- rep("", nrow(mat))
  for (j in seq_len(ncol(mat))) {
    hit <- which(mat[, j])
    if (length(hit)) out[hit] <- paste0(out[hit], ifelse(out[hit] == "", "", ";"),
                                        colnames(mat)[j])
  }
  out[out == ""] <- "."
  out
}

# Empty VEP columns, in the same order as the old left_join added them
.add_empty_vep_cols <- function(df) {
  df %>% mutate(SYMBOL = NA_character_, Gene = NA_character_, Feature = NA_character_,
                HGVSp = NA_character_, MAX_AF = NA_character_, Consequence = NA_character_,
                Existing_variation = NA_character_, BIOTYPE = NA_character_,
                SIFT = NA_character_, PolyPhen = NA_character_, priority_flag = NA_real_)
}

# Rewrite an output file in chunks, filling the VEP columns for variants found
# in `vep_annot`. Lines for variants without an annotation are copied unchanged;
# only annotated lines are split into columns, filled and pasted back together.
# `na_string` should match how the file was first written ("NA" for write_tsv,
# "" for fwrite). `n_skip` = lines before the column header (1 for the #FORMAT
# line of the genotype file).
.fill_vep_cols <- function(file, vep_annot, na_string = "NA", n_skip = 0, chunk_lines = 200000, label = "VEP") {
  if (!file.exists(file)) return(invisible(NULL))
  message("Adding ", label, " annotations to ", file)
  tmp <- paste0(file, ".tmp")
  con <- file(file, "r")
  top_lines <- readLines(con, n = n_skip + 1)
  header <- top_lines[n_skip + 1]
  writeLines(top_lines, tmp)
  out <- file(tmp, "a")
  col_names <- strsplit(header, "\t", fixed = TRUE)[[1]]
  fill_cols <- intersect(setdiff(names(vep_annot), "variant"), col_names)
  # regex that pulls out the `variant` field of each line
  variant_regex <- paste0("^(?:[^\t]*\t){", match("variant", col_names) - 1, "}([^\t]*).*$")
  repeat {
    lines <- readLines(con, n = chunk_lines)
    if (length(lines) == 0) break
    hit <- which(sub(variant_regex, "\\1", lines, perl = TRUE) %in% vep_annot$variant)
    if (length(hit)) {
      chunk <- data.table::fread(text = c(header, lines[hit]), sep = "\t", header = TRUE, colClasses = "character",
                                 na.strings = NULL, quote = "", strip.white = FALSE, data.table = FALSE)
      m <- match(chunk$variant, vep_annot$variant)
      for (col in fill_cols) {
        x <- as.character(vep_annot[[col]][m])
        x[is.na(x)] <- na_string
        chunk[[col]] <- x
      }
      lines[hit] <- do.call(paste, c(chunk, sep = "\t"))
    }
    writeLines(lines, out)
  }
  close(out)
  close(con)
  file.rename(tmp, file)
}

.complex_snv_pos <- function(REF, ALT, complex_snv) {
  out <- rep(NA_real_, length(REF))
  idx <- which(complex_snv)
  out[idx] <- vapply(idx, function(i) {
    suppressWarnings(as.numeric(paste(which(str_split(ALT[i], pattern = "")[[1]] != str_split(REF[i], pattern = "")[[1]]), collapse = ",")))
  }, numeric(1))
  out
}

# Empty co-occurrence columns, filled after the loop for passing variants
.add_empty_cooc_cols <- function(df) {
  df %>% mutate(cooc_n_carriers = NA_integer_, cooc_testable = NA_integer_,
                cooc_partners = NA_integer_, cooc_p = NA_real_)
}

# Per-cell calls for co-occurrence scoring (raw to save memory):
# 2 = confident alt, 1 = confident reference, 0 = uncertain
.cooc_calls <- function(AF, DP, min_af, min_altreads, min_dp_ref) {
  alt_reads <- round(AF * DP)
  calls <- matrix(as.raw(0L), nrow(AF), ncol(AF))
  is_ref <- is.finite(AF) & AF == 0 & is.finite(DP) & DP >= min_dp_ref
  is_alt <- is.finite(alt_reads) & AF >= min_af & alt_reads >= min_altreads
  calls[is_ref] <- as.raw(1L)
  calls[is_alt] <- as.raw(2L)
  calls
}

# Count partners for each query variant against a reference set. Overlap is
# computed over jointly covered cells; nearby variants (and self) are excluded.
.count_partners_cov <- function(carr, cov, query, ref, chrom, pos, jaccard_thresh,
                                overlap_thresh, min_distance, block = NULL) {
  if (is.null(block)) block <- max(50L, floor(5e6 / length(ref)))
  carr_r <- carr[ref, , drop = FALSE] * 1
  cov_r  <- cov[ref, , drop = FALSE] * 1
  out <- integer(length(query))
  for (s in seq(1L, length(query), by = block)) {
    pos_q <- s:min(s + block - 1L, length(query))
    qi <- query[pos_q]
    cq <- carr[qi, , drop = FALSE] * 1
    vq <- cov[qi, , drop = FALSE] * 1
    inter  <- tcrossprod(cq, carr_r)          # alt in both
    a_in_b <- tcrossprod(cq, cov_r)           # query carriers covered in ref variant
    b_in_a <- tcrossprod(vq, carr_r)          # ref carriers covered in query variant
    jac <- inter / (a_in_b + b_in_a - inter)
    ovl <- inter / pmin(a_in_b, b_in_a)
    partner <- inter >= 2 & (jac >= jaccard_thresh | ovl >= overlap_thresh)
    partner[is.na(partner)] <- FALSE
    near <- outer(chrom[qi], chrom[ref], "==") & abs(outer(pos[qi], pos[ref], "-")) < min_distance
    partner[near | outer(qi, ref, "==")] <- FALSE
    out[pos_q] <- rowSums(partner)
  }
  out
}

# Null: keep each variant's carrier count and covered cells, but place carriers
# with probability proportional to each cell's alt rate
.shuffle_carriers_weighted <- function(carr, cov) {
  w <- colSums(carr) / pmax(colSums(cov), 1) + 1e-6
  out <- carr
  out[] <- FALSE
  for (i in seq_len(nrow(carr))) {
    k <- sum(carr[i, ])
    if (k > 0) {
      ok <- which(cov[i, ])
      out[i, ok[sample.int(length(ok), k, prob = w[ok])]] <- TRUE
    }
  }
  out
}

# Partner counts and empirical p-values against the weighted null,
# compared within carrier-count bins
.cooccurrence_scores <- function(calls, info, jaccard_thresh, overlap_thresh, min_distance,
                                 max_reference, n_bins = 20, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  carr <- calls == as.raw(2L)
  cov  <- calls != as.raw(0L)
  rm(calls)
  n <- nrow(carr)
  chrom <- as.integer(factor(info$chromosome))
  pos <- info$start_position
  ref    <- if (n > max_reference) sort(sample.int(n, max_reference)) else seq_len(n)
  null_q <- if (n > max_reference) sort(sample.int(n, max_reference)) else seq_len(n)

  obs <- .count_partners_cov(carr, cov, seq_len(n), ref, chrom, pos,
                             jaccard_thresh, overlap_thresh, min_distance)
  shuf <- .shuffle_carriers_weighted(carr, cov)
  nul <- .count_partners_cov(shuf, cov, null_q, ref, chrom, pos,
                             jaccard_thresh, overlap_thresh, min_distance)

  k_bin <- pmin(rowSums(carr), n_bins)
  cooc_p <- rep(NA_real_, n)
  for (b in unique(k_bin)) {
    null_b <- sort(nul[k_bin[null_q] == b])
    if (length(null_b) == 0) next
    i <- which(k_bin == b)
    n_ge <- length(null_b) - findInterval(obs[i] - 0.5, null_b)
    cooc_p[i] <- (1 + n_ge) / (1 + length(null_b))
  }
  data.frame(variant = info$variant, cooc_partners = obs, cooc_p = signif(cooc_p, 4))
}
