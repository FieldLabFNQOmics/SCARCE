#' prepare_annotation_data.R
#'
#' This function loads annotation data from a tapestri .h5 file and generates a seurat object ready for annotation
#'
#' @param input_file Path to input .dna+protein.h5 file
#' @param sample_ID Sample ID, or NULL to derive from .h5 filname (default)
#' @param cluster_resolutions Numeric vector of resolutions for Seurat::FindClusters (default: c(0.2, 0.35, 0.5, 0.8, 1, 1.2))
#' @param save_directory Output directory to save seurat object as .RDS, or NULL to prevent saving (default)
#' @param qc_plot_directory Directory to save QC plots to
#' @param variance_quantile Threshold for selecting the top quantile of ADT markers based on their variance
#' @param neighbor_metric Distance metric to be used for cell clustering by Seurat::RunUMAP and Seurat::FindNeighbors (eg. euclidean (default), cosine, correlation)
#' @param k_param k.param for Seurat::FindNeighbors, also used for n.neighbors in Seurat::RunUMAP
#' @return Seurat object
#' @export

prepare_annotation_data <- function(
  input_file=NULL,
  data_type="Tapestri_ADT",
  sample_ID=NULL,
  cluster_resolutions=c(0.2, 0.35, 0.5, 0.8, 1, 1.2),
  save_directory=NULL,
  qc_plot_directory=NULL,
  variance_quantile=0.25,
  neighbor_metric="euclidean",
  k_param=50,
  filter=TRUE,
  mito_pattern = "^MT-",
  nfeatures = 2000,
  min_features = 200,
  max_mito_pct = 20,
  convert_ENSEMBL_to_HGNC = FALSE
){
  if(is.null(sample_ID)){
    sample_ID = gsub(".dna[+]protein.h5", "", sort(basename(input_file)))
  }
  message(paste0("Preparing data for ",sample_ID," from: ", input_file))
  
  if (data_type == "Tapestri_ADT") {

    #############
    # Ingestion #
    #############
    protein_count_matrix <- h5read(input_file, "/assays/protein_read_counts/layers/read_counts") 
    rownames(protein_count_matrix) <-  as.character(h5read(input_file, "/assays/protein_read_counts/ca/id"))
    colnames(protein_count_matrix) <- as.character(h5read(input_file, "/assays/protein_read_counts/ra/barcode"))
    protein_count_matrix <- Matrix(protein_count_matrix, sparse=T)
    
    obj <- CreateSeuratObject(counts = protein_count_matrix, assay = "ADT", project = "Tapestri_Protein")
    DefaultAssay(obj) <- "ADT"
    
    #################
    # Normalisation #
    #################
    
    # CLR-normalize across cells
    obj <- NormalizeData(obj, normalization.method = "CLR", margin = 2)
    
    # Variance filtering
    norm_mat <- GetAssayData(obj, assay = "ADT", layer = "data")
    vars <- matrixStats::rowVars(as.matrix(norm_mat))
    thresh <- quantile(vars, variance_quantile)
    
    keep <- vars > thresh & Matrix::rowSums(is.finite(norm_mat)) == ncol(norm_mat)
    if (any(!keep)) {
      message(paste0("Warning: The following features are being removed due to having low variance or non-finite values\n",paste(names(which(keep==F)),collapse = "\n")))
      obj <- subset(obj, features = names(keep)[keep])
    }
    
    # scaled data
    obj <- ScaleData(obj, assay = "ADT", layer = "data")  # writes layer "scale.data"
    
    ###########
    # Run PCA #
    ###########
    # calculate reasonable number of PCs
    npcs <- min(20L, nrow(obj), ncol(obj) - 1L)
    
    obj <- RunPCA(
      object  = obj,
      assay   = "ADT",
      layer   = "scale.data",
      features = names(keep)[keep],
      npcs    = npcs,
      approx  = FALSE,
      verbose = FALSE
    )
    
    ##############
    # Elbow Plot #
    ##############
    # Take only PCs which have a change in variation is > 0.1% from the next -- OTHERWISE, can just use 7
    elb_plot <- ElbowPlot(obj, ndims=30)$data

    pct <- elb_plot$stdev
    max_pca_dim <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
    
    if(is.na(max_pca_dim)) max_pca_dim <- 7
    
    if(!is.null(qc_plot_directory)){
      elb_plot <- ggplot() + 
        geom_point(data = elb_plot %>% filter(dims <= max_pca_dim), mapping = aes(dims, stdev), size = 5, color = 'green') +
        geom_point(data = elb_plot %>% filter(dims > max_pca_dim), mapping = aes(dims, stdev), size = 5, color = 'red') +
        ggtitle(sample_ID) +
        theme_bw()
      ggsave(file.path(qc_plot_directory,paste0(sample_ID,"_pca_elbow_plot.png")),plot = elb_plot)
    }
    
    ########
    # UMAP #
    ########
    set.seed(123)
    obj <- FindNeighbors(obj, reduction = "pca", 
                        dims = 1:max_pca_dim,
                        annoy.metric = neighbor_metric,
                        k.param=k_param)
    nn_name <- grep("_nn$", names(obj@graphs), value = TRUE)[1]
    
    obj <- RunUMAP(obj, 
                  reduction = "pca",
                  dims = 1:max_pca_dim,
                  metric = neighbor_metric,
                  n.neighbors=k_param,
                  min.dist = 0.3)
    
    
    for (res in cluster_resolutions) {
      obj <- FindClusters(obj, resolution = res,random.seed=123)
      if(!is.null(qc_plot_directory)){
        res_col <- paste0("ADT_snn_res.", res)
        DimPlot(obj,
                reduction = "umap",
                group.by  = res_col,
                label     = TRUE,
                repel     = TRUE, label.box = T) +
          ggtitle(paste("Seurat Clusters at Resolution: ",res))
        ggsave(file.path(qc_plot_directory,paste0(sample_ID,"_seurat_clusters_res_",res,"_umap.png")))
        DimPlot(obj,
                reduction = "pca",
                group.by  = res_col,
                label     = TRUE,
                repel     = TRUE, label.box = T) +
          ggtitle(paste("Seurat Clusters at Resolution: ",res))
        ggsave(file.path(qc_plot_directory,paste0(sample_ID,"_seurat_clusters_res_",res,"_pca.png")))
      }
    }
    
    obj@meta.data$sample_ID <- sample_ID
    obj@meta.data <- cbind(obj@meta.data, as.data.frame(obj@reductions$pca@cell.embeddings)[rownames(obj@meta.data),c(1:npcs)])
    obj@meta.data <- cbind(obj@meta.data, as.data.frame(obj@reductions$umap@cell.embeddings)[rownames(obj@meta.data),c(1:2)])
    
    if(!is.null(save_directory)){
      saveRDS(obj,file.path(save_directory,paste0(sample_ID,"_seurat_obj.rds")))
    }
  } else if (data_type == "SEACells_h5ad") {

    #############
    # Ingestion #
    #############

    sce <- zellkonverter::readH5AD(input_file, use_hdf5 = FALSE, reader = "R")
    an  <- SummarizedExperiment::assayNames(sce)
    use <- if ("X" %in% an) "X" else an[1]
    message("Using assay '", use, "' from h5ad (available: ",
            paste(an, collapse = ", "), ")")
    mat <- methods::as(SummarizedExperiment::assay(sce, use), "CsparseMatrix")
    
    if (convert_ENSEMBL_to_HGNC) {
      ids <- rownames(mat)
      is_ens <- grepl("^ENSG[0-9]+(\\.[0-9]+)?$", ids)   # only the still-unresolved ones
      message(sum(is_ens), "/", length(ids), " rownames are Ensembl IDs; the rest are already symbols.")

      if (any(is_ens)) {
        if (!requireNamespace("EnsDb.Hsapiens.v86", quietly = TRUE))
          stop("convert_ENSEMBL_to_HGNC = TRUE requires 'EnsDb.Hsapiens.v86'.")
        sym <- AnnotationDbi::mapIds(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86,
                                    keys = sub("\\.\\d+$", "", ids[is_ens]),
                                    column = "SYMBOL", keytype = "GENEID",
                                    multiVals = "first")
        ids[is_ens] <- ifelse(is.na(sym), ids[is_ens], sym)  # keep original ID where no symbol found
        message("  ...mapped ", sum(!is.na(sym)), " of those to symbols.")
      }
      rownames(mat) <- make.unique(ids)
    }

    obj <- Seurat::CreateSeuratObject(
      counts = mat, assay = "RNA", project = "SEACells_metacell"
    )
    
    Seurat::DefaultAssay(obj) <- "RNA"

    ######
    # QC #
    ######
    if (filter) {
      obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(obj, pattern = mito_pattern)
      keep_cells <- obj$nFeature_RNA >= min_features & obj$percent.mt <= max_mito_pct
      message("QC: keeping ", sum(keep_cells), "/", ncol(obj), " cells.")
      obj <- subset(obj, cells = colnames(obj)[keep_cells])
    }

    ###################################
    # Normalisation + feature select. #
    ###################################
    
    obj <- Seurat::SCTransform(obj, assay = "RNA", verbose = FALSE,
                              variable.features.n = nfeatures)
    var_feats <- Seurat::VariableFeatures(obj)
    scale_layer <- "scale.data"

    ###########
    # Run PCA #
    ###########
    # Clamp npcs to what the data supports (matters for small metacell sets).
    npcs <- min(20, length(var_feats), ncol(obj) - 1L)
    obj <- Seurat::RunPCA(obj, assay = Seurat::DefaultAssay(obj),
                          features = var_feats, npcs = npcs, verbose = FALSE)

    ##############
    # Elbow Plot #
    ##############
    elb_plot <- Seurat::ElbowPlot(obj, ndims = npcs)$data
    pct <- elb_plot$stdev
    max_pca_dim <- sort(which((pct[1:(length(pct) - 1)] - pct[2:length(pct)]) > 0.1),
                        decreasing = TRUE)[1] + 1
    if (is.na(max_pca_dim)) max_pca_dim <- min(30L, npcs)
    max_pca_dim <- min(max_pca_dim, npcs)

    if (!is.null(qc_plot_directory)) {
      p <- ggplot2::ggplot() +
        ggplot2::geom_point(data = subset(elb_plot, dims <= max_pca_dim),
                            ggplot2::aes(dims, stdev), size = 5, colour = "green") +
        ggplot2::geom_point(data = subset(elb_plot, dims > max_pca_dim),
                            ggplot2::aes(dims, stdev), size = 5, colour = "red") +
        ggplot2::ggtitle(sample_ID) + ggplot2::theme_bw()
      ggplot2::ggsave(file.path(qc_plot_directory,
                                paste0(sample_ID, "_pca_elbow_plot.png")), plot = p)
    }

    #####################
    # Neighbours + UMAP #
    #####################
    # Clamp neighbour size below n cells (metacell sets can be tiny).
    k_use <- min(k_param, ncol(obj) - 1L)
    if (k_use < k_param)
      message("Reducing k.param/n.neighbors to ", k_use, " (only ", ncol(obj),
              " observations).")

    set.seed(123)
    obj <- Seurat::FindNeighbors(obj, reduction = "pca", dims = 1:max_pca_dim,
                                annoy.metric = neighbor_metric, k.param = k_use)
    obj <- Seurat::RunUMAP(obj, reduction = "pca", dims = 1:max_pca_dim,
                          metric = neighbor_metric, n.neighbors = k_use,
                          min.dist = 0.3)

    for (res in cluster_resolutions) {
      obj <- Seurat::FindClusters(obj, resolution = res, random.seed = 123)
      if (!is.null(qc_plot_directory)) {
        res_col <- grep(paste0("_snn_res\\.", res, "$"), colnames(obj@meta.data),
                        value = TRUE)[1]
        Seurat::DimPlot(obj, reduction = "umap", group.by = res_col,
                        label = TRUE, repel = TRUE, label.box = TRUE) +
          ggplot2::ggtitle(paste("Seurat Clusters at Resolution:", res))
        ggplot2::ggsave(file.path(qc_plot_directory,
                                  paste0(sample_ID, "_seurat_clusters_res_", res, "_umap.png")))
      }
    }

    ############
    # Metadata #
    ############
    obj@meta.data$sample_ID <- sample_ID
    pca_emb  <- as.data.frame(obj@reductions$pca@cell.embeddings)
    umap_emb <- as.data.frame(obj@reductions$umap@cell.embeddings)
    obj@meta.data <- cbind(obj@meta.data,
                          pca_emb[rownames(obj@meta.data), 1:npcs, drop = FALSE],
                          umap_emb[rownames(obj@meta.data), 1:2, drop = FALSE])

    if (!is.null(save_directory))
      saveRDS(obj, file.path(save_directory, paste0(sample_ID, "_seurat_obj.rds")))
  }
  return(obj)
}
