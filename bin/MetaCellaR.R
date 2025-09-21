# ====== GLOBAL CORE RESTRICTION (25 CORES) ======
# 1. Set R-wide parallel options
#options(mc.cores = 25)  # Controls mclapply/BiocParallel

# 2. Environment variables for system-wide control
#Sys.setenv(
#  OMP_NUM_THREADS = 25,      # OpenMP (uwot, Rcpp, etc.)
#  MKL_NUM_THREADS = 25,      # Intel MKL 
#  NUMEXPR_NUM_THREADS = 25,  # NumExpr-based packages
#  PBLAS_NUM_THREADS = 25,    # Parallel BLAS
#  VECLIB_MAXIMUM_THREADS = 25 # macOS Accelerate Framework
#)

# 3. Force data.table compliance
#if ("data.table" %in% installed.packages()) {
#  data.table::setDTthreads(25)
#}

#cat("===== ACTIVE CORE LIMITS =====\n")
#cat("mc.cores     :", getOption("mc.cores", 1), "\n")            # mclapply
#cat("OMP_THREADS  :", Sys.getenv("OMP_NUM_THREADS"), "\n")       # OpenMP
#cat("MKL_THREADS  :", Sys.getenv("MKL_NUM_THREADS"), "\n")       # Intel MKL
#cat("data.table   :", data.table::getDTthreads(), "threads\n")   # data.table
#if (exists("blas_get_num_procs", where = asNamespace("RhpcBLASctl"))) {
#  cat("BLAS         :", RhpcBLASctl::blas_get_num_procs(), "threads\n")
#}

### To summerize single cells based on k-medoids ###
### The clustering is done per cell type ###
### parameter k is defined based on the number of potential groups that have at least 30 cells inside (k <- ncol(CT_cluster) / 30)

#set.seed(0)
library(cluster)
library(knn.covertree)
library(irlba)
library(data.table)
library(ggplot2)
library(ggdendro)
library(cowplot)
library(parallel)
library(DESeq2)

########################
summary_method <- "kmed_means" #kmed
iter_flag <- F
csv_flag <- F
merge_flag <- T
normalization_flag <- T
########################
####### Begin functions #######
################################
get_outliers <- function(clara_out, outlier_cell_num= 10){
	metacell_counts <- table(clara_out$clustering)
	return(which(metacell_counts < outlier_cell_num))
}
metacell_quality <- function(clara_out, expected_cell_num= 30, slack_ratio= .15){
	cell_num_slack <- expected_cell_num - floor(expected_cell_num * slack_ratio)
	metacell_counts <- table(clara_out$clustering)
	qualities <- metacell_counts > cell_num_slack
	return(qualities)
}
my_knn <- function(train, test, k= 1, threshold= 3 * expected_cells){
	# TODO: Need to assign the ATAC cells prioritizing based on their closest distance to the metacells. Whatever that is leftover (MC_tracking reaching threshold for all MCs) needs to be assigned to an NA cluster
	dist_mat <- matrix(NA, nrow= nrow(test), ncol= nrow(train)) ## Matrix holding the distances of any ATAC point to any RNA metacell
	rownames(dist_mat) <- rownames(test)
	MC_tracking <-seq(ncol(dist_mat)) * 0
	for(i in seq(nrow(test))){
		## Compute the distance between each test point and all the training points
		rep_test <- matrix(rep(test[i, ], nrow(train)), byrow= T, nrow= nrow(train), ncol= ncol(train))
		dist <- sqrt(rowSums((rep_test - train)^2)) ## Euclidean distance
		dist_mat[i, ] <- dist		
	}

	## Take the row-wise min of dist_mat
	all_mins <- apply(dist_mat, 1, FUN= min)
	## sort dist_mat according to the min values obtained above
	dist_mat <- dist_mat[order(all_mins, decreasing= F), ]
	### Assign the labels while keeping track of the label satuaration for any specific metacell
	i <- 1;
	labels <- vector("character", length= nrow(test))
	names(labels) <- rownames(dist_mat)
	while(i <= nrow(dist_mat)){
		hit <- which.min(dist_mat[i, ]);
		if(MC_tracking[hit] < threshold){
			MC_tracking[hit] <- MC_tracking[hit] + 1;
			labels[i] <- rownames(train)[hit];
			i <- i + 1
		}else{
			dist_mat[i, hit] <- Inf;
		}
	}
	if(F){
  	load("mojitoo_counts_metacell_30/kmed_means_clustered.RData")
	  train <- RNA_metacell_umap
	  test <- ATAC_umap;
	  threshold <- 90
  	df <- data.frame(umap1= test[1:15, 1], umap2= test[1:15, 2], assay= "ATAC", label= labels[1:15])
  	df <- rbind(df, data.frame(umap1= train[c(219, 301, 462, 527), 1], umap2= train[c(219, 301, 462, 527), 2], assay= c("RNA", "RNA", "RNA", "RNA"), label= rownames(train)[c(219, 301, 462, 527)]))
	}
	return(labels)

}
#####
cluster_means <- function(clusters){
	temp_cl <- NULL
	for(clst in unique(clusters$clustering)){
		if("numeric" %in% class(clusters$data[clusters$cluster == clst, ])){
			temp_cl <- rbind(temp_cl, clusters$data[clusters$cluster == clst, ])
		}
		else{
			temp_cl <- rbind(temp_cl, apply(clusters$data[clusters$cluster == clst, ], 2, FUN= mean))
		}
	}
	return(temp_cl)
}
#####

invoke_clara <- function(CT_cluster, original_CT_cluster, iter_flag, clusters, RNA_metacell_umap, ct, k){
	class(CT_cluster)
	if(iter_flag){
		for(iter in seq(10)){
			clara_res <- cluster::clara(t(as.matrix(CT_cluster)), k, metric = "euclidean", stand = FALSE, samples= 30, pamLike = FALSE)
			dim(clara_res$medoids)
			dim(CT_cluster)

			clusters[[ct]] <- clara_res
			all_mediods[[ct]] <- rbind(all_mediods[[ct]], clara_res$medoids)
		}
	}
	else{
		clusters[[ct]] <- cluster::clara(t(as.matrix(CT_cluster)), k, metric = "euclidean", stand = FALSE, samples= 30, pamLike = FALSE)
	}
	if(length(assay_slot)){
		RNA_metacell_umap_ct <- NULL
		for(i in unique(clusters[[ct]]$clustering)){
			data_subset <- RNA_umap[which(clusters[[ct]]$clustering == i), ];
			if(is.null(dim(data_subset))){
				RNA_metacell_umap_ct <- rbind(RNA_metacell_umap_ct, data_subset);
			}else{
				RNA_metacell_umap_ct <- rbind(RNA_metacell_umap_ct, colMeans(data_subset))
			}
		}
		rownames(RNA_metacell_umap_ct) <- paste(ct, seq(nrow(RNA_metacell_umap_ct)), sep= "_")
		print(paste("Done clustering", ct))
		RNA_metacell_umap <- rbind(RNA_metacell_umap, RNA_metacell_umap_ct)
	}
	#return(list(clusters= clusters, RNA_metacell_umap= RNA_metacell_umap))
	return(clusters= clusters)
}

#####

invoke_clara_simplified <- function(CT_cluster, k){
	class(CT_cluster)
	cluster_data <- t(as.matrix(CT_cluster))
	clusters <- cluster::clara(t(as.matrix(CT_cluster)), k, metric = "euclidean", stand = FALSE, samples= 30, pamLike = FALSE)
	return(list(clusters= clusters))
}
######
merge_small_mc <- function(clustering, thresh= 30){
	flag <- T
	while(flag){
		metacell_counts <- table(clustering)
		mc_table_sorted <- sort(metacell_counts)
		low_mc_idx <- which(metacell_counts < thresh)
		if(!length(low_mc_idx)){
			return(clustering)
		}
		idx <- order(metacell_counts[low_mc_idx], decreasing= F)
		new_mc <- NULL;
		mc_sum <- metacell_counts[low_mc_idx[idx][1]]
		new_mc <- names(metacell_counts[low_mc_idx[idx][1]])
		max_cluster_id <- max(clustering) + 1
		pointer <- 2
		while(mc_sum < thresh){
			if(length(idx) == 1){
				all_idx <- order(metacell_counts, decreasing= F)
				new_mc <- c(new_mc, names(metacell_counts[all_idx][2]))
				clustering[which(clustering %in% as.numeric(new_mc))] <- max_cluster_id
				flag <- F
				break;
			}
			mc_sum <- mc_sum + metacell_counts[low_mc_idx[idx][pointer]]
			new_mc <- c(new_mc, names(metacell_counts[low_mc_idx[idx][pointer]]))
			pointer <- pointer + 1
			if(pointer > length(idx)){
				break;
			}
		}
		clustering[which(clustering %in% as.numeric(new_mc))] <- max_cluster_id
		if(length(which(table(clustering) < thresh)) == 0){
			flag <- F
		}
	}
	return(clustering)
}
########
####### End of functions #######
################################

#' @import cluster
#' @import knn.covertree
#' @import irlba
#' @import data.table
#' @import ggplot2
#' @import ggdendro
#' @import cowplot
#' @import parallel
#' @import DESeq2
#' @import pheatmap 
#' @import viridis
#' @param file path to the Seurat object or a csv file containing scRNA-seq data
#' @param RNA indicate the slot in the Seurat object where the RNA tags are stored (e.g., assays$RNA@counts)
#' @param ATAC indicate the slot in the Seurat object where the ATAC tags are stored (e.g., assays$peak@counts)
#' @param celltype the slot in the Seurat object where the cell type annotations are stored (e.g., meta.data$Celltypes)
#' @param expected_cells number of expected cells per metacell
#' @param k number of clusters for k-mediods (will be derived from the expected number of cells, if not specified)
#' @param umap TRUE or FALSE indicating if data should be projected to a UMAP space prior to clustering
#' @param output_file path to a directory where the output files should be stored
#' @param assay_slot indicate the slot in the Seurat object where the assay information of the cells are stored
#' @param umap_dim if `umap` is true, the user can accompany it with the `umap_dim` argument allowing to specify the number of dimentions the umap function should reduce the data to. In other words, `umap_dim` is an integer argument indicating the number of UMAP dimensions (default is 20)
#' @param conditions vector of characters containing the conditions that exist in the metadata slot of Seurat for which the metacells should be created separately. For example, if there are multiple replicates and each replicate needs its own metacells, this option can be used to achieve this.
#' @param gtf_file is used to pass the path to a gtf file corresponding to the genome version the RNA reads were aligned to. This is used to obtain the exonic length for the normalization of RNA reads.
#' @return The output will be stored in children nodes (plots, results, and debug) of the `output_file` directory
#' @export

metacellar.run <- function(summary_method = "kmed_means", csv_flag = F,
merge_flag = T, normalization_flag = T, reduction = "umap", file_name, gtf_file = NULL,
umap_flag = F, RNA_count_slot, ATAC_count_slot, celltype_info, assay_slot, output_file = getwd(),
expected_cells = 30, threshold = 3 * expected_cells, umap_dim = 20, k= NULL, metadata= c()){
	summary_method <- "kmed_means" #kmed
	if(!length(metadata) & !csv_flag){
		metadata <- gsub("meta.data$", "", celltype_info)
	}
#	library("Matrix")
	if(!length(file_name)){
		stop("Argument file_name is missing. Please provide the path to your input file using the -file option followed by the path to your file!")
	}else if(csv_flag){# CSV file
		if(!length(celltype_info)){
			stop("For the csv input file the celltype_info argument must be present! Please provide a two column file with headers, linking the cell names (first column) to cell types (second columns).")
		}
		csv_flag <- T
		csv_data <- read.csv(file_name, row.names= 1)
		csv_cells <- read.table(celltype_info, header= T)
		general_data <- csv_data
		print("Done!")
	}
	dir.create(output_file)
	dir.create(paste0(output_file, "/plots"))
	dir.create(paste0(output_file, "/results"))
	dir.create(paste0(output_file, "/debug"))
########################
	if(!csv_flag){
		library(Seurat)
		if(strsplit(file_name, "\\.")[[1]][2] == "rds"){
			print("Reading the Seurat object")
			Sub <- readRDS(file_name)
			Rds_name <- "Sub"
		}else{
			print("Loading the Seurat object")
			Rds_name <- load(file_name)
			Sub <- get(Rds_name) ##To make it generalizable for Suerat object stored with any name
		}
		## Check if the given metadata fields exist in the Seurat object
		metadata_values <- list()
		if(length(metadata)){ ## !!I might need to use a tryCatch statement to properly address this check!!
			for(md in metadata){
				metadata_values[[md]] <- unique(unlist(Sub[[md]]))
			}
		}
		print("Done loading the Seurat object")
		print("Done loading Rds!")
		celltypes <- eval(parse(text= paste0(Rds_name, "@", celltype_info)))
		original_celltypes <- celltypes
		## Identify NA cell annotation and remove them from the object
		na_idx <- is.na(celltypes)
		if(sum(na_idx)){
			## TODO:
			Sub <- subset(Sub, !is.na(celltypes))## I know that it doesn't work as I've already tested it on a Seurat obj... gotta find another way
		}
		print(paste("RNA_count_slot:", RNA_count_slot))
		RNAcounts <- eval(parse(text= paste0(Rds_name, "@", RNA_count_slot))) #Sub@assays$RNA@counts
		# rna_hits <- seq(nrow(Sub@meta.data)) ## in case the input is a Suerat object containing only the RNA assay. This index vector is then later used during the defining of the groups for K-med clustering
		
		gene_names <- rownames(RNAcounts)
#		RNAcounts <- as.matrix(RNAcounts)
		cell_types <- unique(as.character(celltypes))
		if(length(assay_slot)){
			ATACcounts <- eval(parse(text= paste0(Rds_name, "@", ATAC_count_slot))) #Sub@assays$RNA@counts
			assays <- eval(parse(text= paste0(Rds_name, "@", assay_slot)))

			rna_hits <- which(tolower(assays) == "scrna-seq" | tolower(assays) == "scrna" | tolower(assays) == "scrna" | tolower(assays) == "rna")
			ATACcounts <- ATACcounts[, -rna_hits]
			RNAcounts <- RNAcounts[, rna_hits]
			general_data <- RNAcounts
			RNA_celltypes <- celltypes[rna_hits]
			ATAC_celltypes <- celltypes[-rna_hits]
			RNA_umap <- Sub[[reduction]]@cell.embeddings[rna_hits, ]
			ATAC_umap <- Sub[[reduction]]@cell.embeddings[-rna_hits, ]
			print("ATAC_umap")
			print(head(ATAC_umap))
			celltypes <- RNA_celltypes
			cell_types <- unique(celltypes)
			umap_layout <- RNA_umap
		}
	}else{
		cell_types <- unique(csv_cells[, 2])
	}
	if(umap_flag && !length(assay_slot)){
		print("Computing UMAP")
		ptm <- proc.time()
		if(csv_flag){
			ptm_umap <- proc.time(); 
			umap_layout <- uwot::umap(t(as.matrix(csv_data)), pca= 30, pca_center= T, scale= T, n_components= umap_dim)
			rownames(umap_layout) <- colnames(csv_data)
			print(paste("UMAP computation time:", proc.time() - ptm_umap))
			celltypes <- csv_cells[, 2]
		}else{
			ptm_umap <- proc.time();
			umap_layout <- uwot::umap(t(as.matrix(RNAcounts)), pca= 30, pca_center= T, scale= T, n_components= umap_dim)
			print(paste("UMAP computation time:", proc.time() - ptm_umap))
			rownames(umap_layout) <- colnames(RNAcounts)
		}
		RNA_umap <- umap_layout
	}
	clusters <- list()
	all_mediods <- list()
	RNA_metacell_umap <- NULL
	cell2metacell_info <- NULL
	cluster_data <- list()
	print(dim(RNAcounts))

	mc_quality_info <- list()
	mc_outlier_info <- list()
	RNA_metacell_umap_plot <- list()
	RNA_barcodes_ct <- NULL
	mc_distr <- list()
	ks_info <- list()
	pass_num <- 1
	ctr <- 1
	library(pheatmap)
	subsetting_rownames <- NULL
	subsetting_string_all <- NULL
	metadata_comb_names <- NULL
	combs_table <- as.matrix(expand.grid(metadata_values))
	md_combo <- vector("character", length(nrow(combs_table)))
	for(i in seq(nrow(combs_table))){
		subsetting_string <- NULL;
		metadata_names <- NULL
		for(j in seq(ncol(combs_table))){
			subsetting_string <- c(subsetting_string, "Sub[[\"", paste0(colnames(combs_table)[j], "\"]]==", "\"", combs_table[i, j], "\"", ifelse(j == ncol(combs_table), "", "&")));
			metadata_names <- c(metadata_names, paste0(combs_table[i, j], "_"))
		}
		subsetting_string_all <- c(subsetting_string_all, paste(subsetting_string, collapse= ""))
		metadata_comb_names <- c(metadata_comb_names, paste(metadata_names, collapse= ""))
		md_combo_tmp <- subsetting_string_all[i]
		for(idx in seq(length(metadata))){
			md_combo_tmp <- gsub("\"", "", md_combo_tmp)
			md_combo_tmp <- gsub(paste0("Sub\\[\\[", metadata[idx], "\\]\\]=="), "", md_combo_tmp)
		}
		md_combo[i] <- gsub("Sub\\[\\[|\\]\\]", "", gsub("&|==", "_", md_combo_tmp))
	}
	pdf(paste0(output_file, "/plots/", "mc_dendrogram.pdf"))
	for(i in seq(length(subsetting_string_all))){
		hits <- rownames(Sub@meta.data[eval(parse(text= subsetting_string_all[i])), ])
		if(length(assay_slot)){
			hits <- intersect(hits, rownames(Sub@meta.data[tolower(assays) == "rna",]))
		}
		if(!length(hits)){
			next;
		}

		## extract the metadata combo names from subsetting_string_all
		
		if(!csv_flag){
			if(!length(assay_slot)){
				#CT_cluster <- RNAcounts[, rownames(eval(parse(text= paste0(Rds_name, "@meta.data")))[celltypes == ct, ])]
				CT_cluster <- RNAcounts[, hits]
			}
			CT_cluster <- RNAcounts[, hits] #subset the cells based on the current combination of metadata values
			#CT_cluster <- CT_cluster[, rownames(eval(parse(text= paste0(Rds_name, "[[\"", metadata[md], "\"]]"))) == md_vals)]
		}else{
			CT_cluster <- csv_data[, csv_cells[csv_cells[, 2] == ct, 1]] ## this whole csv stuff are falling apart with this new implemntation of allowing to have multiple metadata variables to base the clustering on
		}
		original_CT_cluster <- CT_cluster
		if(umap_flag){
			CT_cluster <- t(umap_layout[hits, ]) # I need to transform it bcuz the data gets transformed in the clara call
			original_CT_cluster <- RNAcounts[, hits]
		}
		print(c("dim(CT_cluster):", dim(CT_cluster)))
		if(length(k)){
			new_k <- as.integer(k)
		}else{
			new_k <- floor(ncol(CT_cluster) / expected_cells)
				if(!new_k){
					new_k <- ncol(CT_cluster)
					print(paste("Setting k to", new_k))
				}else{
					print(paste("Setting k to", new_k))
				}
		}
		if(summary_method == "kmed" || summary_method == "kmed_means"){
			print(paste("k=", new_k, "ncol(CT_cluster)=", ncol(CT_cluster)))
			#########################%%%%%%%%%%%%%%%%%%%^^^^^^^^^^^^^^^%%%%%%%%%%%%%%%%%%%###############
			#########################%%%%%%%%%%%%%%%%%%%^^^^^^^^^^^^^^^%%%%%%%%%%%%%%%%%%%###############
			#########################%%%%%%%%%%%%%%%%%%%^^^^^^^^^^^^^^^%%%%%%%%%%%%%%%%%%%###############
				if(new_k >= ifelse(is.null(ncol(original_CT_cluster)), 1, ncol(original_CT_cluster))){
					print("too small... merging all cells into one metacell")
					clusters[[md_combo[i]]]$clustering <- rep(1, ifelse(is.null(ncol(original_CT_cluster)), 1, ncol(original_CT_cluster)))
					cluster_data[[md_combo[i]]] <- t(as.matrix(original_CT_cluster))
					clusters[[md_combo[i]]]$data <- t(CT_cluster)
					print(paste("dim(cluster_data)", dim(cluster_data[[md_combo[i]]])))
					print(paste("length(clusters$clustering)", length(clusters[[md_combo[i]]]$clustering)))
					if(merge_flag){
						clusters[[md_combo[i]]]$clustering <- merge_small_mc(clusters[[md_combo[i]]]$clustering, thresh= expected_cells)
					}
					print(table(clusters[[md_combo[i]]]$clustering))
					## When original_CT_cluster is a vector, it loses its colname. That's why I had to add the line below to extract the correct annotation, especially for the RNA_barcodes_ct later. However, this fix will only apply to Seurat input data, and not the csv data. N.B. I'm using RNAcounts to infer the names.
					if(length(hits) == 1){
						rownames(cluster_data[[md_combo[i]]]) <- hits # names(RNAcounts[, hits])
					}else{
						rownames(cluster_data[[md_combo[i]]]) <- colnames(RNAcounts[, hits])
					}
					RNA_barcodes_ct <- c(RNA_barcodes_ct, rownames(cluster_data[[md_combo[i]]]))
					cell2metacell_info <- c(cell2metacell_info, paste(md_combo[i], clusters[[md_combo[i]]]$clustering, sep= "_"))

					if(length(assay_slot)){
						RNA_metacell_umap_ct <- NULL
						for(j in unique(clusters[[md_combo[i]]]$clustering)){
							data_subset <- RNA_umap[which(clusters[[md_combo[i]]]$clustering == j), ];
							if(is.null(dim(data_subset))){
								RNA_metacell_umap_ct <- rbind(RNA_metacell_umap_ct, data_subset);
							}else{
								RNA_metacell_umap_ct <- rbind(RNA_metacell_umap_ct, colMeans(data_subset))
							}
						}
						rownames(RNA_metacell_umap_ct) <- paste(md_combo[i], unique(clusters[[md_combo[i]]]$clustering), sep= "_")
						print(paste("Done clustering", md_combo[i]))
						RNA_metacell_umap <- rbind(RNA_metacell_umap, RNA_metacell_umap_ct)
					}
					### merge all cells into one metacell
				}else if(length(new_k) >0 && new_k < ncol(CT_cluster)){
					clara_out <- invoke_clara_simplified(CT_cluster, new_k)
					mc_qual <- metacell_quality(clara_out$clusters, expected_cell_num= expected_cells, slack_ratio= .15)
					outlier_idx <- get_outliers(clara_out$clusters, 10)
					print(c("Pass:", pass_num, "valid metacells:", length(which(mc_qual == T)), "outliers:", length(outlier_idx)))
					mc_quality_info[[md_combo[i]]][[pass_num]] <- length(which(mc_qual == T))
					mc_outlier_info[[md_combo[i]]][[pass_num]] <- length(outlier_idx)
					mc_distr[[md_combo[i]]][[pass_num]] <- clara_out$clusters$clustering
					## Plot heatmap
					ann_col <- data.frame(metacells= factor(clara_out$clusters$clustering))
					rownames(ann_col) <- seq(ncol((CT_cluster)))
					colnames(CT_cluster) <- seq(ncol((CT_cluster)))
					###############
					cluster_data[[md_combo[i]]] <- t(as.matrix(original_CT_cluster))
					## plot the metacell count table with a dendrogram obtained from mean clusters
					cluster_means_res <- cluster_means(clara_out$clusters)
					print(paste("dim(cluster_means_res):", dim(cluster_means_res)))
					tryCatch(
						{
							hc_res <- hclust(d = dist(x = cluster_means_res))
							df <- data.frame(table(clara_out$clusters$clustering));
							df$y <- as.factor(rep(1, nrow(df)))
							df_ordered <- df[hc_res$order, ]
							df_ordered$Var1 <- factor(df_ordered$Var1, levels= df_ordered$Var1)
							p1 <- ggplot2::ggplot(df_ordered, ggplot2::aes(y= Var1, x= y)) + geom_tile(show.legend= F, color= "gray", fill= "white") + theme_classic() + geom_text(ggplot2::aes(label= Freq), size= 1) + theme(axis.text = element_text(size = 3)) + theme(axis.text.y = element_text(size = rel(1), hjust = 1, angle = 0), 
								# margin: top, right, bottom, and left ->  plot.margin = unit(c(1, 0.2, 0.2, -0.5), "cm"), 
								panel.grid.minor = element_blank()) + ggtitle(md_combo[i])
							p2 <- ggdendrogram(data = as.dendrogram(hc_res), rotate= T, segments= T, leaf_labels= F) + geom_text(size= 1)
							print(plot_grid(p1, p2, rel_widths=c(.25, 1), ncol= 2))#align= "h"
						},
						error= function(e){
							message(paste("dim(cluster_means_res):", dim(cluster_means_res)))
							message(e)
						}
					)
					## update k after removing outliars
					ks_info[[md_combo[i]]][[pass_num]] <- new_k
					new_k <- floor(ncol(CT_cluster) / expected_cells)
				
					#clust <- makeCluster(ceiling(detectCores()/2))
					clusters[[md_combo[i]]] <- clara_out$clusters
					print(table(clusters[[md_combo[i]]]$clustering))
					if(merge_flag){
						clusters[[md_combo[i]]]$clustering <- merge_small_mc(clusters[[md_combo[i]]]$clustering, thresh= expected_cells)
					}
					print(table(clusters[[md_combo[i]]]$clustering))
				}
				
				
				if(length(assay_slot)){
					RNA_metacell_umap_ct <- NULL
					for(j in unique(clusters[[md_combo[i]]]$clustering)){
						#data_subset <- RNA_umap[which(clusters[[ct]]$clustering == i), ];
						if(is.null(nrow(clusters[[md_combo[i]]]$data[which(clusters[[md_combo[i]]]$clustering == j), ]))){
							data_subset <- clusters[[md_combo[i]]]$data[, which(clusters[[md_combo[i]]]$clustering == j)];
						}else{
							data_subset <- clusters[[md_combo[i]]]$data[which(clusters[[md_combo[i]]]$clustering == j), ];
						}
						if(is.null(dim(data_subset))){
							RNA_metacell_umap_ct <- rbind(RNA_metacell_umap_ct, data_subset);
						}else{
							RNA_metacell_umap_ct <- rbind(RNA_metacell_umap_ct, colMeans(data_subset))
						}
					}
					rownames(RNA_metacell_umap_ct) <- paste(md_combo[i], unique(clusters[[md_combo[i]]]$clustering), sep= "_")
					print(paste("Done clustering", md_combo[i]))
					RNA_metacell_umap <- rbind(RNA_metacell_umap, RNA_metacell_umap_ct)
					RNA_metacell_umap_plot[[md_combo[i]]] <- RNA_metacell_umap_ct
				}
				RNA_barcodes_ct <- c(RNA_barcodes_ct, rownames(cluster_data[[md_combo[i]]]))
				cell2metacell_info <- c(cell2metacell_info, paste(md_combo[i], clusters[[md_combo[i]]]$clustering, sep= "_"))
				
			#########################%%%%%%%%%%%%%%%%%%%^^^^^^^^^^^^^^^%%%%%%%%%%%%%%%%%%%###############
			#########################%%%%%%%%%%%%%%%%%%%^^^^^^^^^^^^^^^%%%%%%%%%%%%%%%%%%%###############
			#########################%%%%%%%%%%%%%%%%%%%^^^^^^^^^^^^^^^%%%%%%%%%%%%%%%%%%%###############
		}
		else{
			error("Undefined method of summarization. Please pick either kmed, kmeans, or kmed_means!")
		}
	
	}
	dev.off()
	save(RNA_barcodes_ct, cell2metacell_info, ks_info, mc_distr, mc_quality_info, mc_outlier_info, cluster_data, clusters, file= paste0(output_file, "/debug/clusters_debug.Rdata"))
	print("Done clustering!")

	mat <- NULL
	mat_sum <- NULL
	mc_names <- NULL
	for(i in seq(length(clusters))){
		if(summary_method == "kmed"){
			mat <- cbind(mat, t(clusters[[i]]$medoids))
		}else if(summary_method == "kmeans"){
			mat <- cbind(mat, t(clusters[[i]]$centers))
		}
		else{##kmed_means
			temp_cl <- NULL
			temp_cl_sum <- NULL
				if(umap_flag){
					for(clst in unique(clusters[[i]]$clustering)){
						idx <- which(clusters[[i]]$cluster == clst)
						mc_names <- c(mc_names, paste0(names(clusters)[i], "_", clst))
						if("numeric" %in% class(cluster_data[[i]][idx, ])){
							temp_cl <- rbind(temp_cl, cluster_data[[i]][idx, ])
							temp_cl_sum <- rbind(temp_cl_sum, cluster_data[[i]][idx, ])
						}else{
							temp_cl <- rbind(temp_cl, apply(cluster_data[[i]][idx, ], 2, FUN= mean))
							temp_cl_sum <- rbind(temp_cl_sum, apply(cluster_data[[i]][idx, ], 2, FUN= sum))
						}
					}
				}else{
				for(clst in unique(clusters[[i]]$clustering)){
						idx <- which(clusters[[i]]$cluster == clst)
				if("numeric" %in% class(clusters[[i]]$data[idx, ])){
					temp_cl <- rbind(temp_cl, clusters[[i]]$data[idx, ])
					temp_cl_sum <- rbind(temp_cl_sum, clusters[[i]]$data[idx, ])
				}
				else{
							idx <- which(clusters[[i]]$cluster == clst)
					temp_cl <- rbind(temp_cl, apply(clusters[[i]]$data[idx, ], 2, FUN= mean))
					temp_cl_sum <- rbind(temp_cl_sum, apply(clusters[[i]]$data[idx, ], 2, FUN= sum))
				}
				}
				}
			mat <- cbind(mat, t(temp_cl))
			mat_sum <- cbind(mat_sum, t(temp_cl_sum))
		}
	}

	print("done making mat")

	colnames(mat) <- colnames(mat_sum) <- mc_names
	if(!normalization_flag){
		write.csv(mat, paste0(output_file, "/results/RNA_cellSummarized_normalized", ".csv"))
		write.csv(mat_sum, paste0(output_file, "/results/RNA_cellSummarized_normalized_sum", ".csv"))
	}

	print("saving...")
	save(mat, mat_sum, file= paste0(output_file, "/debug/mat.RData"))
	print("done saving!")
	## Normalize the average metacell read counts by sequencing depth to compute the CPM values
	if(normalization_flag){
		dds <- DESeq2::DESeqDataSetFromMatrix(mat_sum, S4Vectors::DataFrame(colnames(mat_sum)), ~1);
		dds <- DESeq2::DESeq(dds);
		write.table(DESeq2::sizeFactors(dds), paste0(output_file, "/debug/DESeq2_sizeFactors.txt"))
		normalized_counts_all <- DESeq2::counts(dds, normalized=TRUE)
		cors <- sapply(seq(nrow(normalized_counts_all)), function(i)cor(as.numeric(mat_sum[i, ]), as.numeric(normalized_counts_all[i, ])))
		inspect_genes <- which(cors > .4 & cors < .6)[seq(5)]
		inspect_genes <- na.omit(inspect_genes)
		print(c("inspect_genes:", inspect_genes))
		pdf(paste0(output_file, "/plots/RNA_DESeq2_inspection.pdf"))
		print(ggplot2::ggplot(data.frame(size_factor= DESeq2::sizeFactors(dds))) + ggplot2::geom_violin(ggplot2::aes(x= "DESeq2", y= size_factor)) + ggplot2::theme_classic())
		print(ggplot2::ggplot(data.frame(correlation = cors)) + ggplot2::geom_violin(ggplot2::aes(x= "normalized vs raw", y= correlation)) + ggplot2::theme_classic())
		plot(log2(1 + rowMeans(normalized_counts_all)), log2(1 + rowMeans(mat_sum)), xlab= "log2(1 + avg. DESeq2 normalized)", ylab= "log2(1 + avg. raw)", pch=20, main= "RNA")
		if(length(inspect_genes)){
			for(i in inspect_genes){
				df <- data.frame(DESeq2= log2(1 + normalized_counts_all[i, ]), raw= log2(1 + mat_sum[i,]))
				pheatmap::pheatmap(df, main= rownames(normalized_counts_all)[i])
			}
		}
		dev.off()
	}
	##############################
	##############################
	final_umap_res <- uwot::umap(t(mat), pca= 30, pca_center= T, scale= T, n_components= umap_dim)
	rownames(final_umap_res) <- colnames(mat)
	colnames(final_umap_res) <- paste0("UMAP", seq(umap_dim))

	#celltypes <- sapply(colnames(mat), function(i) strsplit(i, "_")[[1]][1])
	metadata_ann <- sub("_[^_]+$", "", colnames(mat))
	metadata_tokens <- sapply(metadata_ann, function(i) strsplit(i, "_")[[1]])


	pdf(paste0(output_file, "/plots/umap", ".pdf"))
	for(j in seq(length(metadata))){
		df <- data.frame(UMAP1= final_umap_res[, 1], UMAP2= final_umap_res[, 2], metadata= sapply(metadata_ann, function(i) strsplit(i, "_")[[1]][j]))
		print(ggplot2::ggplot(df, ggplot2::aes(x= UMAP1, y= UMAP2)) + ggplot2::geom_point(ggplot2::aes(color= metadata)) + ggplot2::theme_classic() + ggplot2::geom_text(ggplot2::aes(label= metadata),hjust=0, vjust=0, size= 3, check_overlap = T))
	}
	dev.off()
	##############################
	##############################

	if(length(assay_slot)){
		#kk <- 5L
		#knn_res <- class::knn(train= RNA_metacell_umap, test= ATAC_umap, cl= rownames(RNA_metacell_umap), k= kk)
		save(RNA_metacell_umap, ATAC_umap, threshold, file= paste0(output_file, "/debug/my_knn_debug.Rdata"))
		print("Starting the KNN classification of ATAC cells into metacells:")
		knn_res <- my_knn(train= RNA_metacell_umap, test= ATAC_umap, threshold= threshold)
		print("Finished classifying the ATAC cells")
		#atac2metacell_info <- data.frame(barcode= rownames(ATAC_umap), metacell= knn_res) ## the sorting of the distances in the my_knn function causes the rownames not to match anymore. therefore, I'm commenting this line and changing to the following
		atac2metacell_info <- data.frame(barcode= names(knn_res), metacell= knn_res)
		write.csv(atac2metacell_info, paste0(output_file, "/results/ATAC_cell2metacell_info", ".csv"), row.names= F)
	

		train <- RNA_metacell_umap
		test <- ATAC_umap
		colnames(train) <- NULL
		colnames(test) <- NULL                                                                                                                                   
		df <- data.frame(train, assay= "mcRNA", celltype= sapply(rownames(train), function(i)strsplit(i, "_")[[1]][1]))
		df <- rbind(df, data.frame(test, assay= "scATAC", celltype= "unlabled"))
		colnames(df)[seq(ncol(train))] <- paste0("UMAP", seq(ncol(train)))

		library(viridis)
		pdf(paste0(output_file, "/plots/mcRNA_ATAC_UMAP.pdf"))

		print(ggplot2::ggplot(df, ggplot2::aes(x= UMAP1, y= UMAP2, shape= assay, colour= celltype, label= celltype)) + ggplot2::geom_point(alpha= .5) + ggplot2::ggtitle("UMAP of mcRNAs vs scATACs") + ggplot2::geom_text(ggplot2::aes(label= ifelse(celltype != "unlabled", as.character(celltype), '')), hjust=0, vjust=0, size= 3, check_overlap = T) + viridis::scale_color_viridis(discrete=TRUE) + ggplot2::theme_bw())
		dev.off()
	}

	if(length(gtf_file) & normalization_flag){
	## Compute TPM on normalized metacells ## Following the suggestion here: biostars.org/p/456800/
		library(GenomicFeatures)
		gtf <- as.data.frame(rtracklayer::import(gtf_file))
		gtf_genes <- subset(gtf, type == "gene")

		txdb <- makeTxDbFromGFF(gtf_file, format="gtf")
		exons.list.per.gene <- exonsBy(txdb, by="gene")
		exonic.gene.sizes <- as.data.frame(sum(width(reduce(exons.list.per.gene))))

		gene_df <- merge(gtf_genes[, c("gene_id", "gene_name")], exonic.gene.sizes, by.x= "gene_id", by.y= 0)
		dup_hits <- which(duplicated(gene_df$gene_name) == T)
		for(dpl in unique(gene_df$gene_name[dup_hits])){
			tmp_row <- gene_df[which(gene_df$gene_name == dpl)[1], ]
			tmp_row[, 3] <- floor(mean(gene_df[-which(gene_df$gene_name == dpl), 3])) ## Take the average of exonic size for the duplicated genes
			gene_df <- gene_df[-which(gene_df$gene_name == dpl), ]
			gene_df <- rbind(gene_df, tmp_row)
		}
		colnames(gene_df)[3] <- "exonic_len"

		expr_mat_df <- merge(normalized_counts_all, gene_df[, c(2, 3)], by.x= 0, by.y= "gene_name")
		rownames(expr_mat_df) <- expr_mat_df$Row.names
		expr_mat_df <- expr_mat_df[, -which(colnames(expr_mat_df) == "Row.names")]
		exonic.gene.sizes.expr <- expr_mat_df$exonic_len
		expr_mat_df <- expr_mat_df[, -which(colnames(expr_mat_df) == "exonic_len")]
		norm_len <- expr_mat_df / matrix(rep(exonic.gene.sizes.expr, each= ncol(expr_mat_df)), ncol= ncol(expr_mat_df), byrow= T)
		tpm.mat <- t( t(norm_len) * 1e6 / colSums(norm_len) )
		write.csv(tpm.mat, paste0(output_file, "/results/RNA_cellSummarized_normalized", ".csv"))
		write.csv(normalized_counts_all, paste0(output_file, "/results/cellSummarized_normalized", ".csv"))
	}

	rna2metacell_info <- data.frame(barcode= RNA_barcodes_ct, metacell= cell2metacell_info)
	write.csv(rna2metacell_info, paste0(output_file, "/results/RNA_cell2metacell_info", ".csv"), row.names= F)
	write.csv(mat, paste0(output_file, "/results/cellSummarized", ".csv"))
	write.csv(mat_sum, paste0(output_file, "/results/cellSummarized", "_sum.csv"))
	write.csv(final_umap_res, paste0(output_file, "/results/RNA_metacell_umap", ".csv"))
	if(length(assay_slot)){
		save(atac2metacell_info, ATACcounts, clusters, RNA_metacell_umap, ATAC_umap, mc_names, file= paste0(output_file, "/debug/", "clustered.RData"))

		uniq_mc <- unique(atac2metacell_info$metacell)
		atac_metacell <- NULL;
		atac_metacell_sum <- NULL;
		for(i in seq(length(uniq_mc))){
			hits <- atac2metacell_info$barcode[which(atac2metacell_info$metacell == uniq_mc[i])];
			if(length(hits) > 1){
				atac_metacell <- cbind(atac_metacell, Matrix::rowMeans(ATACcounts[, hits]))
				#atac_metacell <- cbind(atac_metacell, apply(ATACcounts[, hits]), 2, apply= mean)
				atac_metacell_sum <- cbind(atac_metacell_sum, Matrix::rowSums(ATACcounts[, hits]))
			}else{
				atac_metacell <- cbind(atac_metacell, ATACcounts[, hits])
				atac_metacell_sum <- cbind(atac_metacell_sum, ATACcounts[, hits])
			}
		}
		colnames(atac_metacell) <- uniq_mc
		colnames(atac_metacell_sum) <- uniq_mc
		if(normalization_flag){
			#atac_metacell_sum[which(atac_metacell_sum == 0)] <- 1 ## If a cell had zero reads across all peaks, I assign a non-zero value to that total reads to avoid division by zero
			total_ATAC <- colSums(atac_metacell_sum)
			total_ATAC[which(total_ATAC == 0)] <- 1
			#ATACcounts_depthNorm <- ATACcounts / matrix(rep(total_ATAC, times= nrow(ATACcounts)), nrow= ATACcounts, byrow= T)
			ATACcounts_depthNorm <- t(t(atac_metacell_sum) / total_ATAC)
			ATACcounts_depthNorm <- ATACcounts_depthNorm * 1e6
			atac_metacell_sum <- ATACcounts_depthNorm
			#cors <- sapply(seq(nrow(ATACcounts)), function(i)cor(as.nmeric(ATACcounts_raw[i, ]), as.numeric(ATACcounts[i, ])))
			pdf(paste0(output_file, "/plots/ATAC_depthNormalization_inspection.pdf"))
			plot(log2(1 + rowMeans(atac_metacell_sum)), log2(1 + Matrix::rowMeans(ATACcounts)), xlab= "log2(1 + avg depth normalized)", ylab= "log2(1 + raw)", main= "ATAC", pch= 20)
			dev.off()
		}
		write.csv(atac_metacell, paste0(output_file, "/results/cellSummarized_ATAC", ".csv"))
		write.csv(atac_metacell_sum, paste0(output_file, "/results/cellSummarized_ATAC", "_sum.csv"))
	}
	print("Done!")
	Rtnse_plot <- F
	if(Rtnse_plot){
	library(Rtsne)
	Rtsne_whole_res <- Rtsne(as.matrix(RNAcounts), check_duplicates= F)
	pdf(paste0(output_file, "/plots/tSNE", "_Rtsne.pdf"))
	plot(Rtsne_whole_res$Y)
	dev.off()
	}


	if(length(assay_slot)){
	#################################################
	############# BEGIN VISUALIZATION ###############

		addSmallLegend <- function(myPlot, pointSize = 0.5, textSize = 3, spaceLegend = 0.1) {
			myPlot +
				guides(shape = guide_legend(override.aes = list(size = pointSize)),
							color = guide_legend(override.aes = list(size = pointSize))) +
		theme(legend.title = element_text(size = textSize), legend.text  = element_text(size = textSize), legend.key.size = unit(spaceLegend, "lines"))
		}
	####
	###
	##
	## Plot metacells within the scRNA:
	df_mc_code <- NULL;

	for(i in names(RNA_metacell_umap_plot)){
		if(ncol(clusters[[i]]$data) > 1){## in case of having only one data point for the given subset, the dimension of the data gets swapped and I have to catch this exception with this if condition
			df_mc_code <- rbind(df_mc_code, data.frame(umap1= clusters[[i]]$data[, 1], umap2= clusters[[i]]$data[, 2], type= "scRNA", celltype= i));
		}else{
			df_mc_code <- rbind(df_mc_code, data.frame(umap1= clusters[[i]]$data[1, ], umap2= clusters[[i]]$data[2, ], type= "scRNA", celltype= i));
		}
		df_mc_code <- rbind(df_mc_code, data.frame(umap1= RNA_metacell_umap_plot[[i]][, 1], umap2= RNA_metacell_umap_plot[[i]][, 2], type= "mcRNA", celltype= i))
	}
	pdf(paste0(output_file, "/plots/mcRNA_vs_scRNA_UMAP.pdf")); print(ggplot2::ggplot(df_mc_code, ggplot2::aes(x= umap1, y= umap2, colour= celltype, shape= type)) + geom_point());
	for(i in unique(df_mc_code$celltype)){
		print(ggplot2::ggplot(subset(df_mc_code, celltype == i), ggplot2::aes(x= umap1, y= umap2, colour= type, shape= celltype)) + geom_point() + theme_classic() + ggtitle(i))
	}
	dev.off()
	############
	############

		ATAC_info <- atac2metacell_info
		RNA_info <- rna2metacell_info
		df <- data.frame(umap1= Sub[[reduction]]@cell.embeddings[, 1], umap2= Sub[[reduction]]@cell.embeddings[, 2], assays= assays)
		df$Seurat_celltype <- original_celltypes
	#df$Seurat_celltype <- as.character(Sub@meta.data$seurat_clusters)
		df <- cbind(df[c(ATAC_info$barcode, RNA_info$barcode), ], rbind(ATAC_info, RNA_info))
		celltypes <- sapply(as.character(df$metacell), function(i) strsplit(i, "_")[[1]][1])
		uniq_celltypes <- unique(celltypes)
		df$celltype <- celltypes

		pdf(paste0(output_file, "/plots/plot_umap_mc_", expected_cells, ".pdf"));
		print(ggplot2::ggplot(df, ggplot2::aes(x= umap1, y= umap2, shape= assays, colour= Seurat_celltype)) + geom_point(alpha= .6) + theme_classic() + geom_text(ggplot2::aes(label= Seurat_celltype),hjust=0, vjust=0, size= 3, check_overlap = T));
		for(i in seq(length(unique(celltypes)))) {
			myPlot <- ggplot2::ggplot(subset(df, celltype == uniq_celltypes[i]), ggplot2::aes(x= umap1, y= umap2, shape= assays)) + geom_point(ggplot2::aes(colour= metacell)) + ggtitle(uniq_celltypes[i]) + theme_classic() + geom_text(ggplot2::aes(label= metacell),hjust=0, vjust=0, size= 3, check_overlap = T); print(addSmallLegend(myPlot))
		};
		dev.off()

	###############
		df2 <- data.frame(umap1= Sub[[reduction]]@cell.embeddings[, 1], umap2= Sub[[reduction]]@cell.embeddings[, 2], assays= assays, Seurat_celltype= original_celltypes)
		pdf(paste0(output_file, "/plots/plot_umap_seurat_", expected_cells, ".pdf")); for(i in seq(length(unique(celltypes)))) {print(ggplot2::ggplot(subset(df2, Seurat_celltype == uniq_celltypes[i]), ggplot2::aes(x= umap1, y= umap2, shape= Seurat_celltype, colour= assays)) + geom_point(alpha= .6) + ggtitle(uniq_celltypes[i]) + theme_classic())}; dev.off()
	###############

		rna_cnt <- table(RNA_info$metacell)
		atac_cnt <- table(ATAC_info$metacell)

		df_main <- data.frame(RNA_mc_cnt= as.numeric(rna_cnt), row.names= names(rna_cnt))
		df_main$ATAC_mc_cnt <- 0
		df_main[names(atac_cnt), ]$ATAC_mc_cnt <- as.numeric(atac_cnt)
		df_main$celltype <- sapply(rownames(df_main), function(i) strsplit(i, "_")[[1]][1])
		print(paste("expected_cells=", expected_cells))
		pdf(paste0(output_file, "/plots/MC_scatterplots_", expected_cells, ".pdf"))
		print(ggplot2::ggplot(df_main, ggplot2::aes(x= RNA_mc_cnt, y= ATAC_mc_cnt)) + geom_point(colour= "blue") + theme_classic() + ggtitle("all cell types"))
		for(ct in unique(df_main$celltype))
			print(ggplot2::ggplot(subset(df_main, celltype == ct), ggplot2::aes(x= RNA_mc_cnt, y= ATAC_mc_cnt)) + geom_point(colour= "blue") + theme_classic() + ggtitle(ct))

		library(data.table)
		dt_main <- as.data.table(df_main)
		dt <- melt(dt_main, id.vars= "celltype", variable.name= "assay", value.name= "mc_count")


		print(ggplot2::ggplot(dt, ggplot2::aes(x= celltype, y= mc_count)) + geom_boxplot(ggplot2::aes(fill= assay)) + geom_hline(yintercept= expected_cells, linetype="dashed", color = "green") + theme_classic() + geom_text(ggplot2::aes(-1, expected_cells, label = "exp. num. cells", hjust = 0, vjust= 1), color= "green") + theme(axis.text.x = element_text(angle = 45)))
	#######
		rna_rc <- mat_sum
		atac_rc <- atac_metacell_sum 
		rna_rc_sum <- colSums(rna_rc)
		atac_rc_sum <- colSums(atac_rc)
	#######

		rc_dt <- data.table(read_counts= c(rna_rc_sum, atac_rc_sum), assay= c(rep("RNA", length(rna_rc_sum)), rep("ATAC", length(atac_rc_sum))))
		rc_dt_log <- rc_dt
		rc_dt_log$read_counts <- log2(1 + rc_dt_log$read_counts)

		print(ggplot2::ggplot(rc_dt_log, ggplot2::aes(x= assay, y= read_counts)) + geom_boxplot() + theme_classic() + ylab("log2(1 + sum of read counts per metacell)"))

		dev.off()

		if(F){## I have to add the UMAP of ATAC counts and then make the plots... right now it's only for the RNAs.
			pdf(paste0(output_file, "/plots/MC_umaps_", expected_cells, ".pdf"))
			umap_df <- data.frame(UMAP1= umap_res[, 1], UMAP2= umap_res[, 2])
			umap_df_reordered <- umap_df[RNA_info$barcode,]
			umap_df_reordered$MC <- RNA_info$metacell
			ggplot2::ggplot(umap_df_reordered, ggplot2::aes(x= UMAP1, y= UMAP2)) + geom_point(ggplot2::aes(colour= MC), alpha= .5) + theme_classic() + theme(legend.position='none')
			umap_celltypes <- sapply(umap_df_reordered$MC, function(i) strsplit(i, "_")[[1]][1])
			for(ct in unique(umap_celltypes)){
				hits <- which(umap_celltypes == ct)
				print(ggplot2::ggplot(umap_df_reordered[hits, ], ggplot2::aes(x= UMAP1, y= UMAP2)) + geom_point(ggplot2::aes(colour= MC), alpha= .5) + theme_classic() + theme(legend.position='none') + ggtitle(ct))
			}
			dev.off()
		}
	}
}

# #50 cells per metacell
# runtime_50 <- system.time({
# 	metacellar.run(file_name= "/projects/single_cell_stitchit/work/binned_RF/input/paired_PBMC10k/seurat/10xpbmc_10k_scatac_scrna_coembed_Fatemeh.rds", 
# 			RNA_count_slot= "assays$RNA@counts", celltype_info= "meta.data$seurat_annotations", 
# 			output_file= "/projects/single_cell_stitchit/work/binned_RF/input/paired_PBMC10k/metaCellaR_50_normalized", 
# 			umap_flag= T, assay_slot= "meta.data$orig.ident", ATAC_count_slot= "assays$ATAC@counts", expected_cells= 50, 
# 			umap_dim= 20, reduction="umap", gtf_file= "/projects/single_cell_stitchit/work/binned_RF/input/paired_PBMC10k/gencode.v38.annotation.gtf", 
# 			metadata= c("seurat_annotations"), normalization_flag=T)
# })

# print(runtime_50)


# #100 cells per metacell
# runtime_100 <- system.time({
# 	metacellar.run(file_name= "/projects/single_cell_stitchit/work/binned_RF/input/paired_PBMC10k/seurat/10xpbmc_10k_scatac_scrna_coembed_Fatemeh.rds", 
# 			RNA_count_slot= "assays$RNA@counts", celltype_info= "meta.data$seurat_annotations", 
# 			output_file= "/projects/single_cell_stitchit/work/binned_RF/input/paired_PBMC10k/metaCellaR_100_normalized", 
# 			umap_flag= T, assay_slot= "meta.data$orig.ident", ATAC_count_slot= "assays$ATAC@counts", expected_cells= 100, 
# 			umap_dim= 20, reduction="umap", gtf_file= "/projects/single_cell_stitchit/work/binned_RF/input/paired_PBMC10k/gencode.v38.annotation.gtf", 
# 			metadata= c("seurat_annotations"), normalization_flag=T)
# })

# print(runtime_100)

# #150 cells per metacell
# runtime_150 <- system.time({
# 	metacellar.run(file_name= "/projects/single_cell_stitchit/work/binned_RF/input/paired_PBMC10k/seurat/10xpbmc_10k_scatac_scrna_coembed_Fatemeh.rds", 
# 			RNA_count_slot= "assays$RNA@counts", celltype_info= "meta.data$seurat_annotations", 
# 			output_file= "/projects/single_cell_stitchit/work/binned_RF/input/paired_PBMC10k/metaCellaR_150_normalized", 
# 			umap_flag= T, assay_slot= "meta.data$orig.ident", ATAC_count_slot= "assays$ATAC@counts", expected_cells= 150, 
# 			umap_dim= 20, reduction="umap", gtf_file= "/projects/single_cell_stitchit/work/binned_RF/input/paired_PBMC10k/gencode.v38.annotation.gtf", 
# 			metadata= c("seurat_annotations"), normalization_flag=T)
# })

# print(runtime_150)
