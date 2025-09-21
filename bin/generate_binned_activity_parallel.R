# input: cell/metacell-region ATAC count matrix
# output: LSI normalized counts

# steps: convert input count matrix to Seurat object
    #files needed: count matrix, fragment file
    #IDF-TF normalization
    #FindTopFeatures
    #optional: SVD

#fragment estimation from reads if fragment file not present or reads are collapsed across cell barcodes:
    #divide number of reads by two and round up
    #Rsubread featureCounts paired-read mode: quantifies fragments (no 9 bp-shift) from bam file into given regions

#TODO: featureCounts: adapt min overlap? 
#TODO: include 9 bp shift? -> 10x fragment file: start + 4 and end - 5 -> no shift for now
    
library(Seurat)
library(Signac)
library(Rsubread)
library(data.table)
library(parallel)
#library(tidyr)
library(stringi)
#library(EDASeq)
#library(DESeq2)
library(GenomicFeatures)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
#library(stringr)
library(Signac)
library(future) #needed for Signac parallelization (Feature Counts)
library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)
library(ggfortify)

#' Generate binned LSI activity files for ATAC-seq data
#' This function processes ATAC-seq data to generate binned LSI activity files for specified genes
#' in either single-cell or metacell mode. It reads fragment files, RNA count matrices, and GTF files,
#' normalizes the data, and computes LSI activity for each gene across specified chromosomes.
#' The output is saved in the specified output directory with a given file extension.

############### helper functions ####################

# get most 5' tss
get_tss_position <- function(gtf_entry){
    if (gtf_entry$strand == "-"){
        tss_pos <- gtf_entry$end
    } else if (gtf_entry$strand == "+") {
        tss_pos <- gtf_entry$start
    }
    return(tss_pos)
}
####################################################

#TODO: update parameters in documentation
#TODO: check if all parameters are set in the beginning

#' @import Seurat
#' @import Signac
#' @import Rsubread
#' @import data.table
#' @import parallel 
#' @import tidyr
#' @import GenomicRanges
#' @import GenomicFeatures
#' @import BSgenome.Hsapiens.UCSC.hg19
#' @import BSgenome.Hsapiens.UCSC.hg38
#' @import future
#' @import stringi
#' @import dplyr
#' @import ggfortify
#' @import ggplot2
#' @import reshape2
#' @param fragment_file Path to the fragment file
#' @param out_extension Output file extension for the LSI activity files (default: "_lsi_activity.txt.gz")
#' @param normalize_flag Boolean flag to indicate whether to normalize the RNA counts (default: FALSE)
#' @param scaling_factor Scaling factor for normalization (default: 10000)
#' @param mc_flag Boolean flag to indicate whether to run in metacell mode (default: TRUE)
#' @param seed Random seed for reproducibility (default: 123)
#' @param chrom_sizes_file Path to the chromosome sizes file
#' @param cb_mc_file Path to the cell barcode to metacell mapping file
#' @param rna_file Path to the RNA count matrix
#' @param gtf_file Path to the GTF file for gene annotations
#' @param tmp_dir Temporary directory for intermediate files
#' @param out_dir_activity Output path for the binned (LSI-normalized) activity files
#' @param out_dir_fragments Path to the directory where chromosome fragment files will be saved (single-cell or meta-cell)
#' @param test_partition_flag Boolean flag to indicate whether to partition the data into training and testing sets
#' @param test_cells_file Path to the file containing test cells
#' @param test_cells_percent Percentage of cells to be used for testing
#' @param feature_select_flag Boolean flag to indicate whether to select top features based on quantile
#' @param q_top_features Quantile value for feature selection (e.g., 'q85')
#' @param extension_size Size of the extension around target genes (default: 500000)
#' @param bin_size Size of the bins for genome binning (default: 500)
#' @param n_outer_cores Number of outer cores for parallel processing (default: 5)
#' @param n_inner_cores Number of inner cores for parallel processing (default: 5)
#' @param n_threads_feature_counts Number of threads for feature counts (default: 5)
#' @return The output will be stored in children nodes (plots, results, and debug) of the `output_file` directory
#' @export

activity.bin <- function(
    fragment_file,
    cb_mc_file,
    rna_file,
    gtf_file,
    chrom_sizes_file,
    tmp_dir,
    out_dir_activity,
    out_extension = "_lsi_activity.txt.gz", #TODO: better naming scheme
    out_dir_fragments,
    test_partition_flag = TRUE,
    test_cells_file,
    test_cells_percent = 20,
    feature_select_flag = FALSE,
    q_top_features = NULL,
    normalize_flag = FALSE,
    scaling_factor = 10000,
    mc_flag = TRUE, #single-cell or metacell mode
    extension_size = 500000,
    bin_size = 500,
    n_outer_cores = 5,
    n_inner_cores = 5,
    n_threads_feature_counts = 5,
    seed = 123
) {
    # Set the temporary directory
    unlink(tempdir(), recursive = TRUE)
    Sys.setenv(TMPDIR = tmp_dir)
    tempdir(check = TRUE)

    set.seed(seed)

    # Create output directories
    dir.create(out_dir_activity, showWarnings = FALSE, recursive = TRUE)
    dir.create(out_dir_fragments, showWarnings = FALSE, recursive = TRUE)

    #TODO: provide example RNA counts format
    rna_counts <- read.csv(rna_file, header=TRUE)
    rownames(rna_counts) <- rna_counts[,1]
    colnames(rna_counts) <- gsub("\\.", "_", colnames(rna_counts))
    rna_counts <- rna_counts[,-1]

    # Optional CPM normalization with scaling factor 10.000
    if(normalize_flag){
        rna_seurat <- CreateSeuratObject(counts = rna_counts)
        rna_seurat <- NormalizeData(rna_seurat, normalization.method = "RC", scale.factor=scaling_factor)
        rna_counts <- rna_seurat$RNA@data #normalized
    }   
    rna_cells <- colnames(rna_counts) #paired cells: either multi-ome or meta-cells
    gene_list <- rownames(rna_counts)

    if(!test_partition_flag){
        if (!file.exists(test_cells_file)){
            stop("Test cells file does not exist: ", test_cells_file)
        }
        test_cells_df <- read.table(test_cells_file, header=FALSE, sep="\t")
        test_cells <- test_cells_df[,1]
    } else {
        n_test_cells <- floor((test_cells_percent / 100) * length(rna_cells))  # Calculate the number of test cells
        test_cells <- sample(rna_cells, size = n_test_cells, replace = FALSE)  
        train_cells <- rna_cells[which(!rna_cells%in%test_cells)]
        write.table(as.data.frame(test_cells), test_cells_file, row.names=FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
    }      

    #gene names or ensembl gene ids allowed
    ensemble_flag <- grepl("^ENS", gene_list[1])
    gtf <- as.data.frame(rtracklayer::import(gtf_file))
    gtf_genes <- gtf[which(gtf$type=="gene"),] 

    #TODO: remove after testing
    test_chroms <- paste0("chr", 1:5)
    test_genes <- gtf_genes$gene_name[gtf_genes$seqnames %in% test_chroms]
    gene_list <- intersect(gene_list, test_genes)

    if (ensemble_flag){
    gtf_genes$gene_id <- unlist(lapply(gtf_genes$gene_id, function(id) unlist(strsplit(id, split='.', fixed=TRUE))[1]))
    }

    #all genes in gtf
    gtf_genes_sorted_asc <- gtf_genes[order(gtf_genes$seqnames, gtf_genes$start), ]
    gtf_genes_sorted_desc <- gtf_genes[order(gtf_genes$seqnames, gtf_genes$start, decreasing = TRUE), ]

    #get chromosome sizes
    if (file.exists(chrom_sizes_file)) {
        chrom_sizes_df <- read.table(chrom_sizes_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
        chrom_sizes <- setNames(chrom_sizes_df$V2, chrom_sizes_df$V1)
    } else {
        stop(paste("Chromosome sizes file not found at:", chrom_sizes_file))
    }

    #chromosomes of target genes
    if (ensemble_flag){
        chrom_list <- unique(gtf_genes_sorted_asc[which(gtf_genes_sorted_asc$gene_id%in%gene_list),1])
    }else{
        chrom_list <- unique(gtf_genes_sorted_asc[which(gtf_genes_sorted_asc$gene_name%in%gene_list),1])
    }

    if (mc_flag){
        frag_df <- fread(fragment_file, sep="\t", header = FALSE) #,skip=n)
        frag_df_filtered <- frag_df[frag_df$V1 %in% chrom_list, ]
        #hashmap cb - mc 
        cb_mc_df <- read.csv(cb_mc_file, header=TRUE)
        cb_mc_df$barcode <- sub("_2$", "", cb_mc_df$barcode)
        H <- new.env(hash = TRUE, size = nrow(cb_mc_df))
        for (i in seq_len(nrow(cb_mc_df))){
            key <- cb_mc_df$barcode[i] #key = barcode
            H[[key]] <- cb_mc_df$metacell[i] #value = mc
        }    
        #print(ls.str(H))
    }else{      
        n <- as.integer(system(paste("grep -c '^#' ", fragment_file), intern = TRUE))
        frag_df <- fread(fragment_file, sep="\t", header = FALSE ,skip=n)
        frag_df_filtered <- frag_df[frag_df$V1 %in% chrom_list, ]
        #keep only ATAC cells for which RNA cells exist
        multiome_cells <- intersect(rna_cells, frag_df_filtered$V4)
        frag_df_filtered <- frag_df_filtered[frag_df_filtered$V4 %in% multiome_cells, ]

    }

    process_chrom <- function(chrom){
        print(chrom)
        min_gene <- NULL
        max_gene <- NULL
        gene_identifier <- if(ensemble_flag) "gene_id" else "gene_name"   
        # Get the minimum and maximum gene for the current chromosome
        for (gene in gtf_genes_sorted_asc[[gene_identifier]][which(gtf_genes_sorted_asc$seqnames == chrom & gtf_genes_sorted_asc[[gene_identifier]] %in% gene_list)]){
            if (gene%in%gene_list){
                print(gene)
                min_gene <- gene
                break
            }
        }
        for (gene in gtf_genes_sorted_desc[[gene_identifier]][which(gtf_genes_sorted_desc$seqnames == chrom & gtf_genes_sorted_desc[[gene_identifier]] %in% gene_list)]){
            if (gene %in% gene_list){
                print(gene)
                max_gene <- gene
                break
            }
        }    
        gtf_min <- gtf_genes_sorted_asc[which(gtf_genes_sorted_asc[[gene_identifier]] == min_gene),]
        gtf_max <- gtf_genes_sorted_desc[which(gtf_genes_sorted_desc[[gene_identifier]] == max_gene),]
        min_tss <- get_tss_position(gtf_min)
        max_tss <- get_tss_position(gtf_max)

        #generate bins
        upper_bound <- min(max_tss+extension_size, chrom_sizes[chrom])
        lower_bound <- max(0, min_tss-extension_size)

        num_bins <- floor((upper_bound-lower_bound)/bin_size)

        chr <- rep(chrom, num_bins)
        #TODO: include upper bound bin with less bp and bin_size?
        start <- seq(lower_bound, upper_bound - bin_size, bin_size)
        #end <- seq(lower_bound + bin_size, upper_bound, bin_size)
        end <- pmin(start + bin_size - 1, upper_bound)
        strand <- rep('.', num_bins)
        id <- paste(chr, format(start,scientific=FALSE, trim=TRUE), format(end, scientific=FALSE, trim=TRUE), sep="-")
        bins_df <- data.frame(id, chr, format(start, scientific=FALSE, trim=TRUE), format(end, scientific=FALSE, trim=TRUE), strand, stringsAsFactors = FALSE)
        colnames(bins_df) <- c("id", "chr", "start", "end", "strand")

        chr_frag_df <- frag_df_filtered[V1 == chrom]

        if (mc_flag){
                #overwrite cell barcodes in fragment df with metacell ids
                print("generate metacell fragment file...")
                chr_frag_df$mc_id <- unlist(mget(chr_frag_df$V4, envir = H, ifnotfound = NA))
                chr_frag_df_mc <- chr_frag_df[!is.na(chr_frag_df$mc_id), ]

                #save metacell fragment file
                chr_frag_df_mc$V4 <- chr_frag_df_mc$mc_id  
                #chr_frag_df_mc <- chr_frag_df_mc[, which(!colnames(chr_frag_df_mc) == "mc_id")]
                chr_frag_df_mc[, mc_id  := NULL]
                #start
                chr_frag_df_mc$V2 <- format(chr_frag_df_mc$V2,scientific=FALSE, trim=TRUE)
                #end
                chr_frag_df_mc$V3 <- format(chr_frag_df_mc$V3,scientific=FALSE, trim=TRUE)

                #print(head(chr_frag_df_mc))
            
                fragment_file <- paste0(out_dir_fragments,"/",chrom,"_mc_fragments.tsv")
                if (!file.exists(fragment_file)){
                    print("save metacell fragment file....")
                    fwrite(chr_frag_df_mc, fragment_file, row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE)
                    #generate tabix index
                    system(paste("bgzip", fragment_file))
                    system(paste("tabix -p bed", paste0(fragment_file, ".gz")))
                }
            } else {
                    fragment_file <- paste0(out_dir_fragments,"/",chrom,"_sc_fragments.tsv")
                    if (!file.exists(fragment_file)){   
                        print("save chromosome fragment file....")
                        fwrite(chr_frag_df, fragment_file, row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE)
                        #generate tabix index
                        system(paste("bgzip", fragment_file))
                        system(paste("tabix -p bed", paste0(fragment_file, ".gz")))
                    }
            }
            fragments <- CreateFragmentObject(paste0(fragment_file,".gz"))
            bins <- makeGRangesFromDataFrame(bins_df)

            print("generate chromosome count matrix..")
            plan("multicore", workers = n_threads_feature_counts)
            options(future.globals.maxSize = 100000 * 1024^2)
            fragments_mat <- FeatureMatrix(fragments, bins, process_n = 50000, sep = c("-", "-"), verbose = TRUE)

            print("create chromosome Seurat object...")
            chrom_assay <- CreateChromatinAssay(
                counts = fragments_mat,
                sep = c("-", "-"),
                fragments = paste0(fragment_file, ".gz"),
            )

            seurat_atac <- CreateSeuratObject(
                counts = chrom_assay,
                assay = "binned_activity" #,
                #meta.data = metadata
            )

            # Run TF-IDF normalization
            seurat_atac <- RunTFIDF(seurat_atac)

            process_gene <- function(gene) {

                        activity_file <- paste0(out_dir_activity,"/",gene, out_extension)
                        print(activity_file)

                        if(file.exists(activity_file)) {
                            print(paste("Activity file already exists for gene:", gene))
                           return(NULL)
                        }

                        #print(head(gtf_genes_sorted_asc))
                        gene_gtf <-  gtf_genes_sorted_asc[which(gtf_genes_sorted_asc[[gene_identifier]] == gene), ]
                        #print(gene_gtf)
                        tss_pos <- get_tss_position(gene_gtf)
                        #print(tss_pos)
                        tss_query <- IRanges(start = tss_pos, end = tss_pos)
                        #print(tss_query)
                        #find bin containing TSS
                        tss_bin_idx <- subjectHits(findOverlaps(query=tss_query, subject=subject))
                        #print(paste("TSS bin:", tss_bin_idx))
                        min_bin_idx <- max(1, tss_bin_idx-floor(extension_size/bin_size))
                        #print(min_bin_idx)
                        max_bin_idx <- min(tss_bin_idx+floor(extension_size/bin_size),length(chr))
                        #print(max_bin_idx)
                        lsi_activity_df <- as.data.frame(t(seurat_atac$"binned_activity"@data[min_bin_idx:max_bin_idx,]))
                        #print(head(lsi_activity_df))

                        if (feature_select_flag){
                            #perform feature selection only on cells used for training
                            seurat_temp <- subset(seurat_atac, features = colnames(lsi_activity_df))
                            seurat_temp <- subset(seurat_temp, cells = train_cells)
                            seurat_temp <- FindTopFeatures(seurat_temp, min.cutoff = q_top_features)
                            var_regions <- VariableFeatures(seurat_temp)
                            lsi_activity_df <- lsi_activity_df[ ,-which(colSums(lsi_activity_df)==0)]
                            lsi_activity_df <- lsi_activity_df[,colnames(lsi_activity_df)%in%var_regions]
                        }else{
                            #TODO: remove?
                            lsi_activity_df <- lsi_activity_df[ ,-which(colSums(lsi_activity_df)==0)]
                        }
                        if (is.null((colnames(lsi_activity_df)))){
                            print(paste0(gene, ": no variable bins -> gene is skipped"))
                            return(NULL)
                        }

                        colnames(lsi_activity_df) <- sub("^chr", "", colnames(lsi_activity_df))
                        gene_expr <- rna_counts[which(rownames(rna_counts)==gene),]
                        #print(gene_expr)
                        #missing_cbs <- names(gene_expr)[!which(names(gene_expr)%in%rownames(lsi_activity_df))]
                        #print(paste("missing cbs:", missing_cbs))
                        #if (length(missing_cbs)>0){
                        #        zero_cells <- data.frame(matrix(0, nrow = length(missing_cbs), ncol = ncol(lsi_activity_df)),drop=FALSE)
                        #        colnames(zero_cells) <- colnames(lsi_activity_df)
                        #        rownames(zero_cells) <- missing_cbs
                        #        lsi_activity_df <- rbind(lsi_activity_df, zero_cells)
                        #}
                        ##
                        # TODO: should not be necessary
                        common_cells <- intersect(rownames(lsi_activity_df), names(gene_expr))
                        lsi_activity_df <- lsi_activity_df[common_cells, , drop = FALSE]
                        gene_expr <- gene_expr[common_cells]
                        ##
                        #gene_expr <- gene_expr[match(rownames(lsi_activity_df), names(gene_expr))]
                        lsi_activity_df$Expression <- as.numeric(gene_expr[rownames(lsi_activity_df)])
                        print(lsi_activity_df[1:5, 1:5])
                        print(gene_expr)
                        #TODO: better handling instead of transpose gene expression?
                        if (mc_flag){
                            #lsi_activity_df$Expression <- t(gene_expr)
                            lsi_activity_df$Metacell_id <- rownames(lsi_activity_df)
                            lsi_activity_df <- lsi_activity_df[, c("Metacell_id", setdiff(names(lsi_activity_df), "Metacell_id"))]
                        }else{
                            #lsi_activity_df$Expression <- t(gene_expr)
                            lsi_activity_df$Cell_barcode <- rownames(lsi_activity_df)
                            lsi_activity_df <- lsi_activity_df[, c("Cell_barcode", setdiff(names(lsi_activity_df), "Cell_barcode"))]
                        }
                        fwrite(lsi_activity_df, file = activity_file, sep = "\t", row.names = FALSE,  quote = FALSE, compress = "gzip")
                        rm(list=c('gene_expr','lsi_activity_df')) 
                        return(NULL)    
            } 

            chr_genes <- gtf_genes_sorted_asc[[gene_identifier]][which(gtf_genes_sorted_asc$seqnames == chrom & gtf_genes_sorted_asc[[gene_identifier]] %in% gene_list)]
            print(chr_genes)
            subject <- IRanges(start=as.numeric(start), end=as.numeric(end))

            gene_results <- mclapply(chr_genes, function(gene) {
            tryCatch({
                process_gene(gene)
                return(NULL)
            }, error = function(e) {
                return(list(gene = gene, error = conditionMessage(e)))
            })}, mc.cores = n_inner_cores)

            errors <- sapply(gene_results, function(x) !is.null(x$error))
            if (any(errors)) {
                message("Errors processing genes on ", chrom, ":")
                for (err in gene_results[errors]) {
                    message("Gene ", err$gene, ": ", err$error)
                }
            }

        #TODO: free memory
        rm(list=c('chr_frag_df', 'fragments_mat', 'seurat_atac', 'raw_counts_atac', 'raw_counts_atac_filtered_df', 'seurat_temp'))
        gc()
    }

    # Process each chromosome in parallel
    chrom_results <- mclapply(chrom_list, process_chrom, mc.cores = n_outer_cores)
    error_indices <- sapply(chrom_results, function(x) !is.null(x$error))
    if (any(error_indices)) {
        message("Errors encountered in the following chromosomes:")
        for (err in chrom_results[error_indices]) {
            message("Chromosome ", err$chrom, ": ", err$error)
        }
    }
}

####

#function example call
# activity.bin(fragment_file = "/projects/single_cell_stitchit/work/binned_RF/input/paired_PBMC10k/seurat/data/atac/pbmc_granulocyte_sorted_10k_atac_fragments_test.tsv",
#              cb_mc_file = "/projects/single_cell_stitchit/work/binned_RF/input/paired_PBMC10k/metaCellaR_100_normalized/results/ATAC_cell2metacell_info_valid_metacells.csv",
#              rna_file = "/projects/single_cell_stitchit/work/binned_RF/input/paired_PBMC10k/metaCellaR_100_normalized/results/cellSummarized_sum_RNA_valid_metacells.csv",
#              gtf_file = "/projects/single_cell_stitchit/work/binned_RF/input/paired_PBMC10k/gencode.v38.annotation.gtf",
#              chrom_sizes_file = "/projects/single_cell_stitchit/work/binned_RF/input/hg38_chrom_sizes.txt",
#              tmp_dir = "/projects/single_cell_stitchit/work/binned_RF/tmp",
#              out_dir_activity = "/projects/single_cell_stitchit/work/binned_RF/input/paired_PBMC10k/mc100_1MB500bp_binned_lsi_activity_q0_rna_counts",
#              out_dir_fragments = "/projects/single_cell_stitchit/work/binned_RF/input/paired_PBMC10k/chrom_fragments/metacell_100_test",
#              out_extension = "_lsi_activity.txt.gz",
#              test_partition_flag = FALSE,
#              test_cells_file = "/projects/single_cell_stitchit/work/binned_RF/input/paired_PBMC10k/test_meta_cells_100.txt",
#              test_cells_percent = 20,
#              feature_select_flag = FALSE,
#              q_top_features = "q85",
#              normalize_flag = FALSE,
#              scaling_factor = 10000,
#              mc_flag = TRUE,
#              extension_size = 500000,
#              bin_size = 500,
#              n_outer_cores = 5,
#              n_inner_cores = 5,
#              n_threads_feature_counts = 5,
#              seed=123)



#gtf_file <- "/projects/single_cell_stitchit/work/binned_RF/input/paired_PBMC10k/gencode.v38.annotation.gtf"
#frag_file <- "/projects/single_cell_stitchit/work/binned_RF/input/paired_PBMC10k/seurat/data/atac/pbmc_granulocyte_sorted_10k_atac_fragments_test.tsv"
#filtered metacells: at least 200.000 ATAC reads per mc
#cb_mc_file <- "/projects/single_cell_stitchit/work/binned_RF/input/paired_PBMC10k/metaCellaR_100_normalized/results/ATAC_cell2metacell_info_valid_metacells.csv"
#out_dir_fragments <- "/projects/single_cell_stitchit/work/binned_RF/input/paired_PBMC10k/chrom_fragments/metacell_100"
#test_cells_file <- "/projects/single_cell_stitchit/work/binned_RF/input/paired_PBMC10k/test_meta_cells_100.txt"
#rna_file <- "/projects/single_cell_stitchit/work/binned_RF/input/paired_PBMC10k/metaCellaR_100_normalized/results/cellSummarized_sum_RNA_valid_metacells.csv"
#out_path_lsi <- "/projects/single_cell_stitchit/work/binned_RF/input/paired_PBMC10k/mc100_1MB500bp_binned_lsi_activity_q0"