#!/usr/bin/env Rscript
#options(mc.cores = 25)
#Sys.setenv(OMP_NUM_THREADS = 25)

#!/usr/bin/env Rscript

library(tidyverse)

filter_metacells <- function(
    rna.counts.file, 
    atac.counts.file, 
    atac.cell2metacell.file,
    atac.threshold = 200000,
    output.suffix = "_valid_metacells.csv"
) {
  ##### Helper functions ######
  remove_extension <- function(filepath) {
    sub("\\.[^.]*$", "", filepath)  
  }
  clean_names <- function(x) {
  x %>% 
    gsub("\\s", "_", .) %>%    
    gsub("[.-]", "_", .)       
}
  #############################
  rna.counts <- read.csv(rna.counts.file, header = TRUE, check.names = FALSE)
  atac.counts <- read.csv(atac.counts.file, header = TRUE, check.names = FALSE)
  atac.cell2metacell <- read.csv(atac.cell2metacell.file, header = TRUE, check.names = FALSE)
  
  atac.cell2metacell[, 2] <- clean_names(atac.cell2metacell[, 2])
  colnames(atac.counts) <- clean_names(colnames(atac.counts))
  colnames(rna.counts) <- clean_names(colnames(rna.counts))
  
  if (is.character(rna.counts[1, 1])) {
    rownames(rna.counts) <- rna.counts[, 1]
    rna.counts <- rna.counts[, -1, drop = FALSE]
  }
  
  if (is.character(atac.counts[1, 1])) {
    rownames(atac.counts) <- atac.counts[, 1]
    atac.counts <- atac.counts[, -1, drop = FALSE]
  }

  atac.counts <- atac.counts[, colSums(atac.counts) > atac.threshold, drop = FALSE]
  
  mc.ids <- intersect(colnames(rna.counts), colnames(atac.counts))
  rna.counts <- rna.counts[, mc.ids, drop = FALSE]
  atac.counts <- atac.counts[, mc.ids, drop = FALSE]
  atac.cell2metacell <- atac.cell2metacell[atac.cell2metacell[, 2] %in% mc.ids, ]
  
  write.csv(
    atac.cell2metacell,
    paste0(remove_extension(atac.cell2metacell.file), output.suffix),
    row.names = FALSE
  )
  write.csv(
    rna.counts,
    paste0(remove_extension(rna.counts.file), "_RNA", output.suffix),
    row.names = TRUE
  )

  # not needed for pipeline
  write.csv(
    atac.counts,
    paste0(remove_extension(atac.counts.file), "_ATAC", output.suffix),
    row.names = TRUE
  )
  
  return(list(
    rna.counts = rna.counts,
    atac.counts = atac.counts,
    atac.cell2metacell = atac.cell2metacell
  ))
}

# Example run
#x <- filter_metacells(
#  rna.counts.file = "/projects/single_cell_stitchit/work/binned_RF/input/paired_PBMC10k/metaCellaR_200_normalized/results/cellSummarized_sum.csv",
#  atac.counts.file = "/projects/single_cell_stitchit/work/binned_RF/input/paired_PBMC10k/metaCellaR_200_normalized/results/cellSummarized_ATAC_sum.csv",
#  atac.cell2metacell.file = "/projects/single_cell_stitchit/work/binned_RF/input/paired_PBMC10k/metaCellaR_200_normalized/results/ATAC_cell2metacell_info.csv"
#)
