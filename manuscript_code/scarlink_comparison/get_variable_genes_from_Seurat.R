library(Seurat)

rna_obj <- readRDS("/projects/single_cell_stitchit/work/scarlink/pbmc_input/pbmc_scrna_training_cells.rds")

rna_obj <- FindVariableFeatures(rna_obj, selection.method = "vst", nfeatures = 2000)

var_genes_df <- data.frame(variable_genes=VariableFeatures(rna_obj))

write.table(
  var_genes_df, 
  file = "/projects/single_cell_stitchit/work/scarlink/pbmc_input/top_2000_variable_genes.txt", 
  sep = "\t", 
  row.names = FALSE, 
  col.names = TRUE, 
  quote = FALSE      
)


