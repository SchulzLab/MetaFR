library(Seurat)
#options(Seurat.object.assay.version = "v4")

library(ArchR)
library(hdf5r)
library(parallel)
library(ggplot2)

set.seed(1)
addArchRThreads(threads = 25) 
addArchRGenome("hg38")

#test cells
#20% random test set 
test_cells_path <- "/projects/single_cell_stitchit/work/binned_RF/input/paired_PBMC10k/test_cells.txt"
test_cells_df <- read.table(test_cells_path, header=FALSE, sep="\t")
test_cells <- test_cells_df$V1

# use seurat object with celltypes is metadata and CPM (scaling factor = 10.000)
rna.seurat.train <- readRDS("/projects/single_cell_stitchit/work/scarlink/pbmc_input/pbmc_scrna_training_cells_celltype_annotation.rds")
out.scrna.file.train <- "/projects/single_cell_stitchit/work/scarlink/pbmc_input/pbmc_scrna_training_cells_celltype_annotation_cell_prefix_500bp.rds"
rna.seurat.train <- RenameCells(rna.seurat.train, new.names=sub("pbmc_train#", "",colnames(rna.seurat.train)))
rna.seurat.train <- RenameCells(rna.seurat.train, new.names=paste0("pbmc_train_500bp#", colnames(rna.seurat.train)))
saveRDS(rna.seurat.train, file = out.scrna.file.train)

### create arrow file from fragment file ###

#use training cells only
fragment.train.file <- "/projects/single_cell_stitchit/work/binned_RF/input/paired_PBMC10k/seurat/data/atac/pbmc_granulocyte_sorted_10k_atac_fragments_test_train.tsv.gz"

sample.id <- c("pbmc_train_500bp")
names(fragment.train.file) <- sample.id

#tile size = 500
ArrowFiles <- createArrowFiles(
  inputFiles = fragment.train.file,
  sampleNames = names(fragment.train.file),
  minTSS = 4, #Dont set this too high because you can always increase later
  minFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  TileMatParams = list(tileSize = 500, binarize = FALSE),
  offsetPlus = 0, #offset in fragment file
  offsetMinus = 0,
)

# create ArchR project #
scatac.proj <- ArchRProject(
  ArrowFiles = "/projects/single_cell_stitchit/work/scarlink/pbmc_train_500bp.arrow", 
  outputDirectory = "/projects/single_cell_stitchit/work/scarlink/pbmc_input/250kb_500bp_setup/pbmc_scatac_train_500bp",
  copyArrows = TRUE
)

# check if tile matrix is present #
# tm <- getMatrixFromProject(
#  ArchRProj = scatac.proj,
#  useMatrix = 'TileMatrix',
#  binarize = FALSE)

#save ArchR project

scatac.proj <- saveArchRProject(ArchRProj = scatac.proj, outputDirectory = "/projects/single_cell_stitchit/work/scarlink/pbmc_input/250kb_500bp_setup/pbmc_scatac_train_500bp", load = TRUE)

######### create test tile matrix for performance evaluation #########

fragment.test.file <- "/projects/single_cell_stitchit/work/binned_RF/input/paired_PBMC10k/seurat/data/atac/pbmc_granulocyte_sorted_10k_atac_fragments_test_holdout.tsv.gz"
scrna.test.file <- "/projects/single_cell_stitchit/work/scarlink/pbmc_input/pbmc_scrna_test_cells.rds"
out.scrna.file.test <- "/projects/single_cell_stitchit/work/scarlink/pbmc_input/pbmc_scrna_test_cells_cell_prefix_500bp.rds"
rna.seurat.test <- readRDS(scrna.test.file)
rna.seurat.test <- RenameCells(rna.seurat.test, new.names=sub("pbmc_test#", "",colnames(rna.seurat.test)))
rna.seurat.test <- RenameCells(rna.seurat.test, new.names=paste0("pbmc_test_500bp#", colnames(rna.seurat.test)))
saveRDS(rna.seurat.test, file = out.scrna.file.test)
#test archR object

#fragment.test.file <- "/projects/single_cell_stitchit/work/binned_RF/input/paired_PBMC10k/seurat/data/atac/pbmc_granulocyte_sorted_10k_atac_fragments_test_holdout.tsv.gz"
sample.id <- c("pbmc_test_500bp")
names(fragment.test.file) <- sample.id

#tile size
ArrowFiles <- createArrowFiles(
  inputFiles = fragment.test.file,
  sampleNames = names(fragment.test.file),
  minTSS = 4, #Dont set this too high because you can always increase later
  minFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  TileMatParams = list(tileSize = 500, binarize = FALSE),
  offsetPlus = 0, #offset in fragment file
  offsetMinus = 0,
  force = TRUE
)

scatac.proj <- ArchRProject(
  ArrowFiles = "/projects/single_cell_stitchit/work/scarlink/pbmc_test_500bp.arrow", 
  outputDirectory = "/projects/single_cell_stitchit/work/scarlink/pbmc_input/250kb_500bp_setup/pbmc_scatac_test_500bp",
  copyArrows = TRUE
)

scatac.proj <- saveArchRProject(ArchRProj = scatac.proj, outputDirectory = "/projects/single_cell_stitchit/work/scarlink/pbmc_input/250kb_500bp_setup/pbmc_scatac_test_500bp", load = TRUE)

