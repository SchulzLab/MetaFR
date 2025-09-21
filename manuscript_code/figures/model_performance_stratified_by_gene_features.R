library(ggplot2)
library(GenomicFeatures)
library(data.table)
library(dplyr)
library(tidyr)
library(tibble)
library(gridExtra)
library(grid)
#library(MASS)
library(ggpointdensity)
library(viridis)

#model performance
###
#corr_df_mc100_q85_no_scaling <- read.table("/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/lsi_investigation/performance_reports/mc100_1MB100bp_binned_lsi_activity_q85_no_threshold_test_partition_no_scaling_performance_correlation.txt",
#                                           sep="\t",header=TRUE)
#corr_df_aggr_sc_q85_no_scaling_test_all <- read.table("/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/scarlink/sc_1MB100bp_q85_not_scaled_aggregated_sc_to_mc_spearman_test_cells.tsv",
#                                                      sep="\t",header=TRUE)
#rownames(corr_df_aggr_sc_q85_no_scaling_test_all) <- corr_df_aggr_sc_q85_no_scaling_test_all[,1]
#~16888 models
#corr_df_sc_q100_no_scaling <- read.table("/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/lsi_investigation/performance_reports/sc_1MB100bp_binned_lsi_activity_q100_no_threshold_test_partition_no_scaling_performance_correlation.txt",
#                                         sep="\t",header=TRUE)
###
#18937 models
corr_df_sc_q100_no_scaling <- read.table("/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/lsi_investigation/performance_reports/sc_1MB500bp_binned_lsi_activity_q0_no_threshold_performance_correlation.txt",
                                         sep="\t",header=TRUE)
#nrow(corr_df_sc_q100_no_scaling[which(corr_df_sc_q100_no_scaling$Spearman > 0.3),]) = 631
#19341 models
corr_df_mc100_q0_no_scaling <- read.table("/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/lsi_investigation/performance_reports/mc100_1MB500bp_binned_lsi_activity_q0_no_threshold_performance_correlation.txt",
                                          sep="\t",header=TRUE)

### go enrichment gene sets ###
corr_df_mc100_q100_500bp <- read.table("/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/lsi_investigation/performance_reports/mc100_1MB500bp_binned_lsi_activity_q0_no_threshold_performance_correlation.txt",sep="\t",header=TRUE)
corr_df_aggr_sc_q100_500bp <- read.table("/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/lsi_investigation/performance_reports/sc_1MB500bp_q0_aggregated_sc_to_mc_spearman_test_cells.tsv", sep="\t",header=TRUE)

### go enrichment gene sets: high performing genes ###

gene_set_mc100 <- data.frame(rownames(corr_df_mc100_q100_500bp[which(corr_df_mc100_q100_500bp$Spearman>=0.7) ,])) #10368
gene_set_aggr_sc <- data.frame(corr_df_aggr_sc_q100_500bp[which(corr_df_aggr_sc_q100_500bp$Spearman>=0.7),1]) #10430

write.table(gene_set_mc100, file = "/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/go_enrichment/go_gene_set_mc100_1MB500bp_test_spearman_min_07_3959.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gene_set_aggr_sc, file = "/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/go_enrichment/go_gene_set_aggr_sc_1MB500bp_test_spearman_min_07_1750.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# background cleaning
bg <- "/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/go_enrichment/initial_gene_set_no_filtering.tsv"
bg_df <- read.table(bg, sep="\t",header=FALSE)
bg_genes <- bg_df[,1]
bg_clean <- sub("\\.\\d+$", "", bg_genes)
bg_clean <- unique(bg_clean) #32401
writeLines(bg_clean, "/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/go_enrichment/initial_gene_set_no_filtering_no_ext.txt")

query_file <- "/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/go_enrichment/go_gene_set_mc100_1MB500bp_test_spearman_min_05_7246.tsv"
query_genes <- readLines(query_file, warn=FALSE)
query_genes <- trimws(query_genes)
query_genes <- unique(query_genes)
query_genes <- query_genes[query_genes != ""]
writeLines(query_genes, "/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/go_enrichment/initial_gene_set_no_filtering_no_ext.txt")

#corr_df_mc100_q85_no_scaling <- read.table("/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/lsi_investigation/performance_reports/mc100_1MB100bp_binned_lsi_activity_q85_no_threshold_test_partition_no_scaling_performance_correlation.txt",
#                                           sep="\t",header=TRUE)
#corr_df_aggr_sc_q85_no_scaling_test_all <- read.table("/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/scarlink/sc_1MB100bp_q85_not_scaled_aggregated_sc_to_mc_spearman_test_cells.tsv",
#                                                      sep="\t",header=TRUE)

##

### assess influence genomic features on model performance ###

gene_df <- fread("/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/GeneFeatureMat.txt.gz")
#TODO: train?
rna_test_df <- fread("/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/rna/10x_pbmc_rna_counts_test_cells.tsv.gz", sep="\t")
gene_names <- rna_test_df[,1]

#keep only genes where features are available and vice versa
corr_df <- corr_df_sc_q100_no_scaling #single-cell: 18937  
#corr_df <- corr_df_mc100_q0_no_scaling 
common_genes <- intersect(toupper(rownames(corr_df)), toupper(gene_df$'Gene_name')) #18918
# TODO: debug
# genes that are present in performance file but not input genes (from count matrix)? extension
ext_genes <- setdiff(toupper(common_genes), toupper(rna_test_df$V1))
#common_genes_filtered <- setdiff(toupper(rna_test_df$V1), toupper(common_genes)) #TODO: remove when gene extensions are saved
corr_df_filtered <- corr_df[which(toupper(rownames(corr_df))%in%toupper(common_genes)),] 
corr_df_filtered <- corr_df_filtered[-which(toupper(rownames(corr_df_filtered))%in%ext_genes),] #18594 #TODO: remove
corr_df_filtered$`Gene_name` <- toupper(rownames(corr_df_filtered))

gene_df_filtered <- gene_df[which(toupper(gene_df$`Gene_name`)%in%toupper(common_genes)),]
gene_df_filtered$Gene_name_without_ext <- sub("\\..*", "", toupper(gene_df_filtered$`Gene_name`))
gene_df_filtered <- gene_df_filtered[-which(toupper(gene_df_filtered$`Gene_name_without_ext`)%in%ext_genes),] #18611 #19142 #TODO: remove

#TODO: proper gene name storing with extension
#keep only rna counts for genes where model exists
#rna_test_df$V1 <- sub("\\..*", "", toupper(rna_test_df$V1))
rna_test_df_filtered <- rna_test_df[which(toupper(rna_test_df$V1)%in%toupper(common_genes)),] #18594 #19226 -> 18918 when duplicates removed.. which gene is learned??
#rna_test_df_filtered <- rna_test_df[which(rna_test_df$V1%in%rownames(corr_df)),]

rna_mat <- as.matrix(rna_test_df_filtered[, -1])
nonzero_test_rna <- rowSums(rna_mat == 0) #changed to sparsity 
#names(nonzero_test_rna) <- rna_test_df_filtered$V1
rna_df <- data.frame(rna_test_df_filtered$V1, nonzero_test_rna)
colnames(rna_df) <- c('Gene_name', 'nonzero_count')
rna_df$nonzero_fraction <- rna_df$nonzero_count/ncol(rna_mat) #fraction of samples

#corr_df_filtered <- corr_df[-which(is.na(corr_df$Spearman)),]

feature_df <- merge(gene_df_filtered, corr_df_filtered, by = 'Gene_name', all = FALSE) #18339
rna_feature_df <- merge(gene_df_filtered, rna_df,  by = 'Gene_name', all = FALSE)
rna_feature_df <- merge(rna_feature_df, corr_df_filtered,  by = 'Gene_name', all = FALSE) ##18339

#thresholds single-cell 
feature_df <- feature_df %>%
  mutate(corr_class = case_when(
    Spearman >= 0.5 ~ "high",
    Spearman >= 0.1 ~ "medium",
    TRUE ~ "failed"
))

rna_feature_df <- rna_feature_df %>%
  mutate(corr_class = case_when(
    Spearman >= 0.5 ~ "high",
    Spearman >= 0.1 ~ "medium",
    TRUE ~ "failed"
  ))


#thresholds metacell
# feature_df <- feature_df %>%
#   mutate(corr_class = case_when(
#     Spearman >= 0.8 ~ "high",
#     Spearman >= 0.5 ~ "medium",
#     Spearman >= 0 & Spearman < 0.5 ~ "low",
#     TRUE ~ "negative"
# ))

feature_df$corr_class <- factor(feature_df$corr_class, levels = c("high", "medium", "failed"))
rna_feature_df$corr_class <- factor(feature_df$corr_class, levels = c("high", "medium", "failed"))

length(feature_df$corr_class[which(feature_df$corr_class == "high")]) #sc:  159, mc100: 7106
length(feature_df$corr_class[which(feature_df$corr_class == "medium")]) #sc: 1970, mc100: 5536
length(feature_df$corr_class[which(feature_df$corr_class == "failed")]) #sc: 16219, mc100: 6111 

df_nonzero_rna <- rna_feature_df %>% #18339
  dplyr::select('Ensembl ID', corr_class, 'nonzero_fraction') %>%
  dplyr::rename(value = 'nonzero_fraction') %>%
  dplyr::mutate(feature = 'nonzero_fraction')

df_xtss <- feature_df %>% #18339
  dplyr::select('Ensembl ID', corr_class, '#TSS') %>%
  dplyr::rename(value = '#TSS') %>%
  dplyr::mutate(feature = '#TSS')
#zoom
#df_xtss <- df_xtss[df_xtss$value <= 50,]

df_xtranscripts <- feature_df %>%
  dplyr::select('Ensembl ID', corr_class, '#Transcripts') %>%
  dplyr::rename(value = '#Transcripts') %>%
  dplyr::mutate(feature = '#Transcripts')

df_glength <- feature_df %>% #18339
  dplyr::select('Ensembl ID', corr_class, 'Gene_length') %>%
  dplyr::rename(value = 'Gene_length') %>%
  dplyr::mutate(feature = 'Gene_length')
#df_glength <- df_glength[df_glength$value <= 50000,]

df_elength <- feature_df %>% #18339
  dplyr::select('Ensembl ID', corr_class, 'Exons_length') %>%
  dplyr::rename(value = 'Exons_length') %>%
  dplyr::mutate(feature = 'Exons_length')

df_gdensity <- feature_df %>%
  dplyr::select('Ensembl ID', corr_class, 'Gene_density') %>%
  dplyr::rename(value = 'Gene_density') %>%
  dplyr::mutate(feature = 'Gene_density')
#df_gdensity <- df_gdensity[df_gdensity$value <= 100,]

df_ubiquitness <- feature_df %>%
  dplyr::select('Ensembl ID', corr_class, 'Expression_ubiquitousness') %>%
  dplyr::rename(value = 'Expression_ubiquitousness') %>%
  dplyr::mutate(feature = 'Expression_ubiquitousness')

df_plot <- bind_rows(df_xtss, df_xtranscripts, df_glength, df_elength, df_gdensity, df_nonzero_rna)

xtss_failed <- df_xtss$value[df_xtss$corr_class == "failed"]
xtss_medium <- df_xtss$value[df_xtss$corr_class == "medium"]
xtss_high <- df_xtss$value[df_xtss$corr_class == "high"]

glength_failed <- df_glength$value[df_glength$corr_class == "failed"]
glength_medium <- df_glength$value[df_glength$corr_class == "medium"]
glength_high <- df_glength$value[df_glength$corr_class == "high"]

nonzero_failed <- df_nonzero_rna$value[df_nonzero_rna$corr_class == "failed"]
nonzero_medium <- df_nonzero_rna$value[df_nonzero_rna$corr_class == "medium"]
nonzero_high <- df_nonzero_rna$value[df_nonzero_rna$corr_class == "high"]

xtranscripts_failed <- df_xtranscripts$value[df_xtranscripts$corr_class == "failed"]
xtranscripts_medium <- df_xtranscripts$value[df_xtranscripts$corr_class == "medium"]
xtranscripts_high <- df_xtranscripts$value[df_xtranscripts$corr_class == "high"]

elength_failed <- df_elength$value[df_elength$corr_class == "failed"]
elength_medium <- df_elength$value[df_elength$corr_class == "medium"]
elength_high <- df_elength$value[df_elength$corr_class == "high"]

gdensity_failed <- df_gdensity$value[df_gdensity$corr_class == "failed"]
gdensity_medium <- df_gdensity$value[df_gdensity$corr_class == "medium"]
gdensity_high <- df_gdensity$value[df_gdensity$corr_class == "high"]

unpaired.wilcox.fun <- function(gene_set1, gene_set2){
  wilcox <- wilcox.test(gene_set1,gene_set2, paired = FALSE, alternative = "two.sided")
  return(wilcox)
}

#xtss
xtss_failed_medium <- unpaired.wilcox.fun(xtss_failed, xtss_medium)
xtss_medium_high <- unpaired.wilcox.fun(xtss_medium, xtss_high)

#glength
glength_failed_medium <- unpaired.wilcox.fun(glength_failed, glength_medium)
glength_medium_high <- unpaired.wilcox.fun(glength_medium, glength_high)

#nonzero test rna
nonzero_failed_medium <- unpaired.wilcox.fun(nonzero_failed, nonzero_medium)
nonzero_medium_high <- unpaired.wilcox.fun(nonzero_medium, nonzero_high)

#elength
elength_failed_medium <- unpaired.wilcox.fun(elength_failed, elength_medium)
elength_medium_high <- unpaired.wilcox.fun(elength_medium, elength_high)

#xtranscripts
xtranscripts_failed_medium <- unpaired.wilcox.fun(xtranscripts_failed, xtranscripts_medium)
xtranscripts_medium_high <- unpaired.wilcox.fun(xtranscripts_medium, xtranscripts_high)

#gdensity
gdensity_failed_medium <- unpaired.wilcox.fun(gdensity_failed, gdensity_medium)
gdensity_medium_high <- unpaired.wilcox.fun(gdensity_medium, gdensity_high)

#pie chart biotypes genes?

plot_violin_box <- function(data, feature) {
  ggplot(data %>% filter(feature == !!feature), aes(x = corr_class, y = value, fill = corr_class)) +
    geom_violin(trim = FALSE, scale = "width", alpha = 0.7) +
    #geom_boxplot(width = 0.1, outlier.size = 1, outlier.shape = 21, outlier.fill = "white") +
    geom_boxplot(width = 0.1, outlier.shape = NA) +
    labs(#title = paste("Distribution of gene features stratified by Pearson correlation"),
      #subtitle = "high: ≥ 0.8, medium: ≥ 0.5 and < 0.8, low: ≥ 0 and < 0.5, negative: < 0",
      axis.title = element_text(size = 17),  # Adjust axis labels size
      axis.text = element_text(size = 15), 
      x = "correlation class",
      y = feature) +
    theme_minimal() +
    theme(legend.position = "none", axis.text.x = element_text(size = 17), axis.text.y = element_text(size = 17),
          axis.title=element_text(size=20))  +
    scale_fill_brewer(palette = "Set3") #+ ylim(0, 1) #+
  #ggtitle(paste("Feature:", feature))
}

df_xtss <- df_xtss[df_xtss$value <= 20,] 
df_xtranscripts <- df_xtranscripts[df_xtranscripts$value <= 20,]
df_glength <- df_glength[df_glength$value <= 40000,] 
df_elength <- df_elength[df_elength$value <= 10000,] 
df_gdensity <- df_gdensity[df_gdensity$value <= 80,] 
df_plot <- bind_rows(df_xtss, df_xtranscripts, df_glength, df_elength, df_gdensity, df_nonzero_rna)
#features <- c("#TSS", "#Transcripts", "Gene_length", "Exons_length", "Gene_density", "Expression_ubiquitousness", "nonzero_fraction")
#features <- c("#TSS", "#Transcripts", "Gene_density", "Gene_length")
#features <- c('nonzero_fraction')
features <- c("#TSS", "#Transcripts", "Gene_length", "Exons_length", "Gene_density", "nonzero_fraction")
plots <- lapply(features, function(feature) plot_violin_box(df_plot, feature))

#pdf("/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/sc_1MB500bp_corr_classes_01_05.pdf")
pdf("/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/mc100_1MB500bp_corr_classes_01_05.pdf")
plots
dev.off()

#grid_plots <- do.call(grid.arrange, c(plots, ncol = 2))
grid_plots <- grid.arrange(grobs = plots, ncol=2, top = textGrob( "Distribution of gene features stratified by Spearman correlation", gp=gpar(fontsize=13,font=8)))

### model performance ###

# get_density <- function(x, y, ...) {
#   dens <- MASS::kde2d(x, y, ...)
#   ix <- findInterval(x, dens$x)
#   iy <- findInterval(y, dens$y)
#   ii <- cbind(ix, iy)
#   return(dens$z[ii])
# }

# corr_df$density <- get_density(corr_df$train, corr_df$Pearson, n = 100)
# dens <-  MASS::kde2d(corr_df$train, corr_df$Pearson)

#train vs test correlation 
scatter_corr <- ggplot(corr_df, aes(x = train, y = Pearson)) +
  #geom_point(aes(color = density), size = 3) +
  geom_pointdensity() +
  scale_color_viridis_c(name = "density") + 
  geom_abline(intercept = 0.00, slope = 1, linetype = "dashed", color = "black") +
  labs(x = "Pearson train", y = "Pearson test") +
  ggtitle("Model performance: Pearson correlation") +
  theme_minimal(base_size = 12) +
  #coord_equal() +
  theme(plot.title = element_text(face = "bold", size = 15),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 12), aspect.ratio = 1) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01))

#train vs test error 
scatter_err <- ggplot(err_df, aes(x = train_err, y = test_err)) +
  #geom_point(aes(color = density), size = 3) +
  geom_pointdensity() +
  scale_color_viridis_c(name = "density") + 
  geom_abline(intercept = 0.00, slope = 1, linetype = "dashed", color = "black") +
  labs(x = "MSE train", y = "MSE test") +
  ggtitle("Model performance: MSE") +
  theme_minimal(base_size = 12) +
  #coord_equal() +
  theme(plot.title = element_text(face = "bold", size = 15),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 12), aspect.ratio = 1) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1))

#p <- ggplot(err_df, aes(x = train_err, y = test_err)) 
#p <- p + geom_point() + geom_abline(intercept = 0,slope=1, colour='#E41A1C') + ggtitle("Binned RF on simulated single-cell data: train vs. test error") +
#  xlab("train RMSE") + ylab("test RMSE") 

# test MSE and Pearson correlation

performance_df <- merge(corr_df, err_df, by = 'row.names', all = TRUE) 

#boxplot 

boxplot_pearson <- ggplot(performance_df, aes(x = "", y = Pearson)) +
  geom_boxplot(fill = viridis_pal()(1), alpha = 0.5) +
  labs(title = "Test Pearson correlation", x = NULL, y = "test Pearson") +
  theme_minimal()

boxplot_test_err <- ggplot(performance_df, aes(x = "", y = test_err)) +
  geom_boxplot(fill = viridis_pal()(1), alpha = 0.5) +
  labs(title = "Test MSE", x = NULL, y = "test MSE") +
  theme_minimal()

grid.arrange(boxplot_pearson, boxplot_test_err, ncol = 4)

### gene biotypes ###

#biotype_counts_gtf <- table(feature_df$Biotype_gtf)
feature_df <- feature_df[feature_df$Spearman > 0.1,]
biotype_counts <- table(feature_df$Biotype_general)
biotype_counts <- biotype_counts[-7]
biotype_df <- data.frame(biotype = names(biotype_counts), count = as.numeric(biotype_counts))
biotype_df$percentage <- biotype_df$count / sum(biotype_df$count) * 100

custom_labels <- biotype_df %>%
  mutate(label = paste0(
    case_when(
      biotype == "antisense" ~ "antisense",
      biotype == "processed_transcript" ~ "processed transcript",
      biotype == "non-coding RNA" ~ "non-coding RNA",
      biotype == "protein-coding" ~ "protein-coding",
      biotype == "pseudogene" ~ "pseudogene",
      biotype == "sense_intronic" ~ "sense intronic",
      #biotype == "sense_overlapping" ~ "sense overlapping",
      TRUE ~ "Other"
    ),
    " (", round(percentage, 1), "%)"
  )) %>%
  select(biotype, label) %>%
  deframe()

# custom_labels <- c(
#   "antisense" = "antisense",
#   "non-coding RNA" = "non-coding RNA",
#   "processed_transcript" = "processed transcript",
#   "protein-coding" = "protein-coding",
#   "pseudogene" = "pseudogene",
#   "sense_intronic" = "sense intronic",
#   "sense_overlapping" = "sense overlapping"
# )

pie_biotypes <- ggplot(biotype_df, aes(x = "", y = count, fill = biotype)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +  
  scale_fill_manual(values = viridis_pal()(length(custom_labels)), labels = custom_labels) +  # Use custom labels
  #geom_text(aes(label = paste0(round(percentage, 1), "%")), 
  #          position = position_stack(vjust = 0.5), color = "white", size = 4) + 
  #scale_fill_viridis_d(option="magma") +  
  labs(
    title = "Gene model biotypes",
    fill = "Biotypes",
    x = NULL,
    y = NULL
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 15, face = "bold"),
    axis.text = element_blank(),  # Remove axis text
    axis.ticks = element_blank(),  # Remove axis ticks
    panel.grid = element_blank(),  # Remove panel grid lines
    legend.position = "right",   # Position legend to the right
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12)
  )

############################################################

#stratify performance by gene features: two gene sets: aggregtated single-cell or meta-cell performance better

#TODO: 500bp bins, all features
#corr_df_mc100_q85_no_scaling <- read.table("/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/lsi_investigation/performance_reports/mc100_1MB500bp_binned_lsi_activity_q0_no_threshold_performance_correlation.txt",
#                                           sep="\t",header=TRUE)
#corr_df_aggr_sc_q85_no_scaling_test_all <- read.table("/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/scarlink/sc_1MB100bp_q85_not_scaled_aggregated_sc_to_mc_spearman_test_cells.tsv",
#                                                      sep="\t",header=TRUE)
#corr_df_aggr_sc_q85_no_scaling_test_all <- read.table("/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/lsi_investigation/performance_reports/sc_1MB500bp_q0_aggregated_sc_to_mc_spearman_test_cells.tsv", sep="\t",header=TRUE)
corr_df_mc100_q100_500bp <- read.table("/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/lsi_investigation/performance_reports/mc100_1MB500bp_binned_lsi_activity_q0_no_threshold_performance_correlation.txt",sep="\t",header=TRUE)
corr_df_aggr_sc_q100_500bp <- read.table("/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/lsi_investigation/performance_reports/sc_1MB500bp_q0_aggregated_sc_to_mc_spearman_test_cells.tsv", sep="\t",header=TRUE)

gene_df <- fread("/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/GeneFeatureMat.txt.gz")
#TODO: train?
rna_test_df <- fread("/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/rna/10x_pbmc_rna_counts_test_cells.tsv.gz", sep="\t")
gene_names <- rna_test_df[,1]

#corr_df_sc <- corr_df_aggr_sc_q85_no_scaling_test_all
corr_df_sc <- corr_df_aggr_sc_q100_500bp 
rownames(corr_df_sc) <- corr_df_sc[,1]
corr_df_sc <- corr_df_sc[, -1]
colnames(corr_df_sc) <- c("sc_Spearman", "sc_p_value")
#corr_df_mc <- corr_df_mc100_q85_no_scaling
corr_df_mc <- corr_df_mc100_q100_500bp
colnames(corr_df_mc) <- c("mc_train", "mc_Pearson", "mc_Spearman")
common_gene_models <- intersect(rownames(corr_df_mc), rownames(corr_df_sc))
corr_df <- cbind(corr_df_mc[which(rownames(corr_df_mc)%in%common_gene_models),], corr_df_sc[which(rownames(corr_df_sc)%in%common_gene_models),]) #18722
common_genes <- intersect(toupper(rownames(corr_df)), toupper(gene_df$'Gene_name')) #18706 #18665
# TODO: debug
# genes that are present in performance file but not input genes (from count matrix)? extension
#ext_genes <- setdiff(toupper(common_genes), toupper(rna_test_df$V1))
#common_genes_filtered <- setdiff(toupper(rna_test_df$V1), toupper(common_genes)) #TODO: remove when gene extensions are saved
corr_df_filtered <- corr_df[which(toupper(rownames(corr_df))%in%toupper(common_genes)),] 
#corr_df_filtered <- corr_df_filtered[-which(toupper(rownames(corr_df_filtered))%in%ext_genes),] #18594 #TODO: remove
corr_df_filtered$`Gene_name` <- toupper(rownames(corr_df_filtered)) #18706

gene_df_filtered <- gene_df[which(toupper(gene_df$`Gene_name`)%in%toupper(common_genes)),] #18733
#gene_df_filtered$Gene_name_without_ext <- sub("\\..*", "", toupper(gene_df_filtered$`Gene_name`))
#gene_df_filtered <- gene_df_filtered[-which(toupper(gene_df_filtered$`Gene_name_without_ext`)%in%ext_genes),] #18611 #19142 #TODO: remove

#TODO: proper gene name storing with extension
#keep only rna counts for genes where model exists
#rna_test_df$V1 <- sub("\\..*", "", toupper(rna_test_df$V1))
rna_test_df_filtered <- rna_test_df[which(toupper(rna_test_df$V1)%in%toupper(common_genes)),] #18706 #18665
#rna_test_df_filtered <- rna_test_df[which(rna_test_df$V1%in%rownames(corr_df)),]

#single-cell sparsity
rna_mat <- as.matrix(rna_test_df_filtered[, -1])
sparsity_test_rna <- rowSums(rna_mat == 0)
#names(nonzero_test_rna) <- rna_test_df_filtered$V1
rna_df <- data.frame(rna_test_df_filtered$V1, sparsity_test_rna)
colnames(rna_df) <- c('Gene_name', 'sparsity')
rna_df$sparsity <- rna_df$sparsity/ncol(rna_mat) #fraction of samples

gene_df_filtered$Gene_name <- toupper(gene_df_filtered$Gene_name)
rna_df$Gene_name <- toupper(rna_df$Gene_name)
corr_df_filtered$Gene_name <- toupper(corr_df_filtered$Gene_name)

feature_df <- merge(gene_df_filtered, rna_df, by = 'Gene_name', all = FALSE) #18733
feature_df <- merge(feature_df, corr_df_filtered,  by = 'Gene_name', all = FALSE) #18458 #TODO: which genes lost?
#remove duplicates
dup_genes <- feature_df$Gene_name[duplicated(feature_df$Gene_name)]
feature_df_filtered <- feature_df[-which(feature_df$Gene_name %in% dup_genes), ] #18639
#for superior sc performance: only chose gene models with corr > 0.3, min difference between mc and sc correlation = 0.1

#filtering: keep all genes where sc performance (sc_Spearman) has min 0.3 performance 
#sc label: genes where sc Spearman > mc Spearman and min difference between mc and sc correlation = 0.1

#feature_df <- feature_df %>%
#  mutate(corr_class = case_when(
#    mc_Spearman < sc_Spearman &  ~ "sc",
#    TRUE ~ "mc"
#  ))

#feature_df_filtered <- feature_df_filtered %>%
#  filter(sc_Spearman >= 0.3) %>%
#  mutate(
#    corr_class = case_when(
#      sc_Spearman - mc_Spearman >= 0.1 ~ "sc",
#      TRUE ~ "mc"
#    )
#)

feature_df_filtered <- feature_df_filtered %>%
  filter(sc_Spearman >= 0.3 | mc_Spearman >= 0.3) %>%  
  mutate(
    corr_class = case_when(
      sc_Spearman - mc_Spearman >= 0.1 ~ "sc",    
      mc_Spearman - sc_Spearman >= 0.1 ~ "mc",    
      TRUE ~ "neutral"                            
    )
)

sc_genes <- feature_df_filtered[which(corr_class == "sc"), Gene_name]
mc_genes <- feature_df_filtered[which(corr_class == "mc"), Gene_name]
neutral_genes <- feature_df_filtered[which(corr_class == "neutral"), Gene_name]

######### sparsity investigation #########
sc_df <- feature_df_filtered[which(corr_class == "sc"),]
mc_df <- feature_df_filtered[which(corr_class == "mc"),]
all_df <- feature_df_filtered

sc_cor <- cor(sc_df$Exons_length, sc_df$sparsity, method = "spearman")
mc_cor <- cor(mc_df$Exons_length, mc_df$sparsity, method = "spearman")
all_cor_exons <- cor(all_df$Exons_length, all_df$sparsity, method = "spearman")
all_cor_tss <- cor(all_df$'#TSS', all_df$sparsity, method = "spearman")

#sc_df_filtered <- sc_df
sc_df_filtered <- sc_df[which(sc_df$Exons_length <= 10000), ] 
all_df_filtered <- all_df[which(all_df$Exons_length <= 10000), ] 
#sc_df_filtered <- sc_df[which(sc_df$sparsity > 0.9), ] 
sc_df_filtered$Exons_length <- scales::rescale(sc_df_filtered$Exons_length, to = c(0,1))
all_df_filtered$Exons_length <- scales::rescale(all_df_filtered$Exons_length, to = c(0,1))
#all_df_filtered$'#TSS' <- scales::rescale(all_df_filtered$'#TSS', to = c(0,1))

pdf("/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/gene_features_correlation_sparsity_exons_length_number_tss.pdf")
print(ggplot(all_df_filtered, aes(x = Exons_length, y = sparsity)) +
  #geom_density_2d_filled(contour_var = "ndensity", alpha = 0.6) +
  geom_point(color = "blue", alpha = 0.4, size = 1.2) +
  geom_density_2d(color = "black", size = 0.7, alpha = 0.8) +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  annotate("text",
           x = min(all_df_filtered$Exon_length, na.rm = TRUE),
           y = max(all_df_filtered$sparsity, na.rm = TRUE),
           label = paste0("N = ", round(all_cor_exons, 3)),
           hjust = 0, vjust = 1, size = 5) +
  # styling
  theme_minimal(base_size = 14) +
  labs(
    x = "Exon length",
    y = "Sparsity",
  ))
print(ggplot(all_df_filtered, aes(x = factor(`#TSS`), y = sparsity)) +
        #geom_density_2d_filled(contour_var = "ndensity", alpha = 0.6) +
        geom_point(color = "blue", alpha = 0.4, size = 1.2) +
        #geom_jitter(color = "blue", alpha = 0.5, size = 1.5, width = 0.2) +
        #geom_density_2d(color = "black", size = 0.7, alpha = 0.8) +
        geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
        scale_x_discrete(breaks = c(0, 5, 10, 15 ,20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90)) +
        annotate("text",
                 x = min(all_df_filtered$'#TSS', na.rm = TRUE),
                 y = max(all_df_filtered$sparsity, na.rm = TRUE),
                 label = paste0("N = ", round(all_cor_tss, 3)),
                 hjust = 0, vjust = 1, size = 5) +
        # styling
        theme_minimal(base_size = 14) +
        labs(
          x = "Number TSSs",
          y = "Sparsity",
        )) #+ theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
##############################################

### gene sets g:profiler ###
sc_gene_df <- data.frame(superior_genes = sc_genes)
write.table(sc_gene_df, file = "/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/genes_aggr_sc_outperforms_mc_spearman_q0_500bp_min_03_sc_superior_01.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
mc_gene_df <- data.frame(superior_genes = mc_genes)
write.table(mc_gene_df, file = "/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/genes_mc_outperforms_aggr_sc_spearman_q0_500bp_min_03_mc_superior_01.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

### clean mc gene set: unmapped genes
library(biomaRt)
mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

gene_map <- getBM(
  attributes = c("hgnc_symbol"),
  filters = "hgnc_symbol",
  values = mc_genes,
  mart = mart
)

mapped_genes <- unique(gene_map$hgnc_symbol)
unmapped_genes <- setdiff(mc_genes, mapped_genes)
mc_gene_df <- data.frame(superior_genes = mapped_genes)
write.table(mc_gene_df, file = "/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/genes_mc_outperforms_aggr_sc_spearman_q0_500bp_min_03_mc_superior_01_mapped.tsv", sep = "\t", quote = FALSE, row.names = FALSE,col.names = FALSE)
#unmapped_genes = 192  

### 

#remove all neutral genes for analysis
feature_df_filtered <- feature_df_filtered[-which(corr_class == "neutral"),]
feature_df_filtered <- as.data.frame(feature_df_filtered)

df_sparsity <- feature_df_filtered %>% #18339
  dplyr::select('Ensembl ID', corr_class, 'sparsity') %>%
  dplyr::rename(value = 'sparsity') %>%
  dplyr::mutate(feature = 'sparsity')

df_xtss <- feature_df_filtered %>% #18339
  dplyr::select('Ensembl ID', corr_class, '#TSS') %>%
  dplyr::rename(value = '#TSS') %>%
  dplyr::mutate(feature = '#TSS')

df_xtranscripts <- feature_df_filtered %>%
  dplyr::select('Ensembl ID', corr_class, '#Transcripts') %>%
  dplyr::rename(value = '#Transcripts') %>%
  dplyr::mutate(feature = '#Transcripts')

df_glength <- feature_df_filtered %>% #18339
  dplyr::select('Ensembl ID', corr_class, 'Gene_length') %>%
  dplyr::rename(value = 'Gene_length') %>%
  dplyr::mutate(feature = 'Gene_length')

df_elength <- feature_df_filtered %>% #18339
  dplyr::select('Ensembl ID', corr_class, 'Exons_length') %>%
  dplyr::rename(value = 'Exons_length') %>%
  dplyr::mutate(feature = 'Exons_length')

df_gdensity <- feature_df_filtered %>%
  dplyr::select('Ensembl ID', corr_class, 'Gene_density') %>%
  dplyr::rename(value = 'Gene_density') %>%
  dplyr::mutate(feature = 'Gene_density')
#df_gdensity <- df_gdensity[df_gdensity$value <= 100,]

xtss_sc <- df_xtss$value[df_xtss$corr_class == "sc"]
xtss_mc <- df_xtss$value[df_xtss$corr_class == "mc"]

xtranscripts_sc <- df_xtranscripts$value[df_xtranscripts$corr_class == "sc"]
xtranscripts_mc <- df_xtranscripts$value[df_xtranscripts$corr_class == "mc"]

glength_sc <- df_glength$value[df_glength$corr_class == "sc"]
glength_mc <- df_glength$value[df_glength$corr_class == "mc"]

elength_sc <- df_elength$value[df_elength$corr_class == "sc"]
elength_mc <- df_elength$value[df_elength$corr_class == "mc"]

gdensity_sc <- df_gdensity$value[df_gdensity$corr_class == "sc"]
gdensity_mc <- df_gdensity$value[df_gdensity$corr_class == "mc"]

sparsity_sc <- df_sparsity$value[df_sparsity$corr_class == "sc"]
sparsity_mc <- df_sparsity$value[df_sparsity$corr_class == "mc"]

unpaired.wilcox.fun <- function(gene_set1, gene_set2){
  wilcox <- wilcox.test(gene_set1,gene_set2, paired = FALSE, alternative = "two.sided")
  return(wilcox)
}

#xtss
xtss_pval <- unpaired.wilcox.fun(xtss_sc, xtss_mc)
print(xtss_pval)
xtranscripts_pval <- unpaired.wilcox.fun(xtranscripts_sc, xtranscripts_mc)
print(xtranscripts_pval)
glength_pval <- unpaired.wilcox.fun(glength_sc, glength_mc)
print(glength_pval)
elength_pval <- unpaired.wilcox.fun(elength_sc, elength_mc)
print(elength_pval)
gdensity_pval <- unpaired.wilcox.fun(gdensity_sc, gdensity_mc)
print(gdensity_pval)
sparsity_pval <- unpaired.wilcox.fun(sparsity_sc, sparsity_mc)
print(sparsity_pval)

#zoom for plotting
df_xtss <- df_xtss[df_xtss$value <= 20,] #18401
df_xtranscripts <- df_xtranscripts[df_xtranscripts$value <= 20,] #18309
df_glength <- df_glength[df_glength$value <= 40000,] #11780 #length(df_glength$value[df_glength$corr_class == "sc"]): 2847
df_elength <- df_elength[df_elength$value <= 10000,] #18450
df_gdensity <- df_gdensity[df_gdensity$value <= 80,] 

df_plot <- bind_rows(df_xtss, df_xtranscripts, df_glength, df_elength, df_gdensity, df_sparsity)

df_plot <- df_plot %>%
  mutate(method = factor(corr_class, levels = c("sc", "mc"))) 

plot_violin_box <- function(data, feature) {
  ggplot(data %>% filter(feature == !!feature), aes(x = corr_class, y = value, fill = corr_class)) +
    geom_violin(trim = FALSE, scale = "width", alpha = 0.7) +
    geom_boxplot(width = 0.1, outlier.size = 1, outlier.shape = 21, outlier.fill = "white") +
    labs(#title = paste("Distribution of gene features stratified by Pearson correlation"),
      #subtitle = "high: ≥ 0.8, medium: ≥ 0.5 and < 0.8, low: ≥ 0 and < 0.5, negative: < 0",
      axis.title = element_text(size = 17),  # Adjust axis labels size
      axis.text = element_text(size = 15), 
      x = "correlation class",
      y = feature) +
    theme_minimal() +
    theme(legend.position = "none", axis.text.x = element_text(size = 17), axis.text.y = element_text(size = 17),
          axis.title=element_text(size=20))  +
    scale_fill_manual(values = c(
      "sc" = "#66c2a5",
      "mc" = "#8da0cb"
    ))
    #scale_fill_brewer(palette = "Set3") #+ ylim(0, 1) #+
  #ggtitle(paste("Feature:", feature))
}

features <- c("#TSS", "#Transcripts", "Gene_length", "Exons_length", "Gene_density", "sparsity")
plots <- lapply(features, function(feature) plot_violin_box(df_plot, feature))

pdf("/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/aggretated_sc_vs_mc100_1MB500bp_gene_features_01_difference_min_03_corr.pdf")
plots 
dev.off()

##
# investigate genes that have superior sc performance compared to mc100 (aggregated predictions)

pred_df_mc100_q100 <- read.table("/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/10xpmbc_metafr_gex_pred_cpm_all_features_1MB500bp_mc100_sc_superior_genes_100.tsv", sep="\t", header=TRUE)
corr_df_mc100_q0_500bp <-read.table("/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/lsi_investigation/performance_reports/cpm_corrected/mc100_1MB500bp_binned_lsi_activity_q0_no_threshold_performance_correlation.txt", sep="\t", header=TRUE)
pred_df_aggr_sc_q100 <- read.table("/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/sc_1MB500bp_q0_aggregated_sc_to_mc_predictions_test_cells_superior_genes_100_randomly_selected.tsv", sep="\t",header=TRUE) #test single-cells, all meta-cells
corr_df_aggr_sc_q00_500bp <- read.table("/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/sc_1MB500bp_q0_aggregated_sc_to_mc_spearman_test_cells_superior_genes_100_randomly_selected.tsv", sep="\t",header=TRUE) 

#add model performances
rownames(corr_df_aggr_sc_q00_500bp) <- corr_df_aggr_sc_q00_500bp[,1]

#comparison predictions only on test meta-cells
rownames(pred_df_mc100_q100) <- pred_df_mc100_q100[,1]
rownames(pred_df_aggr_sc_q100) <- pred_df_aggr_sc_q100[,1]
mc_ids <- intersect(rownames(pred_df_mc100_q100),rownames(pred_df_aggr_sc_q100))

pred_df_aggr_sc_q100_test <- pred_df_aggr_sc_q100[which(rownames(pred_df_aggr_sc_q100)%in%mc_ids),]

pred_df_mc100_q100_test <- pred_df_mc100_q100[mc_ids, ]
pred_df_aggr_sc_q100_test <- pred_df_aggr_sc_q100[mc_ids, ]

pdf("/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/predictions_aggretated_sc_mc100_on_genes_sc_superior_20_genes_test_mcs.pdf", width = 6, height = 6)
for (gene in colnames(pred_df_mc100_q100_test[2:21])) {
  df_plot <- data.frame(
    mc_pred = pred_df_mc100_q100_test[[gene]],
    sc_pred = pred_df_aggr_sc_q100_test[[gene]]
  )
  
  x_range <- range(df_plot$mc_pred, df_plot$sc_pred, na.rm = TRUE)
  y_range <- x_range
  
  corr_mc100 <- corr_df_mc100_q0_500bp[gene, "Spearman"]
  corr_sc    <- corr_df_aggr_sc_q00_500bp[gene, "Spearman"]
  
  subtitle_txt <- paste0("mc100 corr = ", round(corr_mc100, 3), 
                         ", aggregated sc corr = ", round(corr_sc, 3))
  
  p <- ggplot(df_plot, aes(x = mc_pred, y = sc_pred)) +
    geom_point(color = "blue", alpha = 0.6) +
    geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
    coord_fixed(xlim = x_range, ylim = y_range) +
    labs(
      title = paste("Prediction on test meta-cells for", gene),
      subtitle = subtitle_txt,
      x = "mc100",
      y = "aggregated single-cell"
    ) +
    theme_minimal(base_size = 14)
  
  print(p)  
}

dev.off()
##
#TODO: investigate genes with negative correlation metacell setups

