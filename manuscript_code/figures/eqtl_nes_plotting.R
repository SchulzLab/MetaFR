library(ggplot2)
library(data.table)
library(dplyr)
library(stringr)


eqtl_files_1 <- c("/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/eqtl_analysis/mc50_nes_whole_blood_ext200bp_top100000.tsv",
              "/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/eqtl_analysis/mc150_nes_whole_blood_ext200bp_top100000.tsv",
             "/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/eqtl_analysis/mc100_nes_whole_blood_ext200bp_top100000.tsv",
             "/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/eqtl_analysis/sc_nes_whole_blood_ext200bp_top100000.tsv",
             "/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/eqtl_analysis/scarlink_nes_whole_blood_ext200bp_top100000.tsv")

eqtl_files_2 <- c("/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/eqtl_analysis/scarlink_nes_whole_blood_ext200bp_top100000_mean_no_sampling.tsv",
                  "/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/eqtl_analysis/mc100_nes_whole_blood_ext200bp_top100000_mean_no_sampling.tsv", 
                  "/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/eqtl_analysis/sc_nes_whole_blood_ext200bp_top100000_mean_no_sampling.tsv")


eqtl_files_3 <- c("/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/eqtl_analysis/sc_nes_whole_blood_ext200bp_top100000.tsv",
                  "/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/eqtl_analysis/sc_nes_whole_blood.tsv")

eqtl_files_4 <- c("/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/eqtl_analysis/mc100_nes_whole_blood_ext200bp_top100000.tsv",
                  "/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/eqtl_analysis/mc100_nes_whole_blood.tsv")

eqtl_files_5 <- c("/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/eqtl_analysis/mc150_nes_whole_blood_ext200bp_top100000.tsv",
                  "/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/eqtl_analysis/sc_nes_whole_blood_ext200bp_top100000.tsv")

eqtl_files_6 <- c("/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/eqtl_analysis/sc_nes_whole_blood_ext200bp_top100000_top2000genes.tsv",
                  "/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/eqtl_analysis/mc50_nes_whole_blood_ext200bp_top100000_top2000genes.tsv",
                  "/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/eqtl_analysis/mc100_nes_whole_blood_ext200bp_top100000_top2000genes.tsv")

eqtl_files_7 <- c("/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/eqtl_analysis/mc100_100bp_nes_whole_blood_top100000.tsv",
                  "/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/eqtl_analysis/mc100_250bp_nes_whole_blood_top100000.tsv",
                  "/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/eqtl_analysis/mc100_500bp_nes_whole_blood_top100000.tsv")


eqtl_files_8 <- c("/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/eqtl_analysis/scarlink_nes_whole_blood_ext200bp_top100000_mean_NA.tsv",
                  "/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/eqtl_analysis/sc_nes_whole_blood_ext200bp_top100000_ct_unscpecific.tsv")

eqtl_files_9 <- c("/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/shap_analysis/sc_abs_shap_nes_whole_blood_ext200bp_top100000.tsv",
                  "/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/shap_analysis/sc_abs_zscore_shap_nes_whole_blood_ext200bp_top100000.tsv")

eqtl_files_10 <- c("/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/eqtl_analysis/mc50_nes_whole_blood_ext200bp_top100000.tsv",
                  "/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/eqtl_analysis/mc150_nes_whole_blood_ext200bp_top100000.tsv",
                  "/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/eqtl_analysis/mc100_nes_whole_blood_ext200bp_top100000.tsv",
                  "/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/eqtl_analysis/sc_nes_whole_blood_ext200bp_top100000.tsv")

#TODO: look into negative NES celltypes
eqtl_files_11 <- c("/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/eqtl_analysis/500bp_bins/sc_nes_whole_blood_top100000_abs_shap.tsv",
                   "/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/eqtl_analysis/500bp_bins/mc50_nes_whole_blood_top100000_abs_shap.tsv",
                   "/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/eqtl_analysis/500bp_bins/mc100_nes_whole_blood_top100000_abs_shap.tsv",
                   "/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/eqtl_analysis/500bp_bins/mc150_nes_whole_blood_top100000_abs_shap.tsv")

eqtl_files_12 <- c("/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/eqtl_analysis/500bp_bins/sc_nes_whole_blood_top100000_abs_shap_distance.tsv",
                   "/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/eqtl_analysis/500bp_bins/mc50_nes_whole_blood_top100000_abs_shap_distance.tsv",
                   "/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/eqtl_analysis/500bp_bins/mc100_nes_whole_blood_top100000_abs_shap_distance.tsv",
                   "/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/eqtl_analysis/500bp_bins/mc150_nes_whole_blood_top100000_abs_shap_distance.tsv")

eqtl_files_13 <- c("/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/eqtl_analysis/sc_nes_whole_blood_top100000_abs_shap_distance_robust.tsv",
                   "/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/eqtl_analysis/mc50_nes_whole_blood_top100000_abs_shap_distance_robust.tsv",
                   "/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/eqtl_analysis/mc100_nes_whole_blood_top100000_abs_shap_distance_robust.tsv",
                   "/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/eqtl_analysis/mc150_nes_whole_blood_top100000_abs_shap_distance_robust.tsv")

eqtl_files_14 <- c("/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/eqtl_analysis/sc_nes_whole_blood_top100000_abs_shap_distance_iqr.tsv",
                   "/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/eqtl_analysis/mc50_nes_whole_blood_top100000_abs_shap_distance_iqr.tsv",
                   "/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/eqtl_analysis/mc100_nes_whole_blood_top100000_abs_shap_distance_iqr.tsv",
                   "/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/eqtl_analysis/mc150_nes_whole_blood_top100000_abs_shap_distance_iqr.tsv")

######

eqtl_files_15 <- c("/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/eqtl_analysis/top_1000_variable_genes/sc_nes_whole_blood_top100000_abs_shap_distance_iqr.tsv",
                   "/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/eqtl_analysis/top_1000_variable_genes/mc50_nes_whole_blood_top100000_abs_shap_distance_iqr.tsv",
                   "/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/eqtl_analysis/top_1000_variable_genes/mc100_nes_whole_blood_top100000_abs_shap_distance_iqr.tsv",
                   "/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/eqtl_analysis/top_1000_variable_genes/mc150_nes_whole_blood_top100000_abs_shap_distance_iqr.tsv")

eqtl_files_16 <- c("/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/eqtl_analysis/top_1000_variable_genes/sc_nes_whole_blood_top10000_abs_shap_distance_iqr.tsv",
                   "/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/eqtl_analysis/top_1000_variable_genes/mc50_nes_whole_blood_top10000_abs_shap_distance_iqr.tsv",
                   "/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/eqtl_analysis/top_1000_variable_genes/mc100_nes_whole_blood_top10000_abs_shap_distance_iqr.tsv",
                   "/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/eqtl_analysis/top_1000_variable_genes/mc150_nes_whole_blood_top10000_abs_shap_distance_iqr.tsv")

eqtl_files_17 <- c("/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/eqtl_analysis/top_1000_variable_genes/sc_nes_whole_blood_top10000_abs_shap_robust_zscore.tsv",
                   "/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/eqtl_analysis/top_1000_variable_genes/mc50_nes_whole_blood_top10000_abs_shap_robust_zscore.tsv",
                   "/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/eqtl_analysis/top_1000_variable_genes/mc100_nes_whole_blood_top10000_abs_shap_robust_zscore.tsv",
                   "/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/eqtl_analysis/top_1000_variable_genes/mc150_nes_whole_blood_top10000_abs_shap_robust_zscore.tsv")

eqtl_files_18 <- c("/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/eqtl_analysis/top_1000_variable_genes/sc_nes_whole_blood_top100000_abs_shap_robust_zscore.tsv",
                   "/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/eqtl_analysis/top_1000_variable_genes/mc50_nes_whole_blood_top100000_abs_shap_robust_zscore.tsv",
                   "/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/eqtl_analysis/top_1000_variable_genes/mc100_nes_whole_blood_top100000_abs_shap_robust_zscore.tsv",
                   "/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/eqtl_analysis/top_1000_variable_genes/mc150_nes_whole_blood_top100000_abs_shap_robust_zscore.tsv")

eqtl_df <- data.frame()
for (f in eqtl_files_15) {
  df <- fread(f, select = c("CellType", "Method", "NES", "FDR q-val", "Num_eQTL_hits", "Name"))
  df$Setup <- sub("_nes_whole_blood_top100000_abs_shap_distance_iqr.tsv", "", basename(f))
  eqtl_df <- rbind(eqtl_df, df)      
}

#eqtl_df$Setup <- sub("mc100_.*_(ext200bp).*", "mc100_\\1", eqtl_df$Setup)

training_ct_counts <- read.table("/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/capture_hi-c_analysis/training_celltype_counts.tsv", sep="\t", header=TRUE)
training_ct_counts$celltypes <- gsub(" ", "_", training_ct_counts$celltypes) 

methods <- unique(eqtl_df$Method)
sources <- unique(eqtl_df$Setup)
celltype_counts <- table(eqtl_df$CellType)
common_celltypes <- names(celltype_counts)[celltype_counts == length(methods) * length(sources)]
eqtl_counts <- table(eqtl_df$Num_eQTL_hits)

eqtl_df_filtered <- eqtl_df[eqtl_df$Name == "weighted_score" | eqtl_df$Name == "robust_z" ]
#eqtl_df_filtered <- eqtl_df[eqtl_df$CellType %in% common_celltypes, ] #filtered in NES calculation
eqtl_df_filtered <- eqtl_df[eqtl_df$Name == "weighted_score"]
#neg_nes_df <- eqtl_df_filtered[which(eqtl_df_filtered$NES < 0 & eqtl_df_filtered$Method == "DAP-G"),]

nes_df_ct_counts <- merge(
  eqtl_df_filtered, training_ct_counts, 
  by.x = "CellType", by.y = "celltypes", 
  all = FALSE
)

#length(unique(nes_df_ct_counts_filtered$CellType)) = 13
nes_df_ct_counts_filtered <- nes_df_ct_counts[-which(nes_df_ct_counts$count < 100 | nes_df_ct_counts$Num_eQTL_hits < 800),]
nes_df_ct_counts_filtered <- nes_df_ct_counts_filtered[nes_df_ct_counts_filtered$Name == "weighted_score"]
neg_nes_df <- nes_df_ct_counts_filtered[which(nes_df_ct_counts_filtered$NES < 0 & nes_df_ct_counts_filtered$Method == "DAP-G"),]
#nes_df_ct_counts_filtered <- nes_df_ct_counts[-which(nes_df_ct_counts$count < 100),]
#nes_df_ct_counts_filtered <- nes_df_ct_counts_filtered[-which(nes_df_ct_counts_filtered$CellType == "filtered"),]

plots <- list()

eqtl_df_filtered$Setup <- factor(eqtl_df_filtered$Setup,
                           levels = c("sc", "mc50", "mc100", "mc150"))

nes_df_ct_counts_filtered$Setup <- factor(nes_df_ct_counts_filtered$Setup,
                                          levels = c("sc", "mc50", "mc100", "mc150"))


nes_df_dapg <- nes_df_ct_counts_filtered[nes_df_ct_counts_filtered$Method == "DAP-G", ]

nes_wide <- dcast(
  data = nes_df_dapg,
  formula = CellType ~ Setup,  # One row per CellType
  value.var = "NES"
)

paired.wilcox.fun <- function(mc1, mc2, nes_df){
  wilcox <- wilcox.test(
    nes_df[[mc1]],
    nes_df[[mc2]],
    paired = TRUE,
    alternative = "less"
  )
  return(wilcox)
}

#wilcox_sc_mc50: p-value = 0.0009766
wilcox_sc_mc50 <- paired.wilcox.fun("sc", "mc50", nes_wide)
wilcox_mc50_mc100 <- paired.wilcox.fun("mc50", "mc100", nes_wide)
wilcox_mc100_mc150 <- paired.wilcox.fun("mc100", "mc150", nes_wide)

#eqtl_df$Setup[which(eqtl_df$Setup == "scarlink_nes_whole_blood_ext200bp_top100000_mean_NA.tsv")] <- "scarlink"
#eqtl_df$Setup[which(eqtl_df$Setup == "sc_nes_whole_blood_ext200bp_top100000_ct_unscpecific.tsv")] <- "metafr_sc_100bp"

setup_counts <- table(nes_df_dapg$CellType, nes_df_dapg$Setup)
valid_cts <- rownames(setup_counts)[rowSums(setup_counts > 0) == 4]
nes_df_dapg_filtered <- nes_df_dapg[nes_df_dapg$CellType %in% valid_cts, ]

for (method in methods) {
  #method_eqtls <- eqtl_df_filtered[eqtl_df_filtered$Method == method, ]
  #method_eqtls <- nes_df_ct_counts_filtered[nes_df_ct_counts_filtered$Method == method, ]
  method_eqtls <- nes_df_dapg_filtered[nes_df_dapg_filtered$Method == method, ]
  #method_eqtls <- eqtl_df[eqtl_df$Method == method, ]
  
  p <- ggplot(method_eqtls, aes(x = Setup, y = NES, fill = Setup)) +
    geom_boxplot(outlier.shape = NA, width = 0.7) +
    geom_jitter(width = 0.2, size = 1, alpha = 1) +
    #scale_fill_manual(values = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3","#a6d854")) +
    scale_fill_manual(values = c(
      "sc" = "#66c2a5",
      "mc50" = "#fc8d62",
      "mc100" = "#8da0cb",
      "mc150" = "#e78ac3"
    )) +
    labs(title = paste("NES Distribution for", method),
         x = "MetaFR setup", 
         y = "Normalized Enrichment Score (NES)") +
    theme_bw(base_size = 18, base_family = "Helvetica") +
    theme(
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white"),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.title = element_text(face = "bold"),
      legend.position = "none"
    )
  
  plots[[method]] <- p
}

###investigation negative NES values for single-cell setup

neg_sc_cts <- unique(neg_nes_df$CellType) 

neg_sc_nes_df <- nes_df_ct_counts %>%
  filter(CellType %in% neg_sc_cts)

pdf("/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/eqtl_analysis/iqr_100000_top_1000_genes_number_eqtl_hits.pdf")
print(ggplot(neg_sc_nes_df, aes(x = CellType, y = Num_eQTL_hits, color = Setup, group = Setup)) +
  geom_point(size = 4, alpha = 0.6) +
  #geom_line(aes(group = CellType), color = "grey50", size = 0.8, linetype = "dashed") +
  #geom_text(aes(label = round(NES, 2)), vjust = -1.2, size = 3.5, fontface = "bold") +
  scale_color_manual(values = c(
    "sc" = "#66c2a5",
    "mc50" = "#fc8d62",
    "mc100" = "#8da0cb",
    "mc150" = "#e78ac3"
  )) +
  labs(
    title = "Number of eQTL hits across setups",
    subtitle = "Cell types where sc has negative NES",
    x = "Cell Type",
    y = "Number of eQTL hits",
    color = "Setup"
  ) +
  geom_text(
    data = neg_sc_nes_df %>% filter(Setup == "sc"),
    #aes(label = round(NES, 2)),
    aes(x = CellType, y = Num_eQTL_hits, label = round(NES, 2)),
    vjust = 1.8,
    size = 3.5,
    fontface = "bold",  
    inherit.aes = FALSE) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12),
    legend.title = element_text(face = "bold")
  ) +
  ylim(0, max(neg_sc_nes_df$Num_eQTL_hits) * 1.1))
dev.off()

######################################################################################################

sc_shap_mono14_file <- "/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/shap_analysis/shap_CD14 Mono_interactions_all_genes.tsv"
shap_df_mono_14 <- fread(sc_shap_mono14_file, sep="\t")

p <- ggplot(shap_df_mono_14, aes(x = score)) +
  geom_histogram(binwidth = 0.00002, fill = "steelblue", color = "white") +
  labs(title = "SHAP values Mono14 single-cell MetaFR",
       x = "SHAP",
       y = "Frequency") +
  theme_minimal()




