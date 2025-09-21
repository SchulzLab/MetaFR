library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)

#Manuscript Figure 2: comparison meta-cell and aggregated single-cell performance
#eisvogel_palette <- c("#ed8028","#6ca5c3","#0e28b2")

corr_df_mc100_500bp_q100 <- read.table("/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/lsi_investigation/performance_reports/cpm_corrected/mc100_1MB500bp_binned_lsi_activity_q0_no_threshold_performance_correlation.txt", sep="\t",header=TRUE)
corr_df_aggr_sc_500bp_q100 <- read.table("/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/lsi_investigation/performance_reports/cpm_corrected/sc_1MB500bp_q0_aggregated_sc_to_mc_spearman_test_cells.tsv", sep="\t",header=TRUE)

##### barplot high correlation models ####

#barplot high correlation models
rownames(corr_df_aggr_sc_500bp_q100) <- corr_df_aggr_sc_500bp_q100[,1]
corr_df_aggr_sc_500bp_q100 <- corr_df_aggr_sc_500bp_q100[,-1]
#colnames(corr_df_aggr_sc_500bp_q100) <- c("spearman_sc", "p_value")
corr_df_aggr_sc_500bp_q100 <- corr_df_aggr_sc_500bp_q100[,1, drop=FALSE]

#colnames(corr_df_mc100_500bp_q100) <- c("mc100_train", "mc100_Pearson", "mc100_Spearman")
corr_df_mc100_500bp_q100 <- corr_df_mc100_500bp_q100[,3, drop=FALSE]

corr_df <- bind_rows(
  corr_df_aggr_sc_500bp_q100 %>% mutate(method = "sc"),
  corr_df_mc100_500bp_q100  %>% mutate(method = "mc100")
) %>% 
  filter(Spearman > 0.5) %>%
  arrange(desc(Spearman))

corr_count_df <- corr_df %>%
  group_by(method) %>%
  summarise(
    count = sum(Spearman > 0.5),
    .groups = 'drop'
  ) %>%
  arrange(desc(count))

corr_count_df <- corr_count_df %>%
  mutate(method = factor(method, levels = c("sc", "mc100"))) 

#mc100 = 7217, sc = 4774
pdf("/Users/laurarumpf/Documents/single_cell_stitchit/manuscript/figures/barplot_test_spearman_aggr_sc_mc100_cutoff_05_500bp.pdf", height = 5.5, width = 3.5)
print(ggplot(corr_count_df, aes(x = method, y = count, fill = method)) +
        geom_bar(stat = "identity", width = 0.7) +
        #geom_text( aes(label = count), color = "white", size = 5, position = position_stack(vjust = 0.5)) +
        scale_fill_manual(values = c(
          "sc" = "#66c2a5",
          "mc100" = "#8da0cb"
        )) +
        theme_minimal() +
        theme(
          panel.grid = element_blank(),
          panel.background = element_blank(),
          plot.background = element_blank(),
          axis.line = element_line(color = "black"),
          legend.position = "none",
          axis.ticks = element_line(color = "black"),
        )) 
dev.off()

colnames(corr_df_mc100_500bp_q100) <- c("mc_Spearman")
colnames(corr_df_aggr_sc_500bp_q100) <- c("sc_Spearman")
common_genes <- intersect(rownames(corr_df_aggr_sc_500bp_q100), rownames(corr_df_mc100_500bp_q100)) #18679

corr_df_mc100_500bp_q100_filtered <- corr_df_mc100_500bp_q100[which(rownames(corr_df_mc100_500bp_q100) %in% common_genes), ,drop=FALSE]
corr_df_aggr_sc_500bp_q100_filtered <- corr_df_aggr_sc_500bp_q100[which(rownames(corr_df_aggr_sc_500bp_q100) %in% common_genes), ,drop=FALSE]

corr_df_superior <- cbind(corr_df_mc100_500bp_q100_filtered, corr_df_aggr_sc_500bp_q100_filtered)

corr_df_superior_filtered  <- corr_df_superior  %>%
  filter(sc_Spearman >= 0.3 | mc_Spearman >= 0.3) %>%  
  mutate(
    corr_class = case_when(
      sc_Spearman - mc_Spearman >= 0.1 ~ "sc",    
      mc_Spearman - sc_Spearman >= 0.1 ~ "mc100",    
      TRUE ~ "neutral"                            
  )
)

corr_count_df_superior <- as.data.frame(table(corr_df_superior_filtered$corr_class))
colnames(corr_count_df_superior) <- c("corr_class", "count")
corr_count_df_superior <- corr_count_df_superior[-2,]

corr_count_df_superior <- corr_count_df_superior %>%
  mutate(corr_class = factor(corr_class, levels = c("sc", "mc100"))) 

pdf("/Users/laurarumpf/Documents/single_cell_stitchit/manuscript/figures/barplot_test_spearman_aggr_sc_mc100_500bp_superior_genes_cutoff_03_min_diff_01.pdf", height = 5, width = 4)
print(ggplot(corr_count_df, aes(x = method, y = count, fill = method)) +
        geom_bar(stat = "identity", width = 0.7) +
        #geom_text( aes(label = count), color = "white", size = 5, position = position_stack(vjust = 0.5)) +
        scale_fill_manual(values = c(
          "sc" = "#66c2a5",
          "mc100" = "#8da0cb"
        )) +
        theme_minimal() +
        theme(
          panel.grid = element_blank(),
          panel.background = element_blank(),
          plot.background = element_blank(),
          axis.line = element_line(color = "black"),
          legend.position = "none",
          axis.ticks = element_line(color = "black"),
        )) 
dev.off()

####### mc model performance comparison ######

corr_df_mc100_q0_no_scaling <- read.table("/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/lsi_investigation/performance_reports/cpm_corrected/mc100_1MB500bp_binned_lsi_activity_q0_no_threshold_performance_correlation.txt", sep="\t",header=TRUE)
corr_df_mc50_q0_no_scaling <- read.table("/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/lsi_investigation/performance_reports/cpm_corrected/mc50_1MB500bp_binned_lsi_activity_q0_no_threshold_performance_correlation.txt", sep="\t",header=TRUE)
corr_df_mc150_q0_no_scaling <- read.table("/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/lsi_investigation/performance_reports/cpm_corrected/mc150_1MB500bp_binned_lsi_activity_q0_no_threshold_performance_correlation.txt", sep="\t",header=TRUE)

#remove all models with test correlation NA
corr_df_mc100_q0_no_scaling <- corr_df_mc100_q0_no_scaling[!is.na(corr_df_mc100_q0_no_scaling$Spearman),]  
corr_df_mc50_q0_no_scaling <- corr_df_mc50_q0_no_scaling[!is.na(corr_df_mc50_q0_no_scaling$Spearman),] 
corr_df_mc150_q0_no_scaling <- corr_df_mc150_q0_no_scaling[!is.na(corr_df_mc150_q0_no_scaling$Spearman),] 

#keep models with 0 correlation
corr_df_mc100_q0_no_scaling <- corr_df_mc100_q0_no_scaling[-which(corr_df_mc100_q0_no_scaling$Spearman<0),] 
corr_df_mc50_q0_no_scaling <- corr_df_mc50_q0_no_scaling[-which(corr_df_mc50_q0_no_scaling$Spearman<0),] 
corr_df_mc150_q0_no_scaling <- corr_df_mc150_q0_no_scaling[-which(corr_df_mc150_q0_no_scaling$Spearman<0),] 

common_genes <- Reduce(intersect, list(
  #rownames(corr_df_aggr_sc_q85_no_scaling_test_all), 
  rownames(corr_df_mc100_q0_no_scaling), 
  rownames(corr_df_mc50_q0_no_scaling),
  rownames(corr_df_mc150_q0_no_scaling))) 

keep_common_genes <- function(df){
  df_filtered <- df[common_genes, ,drop=FALSE]
  return(df_filtered)
}

corr_df_mc100_q0_no_scaling_aligned <- keep_common_genes(corr_df_mc100_q0_no_scaling)
corr_df_mc50_q0_no_scaling_aligned <- keep_common_genes(corr_df_mc50_q0_no_scaling)
corr_df_mc150_q0_no_scaling_aligned <- keep_common_genes(corr_df_mc150_q0_no_scaling)

corr_df_mc50_q0_no_scaling_aligned$method <- "mc_50"
corr_df_mc100_q0_no_scaling_aligned$method <- "mc_100"
corr_df_mc150_q0_no_scaling_aligned$method <- "mc_150"

cor_df_all_methods <- data.frame(mc50=corr_df_mc50_q0_no_scaling_aligned$Spearman, 
                                 mc100=corr_df_mc100_q0_no_scaling_aligned$Spearman,
                                 mc150=corr_df_mc150_q0_no_scaling_aligned$Spearman)

cor_df_long <- cor_df_all_methods %>%
  pivot_longer(
    cols = everything(),
    names_to = "method",
    values_to = "corr"
  ) %>%
  mutate(
    method = factor(method, 
                    levels = c("mc50", "mc100", "mc150"))
  )

p.violin <- ggplot(cor_df_long, aes(x = method, y = corr, fill = method)) +
  geom_violin(
    width = 0.9,
    alpha = 0.7,
    trim = TRUE,
    scale = "width"
  ) +
  geom_boxplot(
    width = 0.15,
    alpha = 0.9,
    #outlier.shape = NA,
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = c(
      #"aggregated_sc_100" = "#66c2a5",
      "mc50" = "#fc8d62",
      "mc100" = "#8da0cb",
      "mc150" = "#e78ac3"
    ),
    name = "method"
  ) +
  labs(
    x = "MetaFR setup",
    y = "Test Spearman correlation"
    #caption = "Violin plots show density distribution, boxplots show quartiles"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none",  # Remove legend since x-axis labels are clear
    plot.margin = margin(1, 1, 1, 1, "cm")
  ) 

paired.wilcox.fun <- function(mc1, mc2, corr_df){
  wilcox <- wilcox.test(
    corr_df[[mc1]],
    corr_df[[mc2]],
    paired = TRUE,
    alternative = "less"
  )
  return(wilcox)
}

wilcox_mc50_mc100 <- paired.wilcox.fun("mc50", "mc100", cor_df_all_methods)
wilcox_mc100_mc150 <- paired.wilcox.fun("mc100", "mc150", cor_df_all_methods)

pdf("/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/mc_metafr_setupts_test_spearman_violin_1MB500bp_NA_removed.pdf")
p.violin
dev.off()

####### varying bin sizes #######
corr_df_mc100_100bp <- read.table("/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/lsi_investigation/performance_reports/mc100_1MB100bp_binned_lsi_activity_q100_no_threshold_test_partition_no_scaling_performance_correlation.txt",
                                  sep="\t",header=TRUE)
corr_df_mc100_250bp <- read.table("/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/lsi_investigation/performance_reports/mc100_1MB250bp_binned_lsi_activity_q0_no_threshold_performance_correlation.txt",
                                  sep="\t",header=TRUE)
corr_df_mc100_500bp <- read.table("/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/lsi_investigation/performance_reports/mc100_1MB500bp_binned_lsi_activity_q0_no_threshold_performance_correlation.txt",
                                  sep="\t",header=TRUE)

corr_df_mc100_100bp <- corr_df_mc100_100bp[!is.na(corr_df_mc100_100bp$Spearman),]  
corr_df_mc100_250bp  <- corr_df_mc100_250bp [!is.na(corr_df_mc100_250bp$Spearman),] 
corr_df_mc100_500bp <- corr_df_mc100_500bp[!is.na(corr_df_mc100_500bp$Spearman),] 

#keep models with 0 correlation
corr_df_mc100_100bp <- corr_df_mc100_100bp[-which(corr_df_mc100_100bp$Spearman<0),] 
corr_df_mc100_250bp <- corr_df_mc100_250bp[-which(corr_df_mc100_250bp$Spearman<0),] 
corr_df_mc100_500bp <- corr_df_mc100_500bp[-which(corr_df_mc100_500bp$Spearman<0),] 

common_genes <- Reduce(intersect, list(
  rownames(corr_df_mc100_100bp), 
  rownames(corr_df_mc100_250bp),
  rownames(corr_df_mc100_500bp))) 

keep_common_genes <- function(df){
  df_filtered <- df[common_genes, ,drop=FALSE]
  return(df_filtered)
}

corr_df_mc100_100bp_aligned <- keep_common_genes(corr_df_mc100_100bp)
corr_df_mc100_250bp_aligned <- keep_common_genes(corr_df_mc100_250bp)
corr_df_mc100_500bp_aligned <- keep_common_genes(corr_df_mc100_500bp)

corr_df_mc100_100bp_aligned$method <- "mc100_100bp"
corr_df_mc100_250bp_aligned$method <- "mc100_250bp"
corr_df_mc100_500bp_aligned$method <- "mc100_500bp"

cor_df_all_methods <- data.frame(mc100_100bp=corr_df_mc100_100bp_aligned$Spearman, 
                                 mc100_250bp=corr_df_mc100_250bp_aligned$Spearman, 
                                 mc100_500bp=corr_df_mc100_500bp_aligned$Spearman)

cor_df_long <- cor_df_all_methods %>%
  pivot_longer(
    cols = everything(),
    names_to = "method",
    values_to = "corr"
  ) %>%
  mutate(
    method = factor(method, 
                    levels = c("mc100_100bp", "mc100_250bp", "mc100_500bp"))
  )

b.violin <- ggplot(cor_df_long, aes(x = method, y = corr, fill = method)) +
  geom_violin(
    width = 0.9,
    alpha = 0.7,
    trim = TRUE,
    scale = "width"
  ) +
  geom_boxplot(
    width = 0.15,
    alpha = 0.9,
    #outlier.shape = NA,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = c("#a6cee3", "#1f78b4", "#08519c")
  ) +
  labs(
    x = "MetaFR setup",
    y = "Test Spearman correlation"
    #caption = "Violin plots show density distribution, boxplots show quartiles"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none",  # Remove legend since x-axis labels are clear
    plot.margin = margin(1, 1, 1, 1, "cm")
  ) 

paired.wilcox.fun <- function(mc1, mc2, corr_df){
  wilcox <- wilcox.test(
    corr_df[[mc1]],
    corr_df[[mc2]],
    paired = TRUE,
    alternative = "less"
  )
  return(wilcox)
}

wilcox_bp100_bp250 <- paired.wilcox.fun("mc100_100bp", "mc100_250bp", cor_df_all_methods)
wilcox_bp250_bp500 <- paired.wilcox.fun("mc100_250bp", "mc100_500bp", cor_df_all_methods)

pdf("/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/mc100_metafr_bin_size_test_spearman_violin_NA_filtered.pdf")
b.violin
dev.off()

######## aggr. sc vs mc test correlation ########

corr_df_mc100_500bp_q100 <- read.table("/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/lsi_investigation/performance_reports/cpm_corrected/mc100_1MB500bp_binned_lsi_activity_q0_no_threshold_performance_correlation.txt", sep="\t",header=TRUE)
corr_df_aggr_sc_500bp_q100 <- read.table("/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/lsi_investigation/performance_reports/cpm_corrected/sc_1MB500bp_q0_aggregated_sc_to_mc_spearman_test_cells.tsv", sep="\t",header=TRUE)

filter_df <- function(df) {
  #na_df <- df[which(is.na(df$Spearman)),]
  df_filtered <- df[!is.na(df$Spearman),]
  return(df_filtered)
} 

merge_dfs <- function(df1, df2){
  df1_x <- df1[which(rownames(df1)%in%rownames(df2)),]
  df2_x <- df2[which(rownames(df2)%in%rownames(df1)),]
  df1_2 <- cbind(df1_x, df2_x)
  if (ncol(df1_2) == 6){
    colnames(df1_2) <- c("train1", "Pearson1", "Spearman1", "train2", "Pearson2", "Spearman2")
  }else{
    colnames(df1_2) <- c("train1", "Pearson1", "Spearman1", "Spearman2", "Spearman2_pval")  
  }
  return(df1_2)
}

rownames(corr_df_aggr_sc_500bp_q100) <- corr_df_aggr_sc_500bp_q100[,1]
corr_df_aggr_sc_500bp_q100 <- corr_df_aggr_sc_500bp_q100[,-1]
corr_df_aggr_sc_500bp_q100 <- filter_df(corr_df_aggr_sc_500bp_q100)

corr_df_mc100_500bp_q100 <- filter_df(corr_df_mc100_500bp_q100)
corr_df_aggr_sc_mc100 <- merge_dfs(corr_df_mc100_500bp_q100, corr_df_aggr_sc_500bp_q100)

laufenten_palette <- c("#F6D1A7","#843624")
fangschreckenkrebs_palette <- c("#083508","#8ED755","#F8D918")
fangschreckenkrebs_palette_2 <- c("#083508","#578B7B","#F8D918")
schabrackentapir_palette <- c("#ACB1BF","#62D0E2","#1E2F76")

pdf("/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/scatter_aggregated_single_cell_metacell_spearman_metafr_q0_500bp.pdf")
#print(ggplot(corr_df_aggr_sc_q85_no_scaling_all_mc100_q85_no_scaling, aes(x = Spearman1, y = Spearman2)) + geom_point() + geom_density_2d() + geom_abline(intercept = 0,slope=1, colour='#E41A1C') +
#  xlab("test Spearman meta-cell") + ylab("test Spearman aggr. single-cell"))
print(ggplot(corr_df_aggr_sc_mc100, aes(x = Spearman1, y = Spearman2)) +
        geom_bin2d(bins = 50) +
        # scale_fill_viridis_c(option = "mako", direction = -1) +
        scale_fill_gradientn(colors = schabrackentapir_palette,
                             trans = "log10") +
        geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed", linewidth = 0.6) +
        #coord_cartesian(xlim = c(0,100), ylim = c(0,100)) +    
        coord_fixed(ratio = 1, xlim = c(0,1), ylim = c(0,1), expand = TRUE, clip = "on") +
        labs(
          x = "test Spearman meta-cell",
          y = "test Spearman single-cell",
        ) +
        
        theme_minimal() +
        theme(
          panel.grid = element_blank(),
          panel.background = element_blank(),
          plot.background = element_blank(),
          axis.line = element_line(color = "black"),
          axis.title = element_text(size = 11),
          axis.ticks = element_line(color = "black"),
        ))
dev.off()


########## runtime ###########

# initial gene set metafr: 22281 genes

# mc100: 25 cores
#elapsed: 330.785 for mc generation = ~5.51 min ()
#real:	0m16.019s for mc filtering
#real: 75m11.952s for activity binning
#1h 13m 52s for model learning = ~ 68,112 min 

#real    401m45.394s, user    4069m51.201s, sys     2839m8.431s for activity binning
# 13h 5m 57s for model learning

#single-cell vs ms
runtime_df <- data.frame(method=c("sc", "mc100"), seconds=c(3, 0.42))
runtime_df <- runtime_df %>%
  mutate(method = factor(method, levels = c("sc", "mc100"))) 
#runtime_df <- data.frame(method=c("metafr", "scarlink"), minutes=c(0.05, 4.3))
pdf("/Users/laurarumpf/Documents/single_cell_stitchit/manuscript/figures/barplot_runtime_sc_mc100_500bp.pdf",  height = 5.5, width = 3.5)
print(ggplot(runtime_df, aes(x = method, y = seconds, fill = method)) +
        geom_bar(stat = "identity", width = 0.7) +
        #geom_text( aes(label = count), color = "white", size = 5, position = position_stack(vjust = 0.5)) +
        scale_fill_manual(values = c(
          "sc" = "#66c2a5",
          "mc100" = "#8da0cb"
        )) +
        theme_minimal() +
        theme(
          panel.grid = element_blank(),
          panel.background = element_blank(),
          plot.background = element_blank(),
          axis.line = element_line(color = "black"),
          legend.position = "none",
          axis.ticks = element_line(color = "black"),
        )) 
dev.off()

#varying bin size

  # mc100: 25 cores
  #elapsed: 330.785 for mc generation = ~5.51 min ()
  #real:	0m16.019s for mc filtering
  #real: 339m2.880s for activity binning
  #1h 13m 52s for model learning = ~ 68,112 min 
  runtime_df_mc100 <- data.frame(method=c("100bp", "500bp"), seconds=c(,0.42))
  runtime_df_sc <- data.frame(method=c("100bp", "500bp"), seconds=c(,3))


