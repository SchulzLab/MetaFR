library("gprofiler2")
library(dplyr)
library(ontologyIndex)

#TODO: set background genes 

#go enrichment analysis metafr
sc_file <- "/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/genes_aggr_sc_outperforms_mc_spearman_q0_500bp_min_03_sc_superior_01.tsv"
sc_gene_df <- read.table(sc_file, sep = "\t", header = TRUE)
mc_file = "/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/genes_mc_outperforms_aggr_sc_spearman_q0_500bp_min_03_mc_superior_01.tsv"
mc_gene_df <- read.table(mc_file, sep = "\t", header = TRUE)
bg_file <- "/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/go_enrichment/initial_gene_set_no_filtering_no_ext.txt"
bg_gene_df <- read.table(bg_file, sep = "\t", header = FALSE)

mc_genes <- mc_gene_df$superior_genes
sc_genes <- sc_gene_df$superior_genes
bg_genes <- bg_gene_df$V1

mc_res <- gost(query = mc_genes, organism = "hsapiens", custom_bg = bg_genes, correction_method = "fdr", sources = c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC", "HP", "TF"))
#gostplot(mc_res)
gostplot(mc_res, capped = TRUE, interactive = TRUE)

keep_cols <- c("term_id", "term_name", "source", "p_value",
               "term_size", "intersection_size", "precision", "recall")

mc_res_flat <- mc_res$result[,keep_cols]
write.table(mc_res_flat, file = "/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/go_enrichment/go_enrichment_mc100_spearman_q0_500bp_all_terms.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

#head(mc_res$result) 
top_mc_res <- mc_res$result %>%
  filter(source == "GO:BP") %>%
  arrange(p_value) %>%
  head(100) 

## filter high-level terms accroding to depth in onthology tree

go <- get_ontology("http://purl.obolibrary.org/obo/go.obo", extract_tags = "minimal")

mc_term_depths <- sapply(top_mc_res$term_id, function(x) {
  if (x %in% names(go$ancestors)) {
    length(go$ancestors[[x]])
  } else {
    NA
  }
})

top_mc_res$depth <- mc_term_depths

top_mc_res_specific <- top_mc_res %>%
  filter(depth >= 7)

##
sc_res <- gost(query = sc_genes, organism = "hsapiens", custom_bg = bg_genes,  correction_method = "fdr", sources = c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC", "HP", "TF"))
#gostplot(mc_res)
top_sc_res <- sc_res$result %>%
  #filter(source == "GO:BP") %>%
  arrange(p_value) %>%
  head(100) 
gostplot(sc_res, capped = TRUE, interactive = TRUE)

#check for terms that are enriched in both gene sets
  #number of terms = 0 for BP in top 100 enriched terms
t <- intersect(top_mc_res$term_id, top_sc_res$term_id)

keep_cols_mc <- c("term_id", "term_name", "source", "p_value",
                  "term_size", "intersection_size", "precision", "recall", "depth")
top_mc_res_specific_flat <- top_mc_res_specific[,keep_cols_mc]
write.table(top_mc_res_specific_flat, file = "/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/go_enrichment/go_enrichment_mc100_spearman_q0_500bp_specific_terms_top_100_min_depth_7.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

keep_cols_sc <- c("term_id", "term_name", "source", "p_value",
                  "term_size", "intersection_size", "precision", "recall")
top_sc_res_flat <- top_sc_res[,keep_cols_sc]
write.table(top_sc_res_flat, file = "/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/go_enrichment/go_enrichment_aggr_sc_spearman_q0_500bp_all_terms.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


#### plotting terms of interest as dotplot ####
top_mc_res_specific_flat <- read.table("/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/go_enrichment/go_enrichment_mc100_spearman_q0_500bp_specific_terms_top_100_min_depth_7.tsv",sep="\t", header=TRUE)
top_sc_res_flat <-read.table( "/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/go_enrichment/go_enrichment_sggr_sc_spearman_q0_500bp_all_terms.tsv",sep="\t", header=TRUE)

term_ids <- c("GO:0046649", "GO:0002253", "GO:0002694", "GO:0001817","GO:0002764","GO:0098590", "GO:0046872", "GO:0043169", "GO:0043005", "GO:0022857")

mc_plot_df <- top_mc_res_specific_flat %>%
  mutate(method = "mc100",
         gene_fraction = intersection_size / term_size) 

sc_plot_df <- top_sc_res_flat %>%
  mutate(method = "sc",
         gene_fraction = intersection_size / term_size) 

plot_df <- bind_rows(mc_plot_df, sc_plot_df) %>%
  filter(term_id %in% term_ids)

break_vals <- sort(unique(round(plot_df$gene_fraction,3)))[c(1,5,9)]
pdf("/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/go_enrichment/mc100_aggr_sc_superior_gene_sets_go_enrichment_5_terms_selected.pdf", height=7)
go_p <- ggplot(plot_df, aes(x = method,
                    y = reorder(term_name, -log10(p_value)),
                    size = gene_fraction,
                    color = -log10(p_value))) +
  geom_point() +
  #scale_size_continuous(range = c(3, 10)) +  
  scale_size_continuous(
    #range = c(3,10),  # controls actual dot sizes in both plot and legend
    #breaks = c(0.2132, 0.4877, 0.5431),
    breaks =  c(0.213, 0.229, 0.476, 0.543),
    #breaks = break_vals,
    guide = guide_legend(title = "Gene fraction")
  ) +
  #scale_color_gradient(low = "blue", high = "red") +
  scale_color_gradient(
    low = "red", high = "blue",  # red = low p, blue = high p
    #limits = c(0, 0.05),
    name = "Adjusted p-value"
  ) +
  theme_bw() +
  labs(x = "Gene set",
       y = "GO term",
       size = "Gene fraction",
       color = "-log10(FDR)") +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1)) 
print(go_p)
dev.off()

#############

library(readr)
library(kableExtra)

# Read TSV
go_terms_mc100 <- read_tsv("/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/go_enrichment/go_enrichment_mc100_spearman_q0_500bp_all_terms.tsv")

# Save as PDF table
pdf("/Users/laurarumpf/Documents/single_cell_stitchit/binned_RF/binned_RF/10x_PBMC_multiome/go_enrichment/GO_terms_all_mc100.pdf", width = 8, height = 12) # adjust size
kable(go_terms_mc100, booktabs = TRUE) %>%
  kable_styling(latex_options = c("scale_down")) # auto-scale to page
dev.off()
