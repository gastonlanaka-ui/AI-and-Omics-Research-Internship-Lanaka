

# Breast Cancer DEG GSEA Analysis

# Charger les données
load("C:/Users/dell/Documents/ai biotechno breast cancer project/DEG_2.RData")
#  Load required libraries
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(msigdbr)
library(dplyr)
library(tibble)
library(ggplot2)
library(fgsea)
library(cowplot)

# Load DEG data 
# Assumes deg_results is a data.frame with rownames = gene SYMBOLs
# and a "logFC" column

# Make SYMBOL column from rownames
deg_results <- deg_results %>% rownames_to_column("SYMBOL")

# Remove rows with missing SYMBOL or logFC
deg_results <- deg_results %>% filter(!is.na(SYMBOL) & SYMBOL != "" & !is.na(logFC))

# Map SYMBOL -> EntrezID 
mapped <- bitr(
  deg_results$SYMBOL,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# Keep only one ENTREZID per SYMBOL
mapped_unique <- mapped %>% group_by(SYMBOL) %>% slice(1) %>% ungroup()

# Join back to deg_results
deg_results <- deg_results %>% left_join(mapped_unique, by = "SYMBOL") %>% filter(!is.na(ENTREZID))

cat("Number of genes mapped to ENTREZID:", nrow(deg_results), "\n")

# Prepare ranked gene list
gene_ranks <- deg_results$logFC
names(gene_ranks) <- deg_results$ENTREZID
gene_ranks <- sort(gene_ranks, decreasing = TRUE)

# Load MSigDB gene sets 
# Hallmark (H)
hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, entrez_gene)

# Oncogenic (C6)
oncogenic <- msigdbr(species = "Homo sapiens", category = "C6") %>%
  dplyr::select(gs_name, entrez_gene)

#Run GSEA 
gsea_hallmark <- GSEA(geneList = gene_ranks, TERM2GENE = hallmark, eps = 0)
gsea_oncogenic <- GSEA(geneList = gene_ranks, TERM2GENE = oncogenic, eps = 0)

#  Create plot directory
plot_dir <- "C:/Users/dell/Documents/ai biotechno breast cancer project/enrichment_plots/"
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# Helper function to save plots
save_plot <- function(plot_obj, filename){
  ggsave(filename = file.path(plot_dir, filename),
         plot = plot_obj,
         width = 10, height = 8, dpi = 300)
}

# Plot GSEA results 
# Hallmark
if(inherits(gsea_hallmark, "gseaResult") && nrow(as.data.frame(gsea_hallmark)) > 0){
  save_plot(dotplot(gsea_hallmark, showCategory=20, title="GSEA - Hallmark"), "GSEA_Hallmark_dotplot.png")
  save_plot(gseaplot2(gsea_hallmark, geneSetID = 1, title=gsea_hallmark$Description[1]), "GSEA_Hallmark_top1.png")
}

# Oncogenic
if(inherits(gsea_oncogenic, "gseaResult") && nrow(as.data.frame(gsea_oncogenic)) > 0){
  save_plot(dotplot(gsea_oncogenic, showCategory=20, title="GSEA - Oncogenic"), "GSEA_Oncogenic_dotplot.png")
  save_plot(gseaplot2(gsea_oncogenic, geneSetID = 1, title=gsea_oncogenic$Description[1]), "GSEA_Oncogenic_top1.png")
}

cat("GSEA analysis complete. Plots saved in:", plot_dir, "\n")

