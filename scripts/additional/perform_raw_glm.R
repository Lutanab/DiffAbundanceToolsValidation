library(phyloseq)
library(R6)

project_dir <- "."
source(paste0(project_dir, "/scripts/metrics/getMetrics.R"))

da_tools_dir <- "../DA_tools"
source(paste0(da_tools_dir, "/diff_abundance_tools.R"))

parser <- argparse::ArgumentParser()
parser$add_argument("-p", "--phyloseq", type="character", action="store", help="Filepath to synthetic phyloseq object")
parser$add_argument("-t", "--transform", type="character", action="store", help="Transform of counts before fitting lin reg")
parser$add_argument("-fd", "--da_table_filename", type="character", action="store", help="Filepath for characteristic table for diff abundant taxa")
parser$add_argument("-o", "--output_filepath", type="character", action="store", help="Filepath to output")
parser$add_argument("-pt", "--p_val_treshold", type="numeric", action="store", help="p_val treshold when result is considered to be right")
parser$add_argument("-lr", "--lfc_range", type="numeric", action="store", help="range of delta lfc in which prediction is considered to be true")
args <- parser$parse_args()

sim_ps <- readRDS(args$phyloseq)
da_table <- utils::read.csv(args$da_table_filename, header = TRUE, sep = "\t")

design_formula = ~ group
ref_values <- list("group" = "reference")
da_result <- prepare_raw_glm_feature_table(sim_ps, design_formula, ref_values,
                                           count_transform=args$transform, USE_ADJUSTED_PVAL=TRUE)
da_result$pval_to_stars()

p_val_treshold <- args$p_val_treshold
lfc_range <- args$lfc_range
feature <- "group_target_vs_reference"
count_transform <- paste0("raw_glm_", args$transform)
metrics <- data.frame()

tmp_da_tab <- data.frame()
tmp_nda_tab <- data.frame()
da_taxa <- row.names(da_table)
nda_taxa <- setdiff(da_result$taxa, da_taxa)

tmp_da_tab[da_taxa, "lfc_true"] <- da_table[da_taxa, "lfc"]
tmp_da_tab[da_taxa, "lfc_est"] <- da_result$lfc_table[da_taxa, feature]
tmp_da_tab[da_taxa, "p_val"] <- da_result$pval_table[da_taxa, feature]

tmp_nda_tab[nda_taxa, "lfc_true"] <- 0
tmp_nda_tab[nda_taxa, "lfc_est"] <- da_result$lfc_table[nda_taxa, feature]
tmp_nda_tab[nda_taxa, "p_val"] <- da_result$pval_table[nda_taxa, feature]

# Range metrics calculate
tmp_df <- get_range_metrics_tab(tmp_da_tab, "DA", count_transform, lfc_range=lfc_range, p_val_treshold=p_val_treshold)
rownames(tmp_df) <- nrow(metrics)
metrics <- rbind(metrics, tmp_df)
tmp_df <- get_range_metrics_tab(tmp_nda_tab, "NDA", count_transform, lfc_range=lfc_range, p_val_treshold=p_val_treshold)
rownames(tmp_df) <- nrow(metrics)
metrics <- rbind(metrics, tmp_df)

# Unilateral metrics calculate
tmp_df <- get_unilateral_metrics_tab(tmp_da_tab, "DA", count_transform, p_val_treshold=p_val_treshold)
rownames(tmp_df) <- nrow(metrics)
metrics <- rbind(metrics, tmp_df)

write.table(metrics, file=args$output_filepath, quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE)

print("SUCCESSFUL PERFORMED RAW GLM DA ANALYSIS")