library(ggplot2)
library(cowplot)

project_dir <- "./"
source(paste0(project_dir, "/scripts/metrics/getMetrics.R"))

da_tools_dir <- "../projects/DA_tools"
source(paste0(da_tools_dir, "/diff_abundance_class.R"))
source(paste0(da_tools_dir, "/diff_abundance_tools.R"))


parser <- argparse::ArgumentParser()
# parser$add_argument("-ti", "--title", type="character", action="store", default="", help="Title of plot")
parser$add_argument("-d", "--da_result", type="character", action="store", default="compositional", help="Filepath of diff abundance analysis result")
parser$add_argument("-fd", "--da_table_filename", type="character", action="store", help="Filepath for characteristic table for diff abundant taxa")
parser$add_argument("-o", "--output_filepath", type="character", action="store", help="Filepath to output PCA plot")
parser$add_argument("-pt", "--p_val_treshold", type="numeric", action="store", help="p_val treshold when result is considered to be right")
parser$add_argument("-lr", "--lfc_range", type="numeric", action="store", help="range of delta lfc in which prediction is considered to be true")
args <- parser$parse_args()

if((args$p_val_treshold < 0) | (args$p_val_treshold > 1)){
  stop("p_val_treshold must be between 0 and 1")
}
if(args$lfc_range < 0){
  stop("lfc_range must be positive")
}


da_result <- readRDS(args$da_result)
da_table <- utils::read.csv(args$da_table_filename, header = TRUE, sep = "\t")

p_val_treshold <- args$p_val_treshold
lfc_range <- args$lfc_range
feature <- "group_target_vs_reference"

metrics <- data.frame()

da_taxa <- row.names(da_table)
nda_taxa <- setdiff(da_result$taxa, da_taxa)
for(tool in da_result$tools){
  tmp_da_tab <- data.frame()
  tmp_nda_tab <- data.frame()
  
  tmp_da_tab[da_taxa, "lfc_true"] <- da_table[da_taxa, "lfc"]
  tmp_da_tab[da_taxa, "lfc_est"] <- da_result[[tool]]$lfc_table[da_taxa, feature]
  tmp_da_tab[da_taxa, "p_val"] <- da_result[[tool]]$pval_table[da_taxa, feature]
  
  tmp_nda_tab[nda_taxa, "lfc_true"] <- 0
  tmp_nda_tab[nda_taxa, "lfc_est"] <- da_result[[tool]]$lfc_table[nda_taxa, feature]
  tmp_nda_tab[nda_taxa, "p_val"] <- da_result[[tool]]$pval_table[nda_taxa, feature]
  
  # Range metrics calculate
  tmp_df <- get_range_metrics_tab(tmp_da_tab, "DA", tool, lfc_range=lfc_range, p_val_treshold=p_val_treshold)
  rownames(tmp_df) <- nrow(metrics)
  metrics <- rbind(metrics, tmp_df)
  tmp_df <- get_range_metrics_tab(tmp_nda_tab, "NDA", tool, lfc_range=lfc_range, p_val_treshold=p_val_treshold)
  rownames(tmp_df) <- nrow(metrics)
  metrics <- rbind(metrics, tmp_df)
  
  # Unilateral metrics calculate
  tmp_df <- get_unilateral_metrics_tab(tmp_da_tab, "DA", tool, p_val_treshold=p_val_treshold)
  rownames(tmp_df) <- nrow(metrics)
  metrics <- rbind(metrics, tmp_df)
}

write.table(metrics, file=args$output_filepath, quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE)
print("SUCCESSFULL CALCULATED METRICS")





