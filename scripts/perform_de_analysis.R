library(phyloseq)
library(R6)

da_tools_dir <- "../DA_tools"
source(paste0(da_tools_dir, "/diff_abundance_tools.R"))


parser <- argparse::ArgumentParser()
parser$add_argument("-p", "--phyloseq", type="character", action="store", help="Filepath to synthetic phyloseq object")
parser$add_argument("-o", "--output_filepath", type="character", action="store", help="Filepath to output PCA plot")
args <- parser$parse_args()

sim_ps <- readRDS(args$phyloseq)

design_formula = ~ group
ref_values <- list("group" = "reference")
common_results <- get_multiple_tools_feature_table(sim_ps, design_formula, ref_values)

common_results$DESeq2$pval_to_stars()
common_results$ALDEx2$pval_to_stars()
common_results$ANCOMBC$pval_to_stars()
common_results$glm_raw$pval_to_stars()

saveRDS(common_results, args$output_filepath)

print("SUCCESSFUL PERFORMED DA ANALYSIS")