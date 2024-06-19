library(phyloseq)
library(argparse)
library(R6)

da_tools_dir <- "../DA_tools"
source(paste0(da_tools_dir, "/diff_abundance_class.R"))
source(paste0(da_tools_dir, "/diff_abundance_tools.R"))
source(paste0(da_tools_dir, "/heatmap_tools.R"))


parser <- argparse::ArgumentParser()
parser$add_argument("-ti", "--title", type="character", action="store", default="", help="Title of plot")
# parser$add_argument("-p", "--phyloseq", type="character", action="store", help="Filepath to synthetic phyloseq object")
parser$add_argument("-d", "--da_result", type="character", action="store", default="compositional", help="Filepath of diff abundance analysis result")
parser$add_argument("-fd", "--da_table_filename", type="character", action="store", help="Filepath for characteristic table for diff abundant taxa")
parser$add_argument("-o", "--output_filepath", type="character", action="store", help="Filepath to output PCA plot")
args <- parser$parse_args()


da_result <- readRDS(args$da_result)
da_table <- utils::read.csv(args$da_table_filename, header = TRUE, sep = "\t")
# sim_ps <- readRDS(args$phyloseq)

da_result$tools <- c(da_result$tools, "blank_line", "True_lfc")

tmp_lfc_tab <- data.frame()
tmp_lfc_tab[da_result$taxa,da_result$columns] <- 0
tmp_sig_tab <- data.frame()
tmp_sig_tab[da_result$taxa,da_result$columns] <- ""
da_result$blank_line <- DiffAbundanceComputationResult$from_significance_labels(
  lfc_table=tmp_lfc_tab, significance_labels=tmp_sig_tab
)

tmp_lfc_tab <- data.frame()
tmp_lfc_tab[da_result$taxa, da_result$columns] <- 0
tmp_lfc_tab[row.names(da_table), da_result$columns] <- da_table$mu_tar - da_table$mu_ref
tmp_sig_tab <- data.frame()
tmp_sig_tab[da_result$taxa,da_result$columns] <- "**\n**"
da_result$True_lfc <- DiffAbundanceComputationResult$from_significance_labels(
  lfc_table=tmp_lfc_tab, significance_labels=tmp_sig_tab
)


target_feature <- "group_target_vs_reference"
heatmap_title <- args$title
filename <- args$output_filepath

p <- plot_multiple_tools_on_one_feature_heatmap(sim_ps, da_result, target_feature, heatmap_title, filename, FULL_NAMES=FALSE, REMOVE_IRRELEVANT_TAXA=FALSE,
                                           cluster_taxa=TRUE, cutree_taxa=20)
print("SUCCESSFUL PLOT DA HEATMAP")
