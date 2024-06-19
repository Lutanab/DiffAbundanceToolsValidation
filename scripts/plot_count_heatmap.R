library(phyloseq)
library(microbiome)
library(pheatmap)
library(argparse)

da_tools_dir <- "../DA_tools"
source(paste0(da_tools_dir, "/common_meta_transform_class.R"))
source(paste0(da_tools_dir, "/heatmap_tools.R"))


parser <- argparse::ArgumentParser()
parser$add_argument("-ti", "--title", type="character", action="store", default="", help="Title of plot")
parser$add_argument("-p", "--phyloseq", type="character", action="store", help="Filepath to synthetic phyloseq object")
parser$add_argument("-fd", "--da_table_filename", type="character", action="store", help="Filepath for characteristic table for diff abundant taxa")
parser$add_argument("-o", "--output_filepath", type="character", action="store", help="Filepath to output PCA plot")
args <- parser$parse_args()

sim_ps <- readRDS(args$phyloseq)
da_table <- utils::read.csv(args$da_table_filename, header = TRUE, sep = "\t")

# ==============================================================================
CLUSSTERING_WITHIN_DIVIDING_FEATURE = FALSE

dividing_feature <- 'group'
features_of_interest <- c('test')
# ==============================================================================

sim_ps <- microbiome::transform(sim_ps, "clr")
heatmap_data <- prepare_taxa_abundance_heatmap_table(sim_ps, dividing_feature=dividing_feature, features_of_interest=features_of_interest,
                                                     FULL_FEATURE_NAMES = FALSE, CLUSSTERING_WITHIN_DIVIDING_FEATURE = CLUSSTERING_WITHIN_DIVIDING_FEATURE)
annotation_taxa <- data.frame()
annotation_taxa[taxa(sim_ps), "Status"] <- "Not diff abundant"
annotation_taxa[row.names(da_table), "Status"] <- "Diff abundant"
clustering_distance_rows <- 'euclidean'
# clustering_distance_rows <- 'binary'
# clustering_distance_rows <- 'correlation'
clustering_distance_cols <- 'euclidean'
# clustering_distance_cols <- 'binary'
# clustering_distance_cols <- 'correlation'
filename <- args$output_filepath
breaks <- sum(heatmap_data[['sample_meta']]$group == "reference")

p1 <- pheatmap::pheatmap(heatmap_data[['heat_matrix']], cluster_rows = T, treeheight_row  = 600, treeheight_col = 700,
                         cluster_cols = !CLUSSTERING_WITHIN_DIVIDING_FEATURE,
                         annotation_col = heatmap_data[['sample_meta']][, "group", drop = F],
                         annotation_row = annotation_taxa,
                         annotation_colors = list("Status" = c("Not diff abundant"="green", "Diff abundant"="red")),
                         main=args$title,
                         silent=F,
                         clustering_distance_rows=clustering_distance_rows,
                         clustering_distance_cols=clustering_distance_cols,
                         clustering_method='ward.D2',
                         gaps_col = breaks,
                         cutree_rows = ceiling(ntaxa(sim_ps)/30),
                         cutree_cols = ceiling(nsamples(sim_ps)/30),
                         show_rownames = T,
                         show_colnames = T,
                         filename=filename,
                         cellwidth = 10,
                         cellheight = 10,
                         border_color=NA
)

print("SUCCESSFUL COUNT HEATMAP")