library(microbiome)
library(phyloseq)
library(vegan)
library(argparse)
library(ggplot2)

parser <- argparse::ArgumentParser()
parser$add_argument("-t", "--transform", type="character", action="store", default="compositional", help="Transformations on counts (composition, clr etc)")
parser$add_argument("-ti", "--title", type="character", action="store", default="", help="Title of plot")
parser$add_argument("-d", "--distance", type="character", action="store", default="bray", help="Distance metrics between microbiomes")
parser$add_argument("-p", "--phyloseq", type="character", action="store", help="Filepath to synthetic phyloseq object")
parser$add_argument("-o", "--output_filepath", type="character", action="store", help="Filepath to output PCA plot")
args <- parser$parse_args()

sim_ps <- readRDS(args$phyloseq)

dist = phyloseq::distance(
  microbiome::transform(sim_ps, args$transform),
  method=args$distance)
ordination = phyloseq::ordinate(sim_ps, method="PCoA", distance=dist)

metadf <- data.frame(sample_data(sim_ps))
permanova_full <- vegan::adonis2(dist ~ group ,data = metadf)
R2 <- as.data.frame(permanova_full)["group","R2"]
p_val <- as.data.frame(permanova_full)["group","Pr(>F)"]

plot <- phyloseq::plot_ordination(sim_ps, ordination, color='group') + stat_ellipse() +
  labs(title = args$title,
    subtitle = paste0(
      "Transform: ", args$transform, "; Distance: ", args$distance,
      '\nR^2 = ', round(R2, digits=5), '; p_value = ', round(p_val, digits=5)))
ggsave(args$output_filepath, plot)

print("SUCCESSFUL PLOT PCA")