library(phyloseq)
library(microbiome)

parser <- argparse::ArgumentParser()
parser$add_argument("-p", "--phyloseq", type="character", action="store", help="Filepath to synthetic phyloseq object")
parser$add_argument("-t", "--transform", type="character", action="store", help="Transform of counts before fitting lin reg")
parser$add_argument("-fd", "--da_table_filename", type="character", action="store", help="Filepath for characteristic table for diff abundant taxa")
parser$add_argument("-o", "--output_filepath", type="character", action="store", help="Filepath to output")
parser$add_argument("-lr", "--lfc_range", type="numeric", action="store", help="range of delta lfc in which prediction is considered to be true")
args <- parser$parse_args()

lfc_range <- args$lfc_range
sim_ps <- readRDS(args$phyloseq)
da_table <- utils::read.csv(args$da_table_filename, header = TRUE, sep = "\t")
sim_ps <- microbiome::transform(sim_ps, args$transform)
real_lfc <- colMeans(sim_ps@otu_table[sim_ps@sam_data$group == "target",]) - colMeans(sim_ps@otu_table[sim_ps@sam_data$group == "reference",])

performance_table <- data.frame()
tmp_da_tab <- data.frame()
tmp_nda_tab <- data.frame()
da_taxa <- row.names(da_table)
nda_taxa <- setdiff(names(real_lfc), da_taxa)

tmp_da_tab[da_taxa, "lfc_true"] <- da_table[da_taxa, "lfc"]
tmp_da_tab[da_taxa, "lfc_est"] <- real_lfc[da_taxa]
performance_table[1, "approach"] <- "range_da"
performance_table[1, "correct"] <- (sum(abs(tmp_da_tab$lfc_est - tmp_da_tab$lfc_true) < (lfc_range/2), na.rm=TRUE)) / length(da_taxa)
performance_table[1, "incorrect"] <- 1 - performance_table[1, "correct"]

tmp_nda_tab[nda_taxa, "lfc_true"] <- 0
tmp_nda_tab[nda_taxa, "lfc_est"] <- real_lfc[nda_taxa]
performance_table[2, "approach"] <- "range_nda"
performance_table[2, "correct"] <- (sum(abs(tmp_nda_tab$lfc_est - tmp_nda_tab$lfc_true) < (lfc_range/2), na.rm=TRUE)) / length(nda_taxa)
performance_table[2, "incorrect"] <- 1 - performance_table[2, "correct"]

write.table(performance_table, file=args$output_filepath, quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE)

print("SUCCESSFUL PERFORMED MANUAL RAW GLM DA ANALYSIS")
