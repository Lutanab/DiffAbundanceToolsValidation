library(ggplot2)
library(cowplot)

da_tools_dir <- "../DA_tools"
source(paste0(da_tools_dir, "/diff_abundance_class.R"))
source(paste0(da_tools_dir, "/diff_abundance_tools.R"))

parser <- argparse::ArgumentParser()
parser$add_argument("-ti", "--title", type="character", action="store", default="", help="Title of plot")
parser$add_argument("-d", "--da_result", type="character", action="store", default="compositional", help="Filepath of diff abundance analysis result")
parser$add_argument("-fd", "--da_table_filename", type="character", action="store", help="Filepath for characteristic table for diff abundant taxa")
# parser$add_argument("-fn", "--nda_table_filename", type="character", action="store", help="Filepath for characteristic table for non diff abundant taxa")
parser$add_argument("-o", "--output_filepath", type="character", action="store", help="Filepath to output PCA plot")
parser$add_argument("-cs", "--common_size", type="logical", action="store", default=FALSE, help="Use common xlim and ylim for all plots")
args <- parser$parse_args()


da_result <- readRDS(args$da_result)
da_table <- utils::read.csv(args$da_table_filename, header = TRUE, sep = "\t")


feature <- "group_target_vs_reference"
COMMON_SIZE <- args$common_size
# characteristic_tables <- list()
da_plot_list <- list()
nda_plot_list <- list()
pair_plot_list <- list()

if(COMMON_SIZE){
  min_x <- +Inf
  max_x <- -Inf
  min_y <- +Inf
  max_y <- -Inf
}


for(tool in da_result$tools){
  characteristic_table <- data.frame()
  characteristic_table[da_result$taxa, "Status"] <- "truly_not_diff_expressed"
  characteristic_table[da_result$taxa, "Delta_lfc"] <- da_result[[tool]]$lfc_table[da_result$taxa, feature]
  characteristic_table[da_result$taxa, "minus_log_pval"] <- -log10(da_result[[tool]]$pval_table[da_result$taxa, feature])
  
  characteristic_table[row.names(da_table), "Status"] <- "truly_diff_expressed"
  characteristic_table[row.names(da_table), "Delta_lfc"] <- da_table$lfc - characteristic_table[row.names(da_table), "Delta_lfc"]
  characteristic_table[row.names(da_table), "minus_log_pval"] <- -log10(da_result[[tool]]$pval_table[row.names(da_table), feature])
  
  characteristic_table[(characteristic_table == +Inf) | (characteristic_table == -Inf)] <- NA
  
  g1 <- ggplot2::ggplot(characteristic_table[characteristic_table$Status == "truly_not_diff_expressed",], aes(x=Delta_lfc, y=minus_log_pval), na.rm = TRUE) +
    ggplot2::geom_point(colour="black",pch=21, size=2, na.rm = TRUE) +
    ggplot2::labs(title = paste0(tool,": Not diff_expressed")) +
    ggplot2::theme(plot.title=element_text(size=15))
  g2 <- ggplot2::ggplot(characteristic_table[characteristic_table$Status == "truly_diff_expressed",], aes(x=Delta_lfc, y=minus_log_pval), na.rm = TRUE) +
    ggplot2::geom_point(colour="black",pch=21, size=2, na.rm = TRUE)+
    ggplot2::labs(title = paste0(tool,": Diff_expressed")) +
    ggplot2::theme(plot.title=element_text(size=15))
  
  if(COMMON_SIZE){
    min_x <- min(min_x, min(g1$data$Delta_lfc), min(g2$data$Delta_lfc), na.rm = TRUE)
    min_y <- min(min_y, min(g1$data$minus_log_pval), min(g2$data$minus_log_pval), na.rm = TRUE)
    
    max_x <- max(max_x, max(g1$data$Delta_lfc), max(g2$data$Delta_lfc), na.rm = TRUE)
    max_y <- max(max_y, max(g1$data$minus_log_pval), max(g2$data$minus_log_pval), na.rm = TRUE)
  }
  
  nda_plot_list[[tool]] <- g1
  da_plot_list[[tool]] <- g2
}


for(tool in da_result$tools){
  if(COMMON_SIZE){
    da_plot_list[[tool]] <- da_plot_list[[tool]] + ggplot2::xlim(min_x, max_x) + ylim(min_y, max_y)
    nda_plot_list[[tool]] <- nda_plot_list[[tool]] + ggplot2::xlim(min_x, max_x) + ylim(min_y, max_y)
  }
  
  current_pairplot <- cowplot::plot_grid(da_plot_list[[tool]], nda_plot_list[[tool]], align="v", ncol=2, nrow=1)
  pair_plot_list[[tool]] <- current_pairplot
}


pdf_width <- 22
pdf_heigh_per_row <- 11
title_fontsize <- 8
title_area <- 1.5

n_rows <- length(da_result$tools)
p <- cowplot::plot_grid(plotlist = pair_plot_list, nrow=n_rows, ncol=1)

pdf_height <- n_rows * pdf_heigh_per_row
rel_heights <- title_area / (title_area + pdf_height)
rel_heights <- c(rel_heights, pdf_height / (title_area + pdf_height))
text_plot <- ggplot2::ggplot() + ggplot2::annotate("text", x=1, y=1, size=title_fontsize, label=args$title) + theme_void()
p <- cowplot::plot_grid(text_plot, p, align="v", ncol=1, rel_heights=rel_heights)

ggsave(args$output_filepath, plot=p, width=pdf_width, height=pdf_height, limitsize = FALSE)
print("SUCCESSFUL PLOT RESULT")
