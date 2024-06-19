library(microbiome)
library(phyloseq)
library(MIDASim)
library(argparse)
library(spatstat)

get_kde_sample <- function(
    init_values = ? numeric, sample_size = ? integer
){
  kde_den <- stats::density(init_values, n=5000)
  # tt <- stats::quantile(kde_den, probs = stats::runif(sample_size))
  tt <- quantile(kde_den, probs = stats::runif(sample_size))
  return(as.numeric(tt))
}

parser <- argparse::ArgumentParser()
parser$add_argument("-fso", "--source", type="character", action="store", help="Filepath of source real-data phyloseq Rds for fitting")
parser$add_argument("-tp", "--da_taxa_percentage", type="numeric", action="store", help="Precentage of diff abundant taxa between groups")
parser$add_argument("-ns", "--n_samples", type="integer", action="store", help="Summary number of taxa")
parser$add_argument("-sp", "--target_sample_percentage", type="numeric", action="store", help="Precentage of target samples among all")
parser$add_argument("-pl", "--prevalence_lower_limit", type="numeric", action="store", help="Lower limit of prevalence range, in which diff abundent taxa should be picked")
parser$add_argument("-pu", "--prevalence_upper_limit", type="numeric", action="store", help="Upper limit of prevalence range, in which diff abundent taxa should be picked")
parser$add_argument("-mp", "--mannual_params", type="logical", action="store", help="Should manual mu parameter be used")
parser$add_argument("-fn", "--nda_table_filename", type="character", action="store", help="Filepath for characteristic table for non diff abundant taxa")
parser$add_argument("-fd", "--da_table_filename", type="character", action="store", help="Filepath for characteristic table for diff abundant taxa")
parser$add_argument("-fsy", "--synthetic_phyloseq_filename", type="character", action="store", help="Filepath for generated synthetic data as phyloseq object")
parser$add_argument("-fr", "--report_file", type="character", action="store", help="Filepath for report file")
args <- parser$parse_args()
if((args$da_taxa_percentage < 0) | (args$da_taxa_percentage > 1)){
  stop("--da_taxa_percentage must be between 0 and 1")
}
if((args$target_sample_percentage < 0) | (args$target_sample_percentage > 1)){
  stop("--target_sample_percentage must be between 0 and 1")
}
if(args$prevalence_upper_limit < args$prevalence_lower_limit){
  stop("--prevalence_upper_limit must be more then --prevalence_lower_limit")
}
if((args$prevalence_upper_limit > 1) |  (args$prevalence_lower_limit < 0)){
  stop("--prevalence_upper_limit and --prevalence_upper_limit must be in range 0-1")
}


# Fitting log-gamma model
ps <- readRDS(args$source)
init_count_table <- ps@otu_table@.Data
if(ps@otu_table@taxa_are_rows){
  init_count_table <- t(init_count_table)
}
counts.setup = MIDASim::MIDASim.setup(init_count_table, mode = 'parametric')


# Generating coeficients
n_sam_tar <- as.integer(args$n_samples*args$target_sample_percentage)
n_sam_ref <- args$n_samples-n_sam_tar

taxa_prev <- counts.setup$taxa.1.prop
da_taxa <- names(taxa_prev[(taxa_prev >= args$prevalence_lower_limit) & (taxa_prev <= args$prevalence_upper_limit)])
if(length(da_taxa) == 0){
  stop("No taxa of chosen prevalence fraction")
}
da_taxa <- base::sample(da_taxa, size=ceiling(length(da_taxa)*args$da_taxa_percentage), replace=FALSE)
is.da_taxa <- (counts.setup$taxa.names %in% da_taxa)
n_da_taxa <- sum(is.da_taxa)
nda_taxa <- counts.setup$taxa.names[!is.da_taxa]
theor_lfcs <- stats::runif(n_da_taxa, min=1.5, max=4)
theor_lfcs <- theor_lfcs * base::sample(c(-1, 1), size=n_da_taxa, replace=TRUE, prob=c(0.5, 0.5))

counts.setup.ref <- counts.setup
mu_ref <- counts.setup.ref$mu.est
lib_size_ref <- get_kde_sample(init_values=counts.setup.ref$lib.size, sample_size=n_sam_ref)
counts.modified.ref <- MIDASim.modify(counts.setup.ref, 
                                      lib.size = lib_size_ref)

counts.setup.tar <- counts.setup
mu_tar <- counts.setup.tar$mu.est 
mu_tar[is.da_taxa] <- mu_tar[is.da_taxa] + theor_lfcs
lib_size_tar <- get_kde_sample(init_values=counts.setup.tar$lib.size, sample_size=n_sam_tar)
counts.modified.tar <- MIDASim.modify(counts.setup.tar,
                                      lib.size = lib_size_tar,
                                      gengamma.mu = mu_tar)
if(args$mannual_params){
  counts.modified.tar$mu.est <- mu_tar
}


# Assembling characteristic tables
nda_table <- data.frame()
nda_table[nda_taxa, "mu"] <- mu_ref[!is.da_taxa]
nda_table[nda_taxa, "Q_ref"] <- counts.modified.ref$Q.est[!is.da_taxa]
nda_table[nda_taxa, "sigma_ref"] <- counts.modified.ref$sigma.est[!is.da_taxa]
nda_table[nda_taxa, "Q_tar"] <- counts.modified.tar$Q.est[!is.da_taxa]
nda_table[nda_taxa, "sigma_tar"] <- counts.modified.tar$sigma.est[!is.da_taxa]

da_table <- data.frame()
da_table[da_taxa, "mu_ref"] <- mu_ref[is.da_taxa]
da_table[da_taxa, "Q_ref"] <- counts.modified.ref$Q.est[is.da_taxa]
da_table[da_taxa, "sigma_ref"] <- counts.modified.ref$sigma.est[is.da_taxa]
da_table[da_taxa, "mu_tar"] <- mu_tar[is.da_taxa]
da_table[da_taxa, "Q_tar"] <- counts.modified.tar$Q.est[is.da_taxa]
da_table[da_taxa, "sigma_tar"] <- counts.modified.tar$sigma.est[is.da_taxa]
da_table$lfc <- da_table$mu_tar - da_table$mu_ref


# Perform
sim_data_ref = MIDASim::MIDASim(counts.modified.ref)
sim_data_tar = MIDASim::MIDASim(counts.modified.tar)
sim_data_ref <- sim_data_ref$sim_count
row.names(sim_data_ref) <- paste0("Ref_sample_", seq(1, nrow(sim_data_ref)))
sim_data_tar <- sim_data_tar$sim_count
row.names(sim_data_tar) <- paste0("Tar_sample_", seq(1, nrow(sim_data_tar)))
sim_data <- base::rbind(sim_data_tar, sim_data_ref)

# Assembling phyloseq
count_table <- otu_table(sim_data, taxa_are_rows=FALSE)

sam_data <-  data.frame()
sam_data[row.names(sim_data_tar), "group"] <- "target"
sam_data[row.names(sim_data_ref), "group"] <- "reference"
sam_data[, "test"] <- "test"
sam_data <- sample_data(sam_data)

tt <- data.frame()
tt[colnames(count_table), "Kingdom"] <- "unknown"
tt <- as.matrix(tt)
tt <- tax_table(tt)

sim_ps <- phyloseq::phyloseq(count_table, tt, sam_data)

# Constructing report table
report <- paste0(
  "DATA GENERATION REPORT:\n",
  "Number ovarall taxa: ", ntaxa(ps), ";\n",
  "Number diff abundant taxa: ", n_da_taxa, ";\n",
  "Percentage diff abundant taxa: ", args$da_taxa_percentage, ";\n",
  "Number of samples: ", args$n_samples, ";\n",
  "Target sample percentage: ", args$target_sample_percentage, ";\n",
  "Taxon prevalence lower limit: ", args$prevalence_lower_limit, ";\n",
  "Taxon prevalence upper limit: ", args$prevalence_upper_limit, ";\n",
  "Mannully set params used: ", args$mannual_params, ";"
)

# Saving all results
cat(report ,file=args$report_file,sep="\n")
write.table(nda_table, file=args$nda_table_filename, quote=FALSE, sep='\t')
write.table(da_table, file=args$da_table_filename, quote=FALSE, sep='\t')
saveRDS(sim_ps, args$synthetic_phyloseq_filename)

print("SUCCESSFUL GENERATED DA DATA")
