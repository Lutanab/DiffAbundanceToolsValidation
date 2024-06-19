library(types)
library(phyloseq)

project_dir <- "./"
source(paste0(project_dir, "/scripts/generateData/generators/common.R"))

fit_midasim <- function(source = ? character){
  ps <- readRDS(source)
  init_count_table <- ps@otu_table@.Data
  if(ps@otu_table@taxa_are_rows){
    init_count_table <- t(init_count_table)
  }
  counts.setup = MIDASim::MIDASim.setup(init_count_table, mode = 'parametric')
  return(counts.setup)
}

set_da_params_midasim <- function(
    counts.setup = ? list, 
    da_taxa_percentage = ? numeric, 
    prevalence_lower_limit = ? numeric, 
    prevalence_upper_limit = ? numeric
){
  taxa_prev <- counts.setup$taxa.1.prop
  da_taxa <- names(taxa_prev[(taxa_prev >= prevalence_lower_limit) & (taxa_prev <= prevalence_upper_limit)])
  if(length(da_taxa) == 0){
    stop("No taxa of chosen prevalence fraction")
  }
  da_taxa <- base::sample(da_taxa, size=ceiling(length(da_taxa)*da_taxa_percentage), replace=FALSE)
  is.da_taxa <- (counts.setup$taxa.names %in% da_taxa)
  n_da_taxa <- sum(is.da_taxa)
  nda_taxa <- counts.setup$taxa.names[!is.da_taxa]
  theor_lfcs <- stats::runif(n_da_taxa, min=1.5, max=4)
  theor_lfcs <- theor_lfcs * base::sample(c(-1, 1), size=n_da_taxa, replace=TRUE, prob=c(0.5, 0.5))
  
  mu_nda <- counts.setup$mu.est[!is.da_taxa]
  mu_da_ref <- counts.setup$mu.est[is.da_taxa]
  mu_da_tar <- mu_da_ref + theor_lfcs
  
  # mu_ref <- counts.setup$mu.est
  # mu_tar <- counts.setup$mu.est 
  # mu_tar[is.da_taxa] <- mu_tar[is.da_taxa] + theor_lfcs
  
  nda_table <- data.frame()
  nda_table[nda_taxa, "mu"] <- mu_nda
  nda_table[nda_taxa, "Q_ref"] <- counts.setup$Q.est[!is.da_taxa]
  nda_table[nda_taxa, "sigma_ref"] <- counts.setup$sigma.est[!is.da_taxa]
  nda_table[nda_taxa, "Q_tar"] <- counts.setup$Q.est[!is.da_taxa]
  nda_table[nda_taxa, "sigma_tar"] <- counts.setup$sigma.est[!is.da_taxa]
  
  da_table <- data.frame()
  da_table[da_taxa, "mu_ref"] <- mu_da_ref
  da_table[da_taxa, "Q_ref"] <- counts.setup$Q.est[is.da_taxa]
  da_table[da_taxa, "sigma_ref"] <- counts.setup$sigma.est[is.da_taxa]
  da_table[da_taxa, "mu_tar"] <- mu_da_tar
  da_table[da_taxa, "Q_tar"] <- counts.setup$Q.est[is.da_taxa]
  da_table[da_taxa, "sigma_tar"] <- counts.setup$sigma.est[is.da_taxa]
  da_table$lfc <- da_table$mu_tar - da_table$mu_ref
  
  report <- paste0(
    "DATA GENERATION REPORT:\n",
    "Number ovarall taxa: ", n_da_taxa + length(nda_taxa), ";\n",
    "Number diff abundant taxa: ", n_da_taxa, ";\n",
    "Percentage diff abundant taxa: ", da_taxa_percentage, ";\n",
    "Taxon prevalence lower limit: ", prevalence_lower_limit, ";\n",
    "Taxon prevalence upper limit: ", prevalence_upper_limit, ";"
  )
  
  return(list(
    "da_table" = da_table,
    "nda_table" = nda_table,
    "report" = report
  ))
}

generate_midasim <- function(
    fitted_coefs = ? list, da_table = ? data.frame, n_samples = ? numeric, target_sample_percentage = ? numeric, mannual_params = FALSE ? logical
){
  n_sam_tar <- as.integer(n_samples*target_sample_percentage)
  n_sam_ref <- n_samples-n_sam_tar
  is.da_taxa <- fitted_coeffs$taxa_names %in% row.names(da_table)
  
  fitted_coefs.ref <- fitted_coefs
  mu_ref <- fitted_coefs.ref$mu.est
  lib_size_ref <- get_kde_sample(init_values=fitted_coefs.ref$lib.size, sample_size=n_sam_ref)
  fitted_coefs.ref <- MIDASim::MIDASim.modify(fitted_coefs.ref, 
                                        lib.size = lib_size_ref)
  fitted_coefs.tar <- fitted_coefs
  mu_tar <- fitted_coefs.tar$mu.est 
  mu_tar[is.da_taxa] <- mu_tar[is.da_taxa] + da_table$lfc
  lib_size_tar <- get_kde_sample(init_values=fitted_coefs.tar$lib.size, sample_size=n_sam_tar)
  fitted_coefs.tar <- MIDASim::MIDASim.modify(fitted_coefs.tar,
                                        lib.size = lib_size_tar,
                                        gengamma.mu = mu_tar)
  if(mannual_params){
    fitted_coefs.tar$mu.est <- mu_tar
    fitted_coefs.ref$mu.est <- mu_ref
  }
  
  sim_data_ref = MIDASim::MIDASim(fitted_coefs.ref)
  sim_data_tar = MIDASim::MIDASim(fitted_coefs.tar)
  sim_data_ref <- sim_data_ref$sim_count
  row.names(sim_data_ref) <- paste0("Ref_sample_", seq(1, nrow(sim_data_ref)))
  sim_data_tar <- sim_data_tar$sim_count
  row.names(sim_data_tar) <- paste0("Tar_sample_", seq(1, nrow(sim_data_tar)))
  sim_data <- base::rbind(sim_data_tar, sim_data_ref)
  print("C")
  # Assembling phyloseq
  count_table <- phyloseq::otu_table(sim_data, taxa_are_rows=FALSE)
  
  sam_data <- data.frame()
  sam_data[row.names(sim_data_tar), "group"] <- "target"
  sam_data[row.names(sim_data_ref), "group"] <- "reference"
  sam_data[, "test"] <- "test"
  sam_data <- phyloseq::sample_data(sam_data)
  
  tt <- data.frame()
  tt[colnames(count_table), "Kingdom"] <- "unknown"
  tt <- as.matrix(tt)
  tt <- phyloseq::tax_table(tt)
  
  sim_ps <- phyloseq::phyloseq(count_table, tt, sam_data)
  return(sim_ps)
}

