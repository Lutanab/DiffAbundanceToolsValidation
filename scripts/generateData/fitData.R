library(argparse)

project_dir <- "./"
source(paste0(project_dir, "/scripts/generateData/generators/MIDASim.R"))


parser <- argparse::ArgumentParser()
parser$add_argument("-g", "--data_generator", type="character", action="store", help="Method should be used for data generating")
parser$add_argument("-pf", "--fixed_parameters", type="logical", action="store", help="Shuold generation parameters be fixed within series")
parser$add_argument("-fso", "--source", type="character", action="store", help="Filepath of source real-data phyloseq Rds for fitting")
parser$add_argument("-tp", "--da_taxa_percentage", type="numeric", default=NULL, action="store", help="Precentage of diff abundant taxa between groups")
parser$add_argument("-pl", "--prevalence_lower_limit", default=NULL, type="numeric", action="store", help="Lower limit of prevalence range, in which diff abundent taxa should be picked")
parser$add_argument("-pu", "--prevalence_upper_limit", default=NULL, type="numeric", action="store", help="Upper limit of prevalence range, in which diff abundent taxa should be picked")
parser$add_argument("-fn", "--nda_table_filename", type="character", default=NULL, action="store", help="Filepath for characteristic table for non diff abundant taxa")
parser$add_argument("-fd", "--da_table_filename", type="character", default=NULL, action="store", help="Filepath for characteristic table for diff abundant taxa")
parser$add_argument("-ffc", "--fitted_coeffs_filepath", type="character", action="store", help="Filepath for fitted parameters")
parser$add_argument("-fr", "--report_file", type="character", default=NULL, action="store", help="Filepath for report file")
args <- parser$parse_args()

if(args$fixed_parameters){
  if(any(unlist(lapply(args, is.null)))){
    stop("When fixed generation parameters within a series, all arguments should be specified in fitData")
  }
  
  if((args$da_taxa_percentage < 0) | (args$da_taxa_percentage > 1)){
    stop("--da_taxa_percentage must be between 0 and 1")
  }
  if(args$prevalence_upper_limit < args$prevalence_lower_limit){
    stop("--prevalence_upper_limit must be more then --prevalence_lower_limit")
  }
  if((args$prevalence_upper_limit > 1) |  (args$prevalence_lower_limit < 0)){
    stop("--prevalence_upper_limit and --prevalence_upper_limit must be in range 0-1")
  }
}

fit_data <- list(
  "MIDASim" = fit_midasim
)
set_da_params = list(
  "MIDASim" = set_da_params_midasim
)

fitted_coeffs <- fit_data[[args$data_generator]](args$source)
saveRDS(fitted_coeffs, args$fitted_coeffs_filepath)
if(args$fixed_parameters){
  da_params <- set_da_params[[args$data_generator]](
    counts.setup=fitted_coeffs, da_taxa_percentage=args$da_taxa_percentage, prevalence_lower_limit=args$prevalence_lower_limit,
    prevalence_upper_limit=args$prevalence_upper_limit)

  cat(da_params[["report"]] ,file=args$report_file,sep="\n")
  write.table(da_params[["nda_table"]], file=args$nda_table_filename, quote=FALSE, sep='\t')
  write.table(da_params[["da_table"]], file=args$da_table_filename, quote=FALSE, sep='\t')
}

print("SUCCESSFUL GENERATED DA DATA")

