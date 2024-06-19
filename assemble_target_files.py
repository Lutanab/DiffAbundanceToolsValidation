import json
from snakemake_config_model import GenerationParams
from snakemake_config_model import SnakemakeConfig
from typing import Text
from typing import List


def get_series_dir_name(snakemake_config: SnakemakeConfig) -> Text:
    series_dir_name = f"series_{snakemake_config.series_number}_{snakemake_config.series_name}_fixed_da_taxa_{snakemake_config.fixed_da_taxa.value}" +\
        f"_data_generator_{snakemake_config.data_generator}_source_{snakemake_config.source}"
    if bool(snakemake_config.fixed_da_taxa):
        generation_params = snakemake_config.generation_params_list[0]
        series_dir_name += f"_da_taxa_percentage_{generation_params.da_taxa_percentage}_prevalence_upper_limit_{generation_params.prevalence_upper_limit}" +\
            f"_prevalence_lower_limit_{generation_params.prevalence_lower_limit}"
    return series_dir_name


def get_result_dir_name(generation_params: GenerationParams, series_dir_name: Text, are_da_taxa_fixed: bool) -> Text:
    result_dir = f"results/{series_dir_name}/"
    if not bool(are_da_taxa_fixed):
        result_dir += f"da_taxa_percentage_{generation_params.da_taxa_percentage}_prevalence_upper_limit_{generation_params.prevalence_upper_limit}" +\
            f"_prevalence_lower_limit_{generation_params.prevalence_lower_limit}_"
    result_dir += f"n_samples_{generation_params.n_samples}_tar_samples_percentage_{generation_params.tar_samples_percentage}" +\
        f"_mannual_params_{generation_params.mannual_params.value}"
    return(result_dir)

def assemble_target_files_list(config_filename: Text = "snakemake.config.json") -> List:
    TARGET_LIST = []
    CONFIG = []

    with open(config_filename, 'r') as f:
        CONFIG = json.load(f)
    CONFIG = SnakemakeConfig.model_validate(CONFIG)
    
    are_da_taxa_fixed = bool(CONFIG.fixed_da_taxa)
    series_dir_name = get_series_dir_name(CONFIG)
    for generation_params in CONFIG.generation_params_list:
        result_dir = get_result_dir_name(generation_params, series_dir_name, are_da_taxa_fixed)
        count_heatmap = f"./{result_dir}/count_heatmap.pdf"
        TARGET_LIST.append(count_heatmap)
        da_heatmap = f"./{result_dir}/da_heatmap.pdf"
        TARGET_LIST.append(da_heatmap)
        result_scatter = f"./{result_dir}/scatter_log_pval_x_delta_lfc.pdf"
        TARGET_LIST.append(result_scatter)

        for pca_params in CONFIG.pca_params_list:
            pca_file = f"./{result_dir}/pca_transform_{pca_params.transform.value}_distance_{pca_params.distance.value}.pdf"
            TARGET_LIST.append(pca_file)
        
        for raw_glm_params in CONFIG.raw_glm_params_list:
            for transform in raw_glm_params.transforms:
                glm_result_file = f"./{result_dir}/raw_glm_da_lfc_range_{raw_glm_params.lfc_range}_p_val_treshold_{raw_glm_params.p_val_treshold}"\
                    f"/result_transform_{transform.value}.tsv"
                TARGET_LIST.append(glm_result_file)
        
        for malual_raw_glm_param in CONFIG.manual_raw_glm_params_list:
            for transform in malual_raw_glm_param.transforms:
                manual_glm_result_file = f"./{result_dir}/manual_raw_glm_da_lfc_range_{malual_raw_glm_param.lfc_range}/result_transform_{transform.value}.tsv"
                TARGET_LIST.append(manual_glm_result_file)

        for metrics_params in CONFIG.metrics_params_list:
            metrics_file = f"./{result_dir}/metrics_lfc_range_{metrics_params.lfc_range}_p_val_treshold_{metrics_params.p_val_treshold}.tsv"
            TARGET_LIST.append(metrics_file)
        
    return(TARGET_LIST, are_da_taxa_fixed)


if __name__ == "__main__":
    TARGET_LIST, ARA_DA_TAXA_FIXED = assemble_target_files_list()
    print(ARA_DA_TAXA_FIXED, "\n\n")
    print("\n".join(TARGET_LIST))