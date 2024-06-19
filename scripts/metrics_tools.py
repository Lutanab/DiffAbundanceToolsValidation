import sys
sys.path.insert(1, '../')
import json
import pandas as pd
from assemble_target_files import get_series_dir_name
from assemble_target_files import get_result_dir_name
from snakemake_config_model import GenerationParams
from snakemake_config_model import SnakemakeConfig
from frozendict import frozendict
from typing import Dict
from typing import List
from typing import Text


def get_metrics_list(config_filename: Text = "snakemake.config.json") -> Dict:
    METRICS_FILES = {}
    CONFIG = []
    
    with open(config_filename, 'r') as f:
        CONFIG = json.load(f)
    CONFIG = SnakemakeConfig.model_validate(CONFIG)
    
    are_da_taxa_fixed = bool(CONFIG.fixed_da_taxa)
    series_dir_name = get_series_dir_name(CONFIG)
    series_params_fd = frozendict({
        feature: getattr(CONFIG, feature)
        for feature in ["series_number", "series_name",  "fixed_da_taxa", "data_generator", "source"]
    })
    METRICS_FILES[series_params_fd] = {}
    for generation_params in CONFIG.generation_params_list:
        generation_params_fd = frozendict(generation_params)
        result_dir = get_result_dir_name(generation_params, series_dir_name, are_da_taxa_fixed)

        for metrics_params in CONFIG.metrics_params_list:
            metrics_file = f"./{result_dir}/metrics_lfc_range_{metrics_params.lfc_range}_p_val_treshold_{metrics_params.p_val_treshold}.tsv"
            metrics_params_fd = frozendict(metrics_params)
            if metrics_params_fd not in METRICS_FILES[series_params_fd].keys():
                METRICS_FILES[series_params_fd][metrics_params_fd] = {}
            METRICS_FILES[series_params_fd][metrics_params_fd][generation_params_fd] = metrics_file
        
    return(METRICS_FILES)

def add_data_generation_params(df: pd.DataFrame, data_params: frozendict) -> pd.DataFrame:
    for feature in data_params.keys():
        df[[feature]] = data_params[feature]
    return df

def assemble_metrics_tables(METRICS_FILES: Dict[frozendict, Text]):
    range_da_df = pd.DataFrame()
    range_nda_df = pd.DataFrame()
    uni_df = pd.DataFrame()

    for data_params, file in METRICS_FILES.items():
        current_metrics = pd.read_csv(file, sep = '\t')
        needed_cols = list(current_metrics.columns[2:])
        # Range DA
        tmp_df = current_metrics[(current_metrics["Metrics_type"] == "Range") & (current_metrics["DA_type"] == "DA")]
        tmp_df = tmp_df[needed_cols]
        tmp_df = add_data_generation_params(tmp_df, data_params)
        range_da_df = pd.concat([range_da_df, tmp_df])
        # Range NDA
        tmp_df = current_metrics[(current_metrics["Metrics_type"] == "Range") & (current_metrics["DA_type"] == "NDA")]
        tmp_df = tmp_df[needed_cols]
        tmp_df = add_data_generation_params(tmp_df, data_params)
        range_nda_df = pd.concat([range_nda_df, tmp_df])
        # Unilateral
        tmp_df = current_metrics[(current_metrics["Metrics_type"] == "Unilateral")]
        tmp_df = tmp_df[needed_cols]
        tmp_df = add_data_generation_params(tmp_df, data_params)
        uni_df = pd.concat([uni_df, tmp_df])

    range_da_df = range_da_df.reset_index(drop=True)
    range_nda_df = range_nda_df.reset_index(drop=True)
    uni_df = uni_df.reset_index(drop=True)

    metrics_tables = {"range_da": range_da_df, "range_nda": range_nda_df, "uni": uni_df}
    for approach in metrics_tables.keys():
        metrics_tables[approach] = metrics_tables[approach].rename(columns={
            "Detected_right": "Sure_right",
            "Detected_wrong": "Sure_wrong",
            "Not_detected": "Not_sure"
            })
        metrics_tables[approach]['mannual_params'] = metrics_tables[approach]['mannual_params'].apply(lambda z: z.value)
    return metrics_tables
