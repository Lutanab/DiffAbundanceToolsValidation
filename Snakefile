import json
from assemble_target_files import assemble_target_files_list

TARGET_LIST, FIXED_DA_TAXA = assemble_target_files_list()

# Specify wildcards

series_dir = "series_{series_label}_fixed_da_taxa_{fixed_da_taxa}_data_generator_{data_generator}_source_{source}"
result_dir = "n_samples_{n_samples}_tar_samples_percentage_{tar_samples_percentage}_mannual_params_{mannual_params}"

if FIXED_DA_TAXA:
    series_dir += "_da_taxa_percentage_{da_taxa_percentage}_prevalence_upper_limit_{prevalence_upper_limit}_prevalence_lower_limit_{prevalence_lower_limit}"
    da_table_wc = "./results/" + series_dir + "/da_table.tsv"
    nda_table_wc = "./results/" + series_dir + "/nda_table.tsv"
    report_file_wc = "./results/" + series_dir + "/report.txt"
else:
    result_dir = "da_taxa_percentage_{da_taxa_percentage}_prevalence_upper_limit_{prevalence_upper_limit}_"\
        "prevalence_lower_limit_{prevalence_lower_limit}_" + result_dir
    da_table_wc = "./results/" + series_dir + "/" + result_dir + "/da_table.tsv"
    nda_table_wc = "./results/" + series_dir + "/" + result_dir + "/nda_table.tsv"
    report_file_wc = "./results/" + series_dir + "/" + result_dir + "/report.txt"

source_wc = "./sources/{source}.rds"
ps_syn_wc = "./results/" + series_dir + "/" + result_dir + "/ps_syn.rds"
fitted_coeffs_wc = "./results/" + series_dir + "/fitted_coeffs.rds"

pca_plot_wc = "./results/" + series_dir + "/" + result_dir + "/pca_transform_{transform}_distance_{distance}.pdf"
count_heatmap_wc = "./results/" + series_dir + "/" + result_dir + "/count_heatmap.pdf"
diff_abund_analysis_result_wc = "./results/" + series_dir + "/" + result_dir + "/diff_abund_analysis_result.rds"

raw_glm_diff_abund_result_wc = "./results/" + series_dir + "/" + result_dir +\
    "/raw_glm_da_lfc_range_{lfc_range}_p_val_treshold_{p_val_treshold}/result_transform_{transform}.tsv"

manual_raw_glm_diff_abund_result_wc = "./results/" + series_dir + "/" + result_dir + "/manual_raw_glm_da_lfc_range_{lfc_range}/result_transform_{transform}.tsv"

da_heatmap_wc = "./results/" + series_dir + "/" + result_dir + "/da_heatmap.pdf"
scatter_log_pval_x_delta_lfc_wc = "./results/" + series_dir + "/" + result_dir + "/scatter_log_pval_x_delta_lfc.pdf"
quality_metrics_wc = "./results/" + series_dir + "/" + result_dir + "/metrics_lfc_range_{lfc_range}_p_val_treshold_{p_val_treshold}.tsv"

# Rules
rule all:
    input: TARGET_LIST

if FIXED_DA_TAXA:
    rule fit_data:
        input:
            source = source_wc
        output:
            nda_table = nda_table_wc,
            da_table = da_table_wc,
            fitted_coeffs = fitted_coeffs_wc,
            report_file = report_file_wc
        shell:
            "Rscript ./scripts/generateData/fitData.R -g {wildcards.data_generator} -pf {wildcards.fixed_da_taxa} -fso {input.source}"
            " -tp {wildcards.da_taxa_percentage} -pu {wildcards.prevalence_upper_limit} -pl {wildcards.prevalence_lower_limit} -fn {output.nda_table}"
            " -fd {output.da_table} -ffc {output.fitted_coeffs} -fr {output.report_file}"
    
    rule generate_data:
        input:
            fitted_coeffs = fitted_coeffs_wc
        output:
            ps_syn = ps_syn_wc
        shell:
            "Rscript ./scripts/generateData/generateData.R -g {wildcards.data_generator} -pf {wildcards.fixed_da_taxa} -ffc {input.fitted_coeffs}"
            " -ns {wildcards.n_samples} -sp {wildcards.tar_samples_percentage} -mp {wildcards.mannual_params} -fsy {output.ps_syn}"
else:
    rule fit_data:
        input:
            source = source_wc
        output:
            fitted_coeffs = fitted_coeffs_wc,
        shell:
            "Rscript ./scripts/generateData/fitData.R -g {wildcards.data_generator} -pf {wildcards.fixed_da_taxa} -fso {input.source}"
            " -ffc {output.fitted_coeffs}"
    
    rule generate_data:
        input:
            fitted_coeffs = fitted_coeffs_wc
        output:
            nda_table = nda_table_wc,
            da_table = da_table_wc,
            report_file = report_file_wc,
            ps_syn = ps_syn_wc
        shell:
            "Rscript ./scripts/generateData/generateData.R -g {wildcards.data_generator} -pf {wildcards.fixed_da_taxa} -ffc {input.fitted_coeffs}"
            " -tp {wildcards.da_taxa_percentage} -pu {wildcards.prevalence_upper_limit} -pl {wildcards.prevalence_lower_limit}"
            " -ns {wildcards.n_samples} -sp {wildcards.tar_samples_percentage} -mp {wildcards.mannual_params}"
            " -fn {output.nda_table} -fd {output.da_table} -fr {output.report_file} -fsy {output.ps_syn}"

rule plot_pca:
    input:
        ps_syn = ps_syn_wc
    output:
        pca_plot = pca_plot_wc
    shell:
        "Rscript ./scripts/plot_pca.R -t {wildcards.transform} -d {wildcards.distance} -p {input.ps_syn} -o {output.pca_plot}"
        " -ti 'Source: {wildcards.source}; DA taxa perc: {wildcards.da_taxa_percentage}\\nN samples: {wildcards.n_samples};"
        " Target samples perc: {wildcards.tar_samples_percentage}\\nPrevalence range: "
        "{wildcards.prevalence_lower_limit} - {wildcards.prevalence_upper_limit}; Mannual params used: {wildcards.mannual_params}'"

rule plot_count_heatmap:
    input:
        ps_syn = ps_syn_wc,
        da_table = da_table_wc
    output:
        count_heatmap = count_heatmap_wc
    shell:
        "Rscript ./scripts/plot_count_heatmap.R -p {input.ps_syn} -fd {input.da_table} -o {output.count_heatmap}"
        " -ti 'Source: {wildcards.source}; DA taxa perc: {wildcards.da_taxa_percentage}; N samples: {wildcards.n_samples}; Target samples perc: {wildcards.tar_samples_percentage}; Prevalence range: {wildcards.prevalence_lower_limit} - {wildcards.prevalence_upper_limit}; Mannual params used: {wildcards.mannual_params}'"

rule perform_de_analysis:
    input:
        ps_syn = ps_syn_wc
    output:
        diff_abund_analysis_result = diff_abund_analysis_result_wc
    shell:
        "Rscript ./scripts/perform_de_analysis.R -p {input.ps_syn} -o {output.diff_abund_analysis_result}"

rule perform_glm_da_analysis:
    input:
        ps_syn = ps_syn_wc,
        da_table = da_table_wc
    output:
        result = raw_glm_diff_abund_result_wc
    shell:
        "Rscript ./scripts/additional/perform_raw_glm.R -p {input.ps_syn} -t {wildcards.transform} -fd {input.da_table} -o {output.result}"
        " -lr {wildcards.lfc_range} -pt {wildcards.p_val_treshold}"

rule perform_manual_glm_da_analysis:
    input:
        ps_syn = ps_syn_wc,
        da_table = da_table_wc
    output:
        result = manual_raw_glm_diff_abund_result_wc
    shell:
        "Rscript ./scripts/additional/perform_raw_glm_manually.R -p {input.ps_syn} -t {wildcards.transform} -fd {input.da_table} -o {output.result}"
        " -lr {wildcards.lfc_range}"

rule plot_da_heatmap:
    input:
        diff_abund_analysis_result = diff_abund_analysis_result_wc,
        da_table = da_table_wc
    output:
        da_heatmap = da_heatmap_wc
    shell:
        "Rscript ./scripts/plot_da_heatmap.R -d {input.diff_abund_analysis_result} -fd {input.da_table} -o {output.da_heatmap}"
        " -ti 'Source: {wildcards.source}; DA taxa perc: {wildcards.da_taxa_percentage}\\nN samples: {wildcards.n_samples};"
        " Target samples perc: {wildcards.tar_samples_percentage}\\nPrevalence range:"
        " {wildcards.prevalence_lower_limit} - {wildcards.prevalence_upper_limit}; Mannual params used: {wildcards.mannual_params}'"

rule plot_scatter_log_pval_x_delta_lfc:
    input:
        diff_abund_analysis_result = diff_abund_analysis_result_wc,
        da_table = da_table_wc
    output:
        result = scatter_log_pval_x_delta_lfc_wc
    shell:
        "Rscript ./scripts/scatter_log_pval_x_delta_LFC.R -d {input.diff_abund_analysis_result} -fd {input.da_table} -o {output.result} -cs TRUE"
        " -ti 'Source: {wildcards.source}; DA taxa perc: {wildcards.da_taxa_percentage}; N samples: {wildcards.n_samples};"
        " Target samples perc: {wildcards.tar_samples_percentage}\\nPrevalence range: {wildcards.prevalence_lower_limit} - {wildcards.prevalence_upper_limit};"
        " Mannual params used: {wildcards.mannual_params}'"

rule calculate_quality_metrics:
    input:
        diff_abund_analysis_result = diff_abund_analysis_result_wc,
        da_table = da_table_wc
    output:
        result = quality_metrics_wc
    shell:
        "Rscript ./scripts/calculate_metrics.R -d {input.diff_abund_analysis_result} -fd {input.da_table} -o {output.result}"
        " -lr {wildcards.lfc_range} -pt {wildcards.p_val_treshold}"