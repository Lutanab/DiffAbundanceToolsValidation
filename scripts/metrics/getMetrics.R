get_range_metrics_tab <- function(characteristic_table = ? data.frame, da_label = ? character, tool = ? character,
                                  lfc_range = ? numeric, p_val_treshold = ? numeric){
  curr_taxa <- row.names(characteristic_table)
  tmp_df <- data.frame()
  tmp_df[tool, "Metrics_type"] <- "Range"
  tmp_df[tool, "DA_type"] <- da_label
  tmp_df[tool, "Tool"] <- tool
  tmp_df[tool, "Detected_right"] <- sum(
    (abs(characteristic_table$lfc_est - characteristic_table$lfc_true) < (lfc_range/2)) & (characteristic_table$p_val < p_val_treshold), na.rm=TRUE
  )
  tmp_df[tool, "Detected_wrong"] <- sum(
    (abs(characteristic_table$lfc_est - characteristic_table$lfc_true) >= (lfc_range/2)) & (characteristic_table$p_val < p_val_treshold), na.rm=TRUE
  )
  tmp_df[tool, "Not_detected"] <- length(curr_taxa) - (tmp_df$Detected_right + tmp_df$Detected_wrong)
  tmp_df[tool, c("Detected_right", "Detected_wrong", "Not_detected")] <- tmp_df[tool, c("Detected_right", "Detected_wrong", "Not_detected")] / length(curr_taxa)
  return(tmp_df)
}

get_unilateral_metrics_tab <- function(characteristic_table = ? data.frame, da_label = ? character, tool = ? character, p_val_treshold = ? numeric){
  curr_taxa <- row.names(characteristic_table)
  tmp_df <- data.frame()
  tmp_df[tool, "Metrics_type"] <- "Unilateral"
  tmp_df[tool, "DA_type"] <- da_label
  tmp_df[tool, "Tool"] <- tool
  tmp_df[tool, "Detected_right"] <- sum(((characteristic_table$lfc_est * characteristic_table$lfc_true) > 0) & (characteristic_table$p_val < p_val_treshold), na.rm=TRUE)
  tmp_df[tool, "Detected_wrong"] <- sum(((characteristic_table$lfc_est * characteristic_table$lfc_true) < 0) & (characteristic_table$p_val < p_val_treshold), na.rm=TRUE)
  tmp_df[tool, "Not_detected"] <- length(curr_taxa) - (tmp_df$Detected_right + tmp_df$Detected_wrong)
  tmp_df[tool, c("Detected_right", "Detected_wrong", "Not_detected")] <- tmp_df[tool, c("Detected_right", "Detected_wrong", "Not_detected")] / length(curr_taxa)
  return(tmp_df)
}