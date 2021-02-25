
plot_missing_mechanism_test <- function(orig_data) {
  
  plot_missing_mechanism(orig_data, "plot_missing_mechanism_test")
  
}

plot_pattern_change_test <- function(orig_data) {
  misg_data <- data_with_MV_mask_out(
    orig_data,
    option = "MAR",
    parameter = 1,
    keep = "max"
  )$data_extra_MV
  
  plot_pattern_change(orig_data, misg_data, "plot_pattern_change_test")
}
