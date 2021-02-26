data_with_MV_NRMSE_test <- function(orig_data) {
  misg_data <- data_with_MV_mask_out(orig_data,
                                     option = "MAR",
                                     parameter = 1,
                                     keep = "max")$data_extra_MV
  
  imputed_data <- meanImpute(misg_data)$x
  data_with_MV_NRMSE(orig_data, misg_data, imputed_data)
}

data_with_MV_RMSE_test <- function(log_data) {
  misg_data <- data_with_MV_mask_out(log_data,
                                     option = "MAR",
                                     parameter = 1,
                                     keep = "max")$data_extra_MV
  
  imputed_data <- meanImpute(misg_data)$x
  data_with_MV_RMSE(log_data, misg_data, imputed_data)
}


data_with_MV_SOR_test <- function(orig_data) {
  misg_data <- data_with_MV_mask_out(orig_data,
                                     option = "MAR",
                                     parameter = 1,
                                     keep = "max")$data_extra_MV
  
  imputed_data <- meanImpute(misg_data)$x
  res <- data_with_MV_SOR(orig_data, misg_data, imputed_data)
  # View(res)
}