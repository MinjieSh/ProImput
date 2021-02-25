
data_with_MV_mask_out_test <- function(data_with_MV){
  
  data_with_MV <- matrix(1:100, ncol = 10, nrow = 10) 
  data_with_MV[sample(100, 20)] <- NA
    
  misg_data <- data_with_MV_mask_out(data_with_MV,
                                     option = "MAR",
                                     parameter = 1,
                                     keep = "random")$data_extra_MV
  
  plot_pattern_change(data_with_MV, misg_data, "data_with_MV_mask_out_test")
  
}
