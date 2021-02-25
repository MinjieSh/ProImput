impute_test <- function() {
  data_with_MV <- matrix(1:20,
                         ncol = 4,
                         nrow = 5,
                         byrow = FALSE)
  data_with_MV[sample(20, 4)] <- NA
  
  temp_pred <- impute(data_with_MV, "pwHalfMin (log)")
  pred <- temp_pred$imp
  
  temp_pred <- impute(data_with_MV, "pwMean")
  pred <- temp_pred$imp
  
  data_with_MV <- matrix(
    c(1, 0, 0, 0,
      10, 1, 2, 1.5,
      1.1, 1, 1, 1,
      10, 9, 10.5, 9.8),
    ncol = 4,
    nrow = 4,
    byrow = FALSE
  )
  data_with_MV[c(4, 10)] <- NA
  
  temp_pred <- impute(data_with_MV, "pwKNN", 1)
  pred <- temp_pred$imp
  
  data_with_MV <- matrix(
    c(1, 0, 0, 0,
      10, 1, 2, 1.5,
      1.1, 1, 1, 1,
      10, 9, 10.5, 9.8),
    ncol = 4,
    nrow = 4,
    byrow = TRUE
  )
  data_with_MV[c(4, 10)] <- NA
  
  temp_pred <- impute(data_with_MV, "swKNN", 1)
  pred <- temp_pred$imp
  
  rank <- 3
  temp <- gnr_sim_data(200, 1000, rank, 0.7)
  grnd_truth <- temp$X_wo_n
  misg_data <- grnd_truth
  misg_data[temp$misg_idx] <- NA
  
  temp_pred <- impute(misg_data, "PPCA", parameter = 3)
  pred <- temp_pred$imp
  print(data_with_MV_NRMSE(grnd_truth, misg_data, pred))
  
  temp_pred <- impute(misg_data, "NIPALS", parameter = 3)
  pred <- temp_pred$imp
  print(data_with_MV_NRMSE(grnd_truth, misg_data, pred))
  
  temp_pred <- impute(misg_data, "SVD", parameter = 3)
  pred <- temp_pred$imp
  print(data_with_MV_NRMSE(grnd_truth, misg_data, pred))
  
  temp_pred <- impute(misg_data, "SVT", parameter = 100)
  pred <- temp_pred$imp
  print(data_with_MV_NRMSE(grnd_truth, misg_data, pred))
  
  temp_pred <- impute(misg_data, "FRMF", parameter = 3)
  pred <- temp_pred$imp
  print(data_with_MV_NRMSE(grnd_truth, misg_data, pred))
  
  # grnd_truth <- log2(proteinC)
  misg_data <-
    data_with_MV_mask_out(grnd_truth, option = "MAR", 0.3, keep = "max")$data_extra_MV
  
  temp_pred <-
    impute(
      misg_data,
      "CAM_cmplt",
      parameter = 3,
      dim.rdc = 10,
      cluster.num = 35
    )
  pred <- temp_pred$imp
  print(data_with_MV_NRMSE(grnd_truth, misg_data, pred))
  
  temp_pred <-
    impute(
      misg_data,
      "CAM_NIPALS",
      parameter = 3,
      dim.rdc = 10,
      cluster.num = 35
    )
  pred <- temp_pred$imp
  print(data_with_MV_NRMSE(grnd_truth, misg_data, pred))
  
  temp_pred <-
    impute(
      misg_data,
      "CAM_SVT",
      parameter = 3,
      dim.rdc = 10,
      cluster.num = 35
    )
  pred <- temp_pred$imp
  print(data_with_MV_NRMSE(grnd_truth, misg_data, pred))
  
  temp_pred <-
    impute(
      misg_data,
      "CAM_smpClus",
      parameter = 3,
      dim.rdc = 10,
      cluster.num = 35
    )
  pred <- temp_pred$imp
  print(data_with_MV_NRMSE(grnd_truth, misg_data, pred))
  
  temp_pred <-
    impute(
      misg_data,
      "CAM_smpClus_iter",
      parameter = 3,
      dim.rdc = 10,
      cluster.num = 35
    )
  pred <- temp_pred$imp
  print(data_with_MV_NRMSE(grnd_truth, misg_data, pred))
  
  temp_pred <-
    impute(
      misg_data,
      "FRMF",
      nbr_info_mat = NULL,
      nbr_thld = 0.9,
      parameter = 3
    )
  pred <- temp_pred$imp
  print(data_with_MV_NRMSE(grnd_truth, misg_data, pred))
  
}
