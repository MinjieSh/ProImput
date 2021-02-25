cos_sim_test <- function() {
  mat <- matrix(1:20, ncol = 5, nrow = 4)
  mat[sample(4 * 5, 3)] <- NA
  cosine_similarity <- cos.sim(mat)
}

FRMF_test <- function(grnd_truth) {
  rank <- 3
  
  misg_data <-
     data_with_MV_mask_out(grnd_truth, option = "MAR", 0.3, keep = "max")$data_extra_MV
  
  # No external information incorporated
  res <- FRMF(
    misg_data,
    rank = rank,
    lrn_rate = 0.0005,
    max_step = 10000,
    conv_thld = 0.005,
    lam_A = 0.005,
    lam_S = 0.5,
    nbr_info_mat = NULL,
    nbr_thld = 1,
    coeff_crs_rglr = 0,
    normalization = "std_score"
  )
  print(data_with_MV_NRMSE(grnd_truth, misg_data, res$imp))
  
  # Self information incorporated
  res <- FRMF(
    misg_data,
    rank = rank,
    lrn_rate = 0.0005,
    max_step = 10000,
    conv_thld = 0.005,
    lam_A = 0.005,
    lam_S = 0.5,
    nbr_info_mat = misg_data,
    nbr_thld = 0.9,
    coeff_crs_rglr = 0.5,
    normalization = "std_score"
  )
  print(data_with_MV_NRMSE(grnd_truth, misg_data, res$imp))
  
  # External information incorporated
  res <- FRMF(
    misg_data,
    rank = rank,
    lrn_rate = 0.0005,
    max_step = 10000,
    conv_thld = 0.005,
    lam_A = 0.005,
    lam_S = 0.5,
    nbr_info_mat = grnd_truth,
    nbr_thld = 0.9,
    coeff_crs_rglr = 5,
    normalization = "std_score"
  )
  print(data_with_MV_NRMSE(grnd_truth, misg_data, res$imp))
  
  predicted <- impute(misg_data,
                      "nipals",
                      parameter = rank)$res
  print(data_with_MV_NRMSE(grnd_truth, misg_data, predicted))
  
  predicted <- impute(misg_data,
                      "svdImpute",
                      parameter = rank)$res
  print(data_with_MV_NRMSE(grnd_truth, misg_data, predicted))
}