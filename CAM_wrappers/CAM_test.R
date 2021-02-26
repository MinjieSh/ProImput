gnr_sim_data <- function(smp, prt, rank, misg_rate) {
  A <- matrix(runif(
    n = smp * rank,
    min = 0,
    max = 1
  ), smp, rank)
  
  S <- matrix(runif(n = rank * prt, min = 0, max = 1), rank, prt)
  
  X_wo_n <- abs(A %*% S)
  X_w_n <- abs(A %*% S + matrix(rnorm(smp * prt), smp, prt) / 10)
  
  colnames(X_wo_n) <- 1:ncol(X_wo_n)
  colnames(X_w_n) <- 1:ncol(X_w_n)
  
  misg_idx <- sample(smp * prt, smp  * prt * misg_rate)
  
  return(list(
    A = A,
    S = S,
    X_wo_n = X_wo_n,
    X_w_n = X_w_n,
    misg_idx = misg_idx
  ))
}

CAM_helper_nnls_test <- function() {
  sim_data <- gnr_sim_data(200, 1000, 5, 0.1)
  
  X <- sim_data$X_wo_n
  X[sim_data$misg_idx] <- NA
  
  S <- CAM_helper_nnls(sim_data$A, X)
  
  return(all.equal(S, sim_data$S))
  
}

CAM_helper_test <- function() {
  dataT_true <- t(gnr_sim_data(200, 1000, 10, 0.1)$X_wo_n)
  
  res <- CAM_helper(dataT_true,
                    dim.rdc = 10,
                    cluster.num = 35)
  
  data <- res[[10]]@Aest %*% t(res[[10]]@Sest)
  
  plot(diag(cor(data, t(dataT_true))))
}

CAM_NIPALS_test <- function() {
  rank <- 5
  temp <- gnr_sim_data(100, 1000, rank, 0.1)
  grnd_truth <- temp$X_wo_n
  misg_data <- grnd_truth
  misg_data[temp$misg_idx] <- NA
  
  res <- CAM_NIPALS(
    misg_data,
    dim.rdc = 10,
    cluster.num = 35,
    nPCs = rank,
    normalization = "std_score"
  )
  
  err <- c()
  for (i in 2:10) {
    err <- c(err, data_with_MV_NRMSE(grnd_truth, misg_data, res$res[[i]]))
  }
  plot(err)
  
  print(data_with_MV_NRMSE(grnd_truth, misg_data, res$init_data))
}

CAM_SVT_test <- function() {
  temp <- gnr_sim_data(100, 1000, 5, 0.1)
  grnd_truth <- temp$X_wo_n
  misg_data <- grnd_truth
  misg_data[temp$misg_idx] <- NA
  
  res <- CAM_SVT(
    misg_data,
    dim.rdc = 10,
    cluster.num = 35,
    normalization = "std_score"
  )
  err <- c()
  for (i in 2:10) {
    err <- c(err, data_with_MV_NRMSE(grnd_truth, misg_data, res$res[[i]]))
  }
  plot(err)
  
  print(data_with_MV_NRMSE(grnd_truth, misg_data, res$init_data))
}

CAM_cmplt_test <- function(grnd_truth) {
  misg_data <-
    data_with_MV_mask_out(grnd_truth, option = "MAR", 0.3, keep = "max")$data_extra_MV
  
  res <- CAM_cmplt(
    misg_data,
    dim.rdc = 10,
    cluster.num = 35,
    normalization = "std_score"
  )
  err <- c()
  for (i in 2:10) {
    err <- c(err, data_with_MV_NRMSE(grnd_truth, misg_data, res$res[[i]]))
  }
  plot(err)
  
  predicted <- impute(misg_data,
                      "mean")$res
  print(data_with_MV_NRMSE(grnd_truth, misg_data, predicted))
}

CAM_smpClus_test <- function() {
  rank <- 5
  temp <- gnr_sim_data(200, 1000, rank, 0.7)
  grnd_truth <- temp$X_wo_n
  misg_data <- grnd_truth
  misg_data[temp$misg_idx] <- NA
  
  res <- CAM_smpClus(
    misg_data,
    dim.rdc = 10,
    cluster.num = 35,
    normalization = "none"
  )
  err <- c()
  for (i in 2:10) {
    err <- c(err, data_with_MV_NRMSE(grnd_truth, misg_data, res$res[[i]]))
  }
  plot(err)
  
  predicted <- impute(misg_data,
                      "nipals",
                      parameter = rank)$res
  print(data_with_MV_NRMSE(grnd_truth, misg_data, predicted))
  
}

CAM_smpClus_iter_test <- function() {
  rank <- 3
  temp <- gnr_sim_data(200, 1000, rank, 0.7)
  grnd_truth <- temp$X_wo_n
  misg_data <- grnd_truth
  misg_data[temp$misg_idx] <- NA
  
  res <- CAM_smpClus_iter(
    misg_data,
    dim.rdc = 10,
    cluster.num = 35,
    optimal.source.number = rank,
    normalization = "none",
    cheat_data = grnd_truth,
    iter = 3
  )
  err <- c()
  for (i in 2:10) {
    err <- c(err, data_with_MV_NRMSE(grnd_truth, misg_data, res$res[[i]]))
  }
  plot(err)
  
  predicted <- impute(misg_data,
                      "nipals",
                      parameter = rank)$res
  print(data_with_MV_NRMSE(grnd_truth, misg_data, predicted))
}
