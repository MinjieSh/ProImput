#' CAM imputation, using proteins with no missing values.
#'
#' @param misg_data A matrix in original space.
#' @param dim.rdc Reduce the dimension to dim.rdc.
#' @param cluster.num Number of clusters.
#' @return Imputed result matrix and A matrix (k from 2 to dim.rdc).
CAM_cmplt <- function(misg_data,
                      dim.rdc = 20,
                      cluster.num = 50,
                      normalization = c("total_count", "std_score", "none")) {
  
  colnames(misg_data) <- 1:ncol(misg_data)
  
  misg_idx <- which(is.na(misg_data))
  misg_rate <- apply(misg_data, 2, get_misg_rate)
  cmplt_data <- misg_data[, which(misg_rate == 0)]
  
  if (normalization == "total_count") {
    X_norm <- t(apply(cmplt_data, 1, function(x)
      x / sum(x, na.rm = TRUE)))
    X_norm <- X_norm * 1000000
  }
  
  if (normalization == "std_score") {
    X_norm <- scale(cmplt_data, center = FALSE, scale = TRUE)
  }
  
  if (normalization == "none") {
    X_norm <- cmplt_data
  }
  
  rCAM <- CAM_helper(data = t(X_norm),
                     dim.rdc = dim.rdc,
                     cluster.num = cluster.num)
  
  show("CAM (in original space) done.")
  
  res <- vector('list', dim.rdc)
  S <- vector('list', dim.rdc)
  
  if (normalization == "total_count") {
    rec_fac <- rowSums(cmplt_data, na.rm = TRUE)
    
    misg_data_norm <- misg_data
    for (j in 1:nrow(misg_data)) {
      misg_data_norm[j, ] <- misg_data[j, ] / rec_fac[j]
    }
    
    misg_data_norm <- misg_data_norm * 1000000
  }
  
  if (normalization == "std_score") {
    scale_factor <- rep(NA, length(misg_rate))
    scale_factor[which(misg_rate == 0)] <-
      attr(X_norm, 'scaled:scale')
    
    misg_data_norm <- misg_data
    for (j in 1:ncol(misg_data_norm)) {
      if (is.na(scale_factor[j])) {
        temp <- scale(misg_data_norm[, j], center = FALSE, scale = TRUE)
        scale_factor[j] <- attr(temp, 'scaled:scale')
        misg_data_norm[, j] <- temp
      } else {
        misg_data_norm[, j] <- misg_data_norm[, j] / scale_factor[j]
      }
    }
  }
  
  for (i in 2:dim.rdc) {
    A <- rCAM[[i]]@Aest
    
    S[[i]] <- CAM_helper_nnls(A, misg_data_norm)
    
    X <- A %*% S[[i]]
    
    
    if (normalization == "total_count") {
      X <- X / 1000000
      for (j in 1:nrow(X)) {
        X[j, ] <- X[j, ] * rec_fac[j]
      }
    }
    
    if (normalization == "std_score") {
      for (j in 1:ncol(X)) {
        X[, j] <- X[, j] * scale_factor[j]
      }
    }
    
    res[[i]] <- misg_data
    res[[i]][misg_idx] <- X[misg_idx]
  }
  
  show("Return imputed result in original space.")
  
  return(list(res = res, S = S))
}