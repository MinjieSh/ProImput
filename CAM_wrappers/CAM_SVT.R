#' CAM imputation, initialized by SVT.
#'
#' @param misg_data A matrix in original space.
#' @param dim.rdc Reduce the dimension to dim.rdc.
#' @param cluster.num Number of clusters.
#' @param lambda The penalty on the singular values for SVT.
#' @return Imputed result matrix and A matrix (k from 2 to dim.rdc).
CAM_SVT <- function(misg_data,
                    dim.rdc = 20,
                    cluster.num = 50,
                    lambda = 40000,
                    normalization = c("total_count", "std_score", "none")) {
  colnames(misg_data) <- 1:ncol(misg_data)
  
  temp <- 2 ^ SVTImpute(log2(misg_data),
                             lambda = lambda,
                             max.iters = 500)$x
  
  misg_idx <- which(is.na(misg_data))
  
  init_data <- misg_data
  init_data[misg_idx] <- temp[misg_idx]
  
  
  show("SVT (in log2 space) done.")
  
  if (normalization == "total_count") {
    X_norm <- t(apply(init_data, 1, function(x)
      x / sum(x, na.rm = TRUE)))
    X_norm <- X_norm * 1000000
  }
  
  if (normalization == "std_score") {
    X_norm <- scale(init_data, center = FALSE, scale = TRUE)
  }
  
  if (normalization == "none") {
    X_norm <- init_data
  }
  
  rCAM <- CAM_helper(
    data = t(X_norm),
    dim.rdc = dim.rdc,
    cluster.num = cluster.num
  )
  
  show("CAM (in original space) done.")
  
  res <- vector('list', dim.rdc)
  A <- vector('list', dim.rdc)
  
  if (normalization == "total_count") {
    rec_fac <- rowSums(init_data, na.rm = TRUE)
    
    misg_data_norm <- misg_data
    for (j in 1:nrow(misg_data)) {
      misg_data_norm[j,] <- misg_data[j,] / rec_fac[j]
    }
    
    misg_data_norm <- misg_data_norm * 1000000
  }
  
  if (normalization == "std_score") {
    misg_data_norm <- t(apply(misg_data, 1,
                              function(r)
                                (r / attr(
                                  X_norm, 'scaled:scale'
                                ))))
  }
  
  for (i in 2:dim.rdc) {
    S <- rCAM[[i]]@Sest
    
    A[[i]] <- t(CAM_helper_nnls(S, t(misg_data_norm)))
    
    X <- A[[i]] %*% t(S)
    
    if (normalization == "total_count") {
      X <- X / 1000000
      for (j in 1:nrow(X)) {
        X[j,] <- X[j,] * rec_fac[j]
      }
    }
    
    if (normalization == "std_score") {
      X <- t(apply(X, 1,
                   function(r)
                     (r * attr(
                       X_norm, 'scaled:scale'
                     ))))
      
    }
    
    res[[i]] <- misg_data
    res[[i]][misg_idx] <- X[misg_idx]
    
  }
  
  show("Return imputed result in original space.")
  
  return(list(
    res = res,
    A = A,
    init_data = init_data
  ))
}