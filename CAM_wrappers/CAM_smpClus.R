#' CAM imputation, replace PCA with sample clustering (missing data accommodated).
#'
#' @param misg_data A matrix in original space.
#' @param dim.rdc Reduce the dimension to dim.rdc.
#' @param cluster.num Number of clusters.
#' @return Imputed result matrix and A matrix (optimal.source.number from 2 to dim.rdc).
CAM_smpClus <- function(misg_data,
                        dim.rdc = 20,
                        cluster.num = 50,
                        normalization = c("total_count", "std_score", "none")) {
  
  startTime <- Sys.time()
  
  colnames(misg_data) <- 1:ncol(misg_data)
  
  if (normalization == "total_count") {
    rec_fac <- rowSums(misg_data, na.rm = TRUE)
    
    X_norm <- t(apply(misg_data, 1, function(x)
      x / sum(x, na.rm = TRUE)))
    X_norm <- X_norm * 1000000
  }
  
  if (normalization == "std_score") {
    X_norm <- scale(misg_data, center = FALSE, scale = TRUE)
  }
  
  if (normalization == "none") {
    X_norm <- misg_data
  }
  
  X_norm_T <- t(X_norm)
  misg_idx <- which(is.na(misg_data))
  
  rPrepLL <- CAMPrep3_impu(
    X_norm_T,
    dim.rdc = dim.rdc,
    cluster.num = cluster.num,
    MG.num.thres = 5
  )
  
  rPrep <-
    reClustering(
      rPrepLL,
      cluster.num = cluster.num,
      MG.num.thres = 5,
      quickhull = FALSE
    )
  
  rClusLis <- vector('list', dim.rdc)
  
  for (i in 2:3) {
    rClusLis[[i]] <- CAMMGCluster(i, rPrep)
  }
  
  show('Clustering done')
  
  bInd <-
    c(which(rClusLis[[2]]@corner[1, 1] == colnames(rPrep@centers)),
      which(rClusLis[[2]]@corner[1, 2] == colnames(rPrep@centers)))
  
  cInd <-
    c(which(rClusLis[[2]]@corner[2, 1] == colnames(rPrep@centers)),
      which(rClusLis[[2]]@corner[2, 2] == colnames(rPrep@centers)))
  
  rClusG <-
    CAMMGClusterG(dim.rdc, rPrep, X_norm_T, bInd = bInd, cInd = cInd)
  
  rClusGDim <- vector('list', dim.rdc)
  
  for (i in 2:dim.rdc) {
    rClusGDim[[i]] <- new(
      "CAMMGObj",
      idx = 0,
      corner = rbind(colnames(rPrep@centers)[rClusG@idx[[1]][[i]]],
                     colnames(rPrep@centers)[rClusG@idx[[3]][[i]]]),
      error = matrix(0, 1, 1)
    )
  }
  
  res <- vector('list', dim.rdc)
  A <- vector('list', dim.rdc)
  
  for (i in 2:dim.rdc) {
    rAsest <- CAMASest_impu(rClusGDim[[i]], rPrep, X_norm_T)
    
    A[[i]] <- rAsest@Aest
    
    X <- t(rAsest@Sest %*% t(A[[i]]))
    
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
  
  endTime <- Sys.time()
  show(endTime - startTime)
  
  return(list(res = res, A = A))
}
