#' Iterative CAM imputation, replace PCA with sample clustering (missing data accommodated).
#'
#' @param misg_data A matrix in original space.
#' @param dim.rdc Reduce the dimension to dim.rdc.
#' @param cluster.num Number of clusters.
#' @return Imputed result matrix and A matrix (optimal.source.number from 2 to dim.rdc).
CAM_smpClus_iter <- function(misg_data,
                             dim.rdc = 20,
                             cluster.num = 50,
                             optimal.source.number = 5,
                             normalization = c("total_count", "std_score", "none"),
                             cheat_data = NULL,
                             iter = 2) {
  
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
  
  dataT <- t(X_norm)
  misg_idx <- which(is.na(dataT))
  
  rPrepLL <- CAMPrep3_impu(
    dataT,
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
    CAMMGClusterG(dim.rdc, rPrep, dataT, bInd = bInd, cInd = cInd)
  
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
  
  rAsest <-
    CAMASest_impu(rClusGDim[[optimal.source.number]], rPrep, dataT)
  
  result_t <- rAsest@Sest %*% t(rAsest@Aest)
  
  dataT2 <- dataT
  dataT2[misg_idx] <- result_t[misg_idx]
  
  predicted <- t(dataT2)
  if (normalization == "total_count") {
    predicted <- predicted / 1000000
    for (j in 1:nrow(predicted)) {
      predicted[j,] <- X[j,] * rec_fac[j]
    }
  }
  
  if (normalization == "std_score") {
    predicted <- t(apply(predicted, 1,
                         function(r)
                           (r * attr(
                             X_norm, 'scaled:scale'
                           ))))
    
  }
  if(!is.null(cheat_data))
    show(nrmse_full_data(cheat_data, misg_data, predicted))
  
  for (n in 1:iter) {
    show(paste("Iteration: ", n, ".", sep = ""))
    
    rPrepLL2 <-
      CAMPrep3_impu(
        dataT2,
        dim.rdc = dim.rdc,
        cluster.num = cluster.num,
        MG.num.thres = 5
      )
    
    rPrep2 <-
      reClustering(
        rPrepLL2,
        cluster.num = cluster.num,
        MG.num.thres = 5,
        quickhull = FALSE
      )
    
    rClusLis2 <- vector('list', dim.rdc)
    
    for (i in 2:3) {
      rClusLis2[[i]] <- CAMMGCluster(i, rPrep2)
    }
    
    show('Clustering done')
    
    bInd2 <-
      c(which(rClusLis2[[2]]@corner[1, 1] == colnames(rPrep2@centers)),
        which(rClusLis2[[2]]@corner[1, 2] == colnames(rPrep2@centers)))
    
    cInd2 <-
      c(which(rClusLis2[[2]]@corner[2, 1] == colnames(rPrep2@centers)),
        which(rClusLis2[[2]]@corner[2, 2] == colnames(rPrep2@centers)))
    
    
    rClusG2 <-
      CAMMGClusterG(dim.rdc, rPrep2, dataT2, bInd = bInd2, cInd = cInd2)
    
    rClusGDim2 <- vector('list', dim.rdc)
    
    for (i in 2:dim.rdc) {
      rClusGDim2[[i]] <- new(
        "CAMMGObj",
        idx = 0,
        corner = rbind(colnames(rPrep2@centers)[rClusG2@idx[[1]][[i]]],
                       colnames(rPrep2@centers)[rClusG2@idx[[3]][[i]]]),
        error = matrix(0, 1, 1)
      )
    }
    
    rAsest2 <-
      CAMASest_impu(rClusGDim2[[optimal.source.number]], rPrep2, dataT)
    
    dataT3 <- rAsest2@Sest %*% t(rAsest2@Aest)
    
    dataT2 <- dataT
    dataT2[misg_idx] <- dataT3[misg_idx]
    
    predicted <- t(dataT2)
    
    if (normalization == "total_count") {
      predicted <- predicted / 1000000
      for (j in 1:nrow(X)) {
        predicted[j,] <- X[j,] * rec_fac[j]
      }
    }
    
    if (normalization == "std_score") {
      predicted <- t(apply(predicted, 1,
                           function(r)
                             (r * attr(
                               X_norm, 'scaled:scale'
                             ))))
      
    }
    if(!is.null(cheat_data))
      show(nrmse_full_data(cheat_data, misg_data, predicted))
  }
  
  res <- vector('list', dim.rdc)
  A <- vector('list', dim.rdc)
  
  misg_idx <- which(is.na(misg_data))
  
  for (i in 2:dim.rdc) {
    rAsest <- CAMASest_impu(rClusGDim2[[i]], rPrep2, dataT)
    
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
  
  return(list(res = res, A = A))
}
