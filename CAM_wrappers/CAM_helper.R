#' Non negative least squares, given A and X, solve for S, where X = A %*% S
#'
#' @param A A matrix.
#' @param X A matrix.
#' @return S matrix.
CAM_helper_nnls <- function(A, X) {
  S <- matrix(0, nrow = ncol(A), ncol = ncol(X))
  for (i in 1:ncol(X)) {
    B <- A
    b <- X[, i]
    removeInd <- which(is.na(b))
    if (length(removeInd) != 0) {
      B <- B[-removeInd,]
      b <- b[-removeInd]
    }
    S[, i] <- nnls(B, b)$X
  }
  return(S)
}

#' Wrapper function for CAM
#'
#' @param dataT A matrix.
#' @param dim.rdc Reduce the dimension to dim.rdc.
#' @param cluster.num Number of clusters.
#' @return CAM result object (k from 2 to dim.rdc).
CAM_helper <- function(dataT,
                       dim.rdc = 20,
                       cluster.num = 50) {
  startTime <- Sys.time()
  
  rPrepLL <-
    CAMPrep3(
      dataT,
      dim.rdc = dim.rdc,
      thres.low = 0,
      thres.high = 1,
      cluster.method = "K-Means",
      cluster.num = cluster.num,
      MG.num.thres = 5,
      lof.thres = 0,
      quickhull = FALSE
    )
  
  rPrepL <-
    reClustering(
      rPrepLL,
      cluster.num = cluster.num,
      MG.num.thres = 5,
      quickhull = FALSE
    )
  
  rClusLis <- vector('list', dim.rdc)
  for (i in 2:3) {
    rClusLis[[i]] <- CAMMGCluster(i, rPrepL)
  }
  
  show('Clustering done.')
  
  bInd <-
    c(which(rClusLis[[2]]@corner[1, 1] == colnames(rPrepL@centers)),
      which(rClusLis[[2]]@corner[1, 2] == colnames(rPrepL@centers)))
  
  cInd <-
    c(which(rClusLis[[2]]@corner[2, 1] == colnames(rPrepL@centers)),
      which(rClusLis[[2]]@corner[2, 2] == colnames(rPrepL@centers)))
  
  
  rClusG <-
    CAMMGClusterG(dim.rdc,
                  rPrepL,
                  dataT,
                  bInd = bInd,
                  cInd = cInd)
  
  rClusGDim <- vector('list', dim.rdc)
  for (i in 2:dim.rdc) {
    rClusGDim[[i]] <- new(
      "CAMMGObj",
      idx = 0,
      corner = rbind(colnames(rPrepL@centers)[rClusG@idx[[1]][[i]]],
                     colnames(rPrepL@centers)[rClusG@idx[[3]][[i]]]),
      error = matrix(0, 1, 1)
    )
  }
  
  
  rEstDim <- vector('list', dim.rdc)
  for (i in 2:dim.rdc) {
    rEstDim[[i]] <-
      CAMASestT(rClusGDim[[i]], rPrepL, dataT, 2)
  }
  
  endTime <- Sys.time()
  show(endTime - startTime)
  
  return(rEstDim)
}
