#' Pairwise (Non-missing Element) Cosine Similarity
cos.sim <- function(X) {
  n <- nrow(X)
  idx_comb <- expand.grid(i = 1:n, j = 1:n)
  
  return(matrix(apply(idx_comb, 1, function(r)
    cos.sim.helper(r, X)),
    n,
    n))
}

cos.sim.helper <- function(idx, X) {
  A <- X[idx[[1]],]
  B <- X[idx[[2]],]
  
  misg_idx <- union(which(is.na(A)), which(is.na(B)))
  
  if (length(misg_idx) != 0) {
    A <- A[-misg_idx]
    B <- B[-misg_idx]
  }
  
  if (sqrt(sum(A ^ 2) * sum(B ^ 2)) == 0)
    return(1)
  
  return(sum(A * B) / sqrt(sum(A ^ 2) * sum(B ^ 2)))
}

#' Fused Regularization Matrix Factorization
#' @param X - data (with missing values) to be imputed.
#' @param rank - rank of X (dimension of A and S)
#' @param lrn_rate - learning rate
#' @param max_step - maximum number of steps
#' @param conv_thld - converge threshold
#' @param lam_A - sparsity regularization term coefficient for A
#' @param lam_S  - sparsity regularization term coefficient for S
#' @param nbr_info_mat - additional information to help impute X
#' @param nbr_thld - cosine similarity lower limit for "good neighbors"
#' @param coeff_crs_rglr - cross regularization term coefficient
#' @return An imputed matrix.
FRMF <- function(X,
                 rank = 3,
                 lrn_rate = 0.0005,
                 max_step = 5000,
                 conv_thld = 0.05,
                 lam_A = 0.005,
                 lam_S = 0.005,
                 nbr_info_mat = NULL,
                 nbr_thld = 0.9,
                 coeff_crs_rglr = 0,
                 normalization = c("max_norm", "std_score", "none")) {
  
  if (is.null(nbr_info_mat)) {
    print("No external information incorporated.")
    nbr_cos <- NULL
  } else{
    print("External information incorporated.")
    nbr_cos <- cos.sim(nbr_info_mat)
  }
  
  #### Preprocess X
  misg_idx <- which(is.na(X))
  
  # backup
  X_orig <- X 
  
  if (normalization == "max_norm") {
    col.max <- apply(X, 2, function(x)
      max(x, na.rm = TRUE))
    X <- apply(X, 2, function(x)
      x / max(x, na.rm = TRUE))
  }
  
  if (normalization == "std_score") {
    X_norm <- scale(X, center = TRUE, scale = TRUE)
    X <- X_norm
  }
  
  #### Low-Rank Matrix Factorization
  # row/sample
  smp <- dim(X)[1]
  # column/protein
  prt <- dim(X)[2]
  
  #### Randomly assign initial values to A and S
  A <- matrix(runif(smp * rank, 0, 1), smp, rank)
  S <- matrix(runif(rank * prt, 0, 1), rank, prt)
  
  # Indicator Function
  # 0 - missing
  # 1 - observed
  indic_func <- as.matrix(ifelse(is.na(X), 0, 1))
  X[is.na(X)] <- 0
  
  # initial loss
  loss <- 100000
  
  for (step in 1:max_step) {
    temp_A <- A
    temp_S <- S
    err <- (temp_A %*% temp_S - X) * indic_func
    
    for (i in 1:smp) {
      grad_A <- colSums(err[i, ] * t(temp_S)) + lam_A * temp_A[i, ]
      A[i, ] <- A[i, ] - lrn_rate * grad_A
    }
    
    for (j in 1:prt) {
      grad_S <- colSums(temp_A * err[, j]) + lam_S * temp_S[, j]
      S[, j] <- S[, j] - lrn_rate * grad_S
    }
    
    loss_orig <- loss
    loss <- 0
    
    if (!is.null(nbr_info_mat)) {
      nbrh_idx <- vector('list', smp)
      nbrh_size <- vector('list', smp)
      
      for (i in 1:smp) {
        nbrh_idx[[i]] <- which(nbr_cos[i, ] >= nbr_thld)
        nbrh_idx[[i]] <- nbrh_idx[[i]][nbrh_idx[[i]] != i]
        
        nbrh_size[[i]] <- length(nbrh_idx[[i]])
      }
      
      # update A
      for (i in 1:smp) {
        A_i_rep <-
          matrix(rep(temp_A[i, ], each = nbrh_size[[i]]),
                 nrow = nbrh_size[[i]],
                 ncol = rank)
        cost <- A_i_rep - temp_A[nbrh_idx[[i]], ]
        
        A[i, ] <-
          A[i, ] - lrn_rate * coeff_crs_rglr * colSums(cost)
      }
      
      # update loss
      crs_rglr <- 0
      for (i in 1:smp) {
        A_i_rep <-
          matrix(rep(temp_A[i, ], each = nbrh_size[[i]]),
                 nrow = nbrh_size[[i]],
                 ncol = rank)
        
        diff_norm <- col.norm(t(A_i_rep - A[nbrh_idx[[i]], ]))
        
        crs_rglr <- crs_rglr + sum(diff_norm)
      }
      loss <- loss + coeff_crs_rglr / 2 * crs_rglr
    }
    
    err <- (A %*% S - X) * indic_func
    
    loss <- loss + 1 / 2 * sum(err ^ 2) +
      (lam_A) / 2 * norm(A, type = c("F")) +
      (lam_S) / 2 * norm(S, type = c("F"))
    
    if (abs(loss_orig - loss) < conv_thld) {
      print("Converge.")
      break
    }
    
    if (abs(loss_orig - loss) > 100000) {
      print("Diverge.")
      break
    }
    
  }
  
  print(paste("Steps used:", step, sep = ""))
  
  recons <- A %*% S
  
  if (normalization == "max_norm") {
    for (col in 1:ncol(recons)) {
      recons[, col] <- recons[, col] * col.max[col]
    }
  }
  
  if (normalization == "std_score") {
    recons <- t(apply(recons, 1,
                      function(r)
                        (
                          r * attr(X_norm, 'scaled:scale') +
                            attr(X_norm, 'scaled:center')
                        )))
  }
  
  X_orig[misg_idx] <- recons[misg_idx]
  
  return(list(
    imp = X_orig,
    recons = recons,
    A = A,
    S = S,
    cos = nbr_cos
  ))
}
