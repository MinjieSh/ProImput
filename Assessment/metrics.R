data_with_MV_NRMSE <-
  function(data_with_MV,
           data_extra_MV,
           data_imputed) {
    extra_MV_indices <- setdiff(which(is.na(data_extra_MV)),
                                which(is.na(data_with_MV)))
    
    data_with_MV <- scale(data_with_MV, center = TRUE, scale = TRUE)
    center <- attr(data_with_MV, 'scaled:center')
    scalar <- attr(data_with_MV, 'scaled:scale')
    
    data_imputed <- t(apply(data_imputed, 1,
                            function(r)
                              ((r - center) / scalar)))
    
    error <-
      data_imputed[extra_MV_indices] - data_with_MV[extra_MV_indices]
    
    var <- sd(data_with_MV[extra_MV_indices]) ^ 2
    res <-
      sqrt(sum((error) ^ 2, na.rm = TRUE) / (length(extra_MV_indices) * var))
    
    return(res)
  }



data_with_MV_RMSE <-
  function(data_with_MV,
           data_extra_MV,
           data_imputed) {
    extra_MV_indices <- setdiff(which(is.na(data_extra_MV)),
                                which(is.na(data_with_MV)))
    error <-
      data_imputed[extra_MV_indices] - data_with_MV[extra_MV_indices]
    res <-
      sqrt(sum((error) ^ 2, na.rm = TRUE) / length(extra_MV_indices))
    return(res)
  }



data_with_MV_SOR_helper <- function(data_with_MV,
                                    data_extra_MV,
                                    data_imputed) {
  extra_MV_indices <- setdiff(which(is.na(data_extra_MV)),
                              which(is.na(data_with_MV)))
  if (length(extra_MV_indices) < 2)
    return(0)
  
  error <-
    data_imputed[extra_MV_indices] - data_with_MV[extra_MV_indices]
  var <- sd(data_with_MV[extra_MV_indices]) ^ 2
  res <-
    sqrt(sum((error) ^ 2, na.rm = TRUE) / (length(extra_MV_indices) * var))
  return(res)
}



data_with_MV_SOR <-
  function(data_with_MV,
           data_extra_MV,
           data_imputed) {
    data_with_MV <- scale(data_with_MV, center = TRUE, scale = TRUE)
    center <- attr(data_with_MV, 'scaled:center')
    scalar <- attr(data_with_MV, 'scaled:scale')
    
    data_imputed <- t(apply(data_imputed, 1,
                            function(r)
                              ((r - center) / scalar)))
    
    res <- matrix(0, nrow = ncol(data_with_MV), ncol = 1)
    for (i in 1:ncol(data_with_MV)) {
      res[i,] <-
        data_with_MV_SOR_helper(data_with_MV[, i], data_extra_MV[, i], data_imputed[, i])
    }
    return(res)
  }
