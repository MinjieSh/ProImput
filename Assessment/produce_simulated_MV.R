get_misg_rate <- function(x) {
  return (100 * sum(is.na(x)) / length(x))
}

is_almost_empty_matrix <- function(mat) {
  protein_misg_rate <- apply(mat, 2, get_misg_rate)
  return(which(protein_misg_rate > 100 * (nrow(mat) - 5) / nrow(mat)))
}

is_almost_empty_vector <- function(vec) {
  return(get_misg_rate(vec) > 100 * (length(vec) - 5) / length(vec))
}

data_with_MV_mask_out <- function(data_with_MV,
                                  option = c("MCAR", "MAR", "MNAR"),
                                  parameter,
                                  keep = c("max", "random", "remove (not recommend)", "none")) {
  show("Generating simulation data with set-aside masked values from the full data matrix...")
  
  if (option == "MCAR") {
    show("Option 1 - MCAR.")
    data_extra_MV <- data_with_MV_MCAR(data_with_MV,
                                       extra_MV_percentage = parameter)
  } else if (option == "MAR") {
    show("Option 2 - MAR.")
    data_extra_MV <- data_with_MV_MAR(data_with_MV,
                                      coefficient = parameter)
  } else if (option == "MNAR") {
    show("Option 3 - MNAR.")
    data_extra_MV <- data_with_MV_MNAR(data_with_MV,
                                       quantile = parameter)
  }
  
  modify <- is_almost_empty_matrix(data_extra_MV)
  
  if (length(modify) != 0)
    show("Found almost empty proteins after simulation. Adjusted...")
  
  switch(
    keep,
    
    "max" = {
      show("Keep the 5 largest data points.")
      data_extra_MV[, modify] <- data_with_MV[, modify]
      for (index in modify) {
        threshold <- sort(data_with_MV[, index], decreasing = TRUE)[[5]]
        data_extra_MV[which(data_with_MV[, index] < threshold), index] <-
          NA
      }
    },
    
    
    "random" = {
      show("Keep random 5 data points.")
      data_extra_MV[, modify] <- data_with_MV[, modify]
      for (index in modify) {
        non_misg <- which(!is.na(data_with_MV[, index]))
        if (length(non_misg) >= 5) {
          keep <- sample(non_misg, 5)
          data_extra_MV[-keep, index] <- NA
        }
      }
    },
    
    "remove (not recommend)" = {
      show("Remove those almost empty proteins after simulation.")
      data_extra_MV <- data_extra_MV[, -modify]
      data_with_MV <- data_with_MV[, -modify]
      return(list(data_extra_MV = data_extra_MV, data_with_MV = data_with_MV))
    },
    
    "none" = {
      show("Do nothing.")
    }
  )
  
  extra_MV_indices <- setdiff(which(is.na(data_extra_MV)),
                              which(is.na(data_with_MV)))
  
  extra_MV_percentage <- 100 * length(extra_MV_indices) /
    (ncol(data_with_MV) * nrow(data_with_MV))
  
  show(paste(
    "Total misg rate increases by (%):",
    round(extra_MV_percentage, 3),
    sep = " "
  ))
  
  
  total_MV_percentage <- 100 * sum(is.na(data_extra_MV)) /
    (ncol(data_with_MV) * nrow(data_with_MV))
  show(paste(
    "Current total misg rate is (%):",
    round(total_MV_percentage, 3),
    sep = " "
  ))
  
  return(list(data_extra_MV = data_extra_MV,
              extra_MV_indices = extra_MV_indices))
}

data_with_MV_MCAR <- function(orig_data, extra_MV_percentage) {
  size <- nrow(orig_data) * ncol(orig_data)
  
  replace <- sample(size, round(extra_MV_percentage * size))
  orig_data[replace] <- NA
  
  return(orig_data)
}


data_with_MV_MAR_helper <- function(x, coefficient) {
  if (is_almost_empty_vector(x))
    return(x)
  
  percentage <- (sum(is.na(x)) / length(x)) * coefficient
  
  if (percentage >= 1)
    percentage <- 1
  
  size <- length(x)
  replace <- sample(size, round(percentage * size))
  x[replace] <- NA
  
  return(x)
  
}

data_with_MV_MAR <- function(orig_data, coefficient) {
  result <- apply(orig_data, 2, function(x)
    data_with_MV_MAR_helper(x, coefficient))
  return(result)
}


data_with_MV_MNAR <- function(orig_data, quantile) {
  misg_indices <- which(is.na(orig_data))
  
  nSamples = nrow(orig_data)
  nProt = ncol(orig_data)
  
  threshold_mean <- quantile(orig_data, quantile, na.rm = TRUE)
  
  threshold_matrix = matrix(rnorm(nSamples * nProt, threshold_mean, 0.1),
                            ncol = nProt,
                            nrow = nSamples)
  threshold_matrix[misg_indices] <- NA
  
  indices.MNAR = which(orig_data < threshold_matrix)
  
  orig_data[indices.MNAR] <- NA
  return(orig_data)
}
