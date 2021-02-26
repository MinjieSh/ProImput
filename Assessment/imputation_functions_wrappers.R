impute <- function(log_data = NULL,
                   methods = c(
                     "pwHalfMin_log",
                     "pwMean",
                     "swKNN",
                     "pwKNN",
                     "PPCA",
                     "NIPALS",
                     "svdImpute",
                     "SVT",
                     
                     "FRMF",
                     
                     "CAM_NIPALS",
                     "CAM_SVT",
                     "CAM_cmplt",
                     "CAM_smpClus",
                     "CAM_smpClus_iter"
                     
                   ),
                   parameter = NULL,
                   orig_data = NULL,
                   nPCs = 10,
                   # CAM_NIPALS
                   lambda = 40000,
                   # CAM_SVT
                   dim.rdc = 20,
                   cluster.num = 50,
                   nbr_info_mat = NULL,
                   nbr_thld = 0.9) {
  
  
  
  
  if(!is.null(log_data)){
    log_data <- as.matrix(log_data)
    colnames(log_data) <- 1:ncol(log_data)
    data_t <- t(log_data)
  }
    
  if(!is.null(orig_data)){
    orig_data <- as.matrix(orig_data)
    colnames(orig_data) <- 1:ncol(orig_data)
  }
    
  
  temp_res <- list()
  
  switch(
    methods,
    
    
    "pwHalfMin_log" = {
      show("pwHalfMin (log)")
      imp <- half_min_log(log_data)
    },
    
    
    "PPCA" = {
      show("PPCA")
      imp <- completeObs(pca(
        object = log_data,
        method = "ppca",
        nPcs = parameter,
        center = TRUE,
        scale = "uv"
      ))
    },
    
    "NIPALS" = {
      show("NIPALS")
      imp <- completeObs(pca(
        object = log_data,
        method = "nipals",
        nPcs = parameter,
        center = TRUE,
        scale = "uv"
      ))
    },
    
    "SVD" = {
      show("SVD")
      imp <- completeObs(
        pca(
          object = log_data,
          method = "svdImpute",
          nPcs = parameter,
          center = FALSE,
          scale = "none"
        )
      )
    },
    
    "pwMean" = {
      show("pwMean")
      imp <- meanImpute(log_data)$x
    },
    
    
    "swKNN" = {
      show("Sample-wise kNN")
      imp <- KNN(data_t, k = parameter)
      
    },
    
    "pwKNN" = {
      show("Protein-wise kNN")
      imp <- t(KNN(log_data, k = parameter))
    },
    
    "SVT" = {
      show("SVT")
      temp <- SVTImpute(log_data, lambda = parameter,
                        max.iters = 1500)$x
      imp <- log_data
      imp[which(is.na(log_data))] <-
        temp[which(is.na(log_data))]
      
    },
    
    "FRMF" = {
      show("FRMF")
      temp_res <- FRMF(
        log_data,
        rank = 3,
        lrn_rate = 0.0005,
        max_step = 10000,
        conv_thld = 0.005,
        lam_A = 0.005,
        lam_S = 0.5,
        nbr_info_mat = nbr_info_mat,
        nbr_thld = 0.9,
        coeff_crs_rglr = parameter,
        normalization = "std_score"
      )
      
      imp <- temp_res$imp
    },
    
    "CAM_cmplt" = {
      show("CAM_cmplt")
      temp_res <- CAM_cmplt(
        orig_data,
        dim.rdc = dim.rdc,
        cluster.num = cluster.num,
        normalization = "std_score"
      )
      
      imp <- temp_res$res[[parameter]]
      
    },
    
    "CAM_NIPALS" = {
      show("CAM_NIPALS")
      temp_res <- CAM_NIPALS(
        orig_data,
        dim.rdc = dim.rdc,
        cluster.num = cluster.num,
        nPCs = parameter,
        normalization = "std_score"
      )
      
      imp <- temp_res$res[[parameter]]
    },
    
    
    "CAM_SVT" = {
      show("CAM_SVT")
      temp_res <- CAM_SVT(
        orig_data,
        dim.rdc = dim.rdc,
        cluster.num = cluster.num,
        lambda = lambda,
        normalization = "std_score"
      )
      
      imp <- temp_res$res[[parameter]]
    },
    
    "CAM_smpClus" = {
      show("CAM_smpClus")
      temp_res <- CAM_smpClus(
        orig_data,
        dim.rdc = dim.rdc,
        cluster.num = cluster.num,
        normalization = "none"
      )
      
      imp <- temp_res$res[[parameter]]
    },
    
    "CAM_smpClus_iter" = {
      show("CAM_smpClus_iter")
      temp_res <- CAM_smpClus_iter(
        orig_data,
        dim.rdc = dim.rdc,
        cluster.num = cluster.num,
        optimal.source.number = 20,
        normalization = "none",
        iter = 5
      )
      
      imp <- temp_res$res[[parameter]]
      
    }
    
  )
  
  return(list(imp = imp, more = temp_res))
}

KNN <- function(data_t, k) {
  data_norm_t <- scale(data_t, center = TRUE, scale = TRUE)
  result_norm <- kNNImpute(t(data_norm_t), k, verbose = FALSE)$x
  temp_res <- apply(t(result_norm), 1,
                    function(r)
                      (
                        r * attr(data_norm_t, 'scaled:scale') +
                          attr(data_norm_t, 'scaled:center')
                      ))
  return(temp_res)
}

half_min_log <- function(log_data) {
  imp <- log_data
  colHalfMin <- apply(log_data, 2,
                      function(x)
                        min(x, na.rm = TRUE) - 1)
  index <- which(is.na(log_data), arr.ind = TRUE)
  imp[index] <- colHalfMin[index[, 2]]
  return(imp)
}

half_min <- function(orig_data) {
  imp <- orig_data
  colHalfMin <- apply(orig_data, 2,
                      function(x)
                        min(x, na.rm = TRUE) / 2)
  index <- which(is.na(orig_data), arr.ind = TRUE)
  imp[index] <- colHalfMin[index[, 2]]
  return(imp)
}