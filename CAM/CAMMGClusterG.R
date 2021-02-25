#' MG cluster detection for CAM
#'
#' This function finds corner clusters as MG clusters
#' (clusters containing marker genes).
#' @param K The candidate subpopulation number.
#' @param PrepResult An object of class "\code{\link{CAMPrepObj}}" obtained
#'     from \code{\link{CAMPrep}} function.
#' @param generalNMF If TRUE, the decomposed proportion matrix has no sum-to-one
#'     constraint for each row. Without this constraint, the scale ambiguity of
#'     corner cluster centers will not be removed when computing
#'     reconstruction errors.
#'     The default is FALSE.
#' @param nComb The number of possible combinations of clusters as corner
#'     clusters. Within these possible combinations ranked by margin errors,
#'     we can further select the best one based on reconstruction errors.
#'     The default is 200.
#' @details This function is used internally by \code{\link{CAM}} function to
#' detect clusters containing marker genes,
#' or used when you want to perform CAM step by step.
#'
#' This function provides two solutions. The first is the combination of
#' clusters yielding the minimum sum of margin-of-errors for cluster centers.
#' In the second, nComb possible combinations are selected by ranking sum of
#' margin-of-errors for cluster centers. Then the best one is selected based on
#' reconstruction errors of all data points in original space.
#' @return An object of class "\code{\link{CAMMGObj}}" containing the
#' following components:
#' \item{idx}{Two numbers which are two solutions' ranks by sum of
#'     margin-of-error.}
#' \item{corner}{The indexes of clusters as detected corners. Each row is a
#'     solution.}
#' \item{error}{Two rows. The first row is sum of margin-of-errors for nComb
#'     possible combinations. The second row is reconstruction errors for
#'     nComb possible combinations.}
#' @export
#' @examples
#' #obtain data
#' data(ratMix3)
#' data <- ratMix3$X
#'
#' #preprocess data
#' rPrep <- CAMPrep(data, dim.rdc = 3, thres.low = 0.30, thres.high = 0.95)
#'
#' #Marker gene cluster detection with a fixed K = 3
#' rMGC <- CAMMGCluster(3, rPrep)
CAMMGClusterG <-
    function(K,
             PrepResult,
             generalNMF = FALSE,
             bInd = NULL,
             cInd = NULL) {
        X <- PrepResult@centers
        if (K > nrow(X) || K > ncol(X)) {
            warning("Return NULL for K = ", K, " larger than ", max(dim(X)))
            return(NULL)
        }
        if (K == ncol(X)) {
            return(new(
                "CAMMGObj",
                idx = c(1, 1),
                corner = matrix(as.integer(colnames(X)), 1)[c(1, 1), ]
            ))
        }
        
        scaleRecover <- !generalNMF
        Xall <- PrepResult@Xprep[, !(PrepResult@cluster$cluster %in%
                                         PrepResult@c.outlier)]
        #X[1, ] <- X[1, ] * 10^5
        
        message('1.1')
        #marginErr <- sbfsnnls(X, X)
        marginErrB <- sbfsqpSum(X, X)
        message('1.2')
        marginErrF <- sffsqpSum(X, X, iInd = bInd)
        message('2')
        
        
        # if (scaleRecover) {
        #     W <- PrepResult@W%*%PrepResult@SW
        #     scale <- c(NMF::.fcnnls(X, matrix(W,ncol=1))$coef)
        #     scale[scale<1e-10] <- 0.01/(sqrt(colSums(X^2)))[scale<1e-10]
        #     X<-X%*%diag(scale)
        # }
        
        
        reconsErrB <- sbfsnnls(Xall, X)
        reconsErrF <- sffsnnls(Xall, X, iInd = cInd)
        #reconsErrB <- marginErrB
        #reconsErrF <- marginErrF
        
        message('3')
        # idx1 <- which.min(error1)
        # idx2 <- which.min(error2)
        
        ind1 <- as.integer(colnames(X[, marginErrB$sbfsset[[K]]]))
        ind2 <- as.integer(colnames(X[, marginErrF$sffsset[[K]]]))
        ind3 <- as.integer(colnames(X[, reconsErrB$sbfsset[[K]]]))
        ind4 <- as.integer(colnames(X[, reconsErrF$sffsset[[K]]]))
        
        error1 <- marginErrB$sbfscost
        error2 <- marginErrF$sffscost
        error3 <- reconsErrB$sbfscost
        error4 <- reconsErrF$sffscost
        
        
        return(new(
            "CAMMGObj",
            idx = list(
                marginErrB$sbfsset,
                marginErrF$sffsset,
                reconsErrB$sbfsset,
                reconsErrF$sffsset
            ),
            corner = rbind(ind1, ind2, ind3, ind4),
            error = cbind(error1, error2, error3, error4)
        ))
    }
