qpSum <- function(A, X) {
    #A <- A[,bInd[[2]][1,]]
    
    errorSum <- 0
    for (c in 1:dim(X)[2]) {
        qp <- limSolve::lsei(A, X[, c], matrix(1, 1, dim(A)[2]), 1,
                             diag(dim(A)[2]), rep(0, dim(A)[2]))
        errorSum <- errorSum + qp$solutionNorm
        # errorSum <- errorSum + sqrt(qp$solutionNorm)
    }
    return(errorSum)
}
