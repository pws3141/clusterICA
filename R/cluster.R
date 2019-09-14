# k-means++ algorithm
# Choose one center uniformly at random from among the data points.
# For each data point x, compute D(x), the distance between x and 
# the nearest center that has already been chosen.
# Choose one new data point at random as a new center, using a weighted 
# probability distribution where a point x is chosen with probability proportional to D(x)2.
# Repeat Steps 2 and 3 until k centers have been chosen.
# Now that the initial centers have been chosen, proceed using standard k-means clustering.
clusterProjPlusPlus <- function(X, K) {
    n <- nrow(X)
    DX <- rep(1, n)
    dist <- matrix(0, nrow=n, ncol=K)
    allSamples <- as.numeric() 
    for (k in 1:K) {
        sampleTmp <- sample(n, size=1, prob=DX)
        allSamples <- c(allSamples, sampleTmp)
        pointTmp <- X[sampleTmp, ]
        dist[, k] <- 1 - (X %*% pointTmp)^2
        DX <- apply(dist[,1:k,drop=FALSE], 1, min)
        # overwrite to stop numerical errors so that DX is always positive
        DX[allSamples] <- 0
        # use following to stop "too few positive probabilities" error
        if (all(DX == 0)) DX[-allSamples] <- rep(1, n - k)
        # as sample() only requires relative probabilites
        # transform s.t. max(DX) = 1
        DXlog <- log(DX)
        DXlog <- DXlog - max(DXlog) # make max DXlog equal 0
        DX <- exp(DXlog) # exp(log) to set max to 1
        if(any(is.na(DX))) DX[is.na(DX)] <- 0
        # ensure  no error in sample()
        if (all(DX == 0)) DX[-allSamples] <- rep(1, n - k)
    }
    apply(dist, 1, which.min)
}

.clusterProjKmeans <- function(X, K, iter.max=100, initial, verbose=TRUE) {
    n <- nrow(X)
    if(missing(initial)) {
            if (missing(K)) {
                    emessage <- paste0("missing both 'K' and 'initial'")
                    stop(emessage)
            }
        c <- clusterProjPlusPlus(X, K)
    } else {
            if (missing(K)) {
                    K <- max(initial)
                    if(!(K == as.integer(K))) {
                            emessage <- paste0("max(initial) = ", K, 
                                               " is not integer-valued")
                            stop(emessage)
                    }
            }
            if(!all(sort(unique(initial)) == 1:K)) {
                    emessage <- paste0("'initial' must contain integers 1 to K = ", K)
                    stop(emessage)
            }
            if (!(length(initial) == n)) {
                    emessage <- paste0("'initial' must be of length n = ", n)
                    stop(emessage)
            } # end if
        c <- initial
    } # end else
    j <- 0
    for (i in 1:iter.max) {
        j <- j + 1
        dist <- matrix(0, nrow=n, ncol=K)
        for (k in 1:K) {
            Y <- X[c == k,, drop=FALSE]
            # don't want empty clustering
            if(!(nrow(Y) == 0)) {
                s <- La.svd(Y, nu=0, nv=1)
                centre <- s$vt[1,]
                dist[, k] <- 1 - (X %*% centre)^2
            } # end if
        } # end for
        cOld <- c
        c <- apply(dist, 1, which.min)
        if (all(cOld == c)) {
            break
        } # end if
    } # end for
    if(verbose == TRUE) {
        if(j < iter.max) {
                cat("Converged to ", K, " clusters in ", j, " iterations \n")
        }
        if(j == iter.max) {
                cat("Maximum iterations of ", iter.max, 
                        " used, consider increasing iter.max \n")
        }
    }
    c
}

# The following function is an implementation of the k-Means algorithm
# on projective spaces.  K-means++ is used for the initial cluster assignments.
#' K-means clustering on the projective space
#'
#' Creates K clusters of points on the projective space using the k-means method
#'
#' @param X the data belonging to the projective space
#' @param K the number of clusters required in output
#' @param iter.max the maximum number of iterations
#' @param initial (optional) the initial clustering.
#'                  The argument 'initial' is required to be a vector of
#'                      the same length as number of rows in X.
#'                  Each element of 'initial' is the cluster number that
#'                      the corresponding row of X belongs to.
#'                  If 'initial' is specified, than 'K' is set to be the
#'                      number of clusters in 'initial'
#'                  If no initial is given, then the initial clusters are
#'                      found using the k-means++ method.
#'
#' @return Vector of real numbers from 1 to K representing the cluster
#'     that the corresponding X value belongs to.
#'
#'
#' @author Jochen Voss, \email{Jochen.Voss@@leeds.ac.uk}
#' @seealso clusterProjDivisive
#' @keywords clustering kmeans
#'
#' @examples
#' n1 <- 37; n2 <- 19
#' x1 <- rnorm(n1, 6); y1 <- rnorm(n1, 0); z1 <- rnorm(n1, 0, 0.1)
#' x2 <- rnorm(n2, 8); y2 <- rnorm(n2, 8); z2 <- rnorm(n2, 0, 0.1)
#' X <- rbind(cbind(x1, y1, z1), cbind(x2, y2, z2)) * sample(c(-1, 1), size=n1+n2, replace=TRUE)
#' X <- X / sqrt(rowSums(X^2))
#' (c <- clusterProjKmeans(X, 2))
#'
#' @export
clusterProjKmeans <- function(X, K, iter.max=100, initial) {
    c <- .clusterProjKmeans(X=X, K=K, iter.max=iter.max,
                            initial=initial, verbose=FALSE)
    c
}

# find total sum of squares and
# within cluster sum-of-squares and
# root-mean-squared-error
clusterRMSE <- function(X, c) {
    K <- max(c)
    p <- ncol(X)
    n <- nrow(X)
    # emply list of matrix for splitting into clusters
    clust <- vector(mode = "list", length = K)
    for(i in 1:K) {
        clust[[i]] <- matrix(nrow=length(which(c == i)), ncol = p)
    }
    # split X into clusters
    ticker <- rep(0, K)
    for(i in 1:n) {
        cTmp <- c[i]
        ticker[cTmp] <- ticker[cTmp] + 1
        clust[[cTmp]][ticker[cTmp],] <- X[i,]
    }
    squaredErrorsIndividual <- sapply(clust, function(x) {
                        nTmp <- nrow(x)
                        if(nTmp == 0) {
                            sse <- 0
                        } else {
                            s <- La.svd(x, nu=0, nv=1)
                            sse <- nTmp - s$d[1]^2
                        }
                        mse <- sse / nTmp
                        resTmp <- c(mse, sse)
                    })
    rownames(squaredErrorsIndividual) <- c("mse", "sse") 
    mseIndividual <- squaredErrorsIndividual[1, ]
    sseIndividual <- squaredErrorsIndividual[2, ]
    if(any(sseIndividual < 0)) {
        whichWss <- which(sseIndividual < 0)
        sseIndividual[whichWss] <- 0
        mseIndividual[whichWss] <- 0
    }
    sse <- sum(sseIndividual)
    rmse <- sqrt(sse / n)
    return(list(rmse = rmse, 
                sse = sse, 
                sseIndividual = sseIndividual))
}

# clustering method using hierarchical divisive clustering
# tolerance is a related to change in within cluster sum-of-squares (wcss) with each iteration
# when this change is less than the tolerance then break
#' Divisive (hierarchical) clustering on the projective space
#'
#' Creates clusters of points on the projective space using divisive k-means clustering
#'
#' @param X the data belonging to the projective space
#' @param tol the tolerance that when reached, stops increasing the number of clusters.
#'                  At each step, the (change in wcss) / (original wcss) must be above this tolerance.
#'                  In general, as the tolerance decreases, the number of clusters in the output increases.
#' @param iter.max the maximum number of iterations
#'
#' @return A list with the following components:
#' \itemize{
#'          \item{c} {Vector of real numbers from 1 to K representing the cluster
#'                      that the corresponding X value belongs to.}
#'          \item{sseIndividual} {Vector of the within sum-of-squares for each cluster}
#'          \item{rmse} {The total within-cluster root mean-squared-error for the
#'          obtained cluster} 
#'          \item{rmseSequence} {The change in total within-cluster root mean-squared-error for each
#'                              iteration of the function}
#' }
#' @author Paul Smith, \email{mmpws@@leeds.ac.uk}
#' @seealso clusterProjKmeans
#' @keywords clustering, kmeans
#'
#'
#' @examples
#' n1 <- 37; n2 <- 19
#' x1 <- rnorm(n1, 6); y1 <- rnorm(n1, 0); z1 <- rnorm(n1, 0, 0.1)
#' x2 <- rnorm(n2, 8); y2 <- rnorm(n2, 8); z2 <- rnorm(n2, 0, 0.1)
#' X <- rbind(cbind(x1, y1, z1), cbind(x2, y2, z2)) * sample(c(-1, 1), size=n1+n2, replace=TRUE)
#' X <- X / sqrt(rowSums(X^2))
#' (c <- clusterProjDivisive(X=X, tol=0.1))
#'
#' @export
clusterProjDivisive <- function(X, tol, iter.max=100) {
    stopifnot(tol > 0 && tol <= 1)
    n <- nrow(X)
    i <- 1
    cCurr <- rep(1, n)
    errorRes <- clusterRMSE(X=X, c=cCurr)
    rmse <- errorRes$rmse
    if(rmse < 1e-10) {
        return(list(c=cCurr, 
                    sseIndividual = errorRes$sseIndividual, 
                    rmse = rmse, rmseSequence = NA))
    }
    if (missing(tol)) tol <- 0.1
    sseIndividualMax <- which.max(errorRes$sseIndividual)
    rmseSequence <- rmse
    diffc <- 1
    while(diffc > tol & i < iter.max) {
        i <- i + 1
        # select cluster with largest squared error
        clustMax <- X[cCurr == sseIndividualMax, , drop=FALSE]
        # split this cluster into two using kmeans
        cTmp <- clusterProjKmeans(X = clustMax, K = 2)
        # rearrange cCurr so that we can add c=1 and c=2 for new clust
        whichCurr <- which(cCurr == sseIndividualMax)
        if(sseIndividualMax > 1) { # make sure c = 1 is free
                cCurr[cCurr < sseIndividualMax] <- cCurr[cCurr < sseIndividualMax] + 1
        }
        cCurr[-whichCurr] <- cCurr[-whichCurr] + 1 # make sure c = 2 is free
        # add new clust to beginning, i.e. c = 1 and c = 2
        cCurr[whichCurr] <- cTmp
        # find which cluster has largest squared error
        errorRes <- clusterRMSE(X=X, c=cCurr)
        sseInTmp <- errorRes$sseIndividual
        if (isTRUE(all.equal(sseInTmp, rep(0, times=length(sseInTmp))))) break
        sseIndividualMax <- which.max(errorRes$sseIndividual)
        rmse <- errorRes$rmse
        rmseSequence[i] <- rmse
        diffc <- (rmseSequence[i-1] - rmseSequence[i]) / rmseSequence[1]
    }
    if (i == iter.max) warning("Max iterations reached: consider 
                               increasing kmeans.tol")
    # remove empty clusters
    j <- 1
    while (j <= max(cCurr)) {
        if (length(cCurr[cCurr == j]) == 0) {
            cCurr[cCurr > j] <- cCurr[cCurr > j] - 1
            j <- j - 1
        }
        j <- j + 1
    }
    list(c = cCurr, sseIndividual = errorRes$sseIndividual, rmse = rmse, 
         rmseSequence = rmseSequence)
}
