# k-means++ algorithm
# Choose one center uniformly at random from among the data points.
# For each data point x, compute D(x), the distance between x and the nearest center that has already been chosen.
# Choose one new data point at random as a new center, using a weighted probability distribution where a point x is chosen with probability proportional to D(x)2.
# Repeat Steps 2 and 3 until k centers have been chosen.
# Now that the initial centers have been chosen, proceed using standard k-means clustering.
clusterProjPlusPlus <- function(X, K) {
    n <- nrow(X)
    p <- ncol(X)

    DX <- rep(1/n, n)
    dist <- matrix(0, nrow=n, ncol=K)
    for (k in 1:K) {
        sampleTmp <- sample(n, size=1, prob=DX)
        pointTmp <- X[sampleTmp, ]
        dist[, k] <- 1 - (X %*% pointTmp)^2
        # overwrite to stop numerical errors so that DX is always positive
        dist[sampleTmp, k] <- 0
        DX <- apply(dist[,1:k,drop=FALSE], 1, min)
    }
    apply(dist, 1, which.min)
}

.clusterProjKmeans <- function(X, K, iter.max=100, initial, verbose=TRUE) {
    n <- nrow(X)
    if(missing(initial)) {
        c <- clusterProjPlusPlus(X, K)
        } else {
            # set cluster s.t. they start at 1
            if(!min(initial) == 1) initial <- initial - min(initial) + 1
            if (!(length(initial) == n)) {
                warning("'initial' must be of length n = ", 
                            n, ". Initialising using kmeans++ instead.")
                if(missing(K)) K <- max(initial)
                initial <- clusterProjPlusPlus(X, K)    
            } # end if
            c <- initial
            K <- max(c)
        } # end else
    
    j <- 0
    for (i in 1:iter.max) {
        j <- j + 1
        dist <- matrix(0, nrow=n, ncol=K)
        for (k in 1:K) {
            Y <- X[c == k,, drop=FALSE]
            # need to change this so that the clust mean changes...
            # don't want empty clustering
            if(!(nrow(Y) == 0)) {
                s <- La.svd(Y, nu=0, nv=1)
                centre <- s$vt[1,]
                dist[, k] <- 1 - (X %*% centre)^2
            }

        }
        cOld <- c
        c <- apply(dist, 1, which.min)

        if (all(cOld == c)) {
            break
        }
    }
    if(verbose == TRUE) {
        if(j < iter.max) cat("Converged to ", K, " clusters in ", j, " iterations \n")
        if(j == iter.max) cat("Maximum iterations of ", iter.max, " used, consider increasing iter.max \n")
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
# within cluster sum-of-squares
clusterProjWss <- function(X, c) {
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
    if (K == 1) {

    }
    wssClust <- unlist(sapply(clust, function(x) {
        nTmp <- nrow(x)
        if(nTmp == 0) {
            SEE <- 0
        } else {
            s <- La.svd(x, nu=0, nv=1)
            SSE <- nTmp - s$d[1]^2
        }
        }))
    if(any(wssClust < 0)) {
        whichWss <- which(wssClust < 0)
        wssClust[whichWss] <- 0
    }
    wss <- sum(wssClust)
    return(list(wss=wss, rss=wssClust))
}

# clustering method using hierarchical divisive clustering
# tolerance is a related to change in within sum-of-squares (wss) with each iteration
# when this change is less than the tolerance then break
#' Divisive (hierarchical) clustering on the projective space
#'
#' Creates clusters of points on the projective space using divisive k-means clustering
#'
#' @param X the data belonging to the projective space
#' @param tol the tolerance that when reached, stops increasing the number of clusters. 
#'                  At each step, the (change in wss) / (original wss) must be above this tolerance.
#'                  In general, as the tolerance decreases, the number of clusters in the output increases.
#' @param iter.max the maximum number of iterations
#'
#' @return A list with the following components:
#' \itemize{
#'          \item{c} {Vector of real numbers from 1 to K representing the cluster
#'                      that the corresponding X value belongs to.}
#'          \item{rss} {Vector of the within sum-of-squares for each cluster}
#'          \item{wss} {The total within sum-of-squares for the outputted cluster}
#'          \item{wssAll} {The change in total within sum-of-squares for each 
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

    p <- ncol(X)
    n <- nrow(X)

	i <- 1
	cCurr <- rep(1, n)
	rssAll <- clusterProjWss(X=X, c=cCurr)	
	wss <- rssAll$wss
	if(wss < 10e-16) {
        return(list(c=cCurr, rss=rssAll$rss, wss=wss, wssAll=NA))
    }
    if (missing(tol)) tol <- 0.1

	rssMax <- which.max(rssAll$rss)
	wssAll <- wss
    diffc <- 1
    diffct <- 0
	while(diffc > tol & i < iter.max) {
		i <- i + 1
		# select cluster with largest rss 
		# nb drop=F not needed as clust with only 1 element will not be picked
		clustMax <- X[cCurr == rssMax, , drop=FALSE]

		# split this cluster into two using kmeans
		cTmp <- clusterProjKmeans(X=clustMax, K=2)

		# rearrange cCurr so that we can add c=1 and c=2 for new clust
		whichCurr <- which(cCurr==rssMax)
		if(rssMax > 1) cCurr[cCurr < rssMax] <- cCurr[cCurr < rssMax] + 1
		cCurr[-whichCurr] <- cCurr[-whichCurr] + 1
		# add new clust to beginning
		cCurr[whichCurr] <- cTmp
		
		# find which cluster has largest RSS
		# this is at end of loop so that wss > tol_wss is correct
		rssAll <- clusterProjWss(X=X, c=cCurr)
		rssMax <- which.max(rssAll$rss)
		wss <- rssAll$wss
		wssAll[i] <- wss
        diffc <- (wssAll[i-1] - wssAll[i]) / wssAll[1]
	}
	if (i == iter.max) warning("Max iterations reached: increase iter.max or tol")

	list(c=cCurr, rss=rssAll$rss, wss=wss, wssAll=wssAll)
}