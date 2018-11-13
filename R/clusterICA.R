# ICA function
#' Independent Component Analysis (ICA)
#'
#' Uses random directions, clustering and optimisation to obtain ICA loadings,
#'  using the m-spacing entropy estimation as the objective function to minimise.
#'
#' @param x the data to perform ICA on, ncol(x) = n, nrow(x) = p
#' @param xw (optional) the whitened version of x
#' @param m (optional) the value of m-spacing for calculating approximate entropy, if missing(m), m <- sqrt(n)
#' @param nComp the number of ICA loadings outputted
#' @param p (optional) the size of the whitened matrix, i.e. how many PCA loadings to keep in the whitening step
#' @param randIter the number of random directions to initialise
#' @param randOut the number of the best random directions to keep
#' @param seed (optional) the set.seed number used for initialising the random directions
#' @param kmeanTol the tolerance used in divisive clustering, see clusterProjDivisive
#' @param kmeanIter the maximum number of iterations used in divisive clustering, see clusterProjDivisive
#' @param optimMaxit the maximum number of iterations used in the optimisation step, see optim
#' @param optimMethod the method used in the optimisation step, see optim
#' @param sizeClust (optional) if sizeClust = k > 1, then optimisation is performed on k random directions in each cluster.
#'                      If missing, then optimisation is performed on the best direction in each cluster.
#' @param computeScores if TRUE then scores of the whitened data are outputted
#' @param verbose if TRUE then information is given on the status of the function
#' @return A list with the following components:
#'         \itemize{
#'              \item{xw} {The output from jvcoords::whiten(x)}
#'              \item{IC} {The matrix of the loading vectors for the whitened data}
#'              \item{y} {The matrix of the projections of the whitened data along the loading vectors}
#'              \item{entr} {The m-spacing entropy of each row of y}
#'          }
#'
#' @author Paul Smith, \email{mmpws@@leeds.ac.uk}
# #' @seealso clusterProjKmeans()
#' @keywords independent component analysis, entropy, clustering
#'
#'
#' @examples
#' #---------------------------------------------------
#' #Example 1: un-mixing two stratified independent normals
#' #---------------------------------------------------
#' p <- 2
#' n <- 10000
#' set.seed(1243)
#' x1 <- matrix(rnorm(n*p, mean = 2.5), n, p)
#' x1[,2] <- scale(x1[,2])
#' a1 <- c(0.2,1)
#' a1 <- a1 / sqrt(sum(a1^2))
#' good1 <- cos(20 * x1 %*% a1) >= 0
#' x1Good <- x1[which(good1),]
#' x2 <- matrix(rnorm(n*p, mean = -2.5), n , p)
#' x2[,2] <- scale(x2[,2])
#' good2 <- cos(20 * x2 %*% a1) >= 0
#' x2Good <- x2[which(good2),]
#' xGood <- rbind(x1Good, x2Good)
#' a <- clusterICA(x=xGood, nComp=1, randIter=1000, seed=5)
#' par(mfrow = c(1,3))
#' plot(xGood, main = "Pre-processed data")
#' plot(a$xw$y, main = "Whitened data")
#' plot(density(a$y, bw="sj"), main = "ICA components")
#'
#' #---------------------------------------------------
#' #Example 2: un-mixing two mixed independent uniforms
#' #From fastICA man page
#' #---------------------------------------------------
#' S <- matrix(runif(10000), 5000, 2)
#' A <- matrix(c(1, 1, -1, 3), 2, 2, byrow = TRUE)
#' X <- S %*% A
#' a <- clusterICA(X, p=2, randIter=1000, randOut=50)
#' par(mfrow = c(1, 3))
#' plot(X, main = "Pre-processed data")
#' plot(a$xw$y, main = "Whitened data")
#' plot(a$y, main = "ICA components")
#'
#' #---------------------------------------------------
#' #Example 3: un-mixing iris data
#' #---------------------------------------------------
#' a <- clusterICA(iris[,1:4])
#' plot(a$y, main = "ICA components")
#' pairs(a$y, col=iris$Species)
#' @export
clusterICA <- function(x, xw, m, nComp, p, randIter=5000, randOut=100, seed,
                       kmeanTol=0.1, kmeanIter=100,
                       optimMaxit=1000, optimMethod="Nelder-Mead",
                       sizeClust, computeScores = TRUE, verbose=FALSE) {

    # check if we have whitened data
    # here p is how many PCA loadings we use to do ICA on
    # nComp is how many ICA loadings we want outputted
    if (missing(xw)) {
        xw <- jvcoords::whiten(x, compute.scores=TRUE)
    } else {
        # rescale xw so works with p
        if (!missing(p)) {
            if (!(floor(p) == p)) stop("p must be integer valued")
            if (p > xw$q) stop("p too large for xw")
            xw$loadings <- xw$loadings[,1:p]
            xw$q <- p
            xw$post.mul <- xw$post.mul[1:p]
            xw$y <- xw$y[,1:p]
        }
    }
    z <- xw$y

    n <- nrow(z)
    if(missing(p)) p <- ncol(z)
    if(missing(nComp)) nComp <- p
    if(missing(m)) m <- floor(sqrt(n))

    # some error checking
    if(nComp > (p)) {
        warning("nComp = ", nComp, " must be less than p + 1 = ",
                (p+1), ". Set nComp = p.")
        nComp <- p
    }

    # initiate loadings list
    loadings <- vector(mode = "list", length = nComp)
    entr <- numeric()
    IC <- diag(p)
    loopNum <- nComp
    if (loopNum == p) loopNum <- p - 1
    for(k in 1:loopNum) {
        if (verbose == TRUE) {
            cat("optimising direction", k, "out of", nComp, "\n")
        }
        r <- p - k + 1 # the dimension of the search space

        if (verbose == TRUE) {
            cat("// Finding random starting points", "\n")
        }
        randDir <- randDirs(z=z, IC=IC, k=k, m=m, iter=randIter, out=randOut,
                            seed=seed)
        if (verbose == TRUE) {
            cat("/// Found ", length(randDir$entr), " starting directions", "\n",
                sep="")
        }

        if (verbose == TRUE) {
            cat("/// Sorting these into clusters \n")
        }
        # do we want to save all directions in each cluster,
        # or just the best (pre-optim)
        if(!missing(sizeClust) && sizeClust < 1 && sizeClust > -1) {
            warning("sizeClust must be >= 1. Set sizeClust = 1")
            sizeClust <- 1L
        }
        if(!missing(sizeClust) && sizeClust < -1) {
            warning("sizeClust must be >= 1. Set sizeClust = ",
                    as.integer(-sizeClust))
            sizeClust <- as.integer(-sizeClust)
        }
        if(!missing(sizeClust) && (sizeClust > 1)) {
            bestDirs <- clusterNorm(z = z, IC=IC, k=k, m=m,
                                    dirs=randDir, kmeanTol=kmeanTol,
                                    kmeanIter=kmeanIter, saveAll=TRUE)
        } else {
            bestDirs <- clusterNorm(z = z, IC=IC, k=k, m=m,
                                    dirs=randDir, kmeanTol=kmeanTol,
                                    kmeanIter=kmeanIter)
        }

        if (verbose == TRUE) {
            cat("//// Sorted into ", length(bestDirs), " clusters", "\n", sep="")
        }
        entrPreOptim <- unlist(sapply(bestDirs, function(x) x$entr))
        if (verbose == TRUE) {
            cat("//// Best pre-optim entropy = ", min(entrPreOptim), "\n", sep="")
        }

        # step 2: use local optimisation to find the best solution in the
        # each cluster
        if (verbose == TRUE) {
            cat("//// Optimising ", length(bestDirs), " clusters", "\n", sep="")
        }
        icaLoading <- icaClusters(z=z, IC=IC, k=k, m=m,
                                  bestDirs=bestDirs, maxit = optimMaxit,
                                  optimMethod=optimMethod, sizeClust=sizeClust,
                                  verbose=verbose)
        if(!missing(sizeClust) && (sizeClust > 1)) {
            icaLoading <- icaLoading$best
        }
        if (verbose == TRUE) {
            cat("//// Optimised direction has entropy ",
                icaLoading$dir_entr, "\n", sep="")
        }
        bestEntr <- icaLoading$dir_entr
        bestDir <- icaLoading$dirOptim

        if (verbose == TRUE) {
            cat("///// Householder reflection", "\n")
        }
        # Use a Householder reflection which maps e1 to best.dir to update IC.
        e1 <- c(1, rep(0, r-1))
        # from wiki: take sign of x_k s.t.
        # k is the last col entry of non-zero in UT form A = QR
        signTmp <- sign(bestDir[1])
        v <- bestDir - signTmp * e1
        v <- v / sqrt(sum(v^2))
        P <- diag(r) - 2 * tcrossprod(v)
        IC[,k:p] <- IC[,k:p,drop=FALSE] %*% P
        entr[k] <- bestEntr
    }

    if (nComp == p) {
        bestDir <- 1
        wOrigSpace <- IC %*% c(rep(0, p-1), bestDir)
        zProj <- t(z %*% wOrigSpace)
        bestEntr <- entropy(zProj, m = m)
        entr[p] <- bestEntr
    }

    IC <- IC[, seq_len(nComp), drop=FALSE]

    colnames(IC) <- paste0('IC', seq_len(nComp))
    rownames(IC) <- rownames(xw$loadings)

    res <- list()
    res$xw <- xw
    res$IC <- IC
    if(computeScores == TRUE) {
        y <- z %*% IC
        res$y <- y
    }
    res$entropy <- entr
    class(res) = "clusterICA"
    res
}
