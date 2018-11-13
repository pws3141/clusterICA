# ICA function
#' Independent Component Analysis (ICA)
#'
#' Uses random directions, clustering and optimisation to obtain ICA loadings,
#'  using the m-spacing entropy estimation as the objective function to minimise.
#'
#' @param x the data to perform ICA on, ncol(x) = n, nrow(x) = p
#' @param xw (optional) the whitened version of x
#' @param m (optional) the value of m-spacing for calculating approximate entropy, if missing(m), m <- sqrt(n)
#' @param n.comp the number of ICA loadings outputted
#' @param p (optional) the size of the whitened matrix, i.e. how many PCA loadings to keep in the whitening step
#' @param rand.iter the number of random directions to initialise
#' @param rand.out the number of the best random directions to keep
#' @param seed (optional) the set.seed number used for initialising the random directions
#' @param kmean.tol the tolerance used in divisive clustering, see clusterProjDivisive
#' @param kmean.iter the maximum number of iterations used in divisive clustering, see clusterProjDivisive
#' @param opt.maxit the maximum number of iterations used in the optimisation step, see optim
#' @param opt.method the method used in the optimisation step, see optim
#' @param size.clust (optional) if size.clust = k > 1, then optimisation is performed on k random directions in each cluster.
#'                      If missing, then optimisation is performed on the best direction in each cluster.
#' @param compute.scores if TRUE then scores of the whitened data are outputted
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
#' a <- clusterICA(x=xGood, n.comp=1, rand.iter=1000, seed=5)
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
#' a <- clusterICA(X, p=2, rand.iter=1000, rand.out=50)
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
clusterICA <- function(x, xw, m, n.comp, p, rand.iter=5000, rand.out=100, seed,
		    kmean.tol=0.1, kmean.iter=100,
		    opt.maxit=1000, opt.method="Nelder-Mead",
		    size.clust, compute.scores = TRUE, verbose=FALSE) {
    # check if we have whitened data
    # here p is how many PCA loadings we use to do ICA on
    # n.comp is how many ICA loadings we want outputted
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
    if(missing(n.comp)) n.comp <- p
    if(missing(m)) m <- floor(sqrt(n))

    # some error checking
    if(n.comp > (p)) {
	warning("n.comp = ", n.comp, " must be less than p + 1 = ",
		(p+1), ". Set n.comp = p.")
	n.comp <- p
    }

    # initiate loadings list
    loadings <- vector(mode = "list", length = n.comp)
    entr <- numeric()
    IC <- diag(p)
    loopNum <- n.comp
    if (loopNum == p) loopNum <- p - 1
    for(k in 1:loopNum) {
	if (verbose == TRUE) {
	    cat("optimising direction", k, "out of", n.comp, "\n")
	}
	r <- p - k + 1 # the dimension of the search space

	if (verbose == TRUE) {
	    cat("// Finding random starting points", "\n")
	}
	randDir <- randDirs(z=z, IC=IC, k=k, m=m, iter=rand.iter, out=rand.out,
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
	if(!missing(size.clust) && size.clust < 1 && size.clust > -1) {
		warning("size.clust must be >= 1. Set size.clust = 1")
		size.clust <- 1L
	    }
	if(!missing(size.clust) && size.clust < -1) {
		warning("size.clust must be >= 1. Set size.clust = ",
			    as.integer(-size.clust))
		size.clust <- as.integer(-size.clust)
	    }
	if(!missing(size.clust) && (size.clust > 1)) {
	    best.dirs <- clusterNorm(z = z, IC=IC, k=k, m=m,
				      dirs=randDir, kmean.tol=kmean.tol,
				      kmean.iter=kmean.iter, save.all=TRUE)
	} else {
	    best.dirs <- clusterNorm(z = z, IC=IC, k=k, m=m,
				      dirs=randDir, kmean.tol=kmean.tol,
				      kmean.iter=kmean.iter)
	}

	if (verbose == TRUE) {
	    cat("//// Sorted into ", length(best.dirs), " clusters", "\n", sep="")
	}
	entrPreOptim <- unlist(sapply(best.dirs, function(x) x$entr))
	if (verbose == TRUE) {
	    cat("//// Best pre-optim entropy = ", min(entrPreOptim), "\n", sep="")
	}

	# step 2: use local optimisation to find the best solution in the
	# each cluster
	if (verbose == TRUE) {
	    cat("//// Optimising ", length(best.dirs), " clusters", "\n", sep="")
	}
	icaLoading <- icaClusters(z=z, IC=IC, k=k, m=m,
				    best.dirs=best.dirs, maxit = opt.maxit,
				    opt.method=opt.method, size.clust=size.clust,
				    verbose=verbose)
	if(!missing(size.clust) && (size.clust > 1)) {
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

    if (n.comp == p) {
	bestDir <- 1
	wOrigSpace <- IC %*% c(rep(0, p-1), bestDir)
	zProj <- t(z %*% wOrigSpace)
	bestEntr <- entropy(zProj, m = m)
	entr[p] <- bestEntr
    }

    IC <- IC[, seq_len(n.comp), drop=FALSE]

    colnames(IC) <- paste0('IC', seq_len(n.comp))
    rownames(IC) <- rownames(xw$loadings)

    res <- list()
    res$xw <- xw
    res$IC <- IC
    if(compute.scores == TRUE) {
	y <- z %*% IC
	res$y <- y
    }
    res$entropy <- entr
    class(res) = "clusterICA"
    res
}
