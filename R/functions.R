# The following function is an implementation of the k-Means algorithm
# on projective spaces.  Kmeans++ is used for the initial cluster assignments.
#' Kmeans clustering on the projective space
#'
#' Creates K clusters of points on the projective space using the k-means method
#'
#' @param X the data belonging to the projective space
#' @param K the number of clusters required in output
#' @param maxiter the maximum number of iterations
#' @param initial (optional) the initial clustering. 
#'                  If missing(K), then initial clusters found using the k-means++ method
#' @param verbose gives information of completed clusters
#'
#' @return Vector of real numbers from 1 to K representing the cluster
#'     that the corresponding X value belongs to.
#'
#'
#' @author Jochen Voss, \email{Jochen.Voss@@leeds.ac.uk}
#' @seealso \code{\link{projective.divisive_clust}}
#' @keywords clustering, kmeans
#'
#' @examples
#' n1 <- 37; n2 <- 19
#' x1 <- rnorm(n1, 6); y1 <- rnorm(n1, 0); z1 <- rnorm(n1, 0, 0.1)
#' x2 <- rnorm(n2, 8); y2 <- rnorm(n2, 8); z2 <- rnorm(n2, 0, 0.1)
#' X <- rbind(cbind(x1, y1, z1), cbind(x2, y2, z2)) * sample(c(-1, 1), size=n1+n2, replace=TRUE)
#' X <- X / sqrt(rowSums(X^2))
#' (c <- projective.cluster(X, 2))
#'
#' @export
projective.cluster <- function(X, K, maxiter=100, initial, verbose=TRUE) {
    n <- nrow(X)
    if(missing(initial)) {
        c <- .projective.plusplus(X, K)
        } else {
            if (!is.numeric(initial)) stop("initial must be a vector of clusters")
            c <- initial
            if(missing(K)) K <- max(c)
        }
    
    j <- 0
    for (i in 1:maxiter) {
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
        c.old <- c
        c <- apply(dist, 1, which.min)

        if (all(c.old == c)) {
            break
        }
    }
    if(verbose == TRUE) {
        if(j < maxiter) cat("Converged to ", K, " clusters in ", j, " iterations \n")
        if(j == maxiter) cat("Maximum iterations of ", maxiter, " used, consider increasing maxiter \n")
    }
    c
}


# clustering method using hierarchical divisive clustering
# tolerance is a related to change in wihtin sum-of-squares (wss) with each iteration
# when this change is less than the tolerance then break
#' Divisive (hierarchical) clustering on the projective space
#'
#' Creates clusters of points on the projective space using divisive kmeans clustering
#'
#' @param X the data belonging to the projective space
#' @param tol the tolerance that when reached, stops increasing the number of clusters. 
#'                  At each step, the (change in wss) / (original wss) must be above this tolerance
#' @param maxiter the maximum number of iterations
#'
#' @return A list with the following components:
#'          \list{c} {Vector of real numbers from 1 to K representing the cluster
#'                      that the corresponding X value belongs to.}
#'          \list{rss}{Vector of the within sum-of-squares for each cluster}
#'          \list{wss}{The total within sum-of-squares for the outputted cluster}
#'          \list{wss_all}{The change in total within sum-of-squares for each 
#'                              iteration of the function}
#'
#' @author Paul Smith, \email{mmpws@@leeds.ac.uk}
#' @seealso \code{\link{projective.cluster}}
#' @keywords clustering, kmeans
#'
#' 
#' @examples
#' n1 <- 37; n2 <- 19
#' x1 <- rnorm(n1, 6); y1 <- rnorm(n1, 0); z1 <- rnorm(n1, 0, 0.1)
#' x2 <- rnorm(n2, 8); y2 <- rnorm(n2, 8); z2 <- rnorm(n2, 0, 0.1)
#' X <- rbind(cbind(x1, y1, z1), cbind(x2, y2, z2)) * sample(c(-1, 1), size=n1+n2, replace=TRUE)
#' X <- X / sqrt(rowSums(X^2))
#' (c <- projective.divisive_clust(X=X, tol=0.1))
#'
#' @export
projective.divisive_clust <- function(X, tol, maxiter=100) {
    stopifnot(tol > 0 && tol <= 1)

    p <- ncol(X)
    n <- nrow(X)

	i <- 1
	c_curr <- rep(1, n)
	rss_all <- .projective.wss(X=X, c=c_curr)	
	wss <- rss_all$wss
	if(wss == 0) {
        return(list(c=c_curr, rss=rss_all$rss, wss=wss, wss_all=NA))
    }
    if (missing(tol)) tol <- 0.1

	rss_max <- which.max(rss_all$rss)
	wss_all <- wss
    diffc <- 1
	while(diffc > tol & i < maxiter) {
		i <- i + 1
		# select cluster with largest rss 
		# nb drop=F not needed as clust with only 1 element will not be picked
		clust_max <- X[c_curr == rss_max, , drop=FALSE]

		# split this cluster into two using kmeans
		c_tmp <- projective.cluster(X=clust_max, K=2, verbose=FALSE)

		# rearrange c_curr so that we can add c=1 and c=2 for new clust
		which_curr <- which(c_curr==rss_max)
		if(rss_max > 1) c_curr[c_curr < rss_max] <- c_curr[c_curr < rss_max] + 1
		c_curr[-which_curr] <- c_curr[-which_curr] + 1
		# add new clust to beginning
		c_curr[which_curr] <- c_tmp
		
		# find which cluster has largest RSS
		# this is at end of loop so that wss > tol_wss is correct
		rss_all <- .projective.wss(X=X, c=c_curr)
		rss_max <- which.max(rss_all$rss)
		wss <- rss_all$wss
		wss_all[i] <- wss
        diffc <- (wss_all[i-1] - wss_all[i]) / wss_all[1]
	}
	if (i == maxiter) warning("Max iterations reached: increase maxiter or tol")

    # account for numerical error in if statement
    # TODO: don't bother with this?
	if (any(rss_all$rss < 10e-16)) {
		cat("Sparse or missing clusters. Tolerance increased to tol = ", 
				tol + 0.05, "\n")
		rss_tmp <- rss_all$rss
		while(any(rss_tmp < 10e-16)) {
			tol <- tol + 0.05
			out <- projective.divisive_clust(X=X, tol=tol, maxiter=maxiter)
			rss_tmp <- out$rss_all$rss
			if(!any(rss_tmp < 10e-16)) return(out)
		}
	}
	list(c=c_curr, rss=rss_all$rss, wss=wss, wss_all=wss_all)
}

# ICA function
#' Approximate Independent Component Analysis (ICA) method
#'
#' Uses random directions, clustering and optimisation to obtain approximate ICA loadings,
#'  using the m-spacing entropy approximation as the objective function to minimise
#'
#' @param x the data to perform ICA on, ncol(x) = n, nrow(x) = p
#' @param xw (optional) the whitened version of x
#' @param m (optional) the value of m-spacing for calculating approximate entropy, if missing(m), m <- sqrt(n)
#' @param num_loadings the number of ICA loadings outputted
#' @param p (optional) the size of the whitened matrix, i.e. how many PCA loadings to keep in the whitening step
#' @param rand_iter the number of random directions to initialise
#' @param rand_out the number of the best random directions to keep
#' @param seed (optional) the set.seed number used for initialising the random directions
#' @param kmeans_tol the tolerance used in divisive clustering, see \code{\link[projective.divisive_clust]{tol}}
#' @param kmeans_iter the maximum number of iterations used in divisive clustering, see \code{\link[projective.divisive_clust]{maxiter}}
#' @param optim_maxit the maximum number of iterations used in the optimisation step, see \code{\link[optim]}
#' @param opt_method the method used in the optimisation step, see \code{\link[optim]}
#' @param size_clust (optional) if size_clust = k > 1, then optimisation is performed on k random directions in each cluster. 
#'                      If missing, then optimisation is performed on the best direction in each cluster.
#'
#' @return A list with the following components:
#'          \list{xw} {The output from jvmulti::whiten(x)}
#'          \list{IC}{The matrix of the loading vectors for the whitened data}
#'          \list{y}{The matrix of the projections of the whitened data along the loading vectors}
#'          \list{entr}{The m-spacing entropy of each row of y}
#'
#' @author Paul Smith, \email{mmpws@@leeds.ac.uk}
# #' @seealso \code{\link{projective.cluster}}
#' @keywords independent component analysis, entropy, clustering
#'
#' 
#' @examples
#' #---------------------------------------------------
#' #Example 1: un-mixing two stratified independent normals
#' #---------------------------------------------------
#' p <- 2
#' n <- 10000
#' x1 <- matrix(rnorm(n*p, mean = 2.5), n, p)
#' x1[,2] <- scale(x1[,2])
#' a1 <- c(0,1)
#' a1 <- a1 / sqrt(sum(a1^2))
#' good1 <- cos(30 * x1 %*% a1) >= 0
#' x1_good <- x1[which(good1),]
#' x2 <- matrix(rnorm(n*p, mean = -2.5), n , p)
#' x2[,2] <- scale(x2[,2])
#' x2_good <- x2[which(good1),]
#' x_good <- rbind(x1_good, x2_good)
#' a <- goodICA(x=x_good, num_loadings=1)
#' par(mfrow = c(1,3))
#' plot(x_good, main = "Pre-processed data")
#' plot(a$xw$y, main = "PCA components")
#' plot(density(a$x, bw="sj"), main = "ICA components")
#'
#' #---------------------------------------------------
#' #Example 2: un-mixing two mixed independent uniforms
#' #From fastICA man page
#' #---------------------------------------------------
#' S <- matrix(runif(10000), 5000, 2)
#' A <- matrix(c(1, 1, -1, 3), 2, 2, byrow = TRUE)
#' X <- S %*% A
#' a <- goodICA(X, p=2, rand_iter=1000, rand_out=50)
#' par(mfrow = c(1, 3))
#' plot(X, main = "Pre-processed data")
#' plot(a$xw$y, main = "PCA components")
#' plot(a$x, main = "ICA components")
#' @export
goodICA <- function(x, xw, m, num_loadings, p, rand_iter=5000, rand_out=100, seed, 
                    kmeans_tol=0.1, kmeans_iter=100,
                    optim_maxit=1000, opt_method="Nelder-Mead",
                    size_clust) {
    
    # check if we have whitened data
    # here p is how many PCA loadings we use to do ICA on
    # num_loadings is how many ICA loadings we want outputted
    if (missing(xw)) {
        xw <- jvmulti::whiten(x, compute.scores=TRUE)
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
    if(missing(num_loadings)) num_loadings <- p - 1
    if(missing(m)) m <- floor(sqrt(n))

    # some error checking
    if(num_loadings > p) {
        warning("num_loadings = ", num_loadings, " must be less than p = ", p, 
                    ". Set num_loadings = p - 1.")
        num_loadings <- p - 1
    }

    # initiate loadings list
    loadings <- vector(mode = "list", length = num_loadings)
    entr <- numeric()
    IC <- diag(p)
    for(k in 1:num_loadings) {
        cat("optimising direction", k, "out of", num_loadings, "\n")
        r <- p - k + 1 # the dimension of the search space
        
        cat("// Finding random starting points", "\n")
        rand_dir <- .rand.dirs(z=z, IC=IC, k=k, m=m, iter=rand_iter, out=rand_out,
                              seed=seed)
        cat("/// Found ", length(rand_dir$entr), " starting directions", "\n", 
            sep="")
        
        cat("/// Sorting these into clusters \n")
        # do we want to save all directions in each cluster,
        # or just the best (pre-optim)
        if(!missing(size_clust) && size_clust < 1 && size_clust > -1) {
                warning("size_clust must be >= 1. Set size_clust = 1")
                size_clust <- 1L
            }
        if(!missing(size_clust) && size_clust < -1) {
                warning("size_clust must be >= 1. Set size_clust = ", 
                            as.integer(-size_clust))
                size_clust <- as.integer(-size_clust)
            }
        if(!missing(size_clust) && (size_clust > 1)) {
            best_dirs <- .cluster.norm(z = z, IC=IC, k=k, m=m,
                                      dirs=rand_dir, kmeans_tol=kmeans_tol,
                                      kmeans_iter=kmeans_iter, save.all=TRUE)
        } else {
            best_dirs <- .cluster.norm(z = z, IC=IC, k=k, m=m,
                                      dirs=rand_dir, kmeans_tol=kmeans_tol,
                                      kmeans_iter=kmeans_iter)
        }

        cat("//// Sorted into ", length(best_dirs), " clusters", "\n", sep="")
        entr.pre_optim <- unlist(sapply(best_dirs, function(x) x$entr))
        cat("//// Best pre-optim entropy = ", min(entr.pre_optim), "\n", sep="")
        
        # step 2: use local optimisation to find the best solution in the
        # each cluster
        cat("//// Optimising ", length(best_dirs), " clusters", "\n", sep="")
        ica_loading <- .ica.clusters(z=z, IC=IC, k=k, m=m,
                                    best_dirs=best_dirs, maxit = optim_maxit,
                                    opt_method=opt_method, size_clust=size_clust)
        if(!missing(size_clust) && (size_clust > 1)) {
            ica_loading <- ica_loading$best
        }
        cat("//// Optimised direction has entropy ", 
                ica_loading$dir_entr, "\n", sep="")
        best_entr <- ica_loading$dir_entr
        best_dir <- ica_loading$dir_optim

        cat("///// Householder reflection", "\n")
        # Use a Householder reflection which maps e1 to best.dir to update IC.
        e1 <- c(1, rep(0, r-1))
        # from wiki: take sign of x_k s.t.
        # k is the last col entry of non-zero in UT form A = QR
        sign_tmp <- sign(best_dir[1])
        v <- best_dir - sign_tmp * e1
        v <- v / sqrt(sum(v^2))
        P <- diag(r) - 2 * tcrossprod(v)
        IC[,k:p] <- IC[,k:p,drop=FALSE] %*% P
        entr[k] <- best_entr
    }
    IC <- IC[, seq_len(num_loadings), drop=FALSE]

    colnames(IC) <- paste0('IC', seq_len(num_loadings))

    res <- list(xw=xw, IC=IC, y=z %*% IC, entr=entr)
    class(res) = "goodICA"
    res
}