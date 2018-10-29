
# The following function is an implementation of the k-Means algorithm
# on projective spaces.  PCA-part is used for the initial cluster assignments.
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


# clustering method using heirarchical divisive clustering
# tolerance is a percentage of the total sum of squares
# when total within-sum-of-squares hits this tol then breaks
projective.divisive_clust <- function(X, tol, maxiter=100) {
    stopifnot(tol > 0 && tol <= 1)

    p <- ncol(X)
    n <- nrow(X)

	i <- 1
	c_curr <- rep(1, n)
	rss_all <- .projective.wss(X=X, c=c_curr)	
	wss <- rss_all$wss
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
	list(c=c_curr, rss=rss_all$rss, wss=wss, wss_all=wss_all, tol=tol)
}

# ICA function

goodICA <- function(x, xw, m, num_loadings, p, rand_iter=5000, rand_out=100,
                    kmeans_tol=0.1, kmeans_iter=100,
                    optim_maxit=1000, seed, opt_method="Nelder-Mead",
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

    res <- list(xw=xw, x=z %*% IC, IC=IC, entr=entr, m=m)
    class(res) = "goodICA"
    res
}