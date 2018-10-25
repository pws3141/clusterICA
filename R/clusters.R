# k-means++ algorithm
# Choose one center uniformly at random from among the data points.
# For each data point x, compute D(x), the distance between x and the nearest center that has already been chosen.
# Choose one new data point at random as a new center, using a weighted probability distribution where a point x is chosen with probability proportional to D(x)2.
# Repeat Steps 2 and 3 until k centers have been chosen.
# Now that the initial centers have been chosen, proceed using standard k-means clustering.
proj.plusplus.part <- function(X, K) {
    n <- nrow(X)
    p <- ncol(X)

    DX <- rep(1/n, n)
    dist <- matrix(0, nrow=n, ncol=K)
    for (k in 1:K) {
        sample_tmp <- sample(n, size=1, prob=DX)
        point_tmp <- X[sample_tmp, ]
        dist[, k] <- 1 - (X %*% point_tmp)^2
        # overwrite to stop numerical errors so that DX is always positive
        dist[sample_tmp, k] <- 0
        DX <- apply(dist[,1:k,drop=FALSE], 1, min)
    }
    apply(dist, 1, which.min)
}


# The following function is an implementation of the k-Means algorithm
# on projective spaces.  PCA-part is used for the initial cluster assignments.
proj.cluster <- function(X, K, maxiter=100, initial, verbose=TRUE) {
    n <- nrow(X)
    if(missing(initial)) {
    	c <- proj.plusplus.part(X, K)
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

proj.divisive_clust <- function(X, tol, maxiter=100) {
    stopifnot(tol > 0 && tol <= 1)

    p <- ncol(X)
    n <- nrow(X)

	i <- 0
	c_curr <- rep(1, n)
	wss <- proj.wss(X=X, c=c_curr)$wss
	if (missing(tol)) {
		tol_wss <- 0.4 * wss
		} else {
			tol_wss <- tol *wss
		}
	rss_all <- proj.wss(X=X, c=c_curr)	
	rss_max <- which.max(rss_all$rss)
	wss_all <- wss
	while(wss > tol_wss & i < maxiter) {
		i <- i + 1
		# select cluster with largest rss 
		# nb drop=F not needed as clust with only 1 element will not be picked
		clust_max <- X[c_curr == rss_max, , drop=FALSE]

		# split this cluster into two using kmeans
		c_tmp <- proj.cluster(X=clust_max, K=2, verbose=FALSE)

		# rearrange c_curr so that we can add c=1 and c=2 for new clust
		which_curr <- which(c_curr==rss_max)
		if(rss_max > 1) c_curr[c_curr < rss_max] <- c_curr[c_curr < rss_max] + 1
		c_curr[-which_curr] <- c_curr[-which_curr] + 1
		# add new clust to beginning
		c_curr[which_curr] <- c_tmp
		
		# find which cluster has largest RSS
		# this is at end of loop so that wss > tol_wss is correct
		rss_all <- proj.wss(X=X, c=c_curr)
		rss_max <- which.max(rss_all$rss)
		wss <- rss_all$wss
		wss_all[i] <- wss
	}
	if (i == maxiter) warning("Max iterations reached: increase maxiter or tol")
	#cat(tol, "\n")
	if (any(rss_all$rss == 0)) {
		cat("Sparse or missing clusters. Tolerance increased to tol = ", 
				tol + 0.1, "\n")
		rss_tmp <- rss_all$rss
		while(any(rss_tmp == 0)) {
			tol <- tol + 0.1
			out <- proj.divisive_clust(X=X, tol=tol, maxiter=maxiter)
			rss_tmp <- out$rss_all$rss
			if(!any(rss_tmp == 0)) return(out)
		}
	}
	list(c=c_curr, rss=rss_all$rss, wss=wss, wss_all=wss_all, tol=tol)
}

proj.wss <- function(X, c) {
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
        c_tmp <- c[i]
        ticker[c_tmp] <- ticker[c_tmp] + 1
        clust[[c_tmp]][ticker[c_tmp],] <- X[i,]
    }
    wss_clust <- sapply(clust, function(x) {
        n_tmp <- nrow(x)
        # TODO: make this if statement redundant. Don't want empty clusters
        if(n_tmp == 0) {
            SEE <- 0
        } else {
            s <- La.svd(x, nu=0, nv=1)
            SSE <- n_tmp - s$d[1]^2
        }
    })
    wss <- sum(wss_clust)
    return(list(wss=wss, rss=wss_clust))
}



