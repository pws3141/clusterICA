cluster.norm <- function(z, IC, k, m, dirs, kmeans_tol=0.1,
                         kmeans_iter=100, save.all=FALSE) {
    # convert dirs to listrbose=
    p <- ncol(z)
    n <- nrow(z)
    if(missing(m)) m <- floor(sqrt(n))
    if(missing(IC)) IC <- diag(p)
    if(missing(k)) k <- 1

    #stopifnot(p == ncol(dirs))
    entr <- dirs$entr
    dirs <- dirs$dirs
    dirs.list <- lapply(seq_len(nrow(dirs)), function(i) dirs[i,])

    # list of clusters

    # K-Means Cluster Analysis: Divisive
    c <- projective.divisive_clust(X=dirs, tol=kmeans_tol, maxiter=kmeans_iter)
    clusters <- max(c$c)
    
    # append cluster assignment & put into list
    out_tmp <- vector(mode = "list", length = clusters)
    dirs.cluster_append <- cbind(c$c, entr, dirs)
    for(i in 1:clusters) {
        which.cluster <- which(dirs.cluster_append[,1] == i)
        if (save.all == FALSE) {
            out_tmp[[i]]$entr <- dirs.cluster_append[which.cluster, 2]
            entr_min <- which.min(out_tmp[[i]]$entr)
            out_tmp[[i]]$entr <- out_tmp[[i]]$entr[entr_min]
            out_tmp[[i]]$dirs <- dirs.cluster_append[which.cluster, c(-1, -2)]
            out_tmp[[i]]$dirs <- out_tmp[[i]]$dirs[entr_min,]
        } else {
            out_tmp[[i]]$entr <- dirs.cluster_append[which.cluster, 2]
            out_tmp[[i]]$dirs <- dirs.cluster_append[which.cluster, c(-1, -2)]
        }
    }
    out_tmp
    return(out_tmp)
}


goodICA <- function(x, xw, m, num_loadings, p, rand_iter=5000, rand_out=500,
                    kmeans_tol=0.1, kmeans_iter=100,
                    optim_maxit=1000, seed, opt_method="Nelder-Mead",
                    size_clust) {
    
    # check if we have whitened data
    # here p is how many PCA loadings we use to do ICA on
    # num_loadings is how many ICA loadings we want outputted
    if (missing(xw)) {
        if (missing(p)) {
            if(missing(num_loadings)) stop("Need to specify p or num_loadings")
            cat("Set p = num_loadings + 1. 
                Worth specifying p larger s.t. whitening is not so drastic.")
            p <- num_loadings + 1
        }
        xw <- whiten(x, compute.scores=TRUE)
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
        rand_dir <- rand.dirs(z=z, IC=IC, k=k, m=m, iter=rand_iter, out=rand_out,
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
            best_dirs <- cluster.norm(z = z, IC=IC, k=k, m=m,
                                      dirs=rand_dir, kmeans_tol=kmeans_tol,
                                      kmeans_iter=kmeans_iter, save.all=TRUE)
        } else {
            best_dirs <- cluster.norm(z = z, IC=IC, k=k, m=m,
                                      dirs=rand_dir, kmeans_tol=kmeans_tol,
                                      kmeans_iter=kmeans_iter)
        }

        cat("//// Sorted into ", length(best_dirs), " clusters", "\n", sep="")
        entr.pre_optim <- unlist(sapply(best_dirs, function(x) x$entr))
        cat("//// Best pre-optim entropy = ", min(entr.pre_optim), "\n", sep="")
        
        # step 2: use local optimisation to find the best solution in the
        # each cluster
        cat("//// Optimising ", length(best_dirs), " clusters", "\n", sep="")
        ica_loading <- ica.clusters(z=z, IC=IC, k=k, m=m,
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
