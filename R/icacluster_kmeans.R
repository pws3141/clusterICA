cluster.norm <- function(z, IC, k, m, dirs, clusters,
                         kmeans_iter=100, kmeans_initial,
                         save.all=FALSE) {
    # convert dirs to listrbose=
    p <- ncol(z)
    n <- nrow(z)
    if(missing(m)) m <- floor(sqrt(n))
    if(missing(IC)) IC <- diag(p)
    if(missing(k)) k <- 1
    if(!missing(kmeans_initial) && missing(clusters)) clusters <- max(kmeans_initial)

    #stopifnot(p == ncol(dirs))
    entr <- dirs$entr
    dirs <- dirs$dirs
    dirs.list <- lapply(seq_len(nrow(dirs)), function(i) dirs[i,])

    # draw graph to decide number to clusters
    if(missing(clusters)) {
        wss <- numeric()
        for(i in 1:15) {
            c <- projective.cluster(X=dirs, K=i, maxiter=100, verbose=FALSE)
            wss[i] <- projective.wss(X=dirs, c)
        }
        wss_max <- max(wss)
        plot(1:15, wss, type="b", ylim=c(0,wss_max), xlab="Number of Clusters",
             ylab="Within groups sum of squares")
        clusters <- readline("How many clusters? \n")
        if(!is.null(dev.list())) dev.off()
    }
    clusters <- as.numeric(clusters)
    out_tmp <- vector(mode = "list", length = clusters)
    # list of clusters
    #return(clusters)
    # K-Means Cluster Analysis
    c <- projective.cluster(X = dirs, K = clusters, maxiter=kmeans_iter,
                      initial=kmeans_initial)

    # append cluster assignment & put into list
    dirs.cluster_append <- cbind(c, entr, dirs)
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
                    clusters=NA, kmeans_iter=100,
                    optim_maxit=1000, seed, opt_method="Nelder-Mead",
                    size_clust) {
    # do we have whitened data?
    if (missing(xw)) {
        if (missing(p)) stop("Need to specify p")
        trans <- whiten(x, n.comp=p)
    } else {
        # rescale xw so works with p
        if (!missing(p)) {
            if (!(floor(p) == p)) stop("p must be integer valued")
            if (p > length(xw$inv)) stop("p too large for xw")
            xw$vt <- xw$vt[1:p,]
            xw$inv <- xw$inv[1:p]
            xw$w <- xw$w[,1:p]
        }
        trans <- xw
    }
    z <- trans$w

    n <- nrow(z)
    if(missing(p)) p <- ncol(z)
    if(missing(m)) m <- floor(sqrt(n))

    # some error checking
    if(num_loadings > p) {
        warning("num_loadings = ", num_loadings, " must be less than p = ", p, 
                    ". Set num_loadings = p - 1")
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
        cat("/// Sorting these into clusters", "\n", sep="")
        # do we want to save all directions in each cluster,
        # or just the best (pre-optim)
        if(!missing(size_clust) && (size_clust > 1)) {
            best_dirs <- cluster.norm(z = z, IC=IC, k=k, m=m,
                                      dirs=rand_dir, kmeans_iter=kmeans_iter,
                                      clusters=clusters, save.all=TRUE)
        } else {
            best_dirs <- cluster.norm(z = z, IC=IC, k=k, m=m,
                                      dirs=rand_dir, kmeans_iter=kmeans_iter,
                                      clusters=clusters)
        }

        cat("//// Sorted into ", length(best_dirs), " clusters", "\n", sep="")
        entr.pre_optim <- sapply(best_dirs, function(x) x$entr)
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
        cat("//// Optimised direction has entropy ", ica_loading$dir_entr, "\n", sep="")
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

    res <- list(xw=trans, x=z %*% IC, IC=IC, entr=entr, m=m)
    class(res) = "goodICA"
    res
}
