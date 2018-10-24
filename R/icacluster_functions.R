entropy <- function(x, m, sortt=TRUE) {
    if(is.vector(x)) x <- matrix(x, nrow=1)
    if(ncol(x) == 1 && nrow(x) > 1) {
        # transpose matrx is 1 column
        #cat("assumed n = 1 and  p = ", nrow(x), "\n")
        x <- t(x)
    }
    x <- t(apply(x, 1, function(x) sort(x, method="radix")))
    n <- ncol(x)
    # if x is a vector///
    if(missing(m)) m <- floor(sqrt(n))
    
    d <- x[,(m+1):n, drop=FALSE] - x[,1:(n-m), drop=FALSE]
    apply(log(n * d / m), 1, sum) / n - digamma(m) + log(m)
}


rand.dirs <- function(z, IC, k, m, iter=5000, out, seed) {
    p <- ncol(z)
    n <- nrow(z)
    if(missing(m)) m <- floor(sqrt(n))
    if(missing(k)) k <- 1
    if(missing(IC)) IC <- diag(p)
    r <- p - k + 1 # the dimension of the search space

    if (!missing(seed)) set.seed(seed)
    trials_mat <- matrix(rnorm(r*iter), iter, r)
    trials_mat <- trials_mat / sqrt(rowSums(trials_mat^2))
    trials.orig.space <- trials_mat %*% t(IC[,k:p])
    # switch to columns for each trial so that entr works
    trials.proj <- trials.orig.space %*% t(z)
    entr <- entropy(trials.proj, m=m, sortt=TRUE)

    dir.table <- cbind(entr, trials_mat)
    # arange in order
    dir.table <- dir.table[order(dir.table[,1]),]
    namesW <- paste0('dir', seq_len(iter))
    if(!missing(out)) {
        if(out > iter) {
            warning("out > iter: have set out = iter")
            out <- iter
        }
        dir.table <- dir.table[1:(out),]
        namesW <- paste0('dir', seq_len(out))
    }
    
    entr <- dir.table[,1]
    dirs <- dir.table[,-1]

    rownames(dirs) <- namesW
    colnames(dirs) <- NULL
    output <- list()
    output$entr <- entr
    output$dirs <- dirs
    output
}


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
            c <- proj.cluster(X=dirs, K=i, maxiter=100,
                                plot=FALSE, verbose=FALSE)
            wss[i] <- proj.wss(X=dirs, c)
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
    c <- proj.cluster(X = dirs, K = clusters, maxiter=kmeans_iter,
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

dir.optim <- function(z, IC, k, m, dirs, maxit=1000, 
                        cluster, opt_method="Nelder-Mead") {
    n <- ncol(z)

    opt <- optim(par = dirs,
                 function(w) {
                     w <- w / sqrt(sum(w^2))
                     w.orig.space <- IC %*% c(rep(0, k-1), w)
                     z_proj <- z %*% w.orig.space
                     entropy(z_proj, m = m)
                 }, method = opt_method, control = list(maxit = maxit))
    
    if (opt$convergence == 1) {
        warning("In cluster ", cluster, " optimisation did not converge, consider increasing maxit")
    } else if (opt$convergence != 0) {
        warning("In cluster ", cluster, " optimisation did not converge (error ", opt$convergence, ")")
    }
    
    entr_tmp <- opt$value
    dir_tmp <- opt$par
    dir_tmp <- dir_tmp / sqrt(sum(dir_tmp^2))
    
    output <- list()
    output$entr <- entr_tmp
    output$dirs <- dir_tmp
    output
}

ica.clusters <- function(z, IC, k, m, best_dirs, maxit=1000,
                         opt_method="Nelder-Mead", size_clust) {
    n <- nrow(z)
    p <- ncol(z)
    if(missing(m)) m <- floor(sqrt(n))
    if(missing(IC)) IC <- diag(p)
    if(missing(k)) k <- 1

    clusters <- length(best_dirs)
    cat("////Optimising direction of projection on ", clusters, " clusters \n")

    dir_opt <- matrix(nrow = clusters, ncol = (p  - k + 1 + 1))
    dir_opt_many <- vector(mode="list", length=clusters)
    nn <- numeric()
    for(i in 1:clusters) {
        cat("//// Optimising cluster ", i, "\n")
        dir_tmp <- best_dirs[[i]]
        n_tmp <- length(dir_tmp$entr)
        nn[i] <- n_tmp
        if(n_tmp == 1) {
            dir_opt_tmp <- dir.optim(z = z, IC = IC, dirs = dir_tmp$dirs,
                                     k = k, m = m, maxit = maxit,
                                     cluster=i, opt_method=opt_method)

        } else {
            # randomly choose size_clust dirs to optimise in each cluser
            if(is.numeric(size_clust)) {
                samp <- sample(n_tmp, size = min(size_clust, n_tmp))
            } else {
                samp <- seq_len(n_tmp)
            }
            dir_opt_clust <- lapply(samp, function(j) {
                dirr <- dir_tmp$dirs[j,]
                dir_opt_tmp <- dir.optim(z = z, IC = IC, dirs = dirr,
                                         k = k, m = m, maxit = maxit, cluster=i,
                                         opt_method=opt_method)
            })
            dir_entr_tmp <- sapply(dir_opt_clust, function(x) x$entr)
            dir_dir_tmp <- t(sapply(dir_opt_clust, function(x) x$dir))
            names_tmp <- names(dir_tmp$entr)
            dir_table <- cbind(dir_entr_tmp, dir_dir_tmp)
            ord_tmp <- order(dir_table[,1])
            dir_table <- dir_table[ord_tmp,]
            names_tmp <- names_tmp[ord_tmp]
            dir_entr_tmp <- dir_table[,1]
            names(dir_entr_tmp) <- names_tmp
            dir_dir_tmp <- dir_table[,-1]
            min_entr <- which.min(dir_entr_tmp)

            dir_opt_tmp <- list(entr=dir_entr_tmp[min_entr], 
                                    dirs=dir_dir_tmp[min_entr,])
            dir_opt_many[[i]] <- list(entr=dir_entr_tmp, dirs=dir_dir_tmp)
        }

        dir_opt[i,] <- c(dir_opt_tmp$entr, dir_opt_tmp$dirs)
    }
    cluster_num <- which.min(dir_opt[,1])
    output <- list()
    output$cluster_num <- cluster_num
    output$dir_entr <- dir_opt[cluster_num, 1]
    output$dir_optim <- dir_opt[cluster_num, -1]
    if (any(nn > 1)) {
        return(list(best=output, all=dir_opt_many))
    } else {
        return(output)
    }
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

print.goodICA <- function(y, ...) {
    loadings <- ncol(y$IC)
    length <- nrow(y$IC)
    entr1 <- round(y$entr[1], digits = 5)
    cat("Cluster ICA: ", loadings, " loading(s) found of length ", length,
        ". Best projection has entropy ", entr1, ".\n", sep="")
    invisible(y)
}
