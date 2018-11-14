# produce random directions, and choose the 'out' best directions
# best directions are those that minimise entropy
# The value associated with the less ``important'' whitening loadings
# have more probability of being zero
randDirs <- function(z, IC, k, m, iter=5000, out, seed) {
    p <- ncol(z)
    r <- p - k + 1 # the dimension of the search space

    if (!missing(seed)) set.seed(seed)
    trialsMat <- matrix(rnorm(r*iter), iter, r)
    # lets try with some elements zero
    # seemed to work well when tried a while back
    # probs means that the smaller PC loadings are more
    # likely to be ignored.
    probs <- seq(from=1/r, to=1, length=r)
    trialsMat <- t(apply(trialsMat, 1, function(trials) {
        # always want at least two non-zero elements
        # otherwise would just get the PC loading back
        sampp <- sample(1:r, size=sample(1:(r-2), 1),
                        replace=FALSE, prob=probs)
        trials[sampp] <- 0
        trials
    }))
    trialsMat <- trialsMat / sqrt(rowSums(trialsMat^2))
    trialsOrigSpace <- trialsMat %*% t(IC[,k:p])
    # switch to columns for each trial so that entr works
    trialsProj <- trialsOrigSpace %*% t(z)
    entr <- entropy(trialsProj, m=m)

    dirTable <- cbind(entr, trialsMat)
    # arange in order
    dirTable <- dirTable[order(dirTable[,1]),]
    namesW <- paste0('dir', seq_len(iter))
    if(!missing(out)) {
        if(out > iter) {
            warning("out > iter: have set out = iter")
            out <- iter
        }
        dirTable <- dirTable[1:(out),]
        namesW <- paste0('dir', seq_len(out))
    }

    entr <- dirTable[,1]
    dirs <- dirTable[,-1]

    rownames(dirs) <- namesW
    colnames(dirs) <- NULL
    output <- list()
    output$entr <- entr
    output$dirs <- dirs
    output
}



# put random directions into clusters
# uses divisive kmeans clustering from clusterProjDivisive
clusterNorm <- function(z, IC, k, m, dirs, kmean.tol=0.1,
                        kmean.iter=100, save.all=FALSE, clust.avg=FALSE) {
    p <- ncol(z)

    entr <- dirs$entr
    dirs <- dirs$dirs

    # K-Means Cluster Analysis: Divisive
    c <- clusterProjDivisive(X=dirs, tol=kmean.tol, iter.max=kmean.iter)
    clusters <- max(c$c)

    # append cluster assignment & put into list
    outTmp <- vector(mode = "list", length = clusters)
    dirsClusterAppend <- cbind(c$c, entr, dirs)
    for(i in 1:clusters) {
        whichCluster <- which(dirsClusterAppend[,1] == i)
        if (save.all == FALSE & clust.avg==FALSE) {
            outTmp[[i]]$entr <- dirsClusterAppend[whichCluster, 2]
            entrMin <- which.min(outTmp[[i]]$entr)
            outTmp[[i]]$entr <- outTmp[[i]]$entr[entrMin]
            outTmp[[i]]$dirs <- dirsClusterAppend[whichCluster, c(-1, -2),
                                                  drop=FALSE]
            outTmp[[i]]$dirs <- outTmp[[i]]$dirs[entrMin,]
        } else {
            outTmp[[i]]$entr <- dirsClusterAppend[whichCluster, 2]
            outTmp[[i]]$dirs <- dirsClusterAppend[whichCluster, c(-1, -2)]
            #TODO: Remove clust.avg?
            if (clust.avg == TRUE) {
                s <- La.svd(outTmp[[i]]$dirs, nu=0, nv=1)
                centre <- s$vt[1,]
                outTmp[[i]]$dirs <- centre
                # calc entropy of centre
                centreOrigSpace <- centre %*% t(IC[,k:p])
                centreProj <- centreOrigSpace %*% t(z)
                entr <- entropy(centreProj, m=m)
                outTmp[[i]]$entr <- entr
            }
        }
    }
    outTmp
    return(outTmp)
}

# optimise each direction
# here dir is a single direction (vector)
# cluster arg only used for cat() in clusterICA
dirOptim <- function(z, IC, k, m, dirs, maxit=1000,
                     cluster, opt.method="Nelder-Mead") {
    n <- ncol(z)

    opt <- optim(par = dirs,
                 function(w) {
                     w <- w / sqrt(sum(w^2))
                     wOrigSpace <- IC %*% c(rep(0, k-1), w)
                     zProj <- t(z %*% wOrigSpace)
                     entropy(zProj, m = m)
                 }, method = opt.method, control = list(maxit = maxit))

    if (opt$convergence == 1) {
        warning("In loading", k, ", cluster ", cluster, " optimisation did not converge, consider increasing maxit \n")
    } else if (opt$convergence != 0) {
        warning("In loading", k, ", cluster ", cluster, " optimisation did not converge (error ", opt$convergence, ") \n")
    }

    entrTmp <- opt$value
    dirTmp <- opt$par
    dirTmp <- dirTmp / sqrt(sum(dirTmp^2))

    output <- list()
    output$entr <- entrTmp
    output$dirs <- dirTmp
    output
}

# create a single ICA loading from clustered random projections
# input is from clusterNorm
icaClusters <- function(z, IC, k, m, best.dirs, maxit=1000,
                        opt.method="Nelder-Mead", size.clust,
                        clust.avg=FALSE, verbose=FALSE) {
    n <- nrow(z)
    p <- ncol(z)

    clusters <- length(best.dirs)
    if (verbose == TRUE) {
        cat("////Optimising direction of projection on ",
            clusters, " clusters \n")
    }

    dirOpt <- matrix(nrow = clusters, ncol = (p  - k + 1 + 1))
    dirOptMany <- vector(mode="list", length=clusters)
    nn <- numeric()
    for(i in 1:clusters) {
        if (verbose == TRUE) {
            cat("//// Optimising cluster ", i, "\n")
        }
        dirTmp <- best.dirs[[i]]
        nTmp <- length(dirTmp$entr)
        nn[i] <- nTmp
        if(nTmp == 1) {
            dirOptTmp <- dirOptim(z = z, IC = IC, dirs = dirTmp$dirs,
                                  k = k, m = m, maxit = maxit,
                                  cluster=i, opt.method=opt.method)

        } else {
            # randomly choose size.clust dirs to optimise in each cluser
            if(is.numeric(size.clust)) {
                samp <- sample(nTmp, size = min(size.clust, nTmp))
            } else {
                samp <- seq_len(nTmp)
            }
            dirOpt_clust <- lapply(samp, function(j) {
                dirr <- dirTmp$dirs[j,]
                dirOptTmp <- dirOptim(z = z, IC = IC, dirs = dirr,
                                      k = k, m = m, maxit = maxit, cluster=i,
                                      opt.method=opt.method)
            })
            dirEntrTmp <- sapply(dirOpt_clust, function(x) x$entr)
            dirDirTmp <- t(sapply(dirOpt_clust, function(x) x$dir))
            names_tmp <- names(dirTmp$entr)
            dirTable <- cbind(dirEntrTmp, dirDirTmp)
            ord_tmp <- order(dirTable[,1])
            dirTable <- dirTable[ord_tmp,]
            names_tmp <- names_tmp[ord_tmp]
            dirEntrTmp <- dirTable[,1]
            names(dirEntrTmp) <- names_tmp
            dirDirTmp <- dirTable[,-1]
            entrMin <- which.min(dirEntrTmp)

            dirOptTmp <- list(entr=dirEntrTmp[entrMin],
                              dirs=dirDirTmp[entrMin,])
            dirOptMany[[i]] <- list(entr=dirEntrTmp, dirs=dirDirTmp)
        }

        dirOpt[i,] <- c(dirOptTmp$entr, dirOptTmp$dirs)
    }
    clusterNum <- which.min(dirOpt[,1])
    output <- list()
    output$clusterNum <- clusterNum
    output$dir_entr <- dirOpt[clusterNum, 1]
    output$dirOptim <- dirOpt[clusterNum, -1]
    if (any(nn > 1)) {
        return(list(best=output, all=dirOptMany))
    } else {
        return(output)
    }
}

# for class clusterICA
#' @export
print.clusterICA <- function(x, ...) {
    loadings <- ncol(x$IC)
    length <- nrow(x$IC)
    entr1 <- round(x$entr[1], digits = 5)
    cat("Cluster ICA: ", loadings, " loading(s) found of length ", length,
        ". Best projection has entropy ", entr1, ".\n", sep="")
    invisible(x)
}
