# produce random directions, and choose the 'out' best directions
# best directions are those that minimise entropy
# The value associated with the less ``important'' whitening loadings
# have more probability of being zero
randDirs <- function(z, IC, k, m, iter=5000, out) {
    p <- ncol(IC)

    r <- p - k + 1 # the dimension of the search space

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
    trialsProj <- trialsOrigSpace %*% t(z[,1:p])
    entr <- mSpacingEntropy(trialsProj, m=m)

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
                        kmean.iter, save.all=FALSE, clust.avg=FALSE) {
    p <- ncol(IC)

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
                entr <- mSpacingEntropy(centreProj, m=m)
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
                     cluster, opt.method="BFGS") {
        # if optim method is L-BFGS-B, then use upper bound
        if(opt.method == "L-BFGS-B") {
                opt <- optim(par = dirs,
                             fn=function(w) {
                                     w <- w / sqrt(sum(w^2))
                                     wOrigSpace <- IC %*% c(rep(0, k-1), w)
                                     zProj <- t(z %*% wOrigSpace)
                                     mSpacingEntropy(zProj, m = m)
                             }, 
                             gr=function(w) {
                                     w <- w / sqrt(sum(w^2))
                                     wOrigSpace <- IC %*% c(rep(0, k-1), w)
                                     zProj <- t(z %*% wOrigSpace)
                                     entropyGradOrigSpace <- optimEntropyDeriv(xProj=zProj, x=z, m=m)
                                     #TODO: is this correct?
                                     if(k > 1) {
                                        # use chain rule to obtain w \in \R^r
                                        r <- length(w)
                                        zeroMatrix <- matrix(0, nrow = (k-1), ncol = r)
                                        jacobianG <- z %*% IC %*% rbind(zeroMatrix, diag(r))
                                        entropyGrad <- t(jacobianG) %*% entropyGradOrigSpace
                                     } else {
                                        entropyGradOrigSpace
                                     }
                             },
                             lower=-Inf, upper=(0.5 * (log(2 * pi) + 1)),
                             method = opt.method, control = list(maxit = maxit, trace=0))
        } else {
                opt <- optim(par = dirs,
                             fn=function(w) {
                                     w <- w / sqrt(sum(w^2))
                                     wOrigSpace <- IC %*% c(rep(0, k-1), w)
                                     zProj <- t(z %*% wOrigSpace)
                                     mSpacingEntropy(zProj, m = m)
                             }, 
                             gr=function(w) {
                                     w <- w / sqrt(sum(w^2))
                                     wOrigSpace <- IC %*% c(rep(0, k-1), w)
                                     zProj <- t(z %*% wOrigSpace)
                                     entropyGradOrigSpace <- optimEntropyDeriv(xProj=zProj, x=z, m=m)
                                     #TODO: is this correct?
                                     if(k > 1) {
                                        # use chain rule to obtain \delta w \in \R^r
                                        r <- length(w)
                                        zeroMatrixTop <- matrix(0, nrow = (k-1),
                                                                ncol = r)
                                                                #ncol = (r + k - 1))
                                        #zeroMatrixLeft <- matrix(0, nrow = r, ncol = (k-1))
                                        paddedI <- rbind(zeroMatrixTop, diag(r))
                                                        #cbind(zeroMatrixLeft, diag(r)))
                                        #jacobianG <- z %*% IC %*% paddedI
                                        # with u = IC %*% (rep(0, k-1), w), want du/dw
                                        dudw <- IC %*% paddedI
                                        #entropyGrad <- t(jacobianG) %*% entropyGradOrigSpace
                                        entropyGrad <- entropyGradOrigSpace %*% dudw
                                     } else {
                                        entropyGradOrigSpace
                                     }
                                     #entropyGrad <- optimEntropyDeriv(xProj=zProj, x=z, m=m)
                                     #if(k > 1) {
                                             #entropyGrad[-(1:(k-1))] 
                                     #} else {
                                             #entropyGrad
                                     #}
                             },
                             method = opt.method, control = list(maxit = maxit, trace=0))
        }
        if (opt$convergence == 1) {
                warning("In loading ", k, ", cluster ", cluster, 
                        " optimisation did not converge, consider increasing maxit \n")
        } else if (opt$convergence != 0) {
                warning("In loading ", k, ", cluster ", cluster, 
                        " optimisation did not converge (error ", opt$convergence, ") \n")
        }
        entrTmp <- opt$value
        dirTmp <- opt$par
        dirTmp <- dirTmp / sqrt(sum(dirTmp^2))
        # output
        output <- list()
        output$entr <- entrTmp
        output$dirs <- dirTmp
        output
}

# create a single ICA loading from clustered random projections
# input is from clusterNorm
icaClusters <- function(z, IC, k, m, best.dirs, maxit=1000,
                        opt.method="BFGS",
                        clust.avg=FALSE, verbose=FALSE) {
    n <- nrow(z)
    p <- ncol(IC)

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
        # TODO: remove?
            samp <- seq_len(nTmp)
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
    output$dirEntr <- dirOpt[clusterNum, 1]
    output$dirOptim <- dirOpt[clusterNum, -1]
    if (any(nn > 1)) {
        return(list(best=output, all=dirOptMany))
    } else {
        return(output)
    }
}

householderTransform <- function(IC, bestDir, r, k, p) {
        # Use a Householder reflection which maps e1 to best.dir to update IC.
        e1 <- c(1, rep(0, r-1))
        # take sign of x_k s.t.
        # k is the last col entry of non-zero in UT form A = QR
        signTmp <- sign(bestDir[1])
        v <- bestDir - signTmp * e1
        v <- v / sqrt(sum(v^2))
        P <- diag(r) - 2 * tcrossprod(v)
        IC[,k:p] <- IC[,k:p,drop=FALSE] %*% P
        IC
}

# for class clusterICA
#' @export
print.clusterICA <- function(x, ...) {
    loadings <- ncol(x$r)
    length <- nrow(x$r)
    entr1 <- round(x$entropy[1], digits = 5)
    cat("Cluster ICA: ", loadings, " loading(s) found of length ", length,
        ". Best projection has entropy ", entr1, ".\n", sep="")
    invisible(x)
}
