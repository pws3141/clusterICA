# check if data is whitened or of class "coords"
whitenCheck <- function(x, verbose = FALSE) {
        whitened <- FALSE
        # 1) check whether 'coords' class
        # and if so, if it is actually whitened
        if(class(x) == "coords") {
                xw <- x
                z <- x$y
                zCov <- cov(z)
                zCovZero <- zCov - diag(nrow(zCov))
                covDelta <- 1e-10
                zCovZeroVector <- as.vector(zCovZero)
                if (all(abs(zCovZeroVector) < covDelta)) {
                        whitened <- TRUE
                }
        }
        if(whitened == FALSE) { 
        # 2) class 'coords' but not whitened
                if(class(x) == "coords") {
                        if(verbose == TRUE) {
                                cat("Data 'x' of type 'coords' but not whitened: 
                                    whitening using jvcoords")
                        }
                        xw <- jvcoords::whiten(x$y, compute.scores=TRUE)
                        z <- xw$y
                } else {
        # 3) not class 'coords'
                        if(verbose == TRUE) {
                                cat("Data 'x' not whitened: whitening using jvcoords")
                        }
                        xw <- jvcoords::whiten(x, compute.scores=TRUE)
                        z <- xw$y
                }
        }
        res <- list(xw=xw, z=z)
        res
}

# fastICA initialisation
# find a w using the fastICA objective function
# use this alongside random directions
fastICAInitialisation <- function(z, IC, m, k, norm.sampl) {
        n <- nrow(z)
        nNorm <- length(norm.sampl)
        p <- ncol(z)
        r <- p - k + 1 # the dimension of the search space
        # start with random direction
        w <- rnorm(r)
        w <- w / sqrt(sum(w^2))
        # optim
        opt <- optim(w, function(w) {
                        w <- w / sqrt(sum(w^2))
                        wProj <- IC %*% c(rep(0, k-1), w)
                        xOrigSpace <- z %*% wProj
                        TermOne <- (1 / n) * sum(log(cosh(xOrigSpace)))
                        TermTwo <- (1 / nNorm) * sum(log(cosh(norm.sampl)))
                        output <- (TermOne - TermTwo)^2
                        -output
                        }, method = "BFGS") 
        trial <- opt$par
        trial <- trial / sqrt(sum(trial^2))
        wProj <- IC %*% c(rep(0, k-1), trial)
        xOrigSpace <- z %*% wProj
        # switch to columns for each trial so that entr works
        entropy <- mSpacingEntropy(t(xOrigSpace), m=m)
        res <- list(dir = trial, entropy = entropy)
        res
}

# produce random directions, and choose the 'out' best directions
# best directions are those that minimise entropy
# The value associated with the less ``important'' whitening loadings
# have more probability of being zero
randomSearch <- function(z, IC, k, m, iter, out) {
    p <- ncol(IC)
    r <- p - k + 1 # the dimension of the search space
    trialsMat <- matrix(rnorm(r*iter), iter, r)
    # set some elements to zero
    # probability of zero corresponds to loading number
    # i.e. loadings explaining more variability in the data more
    # likely to be non-zero
    probs <- seq(from=1/r, to=1, length=r)
    trialsMat <- t(apply(trialsMat, 1, function(trials) {
        # need at least two non-zero elements
        # otherwise just get back loadings
        whichIgnore <- sample(1:r, size=sample(1:(r-2), 1),
                                replace=FALSE, prob=probs)
        trials[whichIgnore] <- 0
        trials
    }))
    trialsMat <- trialsMat / sqrt(rowSums(trialsMat^2))
    trialsOrigSpace <- trialsMat %*% t(IC[,k:p])
    # each column corresponds to a trial s.t. 
    # mSpacingEntropy function input is correct
    trialsProj <- trialsOrigSpace %*% t(z[,1:p])
    entr <- mSpacingEntropy(trialsProj, m = m)
    dirTable <- cbind(entr, trialsMat)
    # arange in order
    dirTable <- dirTable[order(dirTable[,1]),]
    namesW <- paste0('dir', seq_len(iter))
    if(!missing(out)) {
        if(out > iter) {
            warning("out > iter: have set out = iter")
            out <- iter
        }
        dirTable <- dirTable[1:out, ]
        namesW <- paste0('dir', seq_len(out))
    }
    entropy <- dirTable[,1]
    dirs <- dirTable[,-1]
    rownames(dirs) <- namesW
    colnames(dirs) <- NULL
    output <- list()
    output$entropy <- entropy
    output$dirs <- dirs
    output
}

# put random directions into clusters
# uses divisive kmeans clustering from clusterProjDivisive
clusterRandomSearch <- function(z, IC, k, m, dirs, kmean.tol,
                        kmean.iter) {
    p <- ncol(IC)
    entr <- dirs$entropy
    dirs <- dirs$dirs
    # K-Means Cluster Analysis: Divisive
    c <- clusterProjDivisive(X=dirs, tol=kmean.tol, iter.max=kmean.iter)
    clusters <- max(c$c)
    # append cluster assignment & put into list
    outTmp <- vector(mode = "list", length = clusters)
    dirsClusterAppend <- cbind(c$c, entropy, dirs)
    for(i in 1:clusters) {
        whichCluster <- which(dirsClusterAppend[,1] == i)
        outTmp[[i]]$entr <- dirsClusterAppend[whichCluster, 2]
        entrMin <- which.min(outTmp[[i]]$entr)
        outTmp[[i]]$entr <- outTmp[[i]]$entr[entrMin]
        outTmp[[i]]$dirs <- dirsClusterAppend[whichCluster, c(-1, -2),
                                          drop=FALSE]
        outTmp[[i]]$dirs <- outTmp[[i]]$dirs[entrMin,]
    }
    outTmp
    return(outTmp)
}

# optimise each direction
# here dir is a single direction (vector)
# cluster arg only used for cat() in clusterICA
.optimiseDirection <- function(z, IC, k, m, dirs, maxit=1000,
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
                        # use chain rule to obtain \delta w \in \R^r
                        r <- length(w)
                        zeroMatrixTop <- matrix(0, nrow = (k-1),
                                                ncol = r)
                        paddedI <- rbind(zeroMatrixTop, diag(r))
                        # with u = IC %*% (rep(0, k-1), w), want du/dw
                        dudw <- IC %*% paddedI
                        entropyGrad <- entropyGradOrigSpace %*% dudw
                        entropyGrad
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
                        if(k > 1) {
                            # use chain rule to obtain \delta w \in \R^r
                            r <- length(w)
                            zeroMatrixTop <- matrix(0, nrow = (k-1),
                                                    ncol = r)
                            paddedI <- rbind(zeroMatrixTop, diag(r))
                            # with u = IC %*% (rep(0, k-1), w), want du/dw
                            dudw <- IC %*% paddedI
                            entropyGrad <- entropyGradOrigSpace %*% dudw
                            entropyGrad
                        } else {
                            entropyGradOrigSpace
                        }
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
# input is from clusterRandomSearch
optimiseAll <- function(z, IC, k, m, clustered.dirs, maxit=1000,
                        opt.method="BFGS", verbose=FALSE) {
    n <- nrow(z)
    p <- ncol(IC)

    clusters <- length(clustered.dirs)
    if (verbose == TRUE) {
        cat("////Optimising direction of projection on ",
            clusters, " clusters \n")
    }

    dirOpt <- matrix(nrow = clusters, ncol = (p  - k + 1 + 1))
    dirOptMany <- vector(mode="list", length=clusters)
    for(i in 1:clusters) {
        if (verbose == TRUE) {
            cat("//// Optimising cluster ", i, "\n")
        }
        dirTmp <- clustered.dirs[[i]]
        nTmp <- length(dirTmp$entr)
        #if (nTmp == 1) {
        dirOptTmp <- .optimiseDirection(z = z, IC = IC, dirs = dirTmp$dirs,
                                k = k, m = m, maxit = maxit,
                                cluster=i, opt.method=opt.method)
        dirOpt[i,] <- c(dirOptTmp$entr, dirOptTmp$dirs)
        #}
    }
    clusterNum <- which.min(dirOpt[,1])
    output <- list()
    output$clusterNum <- clusterNum
    output$dirEntr <- dirOpt[clusterNum, 1]
    output$optimumDirections <- dirOpt[clusterNum, -1]
    return(output)
}

ensureOrder <- function(z, IC, p, m,
                        best.dir, best.entr, entr, 
                        maxit, opt.method, verbose) {
        k <- min(which(best.entr < entr))
        verboseFunction(which.one=5, verbose=verbose, k=k)
        lenBestDir <- length(best.dir)
        r_tmp <- (p - k + 1)
        bestDirOrigSpace <- c(rep(0, times=(r_tmp-lenBestDir)), best.dir)
        trialsOrigSpace <- bestDirOrigSpace %*% t(IC[,k:p])
        # switch to columns for each trial so that entr works
        trialsProj <- trialsOrigSpace %*% t(z[,1:p])
        newEntr <- mSpacingEntropy(trialsProj, m=m)
        r <- p - k + 1 # the dimension of the search space
        dirTmp <- vector("list", length=1)
        dirTmp[[1]]$dirs <- bestDirOrigSpace
        dirTmp[[1]]$entr <- newEntr
        icaLoading <- optimiseAll(z=z, IC=IC, k=k, m=m,
                                  clustered.dirs=dirTmp, maxit = maxit,
                                  opt.method=opt.method, verbose=verbose)
        newDir<- icaLoading$optimumDirection
        newEntr <- icaLoading$dirEntr
        entr <- entr[1:k]
        res <- list(newDir=newDir, newEntr=newEntr, entr=entr, newK=k, newR=r)
        res
}

householderTransform <- function(IC, best.dir, r, k, p) {
        # Use a Householder reflection which maps e1 to best.dir to update IC.
        e1 <- c(1, rep(0, r-1))
        # take sign of x_k s.t.
        # k is the last col entry of non-zero in UT form A = QR
        signTmp <- sign(best.dir[1])
        v <- best.dir - signTmp * e1
        v <- v / sqrt(sum(v^2))
        P <- diag(r) - 2 * tcrossprod(v)
        IC[,k:p] <- IC[,k:p,drop=FALSE] %*% P
        IC
}

verboseFunction <- function(which.one, verbose, k=NA, rand.iter=NA, p.ica=NA, 
                            dir=NA, clustered.dirs=NA, loading=NA) {
        if (verbose == TRUE) {
                if (which.one == 1) {
                        cat("optimising direction", k, "out of", p.ica, "\n")
                        cat("// Finding ", rand.iter, "random starting points", "\n")
                }
                if (which.one == 2) {
                        cat("/// Found ", length(dir$entr), " starting directions", "\n",
                            sep="")
                        cat("/// Sorting these into clusters \n")
                }
                if (which.one == 3) {
                        cat("//// Sorted into ", length(clustered.dirs), " clusters", "\n", sep="")
                        entrPreOptim <- unlist(sapply(clustered.dirs, function(x) x$entr))
                        cat("//// Best pre-optim entropy = ", min(entrPreOptim), "\n", sep="")
                        cat("//// Optimising ", length(clustered.dirs), " clusters", "\n", sep="")
                }
                if (which.one == 4) {
                        cat("//// Optimised direction has entropy ",
                            loading$dirEntr, "\n", sep="")
                }
                if (which.one == 5) {
                        cat("///// Current projection better than ", k, 
                            "th projection", "\n")
                        cat("///// Replacing ", k, "th projection", "\n")
                }
                if (which.one == 6) {
                        cat("///// Householder reflection", "\n")
                }
        }
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
