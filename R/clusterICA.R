# ICA function{{{
#' Independent Component Analysis (ICA)
#'
#' Uses random directions, clustering and optimisation to obtain ICA loadings,
#'  using the m-spacing entropy estimation as the objective function to minimise.
#'
#' @param x the data to perform ICA on, ncol(x) = n, nrow(x) = p
#' @param m (optional) the value of m-spacing for calculating approximate entropy, if missing(m), m <- sqrt(n)
#' @param p.ica the number of ICA loadings outputted
#' @param p.whiten (optional) the size of the whitened matrix, i.e. how many PCA loadings to keep in the whitening step
#' @param rand.iter the number of random directions to initialise
#' @param kmean.tol the tolerance used in divisive clustering, see clusterProjDivisive
#' @param opt.maxit the maximum number of iterations used in the optimisation step, see optim
#' @param opt.method the method used in the optimisation step, see optim
#' @param fast.init if TRUE then the objective function used in the fastICA method
#' is optimised to find an extra direction to cluster 
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
#' a <- clusterICA(x=xGood, p.ica=1, rand.iter=1000)
#' par(mfrow = c(1,3))
#' plot(xGood, main = "Pre-processed data")
#' plot(a$y, main = "Whitened data")
#' plot(density(a$s, bw="sj"), main = "ICA components")
#'
#' #---------------------------------------------------
#' #Example 2: un-mixing two mixed independent uniforms
#' #From fastICA man page
#' #---------------------------------------------------
#' S <- matrix(runif(10000), 5000, 2)
#' A <- matrix(c(1, 1, -1, 3), 2, 2, byrow = TRUE)
#' X <- S %*% A
#' a <- clusterICA(X, p.whiten=2, rand.iter=1000)
#' par(mfrow = c(1, 3))
#' plot(X, main = "Pre-processed data")
#' plot(a$y, main = "Whitened data")
#' plot(a$s, main = "ICA components")
#'
#' #---------------------------------------------------
#' #Example 3: un-mixing iris data
#' #---------------------------------------------------
#' a <- clusterICA(iris[,1:4])
#' plot(a$s, main = "ICA components")
#' pairs(a$s, col=iris$Species)
#' @export
#/*}}}*/

clusterICA <- function(x, m=-1, p.ica, p.whiten, rand.iter,
                        kmean.tol=0.1, opt.maxit=5000, opt.method="BFGS",
                        fast.init=TRUE, compute.scores = TRUE, verbose=FALSE) {
        # check if we have whitened data
        # here p is how many PCA loadings we use to do ICA on
        # p.ica is how many ICA loadings we want outputted
        whitened <- FALSE
        if(class(x) == "coords") {
                xw <- x
                z <- x$y
                zCov <- cov(z)
                zCovZero <- zCov - diag(nrow(zCov))
                covDelta <- 1e-13
                zCovZeroVector <- as.vector(zCovZero)
                if (all(abs(zCovZeroVector) < covDelta)) {
                        whitened <- TRUE
                }
        }
        if(whitened == FALSE && class(x) == "coords") {
                if(verbose == TRUE) {
                        cat("Data 'x' of type 'coords' but not whitened: 
                            whitening using jvcoords")
                }
                xw <- jvcoords::whiten(x$y, compute.scores=TRUE)
                z <- xw$y
        } else {
                if(verbose == TRUE) {
                        cat("Data 'x' not whitened: whitening using jvcoords")
                }
                if(!(class(x) == "coords")) {
                        xw <- jvcoords::whiten(x, compute.scores=TRUE)
                        z <- xw$y
                }
        }
        n <- nrow(z)
        if(missing(p.whiten)) {
                p <- ncol(z)
        } else {
                p <- p.whiten
                z <- z[,1:p]
        }
        if(missing(p.ica)) p.ica <- p
        stopifnot(p.ica <= p)
        if(m == -1) m <- floor(sqrt(n))
        if(missing(rand.iter)) rand.iter <- max(5000, min(35000, 2^(p / 4.5)))
        rand.out <- min(100+p, rand.iter)
        # if fastICA obj function used for initialisation
        if (fast.init == TRUE) normSamp <- rnorm(1e5)
        # initiate loadings list
        loadings <- vector(mode = "list", length = p.ica)
        entr <- numeric()
        IC <- diag(p)
        loopNum <- p.ica
        if (loopNum == p) loopNum <- p - 1

        k <- 1
        while (k <= loopNum) {
                if (verbose == TRUE) {
                        cat("optimising direction", k, "out of", p.ica, "\n")
                }
                r <- p - k + 1 # the dimension of the search space
                if (verbose == TRUE) {
                        cat("// Finding random starting points", "\n")
                }
                randDir <- randDirs(z=z, IC=IC, k=k, m=m, iter=rand.iter, out=rand.out)
                if (fast.init == TRUE) {
                        fastDir <- fastICAInitialisation(z=z, IC=IC, m=m, k=k, norm.sampl=normSamp)
                        randDir$dirs <- rbind(randDir$dirs, fastDir$dir)
                        randDir$entr <- c(randDir$entr, fastDir$entr)
                }
                if (verbose == TRUE) {
                        cat("/// Found ", length(randDir$entr), " starting directions", "\n",
                            sep="")
                        cat("/// Sorting these into clusters \n")
                }
                bestDirs <- clusterNorm(z = z, IC=IC, k=k, m=m,
                                        dirs=randDir, kmean.tol=kmean.tol,
                                        kmean.iter=200)
                if (verbose == TRUE) {
                        cat("//// Sorted into ", length(bestDirs), " clusters", "\n", sep="")
                }
                if (verbose == TRUE) {
                        entrPreOptim <- unlist(sapply(bestDirs, function(x) x$entr))
                        cat("//// Best pre-optim entropy = ", min(entrPreOptim), "\n", sep="")
                }

                # step 2: use local optimisation to find the best solution in the
                # each cluster
                if (verbose == TRUE) {
                        cat("//// Optimising ", length(bestDirs), " clusters", "\n", sep="")
                }
                icaLoading <- icaClusters(z=z, IC=IC, k=k, m=m,
                                          best.dirs=bestDirs, maxit = opt.maxit,
                                          opt.method=opt.method,
                                          verbose=verbose)
                if (verbose == TRUE) {
                        cat("//// Optimised direction has entropy ",
                            icaLoading$dirEntr, "\n", sep="")
                }
                bestEntr <- icaLoading$dirEntr
                bestDir <- icaLoading$dirOptim

                # is this projection better than any previous projections?
                if(any(bestEntr < entr)) {
                        k_tmp <- min(which(bestEntr < entr))
                        lenBestDir <- length(bestDir)
                        r_tmp <- (p - k_tmp + 1)
                        bestDir <- c(rep(0, times=(r_tmp-lenBestDir)), bestDir)
                        trialsOrigSpace <- bestDir %*% t(IC[,k_tmp:p])
                        # switch to columns for each trial so that entr works
                        trialsProj <- trialsOrigSpace %*% t(z[,1:p])
                        newEntr <- mSpacingEntropy(trialsProj, m=m)
                        k <- k_tmp
                        r <- p - k + 1 # the dimension of the search space
                        dirTmp <- vector("list", length=1)
                        dirTmp[[1]]$dirs <- bestDir
                        dirTmp[[1]]$entr <- newEntr
                        icaLoading <- icaClusters(z=z, IC=IC, k=k, m=m,
                                                  best.dirs=dirTmp, maxit = opt.maxit,
                                                  opt.method=opt.method, verbose=verbose)
                        bestDir <- icaLoading$dirOptim
                        bestEntr <- icaLoading$dirEntr
                        entr <- entr[1:k_tmp]
                        if (verbose == TRUE) {
                                cat("///// Current projection better than ", k, "th projection", "\n")
                                cat("///// Replacing ", k, "th projection", "\n")
                        }
                }

                if (verbose == TRUE) cat("///// Householder reflection", "\n")
                entr[k] <- bestEntr
                # Use a Householder reflection which maps e1 to best.dir to update IC.
                IC <- householderTransform(IC=IC, bestDir=bestDir, r=r, k=k, p=p)
                k <- k + 1
        }

        if (p.ica == p) {
                bestDir <- 1
                wOrigSpace <- IC %*% c(rep(0, p-1), bestDir)
                zProj <- t(z %*% wOrigSpace)
                bestEntr <- mSpacingEntropy(zProj, m = m)
                entr[p] <- bestEntr
        }

        IC <- IC[, seq_len(p.ica), drop=FALSE]

        colnames(IC) <- paste0('wIC', seq_len(p.ica))
        # want unmixing matrix for unwhitened data x
        # x with column means removed
        if(class(x) == "coords") {
                X <- jvcoords::fromCoords(xw, xw$y)
                Xt <- t(X) - xw$shift
        } else {
                Xt <- t(x) - xw$shift
        }
        rownames(Xt) <- paste0('Xcentred', seq_len(nrow(Xt)))
        # make into class "coords"
        xw <- jvcoords::appendTrfm(xw, op = "orth", IC)
        # W s.t. S = X %*% t(W)
        # here, X has each column mean subtracted
        W <- t(IC) %*% diag(xw$cmds[[2]][[2]]) %*% t(xw$cmds[[1]][[2]])
        rownames(W) <- paste0('IC', seq_len(nrow(W)))
        colnames(W) <- NULL 
        # R s.t. S = Y %*% R
        R <- IC
        xw$y <- z
        xw$x <- t(Xt)
        xw$w <- W
        xw$r <- R
        if(compute.scores == TRUE) {
                S <- z %*% IC
                colnames(S) <- paste0('S', seq_len(ncol(S)))
                xw$s <- S
        }
        xw$entropy <- entr
        class(xw) = c("clusterICA", "coords")
        xw
}

