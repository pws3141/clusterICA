# ICA function{{{
#' Independent Component Analysis (ICA)
#'
#' Uses random directions, clustering and optimisation to obtain ICA loadings,
#'  using the m-spacing entropy estimation as the objective function to minimise.
#'
#' @param X the data to perform ICA on, ncol(x) = n, nrow(x) = p
#' @param p.ica the number of ICA loadings outputted
#' @param p.whiten (optional) the size of the whitened matrix, i.e. how many PCA loadings to keep in the whitening step
#' @param rand.iter the number of random directions to initialise
#' @param m (optional) the value of m-spacing for calculating approximate entropy, if missing(m), m <- sqrt(n)
#' @param kmean.tol the tolerance used in divisive clustering, see clusterProjDivisive
#' @param opt.maxit the maximum number of iterations used in the optimisation step, see optim
#' @param opt.method the method used in the optimisation step, see optim
#' @param fast.init if TRUE then the objective function used in the fastICA method
#' is optimised to find an extra direction to cluster 
#' @param compute.scores if TRUE then scores of the whitened data are outputted
#' @param verbose if TRUE then information is given on the status of the function
#' @return 
#' An object of class ‘coords’, with the following additional
#' components added:
#' \itemize{
#'  \item{X} {Centered data matrix}
#'  \item{Y} {Whitened data matrix, found using jvcoords::whiten(x)}
#'  \item{S} {Matrix of source signal estimates}
#'  \item{W} {Estimated unmixing matrix, S = X  t(W)}
#'  \item{R} {Orthogonal rotation matrix, S = Y  R}
#'  \item{entr} {The m-spacing entropy of each column of S}
#' }
#' @author Paul Smith, \email{mmpws@@leeds.ac.uk}
# #' @seealso clusterProjKmeans()
#' @keywords independent component analysis, entropy, clustering
#'
#' @importFrom jvcoords whiten appendTrfm fromCoords
#' @importFrom stats optim rnorm cov
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
#' a <- clusterICA(X=xGood, p.ica=1, rand.iter=1000)
#' par(mfrow = c(1,3))
#' plot(xGood, main = "Pre-processed data")
#' plot(a$Y, main = "Whitened data")
#' plot(density(a$S, bw="sj"), main = "ICA components")
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
#' plot(a$Y, main = "Whitened data")
#' plot(a$S, main = "ICA components")
#'
#' #---------------------------------------------------
#' #Example 3: un-mixing iris data
#' #---------------------------------------------------
#' a <- clusterICA(iris[,1:4])
#' plot(a$S, main = "ICA components")
#' pairs(a$S, col=iris$Species)
#' @export
#/*}}}*/

clusterICA <- function(X, p.ica=-1, p.whiten=-1, rand.iter=-1, m=-1, 
                        kmean.tol=0.1, opt.maxit=5000, opt.method="BFGS",
                        fast.init=NA, compute.scores = TRUE, verbose=FALSE) {
        # check whether the data x is whitened
        white <- whitenCheck(x = X, verbose = verbose)
        # set z to be the whitened data
        z <- white$z
        n <- nrow(z)
        # p.whiten is how many whitened loadings do we do ICA on
        if(p.whiten == -1) {
                p <- ncol(z)
        } else {
                p <- p.whiten
                z <- z[,1:p]
        }
        # p.ica is how many ICA loadings we want in output
        if(p.ica == -1) p.ica <- p
        stopifnot(p.ica <= p)
        # m is spacing used in mSpacingEntropy
        if(m == -1) m <- floor(sqrt(n))
        # rand.iter is size of random search
        if(rand.iter == -1) rand.iter <- max(5000, min(35000, 2^(p / 4.5)))
        # rand.out is number of directions saved from random search
        rand.out <- min(100+p, rand.iter)
        # if fastICA obj function used for initialisation
        if (is.na(fast.init)) {
                fast.init <- FALSE
                if (p > 50) fast.init <- TRUE 
        }
        if (fast.init == TRUE) normSamp <- rnorm(1e5)
        # initiate loadings list
        loadings <- vector(mode = "list", length = p.ica)
        entr <- numeric()
        IC <- diag(p)
        loopNum <- p.ica
        if (loopNum == p) loopNum <- p - 1

        k <- 1
        while (k <= loopNum) {
                r <- p - k + 1 # the dimension of the search space
                verboseFunction(which.one=1, verbose=verbose, rand.iter=rand.iter,
                                k=k, p.ica=p.ica)
                # step 1: initialise using random directions
                # save the rand.out best directions
                randomDirections <- randomSearch(z=z, IC=IC, k=k, m=m, 
                                                 iter=rand.iter, out=rand.out)
                if (fast.init == TRUE) {
                        fastDir <- fastICAInitialisation(z=z, IC=IC, m=m, k=k, 
                                                         norm.sampl=normSamp)
                        randomDirections$dirs <- rbind(randomDirections$dirs, 
                                                       fastDir$dir)
                        randomDirections$entropy <- c(randomDirections$entropy, 
                                                        fastDir$entropy)
                }
                verboseFunction(which.one = 2, verbose = verbose, dir = randomDirections)
                # step 2: cluster the rand.out initial directions
                clusteredDirections <- clusterRandomSearch(z = z, IC = IC, k = k, m = m,
                                        dirs = randomDirections, kmean.tol = kmean.tol,
                                        kmean.iter = 200)
                verboseFunction(which.one = 3, verbose=verbose, 
                                clustered.dirs=clusteredDirections)
                # step 3: use local optimisation to find the best solution in the
                # each cluster
                icaLoading <- optimiseAll(z=z, IC=IC, k=k, m=m,
                                          clustered.dirs=clusteredDirections$directions, 
                                          maxit = opt.maxit, opt.method=opt.method,
                                          verbose=verbose)
                verboseFunction(which.one=4, verbose=verbose, loading=icaLoading)
                bestEntr <- icaLoading$optimumEntropy
                bestDir <- icaLoading$optimumDirection
                # step 4: check whether this projection is better than 
                # any previous projections
                if(any(bestEntr < entr)) {
                        v_k <- matrix(c(rep(0, times = p - r), bestDir), ncol = 1)
                        R.v_k <- IC %*% v_k
                        q_k <- R.v_k
                        Y.q_k <- z %*% R.v_k
                        entr.curr <- mSpacingEntropy(t(Y.q_k), m = m)

                        k.replace <- min(which(entr.curr < entr))
                        r.replace <- p - k.replace + 1
                        
                        IC[, c(k.replace, k)] <- IC[, c(k, k.replace)]

                        # bestDir transormed s.t. Y %*% new IC %*% q_k has entr.curr
                        bestDir_transform <- c(bestDir[1], rep(0, times = k - k.replace),
                                               bestDir[-1])
                        # check
                        v_k_transform <- matrix(c(rep(0, times = p - r.replace), 
                                                   bestDir_transform), ncol = 1)
                        Y.q_k_transform <- z %*% IC %*% v_k_transform
                        all.equal(mSpacingEntropy(t(Y.q_k_transform), m = m), entr.curr)

                        verboseFunction(which.one=6, verbose=verbose)
                        # Use a Householder reflection which maps e1 to best.dir to update IC.
                        e1 <- c(1, rep(0, r.replace - 1))

                        # take sign of x_k s.t.
                        # k is the last col entry of non-zero in UT form A = QR
                        signTmp <- sign(bestDir_transform[1])

                        v <- bestDir_transform - signTmp * e1
                        v <- v / sqrt(sum(v^2))
                        P <- diag(r.replace) - 2 * tcrossprod(v)

                        IC[, k.replace:p] <- IC[, k.replace:p, drop=FALSE] %*% P

                        entr <- entr[1:k.replace]
                        entr[k.replace] <- entr.curr

                        k <- k.replace
                        r <- r.replace

                } else {
                        verboseFunction(which.one=6, verbose=verbose)
                        entr[k] <- bestEntr

                        # Use a Householder reflection which maps e1 to best.dir to update IC.
                        IC <- householderTransform(IC=IC, best.dir=bestDir, r=r, k=k, p=p)
                }
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
        xw <- white$xw
        if(class(X) == "coords") {
                X <- jvcoords::fromCoords(xw, xw$y)
                Xt <- t(X) - xw$shift
        } else {
                Xt <- t(X) - xw$shift
        }
        rownames(Xt) <- paste0('Xcentred', seq_len(nrow(Xt)))
        # make into class "coords"
        xw <- jvcoords::appendTrfm(xw, op = "orth", IC)
        # W s.t. S = X %*% t(W)
        # here, X has each column mean subtracted
        W <- t(IC) %*% diag(xw$cmds[[2]][[2]][1:p]) %*%
                t(xw$cmds[[1]][[2]][,1:p,drop=FALSE])
        rownames(W) <- paste0('IC', seq_len(nrow(W)))
        colnames(W) <- NULL 
        # R s.t. S = Y %*% R
        R <- IC
        xw$Y <- z
        xw$X <- t(Xt)
        xw$W <- W
        xw$R <- R
        if(compute.scores == TRUE) {
                S <- z %*% IC
                colnames(S) <- paste0('S', seq_len(ncol(S)))
                xw$S <- S
        }
        xw$entropy <- entr
        class(xw) = c("clusterICA", "coords")
        xw
}

