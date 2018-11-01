# ICA function
#' Independent Component Analysis (ICA)
#'
#' Uses random directions, clustering and optimisation to obtain ICA loadings,
#'  using the m-spacing entropy estimation as the objective function to minimise.
#'
#' @param x the data to perform ICA on, ncol(x) = n, nrow(x) = p
#' @param xw (optional) the whitened version of x
#' @param m (optional) the value of m-spacing for calculating approximate entropy, if missing(m), m <- sqrt(n)
#' @param num_loadings the number of ICA loadings outputted
#' @param p (optional) the size of the whitened matrix, i.e. how many PCA loadings to keep in the whitening step
#' @param rand_iter the number of random directions to initialise
#' @param rand_out the number of the best random directions to keep
#' @param seed (optional) the set.seed number used for initialising the random directions
#' @param kmeans_tol the tolerance used in divisive clustering, see cluster.proj.divisive
#' @param kmeans_iter the maximum number of iterations used in divisive clustering, see cluster.proj.divisive
#' @param optim_maxit the maximum number of iterations used in the optimisation step, see optim
#' @param opt_method the method used in the optimisation step, see optim
#' @param size_clust (optional) if size_clust = k > 1, then optimisation is performed on k random directions in each cluster. 
#'                      If missing, then optimisation is performed on the best direction in each cluster.
#'
#' @return A list with the following components:
#'         \itemize{
#'              \item{xw} {The output from jvmulti::whiten(x)}
#'              \item{IC} {The matrix of the loading vectors for the whitened data}
#'              \item{y} {The matrix of the projections of the whitened data along the loading vectors}
#'              \item{entr} {The m-spacing entropy of each row of y}
#'          }
#'
#' @author Paul Smith, \email{mmpws@@leeds.ac.uk}
# #' @seealso cluster.proj.kmeans()
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
#' x1_good <- x1[which(good1),]
#' x2 <- matrix(rnorm(n*p, mean = -2.5), n , p)
#' x2[,2] <- scale(x2[,2])
#' good2 <- cos(20 * x2 %*% a1) >= 0
#' x2_good <- x2[which(good2),]
#' x_good <- rbind(x1_good, x2_good)
#' a <- clusterICA(x=x_good, num_loadings=1, rand_iter=1000, seed=5)
#' par(mfrow = c(1,3))
#' plot(x_good, main = "Pre-processed data")
#' plot(a$xw$y, main = "Whitened data")
#' plot(density(a$y, bw="sj"), main = "ICA components")
#'
#' #---------------------------------------------------
#' #Example 2: un-mixing two mixed independent uniforms
#' #From fastICA man page
#' #---------------------------------------------------
#' S <- matrix(runif(10000), 5000, 2)
#' A <- matrix(c(1, 1, -1, 3), 2, 2, byrow = TRUE)
#' X <- S %*% A
#' a <- clusterICA(X, p=2, rand_iter=1000, rand_out=50)
#' par(mfrow = c(1, 3))
#' plot(X, main = "Pre-processed data")
#' plot(a$xw$y, main = "Whitened data")
#' plot(a$y, main = "ICA components")
#'
#' #---------------------------------------------------
#' #Example 3: un-mixing iris data
#' #---------------------------------------------------
#' a <- clusterICA(iris[,1:4])
#' plot(a$y, main = "ICA components")
#' pairs(a$y, col=iris$Species)
#' @export
clusterICA <- function(x, xw, m, num_loadings, p, rand_iter=5000, rand_out=100, seed, 
                    kmeans_tol=0.1, kmeans_iter=100,
                    optim_maxit=1000, opt_method="Nelder-Mead",
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
    rownames(IC) <- rownames(xw$loadings)

    res <- list(xw=xw, IC=IC, y=z %*% IC, entr=entr)
    class(res) = "clusterICA"
    res
}