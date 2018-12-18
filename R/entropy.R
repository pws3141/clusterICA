# calculate entropy using the m-spacing method
#' Entropy estimation using m-spacing
#'
#' Calculates entropy using m-spacing.
#'
#'
#' @param x the data, either a vector or matrix.
#'              If x is a matrix, entropy is estimated for each row separately.
#' @param m (optional) the m-spacing. Defaults to m <- sqrt(n)
#'              if missing, where n is length(x) if x is a vector,
#'                  or ncol(x) if a matrix
#'
#' @return Vector of real numbers corresponding to the approximate
#'              entropy for each row of input x.
#'
#'
#' @author Paul Smith \& Jochen Voss, \email{mmpws@@leeds.ac.uk}
#' @references Beirlant, Jan, et al. "Nonparametric entropy estimation: An overview."
#' @keywords entropy
#'
#' @examples
#' X1 <- matrix(rnorm(150), ncol=15, nrow=10)
#' X2 <- matrix(rnorm(150), ncol=10, nrow=15)
#' X3 <- matrix(rnorm(1500), ncol=10, nrow=150)
#' X4 <- matrix(rnorm(1500), ncol=100, nrow=15)
#' X <- list(X1, X2, X3, X4)
#' XiEntr <- vector("list", length = length(X))
#' for (i in 1:length(X)) {
#'     Xi <- X[[i]]
#'     XiEntr_mat <- mSpacingEntropy(Xi)
#'     XiEntr[[i]] <- XiEntr_mat
#' }
#' str(XiEntr)
#'
#' @export
mSpacingEntropy <- function(x, m) {
    if(is.vector(x)) x <- matrix(x, nrow=1)
    if(ncol(x) == 1) stop("require p > 1")

    # change to xt
    xt <- apply(x, 1, function(x) sort(x))
    n <- nrow(xt)
    if(missing(m)) m <- floor(sqrt(n))

    d <- xt[(m+1):n,, drop=FALSE] - xt[1:(n-m),, drop=FALSE]
    apply(d, 2, function(dd) {
        (1/n) * sum(log((n / m) * dd))
    }) - digamma(m) + log(m)
}

#TODO: unused, delete?
optimEntropy <- function(w, x, m) {
	xw <- as.vector(x %*% w)
	res <- mSpacingEntropy(x=xw, m=m)
	res
}
 
optimEntropyDeriv <- function(xProj, x, m) {
	#nb: here we require w to be a vector
	xw <- as.vector(xProj)
	n <- length(xw)
	if(missing(m)) m <- floor(sqrt(n))
	xwOrd <- order(xw, decreasing = FALSE)
	xSort <- x[xwOrd,]
	xwSort <- xw[xwOrd]
	dxw <- xwSort[(m+1):n] - xwSort[1:(n-m)]
    dx <- xSort[(m+1):n,] - xSort[1:(n-m),]
    res <- (1 / n) * colSums(dx / dxw)
	res
}
