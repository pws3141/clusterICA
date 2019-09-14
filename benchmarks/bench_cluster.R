library(clusterICA)

generate.data <- function(n1, n2, d) {
    x1 <- rnorm(n1, d)
    y1 <- rnorm(n1, 0)
    z1 <- rnorm(n1, 0, 0.1)
    x2 <- rnorm(n2, 8)
    y2 <- rnorm(n2, 8)
    z2 <- rnorm(n2, 0, 0.1)
    X <- rbind(cbind(x1, y1, z1), cbind(x2, y2, z2))
    X <- X * sample(c(-1, 1), size=n1+n2, replace=TRUE)
    X / sqrt(rowSums(X^2))
}

test.data <- list(
    `gen(40,10,6)`=function() generate.data(40, 10, 6),
    `gen(80,20,6)`=function() generate.data(80, 20, 6),
    `gen(800,200,6)`=function() generate.data(800, 200, 6)
)

test.method <- function(fn) {
    cat("\n\n*** testing ", deparse(substitute(fn)), "\n\n", sep="")
    times <- numeric(0)
    clusters <- integer(0)
    rmse <- numeric(0)
    for (test.idx in 1:length(test.data)) {
        X <- test.data[[test.idx]]()
        M <- 100
        t0 <- proc.time()["elapsed"]
        rmseTmpVec <- numeric(0)
        numClustersTmp <- numeric(0)
        for (i in 1:M) {
            cl <- fn(X)
            if (is.list(cl)) {
                    clusterTmp <- cl$c
            } else {
                    clusterTmp <- cl
            }
            rmseTmp <- clusterRMSE(X, c = clusterTmp) 
            rmseTmpVec <- c(rmseTmpVec, rmseTmp$rmse)
            numClustersTmp <- c(numClustersTmp, length(unique(clusterTmp)))
        }
        t1 <- proc.time()["elapsed"]
        times <- c(times, (t1 - t0) / M)
        clusters <- c(clusters, mean(numClustersTmp))
        rmse <- c(rmse, mean(rmseTmpVec))
    }

    res <- data.frame(`time [ms]`=1000*times, clusters, rmse, check.names=FALSE)
    row.names(res) <- names(test.data)
    print(res)
}

test.method(function(X) clusterProjKmeans(X, 3))
test.method(function(X) clusterProjDivisive(X, tol = .1))
