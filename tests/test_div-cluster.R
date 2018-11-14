library(clusterICA)

generateData <- function(n1, n2, d) {
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

testData <- list(
    `gen(40,10,6)`=function() generateData(40, 10, 6),
    `gen(80,20,6)`=function() generateData(80, 20, 6),
    `gen(800,200,6)`=function() generateData(800, 200, 6)
)


testMethod <- function(fn) {
    cat("\n\n*** testing ", deparse(substitute(fn)), "\n\n", sep="")
    times <- numeric(0)
    clusters <- integer(0)
    RSS <- numeric(0)
    for (test.idx in 1:length(testData)) {
        X <- testData[[test.idx]]()
        M <- 100
        t0 <- proc.time()["elapsed"]
        for (i in 1:M) {
            cl <- fn(X)
        }
        t1 <- proc.time()["elapsed"]
        times <- c(times, (t1 - t0) / M)
        clusters <- c(clusters, length(unique(cl)))
        RSS <- c(RSS, 0)
    }

    res <- data.frame(`time [ms]`=1000*times, clusters, RSS, check.names=FALSE)
    row.names(res) <- names(testData)
    print(res)
}

testMethod(function(X) clusterProjDivisive(X, tol=0.1))

#######

par(mfrow = c(2,2))
for (i in 1:4) {
    X1 <- generateData(276, 89, 6)
    c <- clusterProjDivisive(X1, tol=0.1)

    diffc <- c$wssAll[-length(c$wssAll)] - c$wssAll[-1]
    plot(diffc / c$wssAll[1])
}
