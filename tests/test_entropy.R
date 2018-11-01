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
    `gen(40,10,6)`=function() t(generate.data(40, 10, 6)),
    `gen(80,20,6)`=function() t(generate.data(80, 20, 6)),
    `gen(800,200,6)`=function() t(generate.data(800, 200, 6))
)

test.method <- function(fn) {
    cat("\n\n*** testing ", deparse(substitute(fn)), "\n\n", sep="")
    times <- numeric(0)
    clusters <- integer(0)
    RSS <- numeric(0)
    for (test.idx in 1:length(test.data)) {
        X <- test.data[[test.idx]]()
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
    row.names(res) <- names(test.data)
    print(res)
}

test.method(function(X) entropy(X))


# check for any numerical errors
set.seed(1234)
X <- matrix(rep(-4:5,by=1,each=100) + rnorm(1000, sd=10e-12), ncol=100, byrow=TRUE)
entropy(X)

# check that matrix = rbind() gives same as indivdual inputs
X1 <- matrix(rnorm(150), ncol=15, nrow=10)
X2 <- matrix(rnorm(150), ncol=10, nrow=15)
X3 <- matrix(rnorm(1500), ncol=10, nrow=150)
X4 <- matrix(rnorm(1500), ncol=100, nrow=15)
X5 <- matrix(rnorm(15000), ncol=100, nrow=150)
X6 <- matrix(rnorm(15000), ncol=1000, nrow=15)

X <- list(X1, X2, X3, X4, X5, X6)
for (i in 1:length(X)) {
    Xi <- X[[i]] 
    Xi.entr_mat <- entropy(Xi)
    for (j in 1:nrow(Xi)) {
        Xi.entrr <- entropy(Xi[j,])
        stopifnot(Xi.entrr == Xi.entr_mat[j])
    }
}


# check: entropy tends towards true when n,m -> \Inf with m/n -> 0
# use m <- sqrt(n)
entrm <- numeric()
iter <- 100
for (i in 1:150) {
    iter[i+1] <- 1.05*iter[i]
    mm <- floor(sqrt(iter[i+1]))
    set.seed(10)
    Xx <- rnorm(iter[i+1])
    entrm[i] <- entropy(Xx, m = 3)
}
plot(iter[-1], entrm, type ="l", log='x', ylim=c(1,1.42))
# true value of entropy
abline(h = 0.5*(log(2*pi)+1), col = "red")

# check: variance of entropy approx as m increases
entrm <- numeric()
var_all <- numeric()
iter <- 1000
len_mm <- 100
mm <- seq(from = log(3)/log(iter), to = 0.9, length = len_mm)

for (j in 1:len_mm) {
    varm <- numeric()
    for (k in 1:50) {
        mm_tmp <- floor(iter^mm[j])
        for (i in 1:50) {
            set.seed(k*(13420 + i*iter))
            Xx <- rnorm(iter)
            entrm[i] <- entropy(Xx, m = mm_tmp)
        }
        varm[k] <- var(entrm)
    }
    var_all[j] <- mean(varm) 
}
# steps in plot relate to floor()
plot(x=mm, var_all, t="l")

# check: density of entropy approx for given m
entrm <- numeric()
iter <- 1000
len_mm <- 100
mm <- seq(from = log(3)/log(iter), to = 0.99, length = len_mm)
varm <- vector("list", length=len_mm)
for (k in 1:len_mm) {
    mm_tmp <- floor(iter^mm[k])
    for (i in 1:100) {
        set.seed(k*(14420 + i*iter))
        Xx <- rgamma(iter, shape = 0.5)
        entrm[i] <- entropy(Xx, m = mm_tmp)
    }
    varm[[k]] <- entrm
}
# steps in plot relate to floor()
plot(density(varm[[1]], bw="sj"))
points(density(varm[[20]], bw="sj"), t="l", col="red")
points(density(varm[[42]], bw="sj"), t="l", col="green")
points(density(varm[[80]], bw="sj"), t="l", col="blue")
xx <- 1:len_mm
yy <- sapply(xx, function(i) var(varm[[i]]))
plot(xx, yy, type="l")


# different m's
# larger m's //seem to// make entropy less accurate
# but smooth out local minima
Xx1 <- rnorm(1000)
Xx2 <- runif(1000)
Xx3 <- rgamma(1000, shape=0.4)
Xx <- rbind(Xx1, Xx2, Xx3)
Xx <- t(apply(Xx, 1, function(x) scale(x, scale=TRUE)))
par(mfrow = c(3,1))
for (i in 1:3) {
Xm <- numeric()
Xxx <- Xx[i,]
    for (j in 2:(length(Xxx) - 3)) {
        Xm[j-1] <- entropy(Xxx, m=j)
        #if (j/100 == floor(j/100)) cat(j)
    }
    plot(x = 2:(length(Xxx) - 3), y = Xm, type = "b")
}






