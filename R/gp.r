##' Constructor function for gpr S3 class
##'
##' Constructor function for gpr S3 class
##' @title Constructor function for gpr S3 class
##' @param pst object to be coerced to S3 gpr class
##' @return S3 gpr class object
##' @author Wajid Jawaid
gpr <- function(pst) {
    if (!is.list(pst$observed)) stop("Observed values not in list being coerced to S3 gpr class")
    if (!is.numeric(pst$observed$x)) stop("Training inputs not numeric or not in list being coerced to S3 gpr class")
    if (!is.numeric(pst$observed$y)) stop("Training outputs not numeric or not in list being coerced to S3 gpr class")
    if (!is.numeric(pst$x)) stop("Test inputs not numeric or present in list being coerced to S3 gpr class")
    if (is.null(pst$mean)) stop("Posterior means not numeric or present in list being coerced to S3 gpr class")
    if (!is.matrix(pst$cv) && !is.numeric(pst$cv)) stop("Posterior covariance not numeric matrix or present in list being coerced to S3 gpr class")
    if (!is.numeric(pst$marginal_likelihood)) stop("Marginal likelihood not numeric or present in list being coerced to S3 gpr class")
    class(pst) <- c("gpr", class(pst))
    return(pst)
}

##' Covariance function
##'
##' Covariance function for - The squared exponential covariance function
##' \eqn{k_{x,y} = \sigma_f^2 exp(-\frac{1}{2l^2}(x - y)^2)}
##' @title Covariance function
##' @param x First value
##' @param y Second value
##' @param sigmaf See equation above
##' @param l See equation above
##' @return Calculated covariance
##' @author Wajid Jawaid
##' @export
##' @examples
##' set.seed(0)
##' library(gpr)
##'
##' d <- .1
##' x <- seq(-5,5,d)
##' K <- covMat(x, x, function(xx,xy) cf(xx,xy, l = 2))
##' L <- cholesky(K)
##'
##' plot(NULL, xlim=c(-5,5), ylim=c(-2.5,2.5), xlab = "time", ylab = "counts")
##' rect(-6, -2, 6, 2, col="#44444411", lty=0)
##' for (i in 1:4) {
##'   y <- sampleFromL(L)
##'   points(x,y, pch=16, xlim=c(-5,5), ylim=c(-2,2), cex=.2, col = i, type = "l")
##' }
cf <- function(x,y, sigmaf = 1, l = 1) {
    z <- sigmaf^2 * exp(- (1/(2 * l^2)) * (sum(x - y))^2)
    z
}

##' Zero mean function
##'
##' Zero mean function
##' @title Zero mean function
##' @param x Value of independent variable
##' @return Returns 0 for all values
##' @author Wajid Jawaid
mf <- function(x) {
    return(0)
}

##' Calculate Covariance matrix
##'
##' Calculate Covariance matrix
##' @title Calculate Covariance matrix
##' @param x Rows of covariance matrix
##' @param y Columns of covariance matrix
##' @param cvFunc Covariance function to be used see \code{\link{cf}}
##' @param sigman See equation above
##' @return Returns a covariance matrix
##' @author Wajid Jawaid
##' @export
##' @examples
##' set.seed(0)
##' library(gpr)
##'
##' d <- .1
##' x <- seq(-5,5,d)
##' K <- covMat(x, x, function(xx,xy) cf(xx,xy, l = 2))
##' L <- cholesky(K)
##'
##' plot(NULL, xlim=c(-5,5), ylim=c(-2.5,2.5), xlab = "time", ylab = "counts")
##' rect(-6, -2, 6, 2, col="#44444411", lty=0)
##' for (i in 1:4) {
##'   y <- sampleFromL(L)
##'   points(x,y, pch=16, xlim=c(-5,5), ylim=c(-2,2), cex=.2, col = i, type = "l")
##' }
covMat <- function(x, y, cvFunc, sigman = 0) {
    cm <- matrix(apply(expand.grid(x, y), 1, function(xa) cvFunc(xa[1], xa[2])), length(x), dimnames = list(x,y))
    diag(cm) <- diag(cm) + sigman^2
    cm
}

##' Sample from Covariance matrix
##'
##' Sample from covariance matrix.
##' To avoid repeated Cholesky decompositions please supply with lower triangular
##' Cholesky decomposition of the covariance matrix rather than the covariance
##' matrix itself.
##' @title Sample from Covariance matrix
##' @param L Cholesky decomposition of Covariance matrix
##' @param mu Means to use for sampling
##' @return Vector of predictions
##' @author Wajid Jawaid
##' @importFrom stats rnorm
##' @export
##' @examples
##' set.seed(0)
##' library(gpr)
##'
##' d <- .1
##' x <- seq(-5,5,d)
##' K <- covMat(x, x, function(xx,xy) cf(xx,xy, l = 2))
##' L <- cholesky(K)
##'
##' plot(NULL, xlim=c(-5,5), ylim=c(-2.5,2.5), xlab = "time", ylab = "counts")
##' rect(-6, -2, 6, 2, col="#44444411", lty=0)
##' for (i in 1:4) {
##'   y <- sampleFromL(L)
##'   points(x,y, pch=16, xlim=c(-5,5), ylim=c(-2,2), cex=.2, col = i, type = "l")
##' }
sampleFromL <- function(L, mu = 0) {
    if (!all(L[upper.tri(L)] == 0)) stop("L must be lower triangular from the cholesky decomposition")
    u <- rnorm(ncol(L), 0, 1)
    y <- mu + L %*% u
}

##' Cholesky decomposition
##'
##' Cholesky decomposition
##' The R function returns an upper triangular matrix. Additionally this adds a little noise
##' to the diagonal
##' @title Cholesky decomposition
##' @param K Covariance matrix - symmetric positive definite
##' @param sigma2 Noise term
##' @param ... Passed to R's built-in function
##' @return Lower triangular matrix
##' @author Wajid Jawaid
##' @export
##' @examples
##' set.seed(0)
##' library(gpr)
##'
##' d <- .1
##' x <- seq(-5,5,d)
##' K <- covMat(x, x, function(xx,xy) cf(xx,xy, l = 2))
##' L <- cholesky(K)
##'
##' plot(NULL, xlim=c(-5,5), ylim=c(-2.5,2.5), xlab = "time", ylab = "counts")
##' rect(-6, -2, 6, 2, col="#44444411", lty=0)
##' for (i in 1:4) {
##'   y <- sampleFromL(L)
##'   points(x,y, pch=16, xlim=c(-5,5), ylim=c(-2,2), cex=.2, col = i, type = "l")
##' }
cholesky <- function(K, sigma2 = 1e-12, ...) {
    L <- t(base::chol((K + diag(sigma2, nrow(K), ncol(K))), ...))
}

##' GP engine - Finds posterior
##'
##' GP engine - Finds posterior
##' @title GP engine - Finds posterior
##' @param x Measured input values
##' @param y Measured output values
##' @param xs Where to make predictions
##' @param cvFunc Which covariance function to use
##' @param mFunc Default 0. Mean function
##' @param sigmaf Passed to \code{\link{cf}}
##' @param l Passed to \code{\link{cf}}
##' @param sigman Passed to \code{\link{cf}}
##' @param verbose To print progress
##' @return Returns a list containing:
##' \describe{
##'   \item{x}{The values where predcitions are made}
##'   \item{mean}{Posterior mean}
##'   \item{cv}{Posterior covariance matrix}
##'   \item{marinal_likelihood}{Marginal likelihood of the model}
##' }
##' @author Wajid Jawaid
##' @export
##' @examples
##' x <- c(-4,-3,-1,0,2)
##' y <- c(-2,0,1,2,-1)
##'
##' pst <- predictGP(x, y, xs = seq(-5,5,.1), cvFunc = cf, sigman = 0.1, l = 1)
##' plot(pst, p.pch = 3, p.cex = 2, xlab = "x", ylab = "y")
predictGP <- function(x, y, xs, cvFunc, mFunc = mf, sigmaf = 1, l = 1, sigman = 0,
                      verbose = FALSE) {
    vcat <- function(xx) if (verbose) cat(xx) 
    vcat("Calculating training covariances ... ")
    kxx <- covMat(x,x, function(xx, xy) cvFunc(xx, xy, sigmaf = sigmaf, l = l), sigman = sigman)
    vcat("done.\nCalculating training to test covariances ... ")
    kxxs <- covMat(x,xs, function(xx, xy) cvFunc(xx, xy, sigmaf = sigmaf, l = l))
    vcat("done.\nCalculating test covariances ... ")
    kxsxs <- covMat(xs,xs, function(xx, xy) cvFunc(xx, xy, sigmaf = sigmaf, l = l))
    vcat("done.\nCalculating means at training values... ")
    mx <- mFunc(x)
    vcat("done.\nCalculating means at test values... ")
    mxs <- mFunc(xs)
    vcat("done.\nCalculating Cholesky ... ")
    L <- cholesky(kxx)
    vcat("done.\nCalculating precision matrix ... ")
    sL <- solve(t(L)) %*% solve(L)
    vcat("done.\nCalculating posterior means ... ")
    mean <- mxs + (t(kxxs) %*% sL %*% (y - mx))
    vcat("done.\nCalculating posterior covariance matrix ... ")
    cv <- kxsxs - t(kxxs) %*% sL %*% kxxs
    vcat("done.\nCalculating marginal likelihood ... ")
    mlik <- -.5 * (t(y) %*% sL %*% y + log(2*pi)) - sum(log(diag(L)))
    vcat("done.\nReturning ... \n")
    pst <- list(observed = list(x = x, y = y), x = xs, mean = mean[,1], cv = cv, marginal_likelihood = mlik[1,1])
    return(gpr(pst))
}

##' Plot Gaussian process regression object
##'
##' Plot Gaussian process regression object
##' Plots gpr class
##' @title Plot Gaussian process regression object
##' @param x gpr class
##' @param xl x-axis lower limit
##' @param xu x-axis upper limit
##' @param yl y-axis lower limit
##' @param yu y-axis upper limit
##' @param ci Default 0.975 (5\% CI). Confidence interval to plot
##' @param ci.col Default "#CCCCCC". Colour of shaded CI interval
##' @param p.pch Default 1. Shape to use for training points
##' @param p.cex Default 0.5. Size of training points
##' @param ci.line.col Default "darkgray". Colour of confidence interval shaded region outline
##' @param ci.line.lty Default 2. Line thickness of confidence interval shaded region outline
##' @param xaxs Deafult "i".
##' @param xlab Passed to plot
##' @param ylab Passed to plot
##' @param ... Parameters passed to plot e.g. main
##' @return Returns a plot
##' @author Wajid Jawaid
##' @importFrom graphics box par plot points polygon
##' @importFrom stats qnorm
##' @rdname plot
##' @method plot gpr
##' @export
plot.gpr <- function(x, xl = NULL, xu = NULL, yl = NULL, yu = NULL, ci = 0.975,
                     ci.col = "#CCCCCC", p.pch = 1, p.cex = 0.5, ci.line.col = "darkgray",
                     ci.line.lty = 2, xaxs = "r", xlab = "x", ylab = "y", ...) {
    ci <- qnorm(0.975) * sqrt(diag(x$cv))
    xl <- c(ifelse(is.null(xl), min(x$observed$x, x$x), xl),
            ifelse(is.null(xu), max(x$observed$x, x$x), xu))
    yl <- c(ifelse(is.null(yl), min(c(x$observed$y, x$mean - ci)), yl),
            ifelse(is.null(yu), max(c(x$observed$y, x$mean + ci)), yu))
    par(oma = c(0,0,0,0))
    plot(NULL, xlim = xl, ylim = yl, xaxs = xaxs, xlab = xlab, ylab = ylab, ...)
    points(x$x, x$mean, type = "p", col = "darkgray", pch = 1, cex = .2)
    polygon(c(x$x, rev(x$x)), c(x$mean + ci, rev(x$mean - ci)), col = ci.col, lty = 0)
    points(x$observed$x, x$observed$y, pch = p.pch, cex = p.cex)
    points(x$x, x$mean, type="l", col = "#FF4444", cex=5)
    points(x$x, x$mean + ci, col=ci.line.col, lty=ci.line.lty, type="l")
    points(x$x, x$mean - ci, col=ci.line.col, lty=ci.line.lty, type="l")
    box()
}
