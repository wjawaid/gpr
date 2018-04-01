##' Covariance function
##'
##' Covariance function for - The squared exponential covariance function
##' \eqn{k_{x,y} = \sigma_f^2 exp(-frac{1}{2l^2}(x - y)^2) + \delta_{xy} \sigma_n^2}
##' @title Covariance function
##' @param x First value
##' @param y Second value
##' @param sigmaf See equation above
##' @param l See equation above
##' @param sigman See equation above
##' @return Calculated covariance
##' @author Wajid Jawaid
##' @export
cf <- function(x,y, sigmaf = 1, l = 1, sigman = 0) {
    z <- sigmaf^2 * exp(- (1/(2 * l^2)) * (sum(x - y))^2)
    if (x == y) z <- z + sigman^2
    z
}

##' Calculate Covariance matrix
##'
##' Calculate Covariance matrix
##' @title Calculate Covariance matrix
##' @param x Rows of covariance matrix
##' @param y Columns of covariance matrix
##' @param cvFunc Covariance function to be used see \code{\link{cf}}
##' @return Returns a covariance matrix
##' @author Wajid Jawaid
##' @export
covMat <- function(x, y, cvFunc) {
    cm <- matrix(apply(expand.grid(x, y), 1, function(xa) cvFunc(xa[1], xa[2])), length(x), dimnames = list(x,y))
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
sampleFromK <- function(L, mu = 0) {
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
##' @param σ2 Noise term
##' @param ... Passed to R's built-in function
##' @return Lower triangular matrix
##' @author Wajid Jawaid
##' @export
cholesky <- function(K, σ2 = 1e-12, ...) {
    L <- t(base:::chol((K + diag(σ2, nrow(K), ncol(K))), ...))
}

##' GP engine - Finds posterior
##'
##' GP engine - Finds posterior
##' @title GP engine - Finds posterior
##' @param x Measured input values
##' @param y Measured output values
##' @param xs Where to make predictions
##' @param cvFunc Which covariance function to use
##' @param sigmaf Passed to \code{\link{cf}}
##' @param l Passed to \code{\link{cf}}
##' @param sigman Passed to \code{\link{cf}}
##' @return Returns a list containing:
##' \describe{
##'   \item{x}{The values where predcitions are made}
##'   \item{mean}{Posterior mean}
##'   \item{cv}{Posterior covariance matrix}
##'   \item{marinal_likelihood}{Marginal likelihood of the model}
##' }
##' @author Wajid Jawaid
predictGP <- function(x, y, xs, cvFunc, sigmaf = 1, l = 1, sigman = 0) {
    kxx <- covMat(x,x, function(xx, xy) cvFunc(xx, xy, sigmaf = sigmaf, sigman = sigman, l = l))
    kxxs <- covMat(x,xs, function(xx, xy) cvFunc(xx, xy, sigmaf = sigmaf, l = l))
    kxsxs <- covMat(xs,xs, function(xx, xy) cvFunc(xx, xy, sigmaf = sigmaf, l = l))
    L <- cholesky(kxx)
    sL <- solve(t(L)) %*% solve(L)
    mean <- t(kxxs) %*% sL %*% y
    cv <- kxsxs - t(kxxs) %*% sL %*% kxxs
    mlik <- -.5 * (t(y) %*% sL %*% y + log(2*pi)) - sum(log(diag(L)))
    return(list(x = xs, mean = mean, cv = cv, marginal_likelihood = mlik))
}
