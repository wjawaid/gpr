d <- .1
x <- seq(-5,5,d)
K <- covMat(x, x, function(xx,xy) cf(xx,xy, l = 2))
L <- cholesky(K)

plot(NULL, xlim=c(-5,5), ylim=c(-2.5,2.5), xlab = "time", ylab = "counts")
rect(-6, -2, 6, 2, col="#44444411", lty=0)
for (i in 1:4) {
    y <- sampleFromK(L)
    points(x,y, pch=16, xlim=c(-5,5), ylim=c(-2,2), cex=.2, col = i, type = "l")
}

x <- c(-4,-3,-1,0,2)
y <- c(-2,0,1,2,-1)

pst <- predictGP(x, y, xs = seq(-5,5,.1), cvFunc = cf, sigman = 0.1, l = 1)
par(oma = c(0,0,0,0))
ci <- qnorm(.975) * sqrt(diag(pst$cv))
plot(pst$x, pst$mean, type = "p", col = "darkgray", xlab = "x", ylab = "y", ylim = range(c(pst$mean + ci, rev(pst$mean - ci))), xaxs = "i", pch = 16, cex = .2)
polygon(c(pst$x, rev(pst$x)), c(pst$mean + ci, rev(pst$mean - ci)), col = "#CCCCCC", lty = 0)
points(pst$x, pst$mean, type="l", col = "#444444", pch=16, cex=.2)
points(x, y, pch=3, cex=2)
points(pst$x, pst$mean + ci, col="darkgray", lty=2, type="l")
points(pst$x, pst$mean - ci, col="darkgray", lty=2, type="l")


L <- cholesky(pst$cv)

for (i in 1:10) {
    yo <- sampleFromK(L, mu = pst$mean)
    points(pst$x,yo, pch=16, xlim=c(-5,5), ylim=c(-2,2), cex=.2, col = i, type = "l")
}
