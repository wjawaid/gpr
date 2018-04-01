\name{predictGP}
\alias{predictGP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
predictGP(x, y, xs, cvFunc, sigmaf = 1, l = 1, sigman = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{xs}{
%%     ~~Describe \code{xs} here~~
}
  \item{cvFunc}{
%%     ~~Describe \code{cvFunc} here~~
}
  \item{sigmaf}{
%%     ~~Describe \code{sigmaf} here~~
}
  \item{l}{
%%     ~~Describe \code{l} here~~
}
  \item{sigman}{
%%     ~~Describe \code{sigman} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, y, xs, cvFunc, sigmaf = 1, l = 1, sigman = 0) 
{
    kxx <- covMat(x, x, function(xx, xy) cvFunc(xx, xy, sigmaf = sigmaf, 
        sigman = sigman, l = l))
    kxxs <- covMat(x, xs, function(xx, xy) cvFunc(xx, xy, sigmaf = sigmaf, 
        l = l))
    kxsxs <- covMat(xs, xs, function(xx, xy) cvFunc(xx, xy, sigmaf = sigmaf, 
        l = l))
    L <- cholesky(kxx)
    sL <- solve(t(L)) \%*\% solve(L)
    mean <- t(kxxs) \%*\% sL \%*\% y
    cv <- kxsxs - t(kxxs) \%*\% sL \%*\% kxxs
    mlik <- -0.5 * (t(y) \%*\% sL \%*\% y + log(2 * pi)) - sum(log(diag(L)))
    return(list(x = xs, mean = mean, cv = cv, marginal_likelihood = mlik))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
