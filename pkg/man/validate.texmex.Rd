\name{validate.texmex}
\alias{validate.texmex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Validate the \code{texmex} package
}
\description{
Using \code{svUnit}, run tests to compare behaviour of \code{texmex} with
expected behaviour
}
\usage{
validate.texmex()
}
%- maybe also 'usage' for other objects documented here.
\details{
The function uses functionality in the \code{svUnit} package to perform many
tests to ensure that \code{texmex} behaves as expected. These include tests
to make sure that \code{texmex} can reproduce published results (often parameter
estimates and model summaries that appear in Coles (2000) or in Heffernan
and Tawn (2004)), logical tests (for
example, to ensure that penalizing the likelihood has the expected effect on
parameter estimates), tests of the structure of objects, tests of proper
vectorization, tests to ensure that simulated posterior means are close to
MAP estimates, and tests to ensure that the code gives the same results as
independently coded calculations (sometimes using code that appears in the
\code{ismev} and \code{evd} packages, sometimes using code written by someone
else and without reference to the \code{texmex} version.

If all is well, the result should indicate that no errors or fails have
occurred.
}
\value{
An object of class \code{svSuiteData}. See the help for \code{svUnit} for
details.
}
\references{
S. Coles, An Introduction to the Statistical Modeling of Extreme Values, Springer, 2000

J. E. Heffernan and J. A. Tawn, A conditional approach for multivariate extreme
values. Journal of the Royal Statistical Society, Series B, volume 66, 2004
}
\author{
Janet E. Heffernan, Harry Southworth
}
\note{
Over 100 tests are performed, some of which run MCMC or bootstrap algorithms.
As a result, the function takes a while to run.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{texmex}}
}
\examples{
val <- validate.texmex()
summary(val)
val
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line