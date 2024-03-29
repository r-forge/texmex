\name{migpd}
\alias{migpd}
\alias{plot.migpd}
\title{ Fit multiple independent generalized Pareto models }
\description{
	Fit multiple independent generalized Pareto models as the first step
  of conditional multivariate extreme values modelling following the
  approach of Heffernan and Tawn, 2004.
}
\usage{
migpd(data, mth, mqu, penalty = "gaussian", maxit = 10000,
      trace = 0, verbose = FALSE, priorParameters = NULL)
\method{plot}{migpd}(x, main=c("Probability plot","Quantile plot","Return level plot","Histogram and density"), xlab=rep(NULL,4), nsim=1000, alpha=.05, ... )
}
\arguments{
  \item{data}{ A matrix or data.frame, each column of which is to be modelled. }
  \item{mth}{ Thresholds above which to fit the models. Only one of
				\code{mth} and \code{mqu} should be supplied. Length one (in which case a common threshold is used) or length equal to the number of coluumns of \code{data} (in which case values correspond to thresholds for each of the columns respectively).}
  \item{mqu}{ Quantiles above which to fit the models. Only one of
				\code{mth} and \code{mqu} should be supplied. Length as for \code{mth} above. }
  \item{penalty}{ How the likelihood should be penalized. Defaults to
		'gaussian'. See documentation for \code{\link{gpd}}.}
  \item{maxit}{ The maximum number of iterations to be used by the
		optimizer.}
  \item{trace}{Whether or not to tell the user how the optimizer is getting on. The argument is passed into \code{\link{optim}} -- see the help for that function.}
  \item{verbose}{ Controls whether or not the function print to screen
		every time it fits a model. Defaults to FALSE.}
  \item{priorParameters}{ Only used if
			\code{penalty = 'gaussian'}. A named list, each element of which contains
			two components: the first should be a vector of length 2
			corresponding to the location of the Gaussian distribution;
			the second should be 2x2 matrix corresponding to the 
			covariance matrix of the distribution. The names should match the
			names of the columns of \code{data}. If not provided, 
			it defaults to independent priors being centred at zero, with variance
			10000 for log(sigma) and 0.25 for xi. See the details section.}
  \item{x}{ Object of class \code{migpd} as returned by function \code{migpd}.}
  \item{main}{Character vector of length four: titles for plots produced by \code{plot} method.}
  \item{xlab}{As \code{main} but for x-axes labels.}
  \item{nsim}{Number of simulations on which to base tolerance envelopes in \code{plot} method.}
  \item{alpha}{Significance level for tolerance and confidence intervals in \code{plot} method.}
  \item{...}{Further arguments to be passed to methods.}
}
\details{
  The parameters in the generalized Pareto distribution are
  estimated for each column of the data in turn, independently
  of all other columns. Note, covariate modelling of GPD parameters is not supported.
  
  Maximum likelihood estimation often fails with generalized
  Pareto distributions because of the likelihood becoming flat
  (see, for example, Hosking et al, 1985).
  Therefore the function allows penalized likelihood estimation,
  which is the same as maximum a posteriori estimation from a 
  Bayesian point of view. 
  
  By default quadratic penalization is used,
  corresponding to using a Gaussian prior. If no genuine
  prior information is available, the following argument can
  be used. If xi = -1, the generalized Pareto distribution 
  corresponds to the uniform distribution, and if xi is 1 or
  greater, the expectation is infinite. Thefore, xi is likely
  to fall in the region (-1, 1). A Gaussian distribution
  centred at zero and with standard deviation 0.5 will have 
  little mass outside of (-1, 1) and so will often be a
  reasonable prior for xi. For log(sigma) a Gaussian
  distribution, centred
  at zero and with standard deviation 100 will often be vague.
  If a Gaussian penalty is specified but no parameters are
  given, the function will assume such indpendent priors.
  
  Note that internally the function works with log(sigma), not
  sigma. The reasons are that quadratic penalization makes more
  sense for log(sigma) than for sigma (because the distribution
  of log(sigma) will be more nearly symmetric), and because
  it was found to stabilize computations.
  
  The associated \code{coef}, \code{print} and \code{summary}
  functions exponentiate the log(sigma) parameter to return results
  on the expected scale. If you are accessesing the parameters
  directly, however, take care to be sure what scale the results
  are on.
  
  Threshold selection can be carried out with the help of functions \code{\link{mrlPlot}} and \code{\link{gpdRangeFit}}.
}
\value{
  An object of class 'migpd'. There are \code{coefficients}, \code{print}, \code{plot}
  and \code{summary} functions available.
}
\references{
J. E. Heffernan and J. A. Tawn, A conditional approach
	for multivariate extreme values, Journal of the Royal Statistical
	society B, 66, 497 -- 546, 2004

J. R. M. Hosking and J. R. Wallis, Parameter and quantile estimation for the genralized Pareto distribution, Technometrics, 29, 339 -- 349, 1987
}
\author{ Harry Southworth }
\seealso{ \code{\link{mex}}, \code{\link{mexDependence}}, \code{\link{bootmex}},
	\code{\link{predict.mex}}, \code{\link{gpdRangeFit}}, \code{\link{mrlPlot}}
		  }
\examples{
mygpd <- migpd(winter, mqu=.7, penalty = "none")
mygpd
summary(mygpd)
plot(mygpd)
}
\keyword{ models }
\keyword{ multivariate }

