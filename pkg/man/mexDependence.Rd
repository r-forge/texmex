\name{mexDependence}
\alias{mexDependence}
\title{ Estimate the dependence parameters in a conditional multivariate extreme values model }
\description{
  Estimate the dependence parameters in a conditional multivariate extreme
  values model using the approach of Heffernan and Tawn, 2004.
}
\usage{
mexDependence(x, which, dqu, margins="laplace", constrain=TRUE, v = 10, maxit=1000000, start=c(.01, .01), marTransform="mixture", nOptim = 1, PlotLikDo=FALSE, PlotLikRange=list(a=c(-1,1),b=c(-3,1)), PlotLikTitle=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ An object of class 'migpd' as returned by \code{migpd}.}
  \item{which}{ The name of the variable on which to condition. This is
		the name of a column of the data that was passed into
		\code{migpd}.}
  \item{dqu}{ See documentation for this argument in \code{\link{mex}}.}
  \item{margins}{The form of margins to which the data are transformed for carrying out dependence estimation.  Defaults to "laplace", with the alternative option being "gumbel".  The choice of margins has an impact on the interpretation of the fitted dependence parameters.  Under Gumbel margins, the estimated parameters a and b describe only positive dependence, while c and d describe negative dependence in this case.  For Laplace margins, only parameters a and b are estimated as these capture both positive and negative dependence.}
  \item{constrain}{Logical value.  Defaults to \code{constrain=TRUE} although this will subsequently be changed to FALSE if \code{margins="gumbel"} for which constrained estimation is not implemented.  If \code{margins="laplace"} and \code{constrain=TRUE} then the dependence parameter space is constrained to allow only combinations of parameters which give the correct stochastic ordering between (positively and negatively) asymptotically dependent variables and variables which are asymptotically independent.}
  \item{v}{Scalar. Tuning parameter used to carry out constrained estimation of dependence structure under \code{constrain=TRUE}. Takes positive values greater than 1; values between 2 and 10 are recommended.}
  \item{maxit}{The maximum number of iterations to be used by the optimizer.
            defaults to \code{maxit = 1000000}.}
  \item{start}{Optional starting value for dependence estimation.  This can be: a vector of length two, with values corresponding to dependence parameters a and b respectively, and in which case \code{start} is used as a starting value for numerical estimation of each of the dependence models to be estimated; a matrix with two rows corresponding to dependence parameters a and b respectively and number of columns equal to the number of dependence models to be estimated (the ordering of the columns will be as in the original data matrix); or a previously estimated object of class \code{'mex'} whose dependence parameter estimates are used as a starting point for estimation.  Note that under \code{constrain=TRUE}, if supplied, \code{start} must lie within the permitted area of the parameter space.} 
  \item{marTransform}{Optional form of transformation to be used for probability integral transform of data from original to Gumbel or Laplace margins.  Takes values \code{marTransform="mixture"} (the default) or \code{marTransform="empirical"}. When \code{marTransform="mixture"}, the rank transform is used below the corresponding GPD fitting threshold used in \code{x}, and the fitted gpd tail model is used above this threshold.  When \code{marTransform="empirical"} the rank transform is used for the entire range of each marginal distribution.}
  \item{nOptim}{Number of times to run optimiser when estimating dependence model parameters. Defaults to 1.  In the case of \code{nOptim > 1} the first call to the optimiser uses the value \code{start} as a starging point, while subsequent calls to the optimiser are started at the parameter value to which the previous call converged.}
\item{PlotLikDo}{Logical value: whether or not to plot the profile likelihood surface for dependence model parameters under constrained estimation.}
\item{PlotLikRange}{This is used to specify a region of the parameter space over which to plot the profile log-likelihood surface.  List of length 2; each item being a vector of length two corresponding to the plotting ranges for dependence parameters a and b respectively. If this argument is not missing, then \code{PlotLikDo} is set equal to TRUE.}
\item{PlotLikTitle}{Used only if \code{PlotLikDo=TRUE}.  Character string.  Optional title added to the profile log-likelihood surface plot.}
  \item{...}{Further arguments to be passed to methods.}
}
\details{
  Estimates the extremal dependence structure of the data in \code{x}.  The precise nature of the estimation depends on the value of \code{margins}.  If \code{margins="laplace"} (the default) then dependence parameters a and b are estimated after transformation of the data to Laplace marginal distributions.  These parameters can describe both positive and negative dependence.   If \code{margins="gumbel"} then the parameters a, b, c and d in the dependence structure described by Heffernan and Tawn (2004) are estimated in the following two steps: first, a and b are estimated; then, if a=0 and b is negative, parameters c and d are estimated (this is the case of negative dependence). Otherwise c and d will be fixed at zero (this is the case of positive dependence). 
  
  If \code{margins="laplace"} then the option of constrained parameter estimation is available by setting argument \code{constrain=TRUE}.  The default is to constrain the values of the parameters (\code{constrain=TRUE}).  This constrained estimation ensures validity of the estimated model, and enforces the consistency of the fitted dependence model with the strength of extremal dependence exhibited by the data.  More details are given in Keef et al. (2011).  The effect of this constraint is to limit the shape of the dependence parameter space so that its boundary is curved rather than following the original box constraints suggested by Heffernan and Tawn (2004).  The constraint brings with it some performance issues for the optimiser used to estimate the dependence parameters, in particular sensitivity to choice of starting value which we describe now.

  The dependence parameter estimates returned by this function can be particularly sensitive to the choice of starting value used for the optimisation.  This is especially true when \code{margins="laplace"} and \code{constrain=TRUE}, in which case the maximum of the objective function can lie on the edge of the (possibly curved) constrained parameter space.  It is therefore up to the user to check that the reported parameter estimates really do correspond to the maximum of the profile lilkelihood surface.  This is easily carried out by using the visual diagnostics invoked by setting \code{PlotLikDo=TRUE} and adjusting the plotting area by using the argument \code{PlotLikRange} to focus on the region containing the surface maximum.  See an example below which illustrates the use of this diagnostic.
}
\value{
An object of class \code{mex} which is a list containing the following three objects:
  \item{margins}{An object of class \code{\link{migpd}}.}
  \item{dependence}{An object of class \code{\link{mexDependence}}.}
  \item{call}{This matches the original function call.}
}
\references{ J. E. Heffernan and J. A. Tawn, A conditional approach
	for multivariate extreme values, Journal of the Royal Statistical
	society B, 66, 497 -- 546, 2004.
  
  Keef, C., Papastathopoulos, I. and Tawn, J.A., (2011) Additional Constraints for the Conditional Multivariate Extremes Model of Heffernan and Tawn, to appear.}
\author{ Harry Southworth, Janet E. Heffernan }

\seealso{ \code{\link{migpd}}, \code{\link{bootmex}}, \code{\link{predict.mex}}, \code{\link{plot.mex}}
		}
\examples{
data(winter)
mygpd <- migpd(winter , mqu=.7, penalty="none")
mexDependence(mygpd , which = "NO", dqu=.7)

# focus on 2-d example with parameter estimates on boundary of constrained parameter space:
NO.NO2 <- migpd(winter[,2:3] , mqu=.7, penalty="none")

# starting value gives estimate far from true max:
mexDependence(NO.NO2, which = "NO",dqu=0.7,start=c(0.01,0.01),
              PlotLikDo=TRUE,PlotLikTitle=c("NO2 | NO"))

# zoom in on plotting region containing maximum:
mexDependence(NO.NO2, which = "NO",dqu=0.7,start=c(0.01,0.01),
              PlotLikDo=TRUE,PlotLikTitle=c("NO2 | NO"), 
              PlotLikRange = list(a=c(0,0.8),b=c(-0.2,0.6)))

# try different starting value:
mexDependence(NO.NO2, which = "NO",dqu=0.7,start=c(0.1,0.1),
              PlotLikDo=TRUE,PlotLikTitle=c("NO2 | NO"), 
              PlotLikRange = list(a=c(0,0.8),b=c(-0.2,0.6)))

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ multivariate }

