\name{bootmex}
\alias{bootmex}
\alias{print.bootmex}
\alias{show.bootmex}
\alias{plot.bootmex}
\title{ Bootstrap a conditional multivariate extreme values model }
\description{
  Bootstrap a conditional multivariate extreme values model following
the method of Heffernan and Tawn, 2004. 
}
\usage{
bootmex(x, R = 100, nPass=3, trace=10)

\method{plot}{bootmex}(x, plots = "gpd", main = "", ...)
\method{print}{bootmex}(x, ...)
\method{show}{bootmex}(x, ...)
}
\arguments{
  \item{x}{ An object of class 'mex' as returned by function \code{\link{mex}}. }
  \item{R}{ The number of bootstrap runs to perform. Defaults to \code{R}=100. }
  \item{nPass}{An integer. Sometimes, particularly with small samples,
          the estimation process fails with
          some bootstrap samples. The function checks which runs
          fail and takes additional bootstrap samples in an attempt to get
          parameter estimates. By default, it has nPass=3 attempts at this
          before giving up.}
  \item{trace}{How often to inform the user of progress. Defaults to \code{trace=10}.}
  \item{plots}{ What type of diagnostic plots to produce.  Defaults to "gpd" in which case gpd parameter estimate plots are produced otherwise plots are made for the dependence parameters. }
  \item{main}{Title for plots.}
  \item{...}{Further arguments to be passed to methods.}
}
\details{
 Details of the bootstrap method are given by Heffernan and Tawn (2004). The procedure is semi-parametric.
 
 Firstly, values of all variables are simulated independently from the parametric Gumbel or Laplace distributions (depending on the choice of \code{margins} in the original call to \code{\link{mex}}). The sample size and data dimension match that of the original data set.  Then an empirical bootstrap sample is generated from the original data after its transformation to the Gumbel/Laplace scale. Again, sample size and structure match the original data set. The empirical bootstrap samples from each margin are then sorted, and then replaced by their corresponding values from the sorted Gumbel/Laplace samples. This procedure preserves the dependence structure of the empirical bootstrap sample while ensuring the marginal properties of the resulting semi-parametric bootstrap sample are those of the parametric Gumbel/Laplace distribution.

The simulated, ordered Laplace/Gumbel sample is then transformed to the scale of the original data by using the Probability Integral Transform. Values beneath the original thresholds for fitting of the GPD tail models are transformed by using the empirical distribution functions and for values above these thresholds, fitted Generalized Pareto Distribution is used.   This completes the sem-parametric bootstrap from the data.  

Parameter estimation is then carried out as follows: The parameters in the generalized Pareto
distributions are estimated by using the bootrap data, these data are then transformed to the Laplace/Gumbel scale using the
orginal threshold, their empirical distribution function and these estimated GPD parameters. The
variables in the dependence structure of these variables are then estimated.

Note that maximum likelihood estimation will often fail for small samples when the generalized
Pareto distribution is being fit. Therefore it will often be useful to use penalized likelihood
estimation. The function \code{bootmex} does whatever was done in the call to \code{migpd} or \code{mex} that generated the object with which it is being called.

Also note that sometimes (again, usually with small data sets) all of the simulated Laplace/Gumbel
random numbers will be beneath the threshold for the conditioning variable. Such samples are
abandoned by \code{bootmex} and a new sample is generated. This probably introduces some
bias into the resulting bootstrap distributions.

The \code{plot} method produces histograms of bootstrap gpd parameters (the default)
or scatterplots of dependence parameters with the point estimates for the
original data shown.

By design, there is no \code{coef} method. The bootstrapping is done to account
for uncertainty. It is not obvious that adjusting the parameters for the
mean bias is the correct thing to do.
}
\value{
 An object of class 'bootmex'. Print and plot functions are available.
}

\references{ J. E. Heffernan and J. A. Tawn, A conditional approach
	for multivariate extreme values, Journal of the Royal Statistical
	society B, 66, 497 -- 546, 2004}
\author{ Harry Southworth }
\seealso{ \code{\link{migpd}} , \code{\link{mexDependence} },
	\code{\link{bootmex}}, \code{\link{predict.mex}}.
		    }
\examples{
# Uncomment the following lines to run example - commented out to keep CRAN robots happy
#mymex <- mex(winter , mqu = .7, dqu = .7, which = "NO")
#myboot <- bootmex(mymex)
#myboot
#plot(myboot,plots="gpd")
#plot(myboot,plots="dependence")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ multivariate }% ..ONLY ONE.. keyword per line
