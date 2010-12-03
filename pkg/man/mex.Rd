\name{mex}
\alias{mex}
\alias{print.mex}
\alias{summary.mex}
\alias{coef.mex}
\alias{plot.mex}
\alias{predict.mex}
\alias{print.predict.mex}
\alias{summary.predict.mex}
\alias{plot.predict.mex}
\title{
Conditional multivariate extreme values modelling
}
\description{
Fit the conditional multivariate extreme value model of Heffernan and Tawn
}
\usage{
mex(data, which, mth, mqu, gth, gqu, penalty = "gaussian", maxit = 10000, trace = 0, verbose = FALSE, priorParameters = NULL)
\method{print}{mex}(x, ...)
\method{summary}{mex}(object, ...)
\method{plot}{mex}(x, quantiles = seq(0.1, by = 0.2, len = 5), col = "grey", ...)
\method{coef}{mex}(object, ...)
\method{predict}
\method{print}{predict.mex}
\method{summary}{predict.mex}
\method{plot}{predict.mex}
}

\arguments{
  \item{data}{A numeric matrix or data.frame, the columns of which are to be modelled.}
  \item{which}{The variable on which to condition.}
  \item{mth}{In \code{mex}, the threshold above which to fit generalized Pareto distributions.
            If this is a vector of length 1, the same threshold will be used for
            each variable. Otherwise, it should be a vector whose length is
            equal to the number of columns in \code{data}.

            In \code{summary.predict.mex}, the thresholds over which to simulate
            data from the fitted multivariate model. If not supplied, it is
            taken to be the thresholds that were used to fit the model on the
            scale of the original data.}
  \item{mqu}{As an alternative to specifying the thresholds via th, you can
            specify the quantile above which to fit generalized Pareto
            distributions. If this is a vector of length 1, the same quantile
            will be used for each variable. Otherwise, it should be a vector
            whose length is equal to the number of columns in data.}
  \item{dth}{The threshold to use on the transformed scale; that is, the
            threshold to use when estimating the dependence structure. The same
            threshold will be used for each variable.}
  \item{dqu}{As an alternative to \code{dth}, specify the quantile at which to
            threshold the transformed data when estimating the dependence
            parameters.}
  \item{penalty}{How to penalize the likelihood when estimating the marginal
            generalized Pareto distributions. Defaults to ``gaussian''. See the
            help file for \code{\link{gpd}} for more information.}
  \item{maxit}{The maximum number of iterations to be used by the optimizer.
            defaults to \code{maxit = 1000}.}
  \item{trace}{Whether or not to inform the user of the progress of the
            optimizer. Defaults to 0, indicating no trace.}
  \item{verbose}{Whether or not to keep the user informed of progress. Defaults
            to \code{verbose = FALSE}.}
  \item{priorParameters}{A named list, each element of which contains two
            components: the first should be a vector of length 2 corresponding
            to the location of the Gaussian distribution; the second should be 2
            by 2 matrix corresponding to the covariance of the distribution.
            This is only used if \code{penalty = 'gaussian'}. The names should match
            the names of the columns of \code{data}. If not provided, it defaults to
            being centred at zero, with variance 10000 for log(sigma) and 0.25
            for xi. See the details section.}
  \item{quantiles}{A vector of quantiles taking values between 0 and 1 specifying 
  the quantiles of the conditional distributions which will be plotted.}
  \item{col}{The color for points on scatterplots of residuals and original data.}
  \item{pqu}{Argument to \code{predict.mex}. The quantile of the conditioning variable over which it will be
thresholded. Defaults to \code{pqu = .99 }. }
  \item{nsim}{Argument to \code{predict.mex}. The number of simulated observations to be generated for each
set of bootstrap parameter estimates. }


}
\details{The parameters in the generalized Pareto distribution are estimated for
         each column of the data in turn, independently of all other columns.
         Then, the conditional multivariate approach of Heffernan and Tawn is
         used to model the dependence between variables.}

         The \code{predict} method for \code{mex} works as follows.
         For each set of parameter estimates and normalized conditional Gumbel
         data, simulated values of the dependent variables are created, given that the
         conditioning variable is above its 100\code{pqu} quantile.

         The returned object has class 'predict.mex'. There are \code{print}, \code{summary} and
         \code{plot} functions available. The simulated values based on the actual data
         and point estimates appear in component \code{object$data$simulated}. The
         simulated data from the bootstrap estimates appears in \code{object$replicates}.


\value{
  \item{margins}{An object of class \code{\link{migpd}}.}
  \item{dependence}{An object of class \code{\link{mexDependence}}.}
}
\references{J. E. Heffernan and J. A. Tawn, A conditional approach for
            multivariate extreme values, Journal of the Royal Statistical
            Society B, 66, 497 – 546, 2004}
\author{
Harry Southworth, Janet E. Heffernan
}

\seealso{\code{\link{migpd}}, \code{\link{mexDependence}}, \code{\link{bootmex}}}
\examples{
w <- mex(winter, mqu=.7)
w
par(mfcol=c(4, 2))
plot(w)
}
\keyword{ models }
\keyword{ multivariate }