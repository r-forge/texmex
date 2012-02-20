\name{predict.gpd}
\alias{predict.gpd}
\alias{predict.link.gpd}
\alias{rl.gpd}
\alias{predict.link.bgpd}
\alias{rl.bgpd}
\alias{predict.link.bootgpd}
\alias{rl.bootgpd}

\title{Predict return levels from Generalized Pareto Distribution models, 
       or obtain the linear predictors.}
\description{Predict return levels from Generalized Pareto Distribution models, 
       or obtain the linear predictors.}

\usage{
predict.gpd(object, newdata = NULL, type = "return level", se.fit = FALSE, ci.fit = FALSE, M = 1000, alpha = 0.05, unique. = TRUE)
predict.link.gpd(object, newdata = NULL, type = "return level", se.fit = FALSE, ci.fit = FALSE, M = 1000, alpha = 0.05, unique. = TRUE, full.cov = FALSE)
rl.gpd(object, M = 1000, newdata = NULL, se.fit = FALSE, ci.fit = FALSE, alpha = 0.050, unique. = TRUE)
predict.bgpd(object, newdata = NULL, type = "return level", M = 1000, se.fit = FALSE, ci.fit = FALSE, alpha = 0.050, unique. = TRUE, all = FALSE, sumfun = NULL)
predict.link.bgpd(object, newdata = NULL, M = 1000, se.fit = FALSE, ci.fit = FALSE, alpha = 0.050, unique. = TRUE, all = FALSE, sumfun = NULL)
rl.bgpd(object, M = 1000, newdata = NULL, se.fit = FALSE, ci.fit = FALSE, alpha = 0.050, unique. = TRUE, all = FALSE, sumfun = NULL)
predict.bootgpd(object, newdata = NULL, type = "return level", M = 1000, se.fit = FALSE, ci.fit = FALSE, alpha = 0.050, unique. = TRUE, all = FALSE, sumfun = NULL)
predict.link.bootgpd(object, newdata = NULL, M = 1000, se.fit = FALSE, ci.fit = FALSE, alpha = 0.050, unique. = TRUE, all = FALSE, sumfun = NULL)
rl.bootgpd(object, M = 1000, newdata = NULL, se.fit = FALSE, ci.fit = FALSE, alpha = 0.050, unique. = TRUE, all = FALSE, sumfun = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An object of class \code{gpd}, \code{bgpd} or \code{bootgpd}.
}
  \item{newdata}{
The new data that you want to make the prediction for. Defaults in
\code{newdata = NULL} in which case the data used in fitting the model
will be used.}

  \item{type}{
The type of prediction, either ``return level'' (or ``rl'') or
``link'' (or ``lp''). Defaults to \code{type = "return level")}. When
a return level is wanted, the user can specify how high the level
is via the \code{M} argument. If \code{type = ``link''} the linear
predictor(s) for \code{phi} and/or \code{xi} are returned. 
}
  \item{se.fit}{
Whether or not to return the standard error of the predicted value.
Defaults to \code{se.fit = FALSE} and is not implemented for
\code{predict.bgpd} or \code{predict.bootgpd}.
}
  \item{ci.fit}{
Whether or not to return a confidence interval for the predicted
value. Defaults to \code{ci.fit = FALSE}. If set to \code{TRUE}
then the confidence interval is a simple symmetric confidence interval
based on the estimated approximate standard error when \code{predict.gpd}
is used. For the \code{bgpd} and \code{bootgpd} methods, the confidence
interval represents quantiles of the simulated distribution of the
parameters.
}
  \item{M}{
The return level. Defaults to \code{M = 1000}. If a vector is passed,
a 3-dimensional array is returned, with the third dimension representing
the values of the vector \code{M}.
}
  \item{alpha}{
If \code{ci.fit = TRUE}, a (1 - alpha/2)\% confidence interval is returned.
Defaults to \code{alpha = 0.050}.
}
  \item{unique.}{
If \code{unique = TRUE}, predictions for only the unique values of
the linear predictors are returned, rather than for every row of
\code{newdata}. Defaults to \code{unique. = TRUE}.
}
   \item{all}{
For the \code{bgpd} and \code{bootgpd}, if \code{all = TRUE}, the
predictions are returned for every simulated parameter vector. Otherwise,
only a summary of the posterior/bootstrap distribution is returned.
Defaults to \code{all = FALSE}.
}
    \item{full.cov}{
Should the full covariance matrix be returned as part of a \code{list}
object. This is used internally and not intended for direct use.
Defaults to \code{full.cov = FALSE}
    }
    \item{sumfun}{
For the \code{bgpd} and \code{bootgpd} methods, a summary function
can be passed in. If \code{sumfun = FALSE}, the default, the
summary function used returns the mean and median, and the quantiles
implied by \code{alpha}.
    }

}
\details{
By default, return levels predicted from the unique values of the
linear predictors are returned. For the \code{bootgpd} method,
confidence intervals are simple percentile intervals. The \code{bootgpd}
method is just a simple wrapper for the \code{bgpd} method.
}
\value{
An array, with 2 dimensions for a single value of M, or 3 dimensions
for multiple values of M.}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Harry Southworth and Janet E. Heffernan
}
\note{
At present, the confidence intervals returned for an object of class
\code{gpd} are simple confidence intervals based on assumptions
of normality that are likely to be far from the truth in many cases.
A better approach would be to use profile likelihood, and we intend
to implement this method at a future date.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
}
\keyword{ methods }
