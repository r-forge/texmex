\name{extremalIndex}
\alias{extremalIndex}
\alias{extremalIndexRangeFit}
\alias{declust}
\alias{declust.default}
\alias{declust.extremalIndex}
\alias{plot.declustered}
\alias{bootExtremalIndex}
\alias{gpd.declustered}
\title{
Extremal index estimation and automatic declustering
}
\description{
Given a threshold which defines excesses above that threshold, estimate the extremal index of a dependent sequence by using the method of Ferro and Segers, 2003.  The extremal index estimate can then be used to carry out automatic declustering of the sequence to identify independent clusters and estimate the GPD for cluster maxima.  Graphical diagnostics of model fit are available.
}
\usage{
extremalIndex(y, data = NULL, threshold)

extremalIndexRangeFit(y, data = NULL, umin = quantile(y,.5), 
   umax = quantile(y, 0.95), nint = 10, nboot = 100, alpha = .05, 
   xlab = "Threshold", addNexcesses = TRUE, estGPD=TRUE, verbose = TRUE, 
   trace = 10, ...)
   
bootExtremalIndex(x)

declust(y, r=NULL, data = NULL, ...)

\method{declust}{extremalIndex}(y, r=NULL,...)

\method{plot}{declustered}(x, ylab = "Data",...)

\method{gpd}{declustered}(y, ...)
}
\arguments{
  \item{y}{ Argument to function \code{extremalIndex}: either a numeric vector or the name of a variable in \code{data}.}
  \item{data}{ A data frame containing \code{y} and any covariates. }
  \item{threshold}{ The threshold for \code{y}, exceedances above which will be used to   estimate the extremal index and carry out automatic declustering.}
  \item{x}{ Objects passed to methods. }
  \item{r}{ Positivie integer: run length to be used under "runs" declustering. If specified then so-called "runs" declustering will be carried out, otherwise defaults to NULL in which case the automatic "intervals" declustering method of Ferro and Segers is used.}
  \item{umin}{ The minimum threshold above which to esimate the parameters. }
  \item{umax}{ The maximum threshold above which to esimate the parameters. }
  \item{nint}{ The number of thresholds at which to perform the estimation. }
  \item{nboot}{ Number of bootstrap samples to simulate at each threshold for estimation.}
  \item{alpha}{ (1 - alpha)\% confidence intervals will be plotted with the 
               point estimates. Defaults to \code{alpha = 0.05}.}
  \item{xlab}{ Label for the x-axis. }
  \item{ylab}{ Label for the y-axis. }
  \item{addNexcesses}{ Whether to annotate the top axis of plots with the number of excesses above the corresponding threhsold. Defaults to \code{TRUE}. }
  \item{estGPD}{ Whether to estimate GPD parameters at each choice of thereshold -- defaults to \code{TRUE} in which case the GPD parameters are estimated.}
  \item{verbose}{Whether to report on progress in RangeFit calculations. Defaults to \code{TRUE}.}
  \item{trace}{How frequently to report bootstrap progress in EangeFit calculations.  Defaults to 10.} 
  \item{...}{Further arguments to methods.}
}
\details{
The function \code{extremalIndex} estimates the extremal index of a dependent series of observations above a given threshold \code{threshold}, returning an object of class "extremalIndex".  Plot, print and show methods are available for this class. A graphical diagnostic akin to Figure 1 in Ferror and Segers (2003) is produced by the \code{plot} method for this class.  This plot is used to test the model assumption underpinning the estimation, with good fit being indicated by interexceedance times which correspond to inter-cluster times lying close to the diagonal line indicated.

In addition to good model fit, an appropriate choice of threshold is one above which the estimated extremal index is stable over further, higher thresholds (up to estimation uncertainty).  This can be assessed by using the function \code{extremalIndexRangeFit}, which examines a range of threshold values.  At each threshold, the extremal index is estimated; that estimate is used to decluster the series and the parameters of the GPD are optionally estimated for the resulting declustered series.  Uncertainty in the estimation of the extremal index and GPD parameters is assessed by using a bootstrap scheme which accounts for uncertainty in the extremal index estimation, and the corresponding uncertainty in the declustering of the series.

The function \code{declust} returns an object of class "declustered", identifying independent clusters in the original series. Print, plot and show methodsa are available for this class. The GPD model can be fitted to objects of this class, including the use of covariates in the linear predictors for the parameters of the GPD.  See examples below.
}
\value{
The function \code{extremalIndex} returns a list of class "extremalIndex":
\item{EIintervals}{Estimate of the extremal index by using the intervals estimator of Ferro and Segers.}
\item{threshold}{threshold for declustering and estimation}
\item{TotalN}{length of original data series}
\item{nExceed}{number of exceedances of \code{threshold} in original series.}
\item{thExceedanceProb}{probablity of threshold exceedance in original series.}
\item{call}{the original function call }
\item{interExceedTimes}{times between threshold exceedances}
\item{thExceedances}{observation from the original series which are above \code{threshold}}
\item{exceedanceTimes}{times of occurrance of threshold exceedances}
\item{y}{original dependent series}
\item{data}{data frame or NULL}

The function \code{declust} returns a list of type "declustered":

\item{clusters}{integer labels assigning threshold exceedances to clusters}
\item{sizes}{number of exceedances in each cluster}
\item{clusterMaxima}{vector made up of the largest observation from each distinct cluster.  In the case of ties, the first value is taken.}
\item{isClusterMax}{logical; length equal to number of threshold exceedances, value is \code{TRUE} for threshold exceedances which correspond to cluster maxima}
\item{y}{see entry for object of class "extremalIndex" above}
\item{data}{see entry for object of class "extremalIndex" above}
\item{threshold}{see entry for object of class "extremalIndex" above}
\item{EIintervals}{see entry for object of class "extremalIndex" above}
\item{call}{see entry for object of class "extremalIndex" above}
\item{InterExceedTimes}{times between threshold exceedances, length is one less than the number of threshold exceedances}
\item{InterCluster}{logical: indicates inter exceedance times larger than \code{r} the run length used for declustering}
\item{thExceedances}{see entry for object of class "extremalIndex" above}
\item{exceedanceTimes}{see entry for object of class "extremalIndex" above}
\item{r}{run length used for declustering}
\item{nClusters}{Number of indenendent clusters identified}
\item{method}{Method used for declustering (either "intervals" or "runs")}

The function \code{bootExtremalIndex} return a single vector corersponding to a bootstrap sample from the original series: observations are censored at \code{threshold} so that values below this threshold are indicated by the value -1.

The method \code{gpd} for class "declustered" returns an object of type "gpd" or "bgpd" depending on the precise function call - see documentation for \code{\link{gpd}}.
}
\references{
Ferro, C.A.T. and Segers, J., (2003) "Inference for clusters of Extreme Values", JRSS B 65, Part 2, pp 545--556.
}
\author{
Janet E. Heffernan
}
\seealso{
\code{\link{gpd}}
}
\examples{
par(mfrow=c(2,2)); 
extremalIndexRangeFit(summer$O3,nboot=10)
ei <- extremalIndex(summer$O3,threshold=45)
plot(ei)
d <- declust(ei)
plot(d)
gpd(d)

## fitting with covariates:

so2 <- extremalIndex(SO2,data=winter,threshold=15)
plot(so2)
so2 <- extremalIndex(SO2,data=winter,threshold=20)
plot(so2) ## fits better

so2.d <- declust(so2)
par(mfrow=c(1,1)); plot(so2.d)
so2.d.gpd <- gpd(so2.d) # AIC 661.1

gpd(so2.d,phi=~NO)
gpd(so2.d,phi=~NO2)
gpd(so2.d,phi=~O3) # better AIC 651.9
gpd(so2.d,phi=~PM10)

so2.d.gpd.o3 <- gpd(so2.d,phi=~O3)

par(mfrow=c(2,2)); plot(so2.d.gpd.o3)
}

