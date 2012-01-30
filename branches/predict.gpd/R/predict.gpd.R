# Author: Harry Southworth
# Date: 2011-11-25
## Purpose: Create a predict method for objects of class gpd and bgpd that
##          returns parameters, return levels or (maybe) return periods, 
##          depending on arguments given.
#
# predict.gpd
# predict.bgpd
# predict.bootgpd
# predict.link.gpd
# rl
# rl.gpd
# predict.link.bgpd
# rl.bgpd
# predict.link.bootgpd
# rl.bootgpd

################################################################################
## gpd

predict.gpd <-
    # Get predictions for a gpd object. These can either be the linear predictors
    # or return levels.
function(object, newdata=NULL, type="return level", se.fit=FALSE,
         ci.fit=FALSE, M=1000, alpha=.050, unique.=TRUE){
    theCall <- match.call()
    
    res <- switch(type,
                  "rl"=, "return level" = rl.gpd(object, M, newdata,
                                                 se.fit=se.fit, ci.fit=ci.fit,
                                                 alpha=alpha, unique.=unique.),
                  "lp" =,"link" = predict.link.gpd(object, newdata, se.fit,
                                                   ci.fit, alpha, unique.=unique.)
                  )
    res
}

## Linear predictor functions for GPD

predict.link.gpd <- function(object, newdata=NULL, se.fit=FALSE, ci.fit=FALSE,
                             alpha=.050, unique.=TRUE, full.cov=FALSE){

    if (!is.null(newdata)){
        xi.fo <- object$call$xi
        phi.fo <- object$call$phi
        X.xi <-  if (!is.null(xi.fo)) { model.matrix(as.formula(xi.fo),  newdata) } else { matrix(1, nrow(newdata)) }
        X.phi <- if (!is.null(phi.fo)){ model.matrix(as.formula(phi.fo), newdata) } else { matrix(1, nrow(newdata)) }
    } else {
        X.xi <- object$X.xi
        X.phi <- object$X.phi
    }

    if (unique.){
        u <- (1 - duplicated(X.phi)) + (1 - duplicated(X.xi)) > 0
        X.xi  <- if(is.matrix(X.xi[u,]))  X.xi[u, ]  else if(ncol(X.xi) == 1)  cbind(X.xi[u,])  else t(cbind(X.xi[u,]))
        X.phi <- if(is.matrix(X.phi[u,])) X.phi[u, ] else if(ncol(X.phi) == 1) cbind(X.phi[u,]) else t(cbind(X.phi[u,]))
    }

    whichPhi <- 1:ncol(X.phi)
    whichXi  <- (ncol(X.phi) + 1):length(object$coefficients)
    phi <- c(object$coefficients[whichPhi] %*% t(X.phi))
    xi  <- c(object$coefficients[whichXi]  %*% t(X.xi))

    res <- cbind(phi, xi)

    if(ci.fit | se.fit | full.cov){
       phi.cov <- as.matrix(object$cov[whichPhi, whichPhi])
       xi.cov  <- as.matrix(object$cov[whichXi,  whichXi])
       
       if(ci.fit | se.fit){
          phi.se <- sqrt(rowSums((X.phi %*% phi.cov) * X.phi))
          xi.se <-  sqrt(rowSums((X.xi  %*% xi.cov)  * X.xi))
       }
    }
    
    if (ci.fit){
        z <- qnorm(1 - alpha/2)

        phi.lo <- phi - phi.se*z
        phi.hi <- phi + phi.se*z
        xi.lo <- xi - xi.se*z
        xi.hi <- xi + xi.se*z

        res <- cbind(res, phi.lo, phi.hi, xi.lo, xi.hi)
    } # Close if(ci.fit

    if (se.fit){
        res <- cbind(res, phi.se, xi.se)
    } # Close if(se.fit

    if(dim(X.phi)[2] > 1){
      res <- cbind(res,X.phi)
    }
    if(dim(X.xi)[2] > 1){
      res <- cbind(res,X.xi)
    }

    if (full.cov){ # Covariance of (phi, xi) for each unique (phi, xi) pair
        covar <- rep(0, nrow(X.xi))
        for(k in 1:length(covar)){
            for (i in 1:ncol(X.phi)){
                for (j in 1:ncol(X.xi)){
                    covar[k] <- covar[k] + X.phi[k, i] * X.xi[k, j] * object$cov[i, ncol(X.phi) + j] 
                } # Close for j
            } # Close for i
        } # Close for k
        phi.var <- rowSums((X.phi %*% phi.cov) * X.phi)
        xi.var <- rowSums((X.xi %*% xi.cov) * X.xi)
        
        res <- list(link=res, cov=list(phi.var=phi.var, xi.var=xi.var, covariances=covar))
    }
    oldClass(res) <- "predict.link.gpd"
    res
}

## Return level functions for GPD

## Reversing arguments M and newdata for anyone who wants to call these functions
## directly

## Will want to get return levels when using GEV rather than GPD, so make
## rl generic

rl <- function(object, M, newdata, se.fit=FALSE, ci.fit=FALSE, alpha=.050, unique.=TRUE, ...){
    UseMethod("rl")
}

gpd.delta <- function(a, m){
   # This is not exact if a prior (penalty) function is used, but
   # the CI is approximate anyway.
        
    out <- matrix(0, nrow=3, ncol=length(m))
        
    if (a[3] == 0){ # exponential case
        out[1,] <- exp(a[2]) / a[1]
        out[2,] <- exp(a[2]) * log(m * a[1])
    } else {
        out[1,] <- exp(a[2]) * m^a[3] * a[1]^(a[3] - 1)
        out[2,] <- exp(a[2]) / a[3] * ((m*a[1])^a[3] - 1) 
        out[3,] <- -exp(a[2]) / (a[3]*a[3]) * ( (m * a[1] )^a[3] - 1 ) +
                   exp(a[2]) / a[3] * (m * a[1])^a[3] * log(m * a[1])
    } 
    
   out
} 

rl.gpd <- function(object, M=1000, newdata=NULL, se.fit=FALSE, ci.fit=FALSE,
                   alpha=.050, unique.=TRUE){
    co <- predict.link.gpd(object, newdata=newdata, unique.=unique., full.cov=TRUE)
    covs <- co[[2]] # List list(phi.var=phi.var, xi.var=xi.var, covariances=covar)
    co <- co[[1]]
 
    gpdrl <- function(u, theta, phi, xi, m){
        res <- u + exp(phi) / xi *((m * theta)^xi -1) 
        cbind(RL=res)
    }
 
    res <- lapply(M, gpdrl,
                  u=object$threshold, theta=object$rate, phi=co[,1], xi=co[,2])

    getse <- function(o, co, M){
        dxm <- lapply(split(co, 1:nrow(co)), gpd.delta, m=M)

        V <- lapply(1:nrow(co),
                    function(i, x, rate, n){
                        # Construct covariance matrix
                        cov <- matrix(c(x[[1]][i], rep(x[[3]][i], 2), x[[2]][i]), ncol=2)
                        matrix(c(rate * (1 - rate) / n, 0, 0,
                               0, cov[1,], 
                               0, cov[2,]), ncol=3)
                    }, rate = o$rate, n = length(o$y) / o$rate, x=covs)

        # Get (4.15) of Coles, page 82, adjusted for phi = log(sigma)
        se <- sapply(1:length(V),
                     function(i, dxm, V){
                        V <- V[[i]]; dxm <- c(dxm[[i]])
                        sqrt(mahalanobis(dxm, center=c(0, 0, 0), cov=V, inverted=TRUE))
                     }, dxm=dxm, V=V)
        se
    }

    if (ci.fit){ # need to update plotrl.gpd too once profile lik confidence intervals implemented here
        co <- cbind(rep(object$rate, nrow(co)), co)
        ci.fun <- function(i, object, co, M, res, alpha){
            wh <- res[[i]];
            se <- getse(object, co, M[i])
            lo <- wh - qnorm(1 - alpha/2)*se
            hi <- wh + qnorm(1 - alpha/2)*se
            wh <- cbind(wh, lo=lo, hi=hi)

            colnames(wh) <- c("RL", paste(100*alpha/2, "%", sep = ""),
                              paste(100*(1 - alpha/2), "%", sep = ""))
            #rownames(wh) <- NULL
            wh
        } # ci.fun
        res <- lapply(1:length(M), ci.fun, object=object, co=co, M=M, res=res, alpha=alpha)
    } # Close if (ci.fit

    if (se.fit){
        if (!ci.fit) { co <- cbind(rep(object$rate, nrow(co)), co) }
        se.fun <- function(i, object, co, M, res, alpha){
            wh <- res[[i]];
            se <- getse(object, co, M[i])
            wh <- cbind(RL=wh, se.fit=se)
            #rownames(wh) <- NULL
            wh
        } # ci.fun
        res <- lapply(1:length(M), se.fun, object=object, co=co, M=M, res=res, alpha=alpha)
    }
    
    names(res) <- paste("M.", M, sep = "")
    oldClass(res) <- "rl.gpd"
    res
}

################################################################################
## bgpd

predict.bgpd <- function(object, newdata=NULL, type="return level", M=1000,
                         se.fit=FALSE, ci.fit=FALSE, alpha=.050, unique.=TRUE,
                         all=FALSE, sumfun=NULL){
    theCall <- match.call()
        
    res <- switch(type,
                  "rl" = , "return level" = rl.bgpd(object, M=M, newdata=newdata,
                                                    se.fit=se.fit, ci.fit=ci.fit,
                                                    alpha=alpha, unique.=unique., all=all,
                                                    sumfun=sumfun),
                  "lp" = , "link" = predict.link.bgpd(object, newdata=newdata,
                                                      se.fit=se.fit, ci.fit=ci.fit,
                                                      alpha=alpha, unique.=unique., all=all,
                                                      sumfun=sumfun)
                  )
    res
}

predict.link.bgpd <- function(object, newdata=NULL, se.fit=FALSE, ci.fit=FALSE,
                              alpha=.050, unique.=TRUE, all=FALSE, sumfun=NULL){
    if (!is.null(newdata)){
        xi.fo <- object$map$call$xi
        phi.fo <- object$map$call$phi
        X.xi <-  if (!is.null(xi.fo))  model.matrix(as.formula(xi.fo),  newdata) else matrix(1, nrow(newdata)) 
        X.phi <- if (!is.null(phi.fo)) model.matrix(as.formula(phi.fo), newdata) else matrix(1, nrow(newdata))
    } else {
        X.xi <- object$X.xi
        X.phi <- object$X.phi
    }

    if (unique.){
        u <- (1 - duplicated(X.phi)) + (1 - duplicated(X.xi)) > 0
        X.xi  <- if(is.matrix(X.xi[u,]))  X.xi[u, ]  else if(ncol(X.xi) == 1)  cbind(X.xi[u,])  else t(cbind(X.xi[u,]))
        X.phi <- if(is.matrix(X.phi[u,])) X.phi[u, ] else if(ncol(X.phi) == 1) cbind(X.phi[u,]) else t(cbind(X.phi[u,]))
    }

    phi <- cbind(object$param[, 1:ncol(X.phi)])
    xi <- cbind(object$param[, (ncol(X.phi) + 1):(ncol(X.phi) + ncol(X.xi))])

    # Get point estimates (means)
    res <- lapply(1:nrow(X.xi),
                  function(i, phi, xi, Xphi, Xxi){
                      phi <- rowSums(t(t(phi) * c(Xphi[i, ])))
                      xi <- rowSums(t(t(xi) * c(Xxi[i, ])))
                      cbind(phi=phi, xi=xi)
                    }, phi=phi, xi=xi, Xphi=X.phi, Xxi=X.xi)

    ############################################################################
    ## Hard part should be done now. Just need to summarize

    if (ci.fit){
        if (is.null(sumfun)){
            sumfun <- function(x){
                c(quantile(x, prob=c(alpha/2, .50, 1 - alpha/2)), mean(x))
            }
            neednames <- TRUE
        } else { 
            neednames <- FALSE 
        }

        res <- t(sapply(res, function(x, fun){ apply(x, 2, sumfun) }, fun=sumfun))
        
        if (neednames){
            nms <- c(paste(100*alpha/2, "%", sep = ""),
                    "50%", paste(100*(1-alpha/2), "%", sep = ""),
                    "Mean")

            colnames(res) <- c(paste("phi:", nms), paste("xi:", nms))
        }
    }
    else if (se.fit){ warning("se.fit not implemented - ignoring") }
    else if (all){ res <- res }
    else { # Just point estimates
        res <- t(sapply(res, function(x){ apply(x, 2, mean) }))
    }
    oldClass(res) <- "predict.link.bgpd"
    res
}

rl.bgpd <- function(object, M, newdata=NULL, unique.=unique., se.fit=FALSE,
                    ci.fit=FALSE, all=FALSE, alpha=.050, sumfun=NULL){
    co <- predict.link.bgpd(object, newdata=newdata, unique.=unique., all=TRUE, sumfun=NULL)

    bgpdrl <- function(o, u, theta, m){
        res <- u + exp(o[, 1]) / o[, 2] *((m * theta)^o[, 2] -1) 
        cbind(RL=res)
    }
    
    # co is (probably) a list with one element for each unique item in
    # new data. Need to loop over vector M and the elements of co
    
    getrl <- function(m, co, u, theta, ci.fit, alpha, all){
        res <- sapply(co, bgpdrl, u=u, theta=theta, m=m)

        if (ci.fit){
            if (is.null(sumfun)){
                sumfun <- function(x){
                    c(quantile(x, prob=c(alpha/2, .50, 1 - alpha/2)), mean(x))
                }
                neednames <- TRUE
            } else { 
               neednames <- FALSE 
            }
            res <- t(apply(res, 2, sumfun))
            if (neednames){
                colnames(res) <- c(paste(100*alpha/2, "%", sep = ""),
                                   "50%", paste(100*(1-alpha/2), "%", sep = ""),
                                   "Mean")
            }

        } # Close if (ci.fit
        else if (se.fit){
            warning("se.fit not implemented")
            res <- apply(res, 2, mean)
        }
        else if (!all){ res <- apply(res, 2, mean) }
        res
    }
    
    res <- lapply(M, getrl, co=co, u=object$threshold, theta=object$map$rate, ci.fit=ci.fit, alpha=alpha, all=all)
    names(res) <- paste("M.", M, sep = "")
    oldClass(res) <- "rl.bgpd"
    res
}

################################################################################
## bootgpd

predict.bootgpd <- function(object, newdata=NULL, type="return level",
                            unique.=TRUE, ci.fit=FALSE, se.fit=FALSE, M=1000,
                            all=FALSE, alpha=.050){
    theCall <- match.call()

    res <- switch(type,
                  "rl" = , "return level" = rl.bootgpd(object, M, newdata=newdata,
                                                       unique.=TRUE,
                                                       se.fit=se.fit, ci.fit=ci.fit,
                                                       all=all, alpha=alpha),
                  "lp" = , "link" = predict.link.bootgpd(object, newdata=newdata,
                                                         unique.=TRUE,
                                                         se.fit=se.fit, ci.fit=ci.fit,
                                                         all=all, alpha=alpha)
                  )
    res

}

namesBoot2bgpd <- function(object){
    names(object) <- c("call", "param", "original", "map")
    object$X.phi <- object$map$X.phi
    object$X.xi <- object$map$X.xi
    object$threshold <- object$map$threshold
    object
}

predict.link.bootgpd <- function(object, newdata=NULL, se.fit=FALSE, ci.fit=FALSE,
                                 unique.=TRUE, all=FALSE, alpha=.050){
    # This should just be the same as for a bgpd object, but some
    # names and stuff are different.
  object <- namesBoot2bgpd(object)
  res <- predict.link.bgpd(object, newdata=newdata, se.fit=se.fit, ci.fit=ci.fit)
  oldClass(res) <- "predict.link.bootgpd"
  res
}

rl.bootgpd <- function(object, M, newdata=NULL, se.fit=FALSE, ci.fit=FALSE, all=FALSE,
                       unique.=TRUE, alpha=alpha){
    # This should just be the same as for a bgpd object, but some
    # names and stuff are different.
  object <- namesBoot2bgpd(object)
  res <- rl.bgpd(object, M=M, newdata=newdata, se.fit=se.fit, ci.fit=ci.fit, all=all, unique.=unique.,alpha=alpha)
  oldClass(res) <- "rl.bootgpd"
  res
}



