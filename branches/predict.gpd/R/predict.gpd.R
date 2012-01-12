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

## TODO: Need to add rownames (or something) resulting from using unique() on
##       the design matrices.

################################################################################
## gpd

predict.gpd <-
    # Get predictions for a gpd object. These can either be the linear predictors
    # or return levels.
function(object, newdata=NULL, type="return level", se.fit=FALSE,
         ci.fit=FALSE, M=1000, alpha=.050, unique.=TRUE){
    theCall <- match.call()
    
#    type <- match.arg(type)
    
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
        X.xi <- if (!is.null(xi.fo)){ model.matrix(as.formula(xi.fo), newdata) }
                else { matrix(1, nrow(newdata)) }
        X.phi <- if (!is.null(phi.fo)){ model.matrix(as.formula(object$call$phi), newdata) }
                 else { matrix(1, nrow(newdata)) }
    }

    else {
        X.xi <- object$X.xi
        X.phi <- object$X.phi
    }

    if (unique.){
        u <- (1 - duplicated(X.phi)) + (1 - duplicated(X.xi)) > 0
        X.xi <- cbind(X.xi[u, ])
        X.phi <- cbind(X.phi[u, ])
    }

    phi <- c(object$coefficients[1:ncol(X.phi)] %*% t(X.phi))
    xi <- c(object$coefficients[(ncol(X.phi) + 1):length(object$coefficients)] %*% t(X.xi))

    res <- cbind(phi, xi)

    if (ci.fit){
        phi.cov <- as.matrix(object$cov[1:ncol(X.phi), 1:ncol(X.phi)])
        xi.cov <- as.matrix(object$cov[(ncol(X.phi) + 1):length(object$coefficients), (ncol(X.phi) + 1):length(object$se)])
        phi.se <- sqrt(rowSums((X.phi %*% phi.cov) * object$X.phi))
        xi.se <- sqrt(rowSums((X.xi %*% xi.cov) * X.xi))

        z <- qnorm(1 - alpha/2)

        phi.lo <- phi - phi.se*z
        phi.hi <- phi + phi.se*z
        xi.lo <- xi - xi.se*z
        xi.hi <- xi + xi.se*z

        res <- cbind(res, phi.lo, phi.hi, xi.lo, xi.hi)
    } # Close if(ci.fit

    if (se.fit){
        if (!ci.fit){ # Because if ci.fit, phi.se and xi.se already exist
            phi.cov <- as.matrix(object$cov[1:ncol(X.phi), 1:ncol(X.phi)])
            xi.cov <- as.matrix(object$cov[(ncol(X.phi) + 1):length(object$coefficients), (ncol(X.phi) + 1):length(object$se)])
            phi.se <- sqrt(rowSums((X.phi %*% phi.cov) * X.phi))
            xi.se <- sqrt(rowSums((X.xi %*% xi.cov) * X.xi))
        }
        res <- cbind(res, phi.se, xi.se)
    } # Close if(se.fit
    
    if (full.cov){ # Covariance of (phi, xi) for each unique (phi, xi) pair
        phi.cov <- as.matrix(object$cov[1:ncol(X.phi), 1:ncol(X.phi)])
        xi.cov <- as.matrix(object$cov[(ncol(X.phi) + 1):length(object$coefficients), (ncol(X.phi) + 1):length(object$se)])

        covar <- rep(0, nrow(X.xi))
        for(k in 1:length(covar)){
            for (i in 1:ncol(X.phi)){
                for (j in 1:ncol(X.xi)){
                    covar[k] <- covar[k] + X.phi[k, i] * X.xi[k, j] * object$cov[i, j]
                } # Close for j
            } # Close for i
        } # Close for k
        phi.var <- rowSums((X.phi %*% phi.cov) * X.phi)
        xi.var <- rowSums((X.xi %*% xi.cov) * X.xi)
        
        res <- list(link=res, cov=list(phi.var=phi.var, xi.var=xi.var, covariances=covar))
    }
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
        # the CI is approxima#te anyway.
        
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
    covs <- co[[2]] # List of covariance matrices
    co <- co[[1]]
    
    res <- object$threshold + exp(co[,1]) / co[,2] * ((M * object$rate)^co[,2] -1)
    res <- cbind(RL=res)

    getse <- function(o, co, M){
#        dxm <- t(apply(co, 1, gpd.delta, m=M))
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

    if (ci.fit){
        co <- cbind(rep(object$rate, nrow(co)), co)
        se <- getse(object, co, M)
        lo <- res - qnorm(1 - alpha/2)*se
        hi <- res + qnorm(1 - alpha/2)*se

        res <- cbind(res, lo=lo, hi=hi)
        colnames(res) <- c("RL", paste(100*alpha/2, "%", sep = ""),
                                 paste(100*(1 - alpha/2), "%", sep = ""))
        rownames(res) <- NULL
    } # Close if (ci.fit

    if (se.fit){
        if (!ci.fit){
            co <- cbind(rep(object$rate, nrow(co)), co)
            se <- getse(object, co, M)
        }
        cnr <- colnames(res)
        res <- cbind(res, se=se)
    }

    res
}

################################################################################
## bgpd

predict.bgpd <- function(object, newdata=NULL, type="return level", M=1000, alpha=.050, unique.=TRUE){
    theCall <- match.call()
        
    res <- switch(type,
                  "rl" = , "return level" = rl.bgpd(object, M),
                  "lp" = , "link" = predict.link.bgpd(object, newdata)
                  )
    res
}

predict.link.bgpd <- function(object, newdata, se.fit, ci.fit, apha=.050, unique.=TRUE){
    if (!is.null(newdata)){
        xi.fo <- object$call$xi
        phi.fo <- object$call$phi
        X.xi <- if (!is.null(xi.fo)){ model.matrix(as.formula(xi.fo), newdata) }
                else { matrix(1, nrow(newdata)) }
        X.phi <- if (!is.null(phi.fo)){ model.matrix(as.formula(object$call$phi), newdata) }
                 else { matrix(1, nrow(newdata)) }
    }

    else {
        X.xi <- object$X.xi
        X.phi <- object$X.phi
    }

    if (unique.){
        u <- (1 - duplicated(X.phi)) + (1 - duplicated(X.xi)) > 0
        X.xi <- cbind(X.xi[u, ])
        X.phi <- cbind(X.phi[u, ])
    }

    phi <- cbind(object$param[, 1:ncol(X.phi)])
    xi <- cbind(object$param[, (ncol(X.phi) + 1):(ncol(X.phi) + ncol(X.xi))])
    
    res <- lapply(1:nrow(X.xi),
                  function(i, phi, xi, Xphi, Xxi){
#browser()
                      phi <- rowSums(t(t(phi) * c(Xphi[i, ])))
                      xi <- rowSums(t(t(xi) * c(Xxi[i, ])))
                      cbind(phi=phi, xi=xi)
                    }, phi=phi, xi=xi, Xphi=X.phi, Xxi=X.xi)

    res <- t(sapply(res, function(x){ apply(x, 2, mean) }))
    res
}

rl.bgpd <- function(object, M){


}

################################################################################
## bootgpd

predict.bootgpd <- function(object, type=c("return level", "link"), M=1000){
    theCall <- match.call()

    type <- match.arg(type)

    res <- switch(type,
                  "return level" = rl.bgpd(object, M),
                  "parameters" = cobgpd(object, newdata)
                  )
    res <- list(rl = res, call = theCall)
    oldClass(res) <- "returnLevel"
    res

}
predict.link.bootgpd <- function(object, newdata, se.fit, ci.fit){


}

rl.bootgpd <- function(object, M){


}



