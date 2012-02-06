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
        u <- !duplicated(cbind(X.phi,X.xi))
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

    res <- addCov(res,X.phi)
    res <- addCov(res,X.xi)
    
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
    covs <- co[[2]] # list(phi.var=phi.var, xi.var=xi.var, covariances=covar)
    co <- co[[1]]
    X <- co[,-(1:2)]
 
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
    
    co <- cbind(rep(object$rate, nrow(co)), co)
    
    if (ci.fit){ # need to update plotrl.gpd too once profile lik confidence intervals implemented here
        ci.fun <- function(i, object, co, M, res, alpha){
            wh <- res[[i]];
            se <- getse(object, co, M[i])
            lo <- wh - qnorm(1 - alpha/2)*se
            hi <- wh + qnorm(1 - alpha/2)*se
            wh <- cbind(wh, lo=lo, hi=hi)

            colnames(wh) <- c("RL", paste(100*alpha/2, "%", sep = ""),
                              paste(100*(1 - alpha/2), "%", sep = ""))
            wh
        } # ci.fun
        res <- lapply(1:length(M), ci.fun, object=object, co=co, M=M, res=res, alpha=alpha)
    } # Close if (ci.fit

    if (se.fit){
        se.fun <- function(i, object, co, M, res, alpha){
            wh <- res[[i]]
            se <- getse(object, co, M[i])
            wh <- cbind(RL=wh, se.fit=se)
            wh
        } # ci.fun
        res <- lapply(1:length(M), se.fun, object=object, co=co, M=M, res=res, alpha=alpha)
    }
    
    cov.fun <- function(i,res){
      wh <- res[[i]]
      wh <- addCov(wh,X)
      wh
    }
    res <- lapply(1:length(M), cov.fun,res=res)
    
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
        u <- !duplicated(cbind(X.phi,X.xi))
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
    } else if (se.fit){ warning("se.fit not implemented - ignoring") }
    else if (all){ res <- res }
    else { # Just point estimates
        res <- t(sapply(res, function(x){ apply(x, 2, mean) }))
    }
    
    if(!all){
      res <- addCov(res,X.phi)
      res <- addCov(res,X.xi)
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
    
    # co is a list with one element for each unique item in
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

addCov <- function(res,X){ # used in predict.link.* to add covariates to columns reported in output
  if(!is.null(dim(X))){
    if(dim(X)[2] > 1){
       cov <- X[,colnames(X) != "(Intercept)"]
       res <- cbind(res,cov)
       if(is.vector(cov)) colnames(res)[dim(res)[2]] <- colnames(X)[colnames(X) != "(Intercept)"]
    }
  }
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
  res <- predict.link.bgpd(object, newdata=newdata, se.fit=se.fit, ci.fit=ci.fit, all=all, unique.=unique.,alpha=alpha)
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

################################################################################
## Method functions

print.rl.gpd <- function(object, digits){
    nms <- names(object)
    newnms <- paste("M =", substring(nms, 3), "predicted return level:\n")
    lapply(1:length(object), function(i, o, title){
                                 cat(title[i])
                                 print(o[[i]])
                                 cat("\n")
                                 NULL}, o=object, title=newnms)
    invisible(object)
}

show.rl.gpd <- summary.rl.gpd <- print.rl.gpd








################################################################################
## test.predict.gpd()

test.predict.gpd <- function(){
# no covariates
  u <- 14
  r.fit <- gpd(rain,th=u)
  co <- coef(r.fit)
  qgpd(0.9,exp(co[1]),co[2],u=u)

  checkEqualsNumeric(target=u,current = predict(r.fit,M=1/r.fit$rate)[[1]],msg="predict.gpd: retrieve threshold")

  t.fit <- r.fit
  t.fit$rate <- 1
  prob <- c(0.5,0.9,0.95,0.99,0.999)
  for(p in prob){
    checkEqualsNumeric(target = qgpd(p,exp(co[1]),co[2],u=u),
                       current = unlist(predict(t.fit,M=1/(1-p))),msg="predict.gpd: est ret levels no covariates")
  }

# with covariates

  n <- 1000
  M <- 1000
  X <- data.frame(a = rnorm(n),b = runif(n,-0.3,0.3))
  Y <- rgpd(n,exp(X[,1]),X[,2])
  X$Y <- Y
  fit <- gpd(Y,data=X,phi=~a,xi=~b,th=0)
  co <- coef(fit)
  sig <- exp(cbind(rep(1,n),X[,1]) %*% co[1:2])
  xi <- cbind(rep(1,n),X[,2]) %*% co[3:4]

  checkEqualsNumeric(target=qgpd(1-1/M,sig,xi,u=0),
                     current = predict(fit,M=M)[[1]][,1],msg="predict.gpd: ret level estimation with covariates")
 
# check multiple M
  M <- c(10,50,100,500,1000)

  target <- sapply(M,function(m,sig,xi,u) qgpd(1-1/m,sig,xi,u=u),sig=sig,xi=xi,u=0)
  current <- predict(fit,M=M)

  for(i in 1:length(M)){
    checkEqualsNumeric(target[,i],current[[i]][,1],msg="predict.gpd: ret level estimation multiple M")
  }

# new data
  nx <- 20
  M <- 1000
  newX <- data.frame(a=runif(nx,0,5),b=runif(nx,-0.1,0.5))

  sig <- exp(co[1] + newX[[1]] * co[2])
  xi <- co[3] + newX[[2]] * co[4]

  checkEqualsNumeric(target=qgpd(1-1/M,sig=sig,xi=xi,u=0),current=predict(fit,M=M,newdata=newX)[[1]][,1],msg="predict.gpd: ret level ests with new data")

  checkEqualsNumeric(c(nx,5),dim(predict(fit,ci=TRUE,newdata=newX)[[1]]), msg="predict.gpd: dimension or return object for ci calc")
  checkEqualsNumeric(c(nx,4),dim(predict(fit,se=TRUE,newdata=newX)[[1]]), msg="predict.gpd: dimension or return object for se calc") 
  checkEqualsNumeric(c(nx,6),dim(predict(fit,se=TRUE,ci=TRUE,newdata=newX)[[1]]), msg="predict.gpd: dimension or return object for se and ci calc")

  checkEquals(c("RL","2.5%","97.5%","se.fit","a","b"), colnames(predict(fit,se=TRUE,ci=TRUE)[[1]]), msg="predict.gpd: colnames of return obejct for se and ci calc, default alpha")
  checkEquals(c("RL","5%","95%","se.fit","a","b"), colnames(predict(fit,se=TRUE,ci=TRUE,alpha=0.1)[[1]]), msg="predict.gpd: colnames of return obejct for se and ci calc, alpha=0.1")

# alpha
  alpha <- c(0.01,0.05,0.1,0.2,0.5,0.9,0.99)
  z <- matrix(qnorm(c(alpha/2,1-alpha/2)),ncol=2)

  for(a in 1:length(alpha)){
    p <- predict(fit,alpha=alpha[a],ci=TRUE,se=TRUE)[[1]]
    checkEquals(current = colnames(p)[2:3],target = c(paste(100*alpha[a]/2,"%",sep=""),paste(100*(1-alpha[a]/2),"%",sep="")),msg="predict.gpd: labelling of confidence intervals")
    checkEqualsNumeric(target = p[,1] + z[a,1]*p[,4],current = p[,2], msg="predict.gpd: ret level Conf Interval calc for different alpha")
    checkEqualsNumeric(target = p[,1] + z[a,2]*p[,4],current = p[,3], msg="predict.gpd: ret level Conf Interval calc for different alpha")
  }

# linear predictors

  checkEqualsNumeric(target = cbind(co[1] + newX[[1]] * co[2],
                                    co[3] + newX[[2]] * co[4]),
                     current = predict(fit,newX,type="lp")[,1:2], msg="predict.gpd: linear predictors")

  checkEqualsNumeric(target = c(nx,6),dim(predict(fit,newX,se=TRUE,type="lp")), msg="predict.gpd: dimension of return object, linear predictor, se calc")
  checkEqualsNumeric(target = c(nx,8),dim(predict(fit,newX,ci=TRUE,type="lp")), msg="predict.gpd: dimension of return object, linear predictor, ci calc")
  
  name <- c("phi", "xi", "phi.lo", "phi.hi", "xi.lo", "xi.hi")
  checkEquals(target = name, current = colnames(predict(fit,newX,ci=TRUE,type="lp"))[1:6],msg="predict.gpd: colnames for linear predictor return object")
  checkEquals(target = name, current = colnames(predict(fit,newX,ci=TRUE,se=TRUE,type="lp"))[1:6],msg="predict.gpd: colnames for linear predictor return object")

# unique

  newX <- data.frame(a=c(0,0,0,1,1,1,2,2,2,3,3,3,4,4,4),b=c(-.1,.1,.1,-.1,.1,.1,-.1,.1,.1,-.1,.1,.1,-.1,.1,.1)) # checks for duplicates in one and both covariates.
  U <- !duplicated(newX)
  checkEqualsNumeric(current = predict(fit,newX,type="lp"),
                     target = predict(fit,newX,unique.=FALSE,type="lp")[U,], msg="predict.gpd: functioning of argument unique, for linear predictor")
  checkEqualsNumeric(current = predict(fit,newX)[[1]],
                     target =  predict(fit,newX,unique.=FALSE)[[1]][U,], msg="predict.gpd: functioning of argument unique, for return levels")

# check standard errors - this takes a while since using bootstrap
    
  M <- c(10,100,500,1000,2000)
  newX <- data.frame("a"=c(1,2),"b"=c(-0.1,0.1))
  fit.p <- predict(fit, newdata=newX,se=TRUE,M=M)
  fit.seest <- unlist(lapply(fit.p,function(x) x[,2]))

  fit.b <- bootgpd(fit,R=1000)
  fit.bp <- predict(fit.b,newdata=newX,all=TRUE,M=M)
  fit.seb <- lapply(fit.bp,function(X) apply(X,2,sd))
  fit.seboot <- unlist(fit.seb)

  checkEqualsNumeric(rep(0,length(fit.seboot)), (fit.seboot -  fit.seest) / fit.seest, tol=0.1,msg="predict.gpd: standard error estimate compared with bootstrap standard errors")
}

test.predict.bgpd <- function(){
# no covariates
  u <- 14
  r.fit <- gpd(rain,th=u,method="sim",trace=20000)

  checkEqualsNumeric(target=u,current=predict(r.fit,M=1/r.fit$map$rate)[[1]], msg="predict.bgpd: retreive threshold")

  t.fit <- r.fit
  t.fit$map$rate <- 1
  p <- c(0.5,0.9,0.95,0.99,0.999)
  checkEqualsNumeric(target = sapply(p,function(p)mean(qgpd(p,exp(t.fit$param[,1]),t.fit$param[,2],u))),
                     current = unlist(predict(t.fit,M=1/(1-p))),msg="predict.bgpd: ret level estimation")

# with covariates

  n <- 100
  M <- 1000
  X <- data.frame(a = rnorm(n),b = runif(n,-0.3,0.3))
  Y <- rgpd(n,exp(X[,1]),X[,2])
  X$Y <- Y
  fit <- gpd(Y,data=X,phi=~a,xi=~b,th=0,method="sim",trace=20000)

  sig <- apply(fit$param,1,function(v)exp(cbind(rep(1,n),X[,1]) %*% v[1:2]))
  xi <-  apply(fit$param,1,function(v)    cbind(rep(1,n),X[,2]) %*% v[3:4]) 
  
  qbgpd <- function(p,sig,xi,u)sapply(1:dim(sig)[1],function(i) mean(qgpd(p,sig[i,],xi[i,],u)))

  checkEqualsNumeric(target = qbgpd(1-1/M,sig,xi,0),
                     current = predict(fit,M=M)[[1]],msg="predict.bgpd: ret level estimation with covariates")

# check multiple M
  M <- c(10,50,100,500,1000)

  temp1 <- sapply(M,function(m)qbgpd(1-1/m,sig,xi,u=0))
  temp2 <- predict(fit,M=M)

  checkEqualsNumeric(target = temp1[,1],current = temp2[[1]],msg="predict.bgpd multiple M")
  checkEqualsNumeric(target = temp1[,2],current = temp2[[2]],msg="predict.bgpd multiple M")
  checkEqualsNumeric(target = temp1[,3],current = temp2[[3]],msg="predict.bgpd multiple M")
  checkEqualsNumeric(target = temp1[,4],current = temp2[[4]],msg="predict.bgpd multiple M")
  checkEqualsNumeric(target = temp1[,5],current = temp2[[5]],msg="predict.bgpd multiple M")

# new data
  nx <- 20
  M <- 1000
  newX <- data.frame(a=runif(nx,0,5),b=runif(nx,-0.1,0.5))
  predict(fit,newX)

  sig <- exp(cbind(rep(1,nx),newX[,1]) %*% t(fit$param[,1:2]))
  xi <-      cbind(rep(1,nx),newX[,2]) %*% t(fit$param[,3:4])

  checkEqualsNumeric(target = qbgpd(1-1/M,sig=sig,xi=xi,u=0),
                     current = predict(fit,M=M,newdata=newX)[[1]],msg="predict.bgpd : ret level estimation with new data")

  checkEqualsNumeric(target = c(n,4), current = dim(predict(fit,ci=TRUE)[[1]]), msg="predict.bgpd: dimension of output with ci calculation")
  checkEqualsNumeric(target = c(n,4), current = dim(predict(fit,se=TRUE,ci=TRUE)[[1]]), msg="predict.bgpd: dimension of output with ci calculation")

  checkEquals(target = c("2.5%","50%","97.5%","Mean"), colnames(predict(fit,ci=TRUE)[[1]]), msg="predict.bgpd: colnames of ret level ests with CI estimation")
  checkEquals(target = c("5%","50%","95%","Mean"), colnames(predict(fit,alpha=0.1,ci=TRUE)[[1]]), msg="predict.bgpd: colnames of ret level ests with CI estimation, alpha=0.1")

# check linear predictors

  checkEqualsNumeric(current = predict(fit,newX,type="lp")[,1:2],
                    target = cbind(apply(cbind(rep(1,nx),newX[,1]) %*% t(fit$param[,1:2]),1,mean),
                                   apply(cbind(rep(1,nx),newX[,2]) %*% t(fit$param[,3:4]),1,mean)),
                    msg = "predict.bgpd: linear predictor estimates")

  checkEqualsNumeric(target = c(nx,10), dim(predict(fit,newX,ci=TRUE,type="lp")), msg="predict.bgpd: dimension of pinear predictor return object")

  cnames <- c("phi: 2.5%", "phi: 50%", "phi: 97.5%", "phi: Mean", "xi: 2.5%", "xi: 50%", "xi: 97.5%", "xi: Mean")# this specific format assumed by plot.rl.bgpd 

  checkEquals(current = cnames, target = colnames(predict(fit,newX,ci=TRUE,type="lp"))[1:8], msg="predict.bgpd: col names of lin predictors with CI calcs")
  checkEquals(current = cnames, target = colnames(predict(fit,newX,ci=TRUE,se=TRUE,type="lp"))[1:8], msg="predict.bgpd: col names of lin predictors with CI+SE calcs")

# unique
  newX <- data.frame(a=c(0,0,0,1,1,1,2,2,2,3,3,3,4,4,4),b=c(-.1,.1,.1,-.1,.1,.1,-.1,.1,.1,-.1,.1,.1,-.1,.1,.1))

  checkEqualsNumeric(current = predict(fit,newX)[[1]], target = unique(predict(fit,newX,unique.=FALSE)[[1]]),msg="predict.bgpd: unique functioning for ret level ests")
  checkEqualsNumeric(current = predict(fit,newX,type="lp")[,], target = unique(predict(fit,newX,unique.=FALSE,type="lp")[,]),msg="predict.bgpd: unique functioning for lin pred ests")

}

