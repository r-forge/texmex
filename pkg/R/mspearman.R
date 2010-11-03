#---------------------------------------------------------------------------
#Function MCS() evaluates the MCS wrt to Schmid et al. 2006
#MCS() uses MCSlower() and MCSupper() according to the 
#3rd argument user feeds. lower and upper correspond to 
#the upper and lower extremes. First and second args are
#X (matrix) and p (vector of probabilities).
#Matrix X has to be i*j where i=1,...,d and j=1,...,n, and d,n
#denote dimension and length, respectively.	
#---------------------------------------------------------------------------

MCSlower <- function(X,p)
  {
    d <- dim(X)[1]
    n <- dim(X)[2]
    U <- t(apply(X,1,edf)) #transpose cause apply transposes g(X), g:edf
    res1         <- p-U
    res1[res1<0] <- 0
    res2         <- apply(res1,2,prod)
    res3         <- sum(res2)
    res4         <- ( (1/n)*res3-((p^2)/2)^(d) )/
      ( (p^(d+1))/(d+1) -((p^2)/2)^(d))
    return(res4)
  }


MCSupper <- function(X,p)
  {
    #Matrix X has to be i*j where i=1,...,d and j=1,...,n
    d    <- dim(X)[1]
    n    <- dim(X)[2]
    U    <- t(apply(X,1,edf)) #transpose cause apply transposes g(X), g:edf
    res1 <- ifelse( U<p, (1-p), 1-U )
    res2 <- apply(res1,2,prod)
    res3 <- sum(res2)
    res4 <-  ( (1/n)*res3-((1/2)-(p^2)/2)^(d) )
    res5 <-  (((-1)^(d))/(d+1))*((p-1)^(d))*(1+d*p)-((1/2)-(p^2)/2)^(d) 
    res6 <- res4/res5
    return(res6)
  }


MCS <- function(X,p=seq(.1, .9, by=.1),method="upper") {
    theCall <- match.call()
    method <- ifelse(method == "upper", "MCSupper",
                 ifelse(method == "lower", "MCSlower", "error")
              ) # Close ifelse
    if (method == "error"){
        stop("method should be \'upper\' or \'lower\'")
    }
    
    X <- t(X) # Yiannis's original code had variables as rows
    n    <- length(p)
    res1 <- vector('list',n)
    res2 <- NULL
    for(i in 1:n){
        res1[[i]] <- X
      }
    res2 <- mapply(method,res1,p)
    res <- list(mcs=res2, p=p, method=method, call=theCall)
    oldClass(res) <- "MCS"
    res
  }

plot.MCS <- function(x, xlab="p", ylab= if(x$method=="MCSupper"){ "Upper MCS" } else{"Lower MCS"}, ...){
   plot(x$p, x$mcs, type="l", xlab=xlab, ylab=ylab, ...)
   invisible()
}

print.MCS <- function(x, ...){
    print(x$call)
    wh <- ifelse(x$method == "MCSupper", "Upper", "Lower")
    cat("Multivariate conditional Spearman's rho:\n",
        wh, " tail method.\n\n", sep = "")
    res <- x$mcs
    names(res) <- x$p
    print(res)
    invisible(res)
}
summary.MCS <- show.MCS <- print.MCS

#------------------------------------------------
#Bootstrap
#------------------------------------------------

# XXX HS Rewrite to make it take an object returned
# by MCS. (Maybe.)

bootMCS <- function(X,p=seq(.1, .9, by=.1),method="upper", B=100, trace=10) {
   theCall <- match.call()
   bfun <- function(i, data, p, method){
       if (i %% trace == 0){ cat("Replicate", i, "\n") }
       d <- data[sample(1:nrow(data), replace=TRUE),]
       MCS(d, p, method)$mcs
   }

   res <- sapply(1:B, bfun, data=X, p=p, method=method)
   res <- list(replicates=res, p=p, method=method, B=B, call=theCall)
   oldClass(res) <- "bootMCS"
   invisible(res)
}

plot.bootMCS <- function(x, xlab="p", ylab= if(x$method=="upper"){ "Upper MCS" } else{"Lower MCS"},alpha=.05, ...){
   m <- rowMeans(x$replicates)
   ci <- apply(x$replicates, 1, quantile, prob=c(1-alpha/2, alpha/2))
   plot(x$p, m, type="l", ylim=range(ci),
        xlab=xlab, ylab=ylab,
        sub=paste(100*(1-alpha), "% interval. ", x$B, " bootstrap samples were performed", sep=""))
   lines(x$p, ci[1,], lty=2)
   lines(x$p, ci[2,], lty=2)
   invisible(ci)
}

print.bootMCS <- function(x, ...){
    print(x$call)
    wh <- ifelse(x$method == "upper", "Upper", "Lower")
    cat("Multivariate conditional Spearman's rho:\n",
        wh, " tail method.\n", x$B, " bootstrap samples were performed.\n\n",
        sep = "")
     
    m <- rowMeans(x$replicates)
    names(m) <- x$p
    print(m)
    invisible(m)
}
show.bootMCS <- print.bootMCS

summary.bootMCS <- function(x, alpha=.05){
    wh <- ifelse(x$method == "upper", "Upper", "Lower")
    cat("Multivariate conditional Spearman's rho:\n",
        wh, " tail method.\n", x$B, " bootstrap samples were performed.\n\n",
        sep = "")
    m <- rowMeans(x$replicates)
    ci <- apply(x$replicates, 1, quantile, prob=c(alpha/2, 1 - alpha/2))
    res <- cbind(m, t(ci))
    dimnames(res) <- list(x$p, c("Mean", rownames(ci)))
    res
}