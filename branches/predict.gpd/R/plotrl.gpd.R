plotrl.gpd <- # intended as a diagnostic for a gpd fitted with no covariates. Called by plot.gpd
function(object, alpha = .050,
         xlab, ylab, main,
         pch= 1, ptcol =2 , cex=.75, linecol = 4 ,
         cicol = 0, polycol = 15, smooth = TRUE,RetPeriodRange=NULL ){

    if(dim(object$X.phi)[2] > 1 | dim(object$X.xi)[2] > 1){
      stop("use plot method for object returned by predict.gpd to see rl plots if covariates in model")
    }
    if (missing(xlab) || is.null(xlab)) { xlab <- "Return period" }
    if (missing(ylab) || is.null(ylab)) { ylab <- "Return level" }
    if (missing(main) || is.null(main)) { main <- "Return Level Plot" }

    xdat <- object$y
    n <- length(xdat) / object$rate # Number of obs prior to thresholding

    if(!is.null(RetPeriodRange)){
      ran <- log10(RetPeriodRange)
      jj <- seq(ran[1],ran[2],by=0.1)
    } else {
      jj <- seq(-1, max(3.75,log10(n)),by=0.1)
    }

    m <- unique( c(1/object$rate, 10^jj) )
    xm <- matrix(unlist(rl(object,M=m,ci=TRUE,alpha=alpha)),ncol=3,byrow=TRUE)
    colnames(xm) <- c("RL","ci.l","ci.u")

    U <- object$threshold - abs(object$threshold/100)
    plotX <- xm[,"RL"] > U
    
    xrange <- range(m)
    yrange <- range(c(xdat, range(xm[plotX,c(1,3)])))
    
    plotRLgpd(m[plotX],xm[plotX,],polycol,cicol,linecol,ptcol,n,xdat,pch,
              smooth,xlab,ylab,main,xrange=xrange,yrange=yrange)
              
    invisible(list(m=m, xm=xm))
}


plot.rl.gpd <- function(object, # method for rl.(boot or b)gpd object, which may have covariates.  Plots return level for each unique row in design matrix
         xlab, ylab, main,
         pch= 1, ptcol =2 , cex=.75, linecol = 4 ,
         cicol = 0, polycol = 15, smooth = TRUE, sameAxes=TRUE, type="median" ){
    if (missing(xlab) || is.null(xlab)) { xlab <- "Return period" }
    if (missing(ylab) || is.null(ylab)) { ylab <- "Return level" }
    if (missing(main) || is.null(main)) { main <- "Return Level Plot" }

    nm <- length(names(object))
    nd <- dim(object[[1]])[2]
    ncov <- length(unlist(object)) / (nm * nd)
    ValNames <- colnames(object[[1]])

    if(length(ValNames) == 1 |  !any(ValNames == "2.5%")){
      stop("Please use ci.fit=TRUE in call to predict, to calculate confidence intervals")
    }
      
    Array <- array(unlist(object),c(ncov,nd,nm),dimnames=list(NULL,ValNames,names(object)))

    m <- as.numeric(substring(dimnames(Array)[[3]],first=3))
   
    if(length(main) == 1){
      Main <- rep(main,ncov)
    } else {
      Main <- main
    }
    
    if( class(object) == "rl.bgpd" | class(object) == "rl.bootgpd"){
      if(casefold(type) == "median"){
        Array <- Array[,c(2,1,3),]
      } else if(casefold(type) == "mean") {
        Array <- Array[,c(4,1,3),]
      } else {
        stop("type must be \"mean\" or \"median\" ")
      }
    }

    if(sameAxes){
       yrange <- range(Array)
    }    
  
    for(i in 1:ncov){
      xm <- t(Array[i,,])
      
      if(!sameAxes){ 
        yrange <- range(xm)
      }
      plotRLgpd(m,xm,polycol = polycol,cicol=cicol,linecol=linecol,ptcol=ptcol,pch=pch,
                smooth=smooth,xlab=xlab,ylab=ylab,main=Main[i],xrange=range(m),yrange=yrange)
    }
    
    invisible(list(m=m,xm=Array))
}

plot.rl.bootgpd <- plot.rl.bgpd <- plot.rl.gpd

plotRLgpd <- function(m,xm,polycol,cicol,linecol,ptcol,n,xdat,pch,smooth,xlab,ylab,main,xrange,yrange){ 
# worker function - called by plotrl.gpd, plot.rl.gpd, plot.rl.bgpd

    o <- order(m) # in case the return period are not in ascending order.
    m <- m[o]
    xm <- xm[o,]
    
    plot(m, xm[,1], log = "x", type = "n",
         xlim=xrange, ylim=yrange, xlab = xlab, ylab = ylab, main = main)

      if (smooth) {
        splo <- spline(log(m), xm[,2] , 200)
        sphi <- spline(log(m), xm[,3] , 200)
        if ( polycol != 0 ) {
            polygon( exp(c( splo$x, rev( sphi$x ) )),
	                       c( splo$y, rev( sphi$y ) ),
                     col = polycol , border=FALSE)
         } # Close if (polycol
         lines( exp(splo$x), splo$y, col = cicol )
         lines( exp(sphi$x), sphi$y, col = cicol )
      } else{
        if (polycol != 0){
            polygon(c( m,        rev( m)),
                    c(xm[,2],rev(xm[,3])),
                    col=polycol, border = FALSE) # Close polygon
        } else {
            lines(m, xm[,2], col = cicol)
            lines(m, xm[,3], col = cicol)
        }
      }      
	
      lines(m, xm[,1], col = linecol[ 1 ] )

    # Add observed data to the plot
    if(!missing(xdat) & !missing(n)){
      ly <- length(xdat)
      points(1 / (1 - ((n - ly + 1):n) / (n + 1)), sort(xdat), pch=pch, col=ptcol)
      box()
    }
 }
    
