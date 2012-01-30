plotrl.gpd <-
function(object, alpha = .050,
         xlab, ylab, main,
         pch= 1, ptcol =2 , cex=.75, linecol = 4 ,
         cicol = 0, polycol = 15, smooth = TRUE ){

    if(dim(object$X.phi)[2] > 1 | dim(object$X.xi)[2] > 1){
      stop("use plot method for object returned by predict.gpd to see rl plots if covariates in model")
    }
    if (missing(xlab) || is.null(xlab)) { xlab <- "Return period" }
    if (missing(ylab) || is.null(ylab)) { ylab <- "Return level" }
    if (missing(main) || is.null(main)) { main <- "Return Level Plot" }

    xdat <- object$y
    n <- length(xdat) / object$rate # Number of obs prior to thresholding

    jj <- seq(-1, max(3.75,log10(n)),by=0.1)

    m <- unique( c(1/object$rate, 10^jj) )
    xm <- matrix(unlist(rl(object,M=m,ci=TRUE,alpha=alpha)),ncol=3,byrow=TRUE)
    colnames(xm) <- c("RL","ci.l","ci.u")

    U <- object$threshold - abs(object$threshold/100)
    plotX <- xm[,"RL"] > U
    
    plotRLgpd(m[plotX],xm[plotX,],polycol,cicol,linecol,ptcol,n,xdat,pch,
              smooth,xlab,ylab,main,xrange=range(m),yrange=range(c(xdat, range(xm[plotX,c(1,3)]))))
              
    invisible(list(m=m, xm=xm))
}

plot.rl.gpd <- function(object,
         xlab, ylab, main,
         pch= 1, ptcol =2 , cex=.75, linecol = 4 ,
         cicol = 0, polycol = 15, smooth = TRUE, sameAxes=TRUE ){
    if (missing(xlab) || is.null(xlab)) { xlab <- "Return period" }
    if (missing(ylab) || is.null(ylab)) { ylab <- "Return level" }
    if (missing(main) || is.null(main)) { main <- "Return Level Plot" }

    nm <- length(names(object))
    nd <- dim(object[[1]])[2]
    ncov <- length(unlist(object)) / (nm * nd)
    if(nd == 3){
      ValNames <- c("RL","ci.l","ci.u")
    } else {
      stop("Please use ci.fit=TRUE in call to predict, to calculate confidence intervals")
    }
      
    Array <- array(unlist(object),c(ncov,nd,nm),dimnames=list(NULL,ValNames,names(object)))

    m <- as.numeric(substring(dimnames(Array)[[3]],first=3))
   
    if(sameAxes){
       yrange <- range(Array)
    }
    if(length(main) == 1){
      Main <- rep(main,ncov)
    } else {
      Main <- main
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

plotRLgpd <- function(m,xm,polycol,cicol,linecol,ptcol,n,xdat,pch,smooth,xlab,ylab,main,xrange,yrange){

    plot(m, xm[,"RL"], log = "x", type = "n",
         xlim=xrange, ylim=yrange, xlab = xlab, ylab = ylab, main = main)

      if (smooth) {
        splo <- spline(log(m), xm[,"ci.l"] , 200)
        sphi <- spline(log(m), xm[,"ci.u"] , 200)
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
                    c(xm[,"ci.l"],rev(xm[,"ci.u"])),
                    col=polycol, border = FALSE) # Close polygon
        } else {
            lines(m, xm[,"ci.u"], col = cicol)
            lines(m, xm[,"ci.l"], col = cicol)
        }
      }      
	
      lines(m, xm[,"RL"], col = linecol[ 1 ] )

    # Add observed data to the plot
    if(!missing(xdat) & !missing(n)){
      ly <- length(xdat)
      points(1 / (1 - ((n - ly + 1):n) / (n + 1)), sort(xdat), pch=pch, col=ptcol)
      box()
    }
 }
    
