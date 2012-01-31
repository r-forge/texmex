plot.predict.link.gpd <- function(object, main=NULL,
         pch= 1, ptcol =2 , cex=.75, linecol = 4 ,
         cicol = 1, polycol = 15){

  if( !any(colnames(object) == "phi.lo") ){
    stop("Please use ci.fit=TRUE in call to predict, to calculate confidence intervals")
  }

  Ests <- list(phi=object[,c(1,3,4)],xi=object[,c(2,5,6)])
  Names <- c("phi","xi")
  cn <- colnames(object)
  which <- cn != "phi"    & cn != "xi" & 
           cn != "phi.lo" & cn != "phi.hi" & 
           cn != "xi.lo"  & cn != "xi.hi" & 
           cn != "phi.se" & cn != "xi.se"  
           
  X <- object[,which]

  for(i in 1:2){
    for(j in 1:dim(X)[2]){
      if(length(unique(Ests[[i]][,1])) > 1){
        if(length(unique(X[,j])) > 1){
          ord <- order(X[,j])
          x <- X[ord,j]
          y <- Ests[[i]][ord,]
          plot(x, y[,1],type="n",ylab=Names[i],xlab=colnames(X)[j],main=main,ylim=range(y))
          
          if (polycol != 0){
            polygon(c( x,        rev(x)),
                    c(y[,2],rev(y[,3])),
                    col=polycol, border = FALSE) # Close polygon
          } else {
            lines(x, y[,2], col = cicol)
            lines(x, y[,3], col = cicol)
          }
	
          lines(x, y[,1], col = linecol[ 1 ] )
        }
      }
    }
  }
  invisible()
}

plot.predict.link.bgpd <- function(object,type="median",...){

# re-format to same column structure as predict.link.gpd object

  if( casefold(type) == "median"){
    object <- object[,c(2,6,1,3,5,7,9:dim(object)[2])]
  } else if(casefold(type) == "mean") {
    object <- object[,c(4,8,1,3,5,7,9:dim(object)[2])]
  } else {
    stop("type must be \"mean\" or \"median\" ")
  }
  
  colnames(object)[1:6] <-  c("phi", "xi", "phi.lo", "phi.hi", "xi.lo", "xi.hi")

  plot.predict.link.gpd(object,...)
}

plot.predict.link.bootgpd <- plot.predict.link.bgpd

