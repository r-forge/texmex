`gpdRangeFit` <-
function (data, umin=quantile(data, .05), umax=quantile(data, .95),
          nint = 10, show = FALSE,
          penalty="gaussian", priorParameters=NULL,
          xlab="Threshold", ylab=NULL,
		  main=NULL, ... ) {
	if ( missing( ylab ) ){
		ylab = c( "log(scale)", "shape" )
	}
	else if ( length( ylab ) != 2 ){
		stop( "length of ylab should be 2" )
	}

	if ( !missing( main ) && length( main ) != 2 ){
		stop( "length of main should be 2" )
	}

    m <- s <- up <- ul <- matrix(0, nrow = nint, ncol = 2)
    u <- seq(umin, umax, length = nint)
    for (i in 1:nint) {
        z <- gpd(data, th=u[i], penalty=penalty, priorParameters=priorParameters)
        m[i, ] <- z$mle
        m[i, 1] <- m[i, 1] - m[i, 2] * u[i]
        d <- matrix(c(1, -u[i]), ncol = 1)
        v <- t(d) %*% z$cov %*% d
        s[i, ] <- sqrt( diag( z$cov ) )
        s[i, 1] <- sqrt(v)
        up[i, ] <- m[i, ] + 1.96 * s[i, ]
        ul[i, ] <- m[i, ] - 1.96 * s[i, ]
    }
    names <- c("Modified Scale", "Shape")
    for (i in 1:2) {
        um <- max(up[, i])
        ud <- min(ul[, i])
        plot(u, m[, i], ylim = c(ud, um), type = "b",
			xlab=xlab, ylab=ylab[i], main=main[i], ...)
        for (j in 1:nint) lines(c(u[j], u[j]), c(ul[j, i], up[j, i]))
    }
    invisible()
}
