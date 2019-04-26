library(BiocParallel) 

sva.class2Model <- function(classes) {
  return(model.matrix(~factor(classes)))
}


modefunc <- function(x) {
  return(as.numeric(names(sort(-table(x)))[1]))
}

mono <- function(lfdr){
    .Call("monotone", as.numeric(lfdr), PACKAGE="sva")
}

edge.lfdr <- function(p, trunc=TRUE, monotone=TRUE, transf=c("probit", "logit"), adj=1.5, eps=10^-8, lambda=0.8, ...) {
    pi0 <- mean(p >= lambda)/(1 - lambda)
    pi0 <- min(pi0, 1)
    
    n <- length(p)
    transf <- match.arg(transf)
    
    if(transf=="probit") {
        p <- pmax(p, eps)
        p <- pmin(p, 1-eps)
        x <- qnorm(p)
        myd <- density(x, adjust=adj)
        mys <- smooth.spline(x=myd$x, y=myd$y)
        y <- predict(mys, x)$y
        lfdr <- pi0*dnorm(x)/y
    }
    
    if(transf=="logit") {
        x <- log((p+eps)/(1-p+eps))
        myd <- density(x, adjust=adj)
        mys <- smooth.spline(x=myd$x, y=myd$y)
        y <- predict(mys, x)$y
        dx <- exp(x) / (1+exp(x))^2
        lfdr <- pi0 * dx/y
    }
    
    if(trunc) {
        lfdr[lfdr > 1] <- 1
    }
    if(monotone) {	
        lfdr <- lfdr[order(p)]
        lfdr <- mono(lfdr)
        lfdr <- lfdr[rank(p)]
    }
    
    return(lfdr)
}



# Trims the data of extra columns, note your array names cannot be named 'X' or start with 'X.'
trim.dat <- function(dat){
    tmp <- strsplit(colnames(dat),'\\.')
    tr <- NULL
    for (i in 1:length(tmp)) {
        tr <- c(tr, tmp[[i]][1] != 'X')
    }
    tr
}

# Following four find empirical hyper-prior values
aprior <- function(gamma.hat) {
    m <- mean(gamma.hat)
    s2 <- var(gamma.hat)
    (2*s2 + m^2) / s2
}

bprior <- function(gamma.hat){
    m <- mean(gamma.hat)
    s2 <- var(gamma.hat)
    (m*s2 + m^3) / s2
}

postmean <- function(g.hat,g.bar,n,d.star,t2){
    (t2*n*g.hat + d.star*g.bar) / (t2*n + d.star)
}

postvar <- function(sum2,n,a,b){
    (.5*sum2 + b) / (n/2 + a - 1)
}

# Inverse gamma distribution density function. (Note: does not do any bounds checking on arguments)
dinvgamma <- function (x, shape, rate = 1/scale, scale = 1) {
    # PDF taken from https://en.wikipedia.org/wiki/Inverse-gamma_distribution
    # Note: alpha = shape, beta = rate
    stopifnot(shape > 0)
    stopifnot(rate > 0)
    ifelse(x <= 0, 0, ((rate ^ shape) / gamma(shape)) * x ^ (-shape - 1) * exp(-rate/x))
}

# Pass in entire data set, the design matrix for the entire data, the batch means, the batch variances, priors (m, t2, a, b), columns of the data  matrix for the batch. Uses the EM to find the parametric batch adjustments

it.sol  <- function(sdat,g.hat,d.hat,g.bar,t2,a,b,conv=.0001){
    n <- rowSums(!is.na(sdat))
    g.old <- g.hat
    d.old <- d.hat
    change <- 1
    count <- 0
    while(change>conv){
        g.new <- postmean(g.hat, g.bar, n, d.old, t2)
        sum2 <- rowSums((sdat - g.new %*% t(rep(1,ncol(sdat))))^2, na.rm=TRUE)
        d.new <- postvar(sum2, n, a, b)
        change <- max(abs(g.new-g.old) / g.old, abs(d.new-d.old) / d.old)
        g.old <- g.new
        d.old <- d.new
        count <- count+1
    }
    ## cat("This batch took", count, "iterations until convergence\n")
    adjust <- rbind(g.new, d.new)
    rownames(adjust) <- c("g.star","d.star")
    adjust
}

## likelihood function used below
L <- function(x,g.hat,d.hat){
    prod(dnorm(x, g.hat, sqrt(d.hat)))
}

## Monte Carlo integration functions
int.eprior <- function(sdat, g.hat, d.hat){
    g.star <- d.star <- NULL
    r <- nrow(sdat)
    for(i in 1:r){
        g <- g.hat[-i]
        d <- d.hat[-i]		
        x <- sdat[i,!is.na(sdat[i,])]
        n <- length(x)
        j <- numeric(n)+1
        dat <- matrix(as.numeric(x), length(g), n, byrow=TRUE)
        resid2 <- (dat-g)^2
        sum2 <- resid2 %*% j
        LH <- 1/(2*pi*d)^(n/2)*exp(-sum2/(2*d))
        LH[LH=="NaN"]=0
        g.star <- c(g.star, sum(g*LH)/sum(LH))
        d.star <- c(d.star, sum(d*LH)/sum(LH))
        ## if(i%%1000==0){cat(i,'\n')}
    }
    adjust <- rbind(g.star,d.star)
    rownames(adjust) <- c("g.star","d.star")
    adjust	
} 

## fits the L/S model in the presence of missing data values

Beta.NA <- function(y,X){
    des <- X[!is.na(y),]
    y1 <- y[!is.na(y)]
    B <- solve(crossprod(des), crossprod(des, y1))
    B
}
