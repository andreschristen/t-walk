
GetAutoCorr <- function( info, par=0, from=1, to=info$Tr, lag=30*n) {
	
	if (par>0)
		cor( info$output[from:(to-lag), par], info$output[(from+lag):to, par])
	else
		cor( info$Us[from:(to-lag)], info$Us[(from+lag):to]);
}
	
	

GetAutoCov <- function( dt, lags)
{
	n <- length(dt);
	
	aut <- rep( 0.0, length(lags));
	mu <- mean(dt);
	
	for (i in 1:length(lags))
	{
		lg <- lags[i]
		aut[i] <- sum( (dt[1:(n-lg)]-mu)*(dt[(lg+1):n]-mu) )/n;
		#cat( i, lg, aut[i], "\n")
	}
	
	aut;
}




PlotAutoCorrTime <- function( info, rholimit= 0.1,
	maxlag=0, lagslen=min(maxlag, 50), par=0, from=1, to=info$Tr) {

	##-lag/log(GetAutoCorr( info, lag, par=par, from=from, to=to));
	
	if (par>0) {
		if (info$n > 1)
			dt <-  info$output[from:to, par]
		else
			dt <-  info$output[from:to]
	}
	else
		dt <-  info$Us[from:to];


	####  Automatically find a maxlag for IAT calculations
	if (maxlag == 0) { 
		s2 <- var(dt);   
		rho <- GetAutoCov( dt, lags=1)/s2; ### lag one autocorrelation
	##cat("maxlag=", maxlag, "rho", abs(rho), "\n");
		
		#### if autocorrelation is like exp(- lag/lam) then
		lam <- -1.0/log(abs(rho));
		### Our initial guess for maxlag is twice lam (eg. twice the mean life)
		maxlag <- floor(2.0*lam)+1;
				
		### We take 1% of lam to jump forward and look for the		### rholimit threshold		jmp <- ceiling(-0.01/log(abs(rho)))+1;
			
		while ((abs(rho) > rholimit) && (maxlag < length(dt)/2)) {
			rho <- GetAutoCov( dt, lags=maxlag)/s2;
			maxlag <- maxlag + jmp;
	##cat("maxlag=", maxlag, "rho", abs(rho), "\n");
		}
		
		maxlag <- floor(1.1*maxlag);  #10% more	
	}

	if (maxlag >= length(dt)/2) { ###not enough data
		cat("Not enough data, maxlag=", maxlag, "\n");
	NaN;
	}
	else {
	
	maxlag <- max( 2, maxlag);
	
	cat("maxlag=", maxlag, "rho", abs(rho), "\n");
	
	### We always take a maximum of 50 autocovariances
	lagslen <- min(maxlag, 50);
	### from 0 to maxlag
	lags <- floor(seq( 0, maxlag, length=lagslen));
	
	aut <- GetAutoCov( dt, lags=lags);
	aut <- aut/aut[1];  ##correlations, aut[1] = variance
	
	##intaut <- 1+2*cumsum(aut[2:lagslen]);
	##intaut <- -1*(2:lagslen)/log(aut[2:lagslen]);
	intaut <- 1 +
	   2*(cumsum(aut[1:(lagslen-1)]*diff(lags) + diff(aut)*(diff(lags)-1)/2) -1);
	
	plot( lags[2:lagslen], intaut, type="l",
		xlim=c(0,lags[lagslen]), ylim=c(0,intaut[lagslen-1]),
		xlab="Lag", ylab="Est. int. Autocorrelation time");
		
	lines( lags, intaut[lagslen-1]*aut, col="red"); ##also plot the autocorelations
	
	axis( 4, seq(0,intaut[lagslen-1], length=5), paste(seq(0,1,length=5)), col="red");
	
	##cat( "tau/n=", (-1*(lags[lagslen])/log(aut[lagslen]))/n , "\n");
	intaut[lagslen-1];
	}
}
	

##### A much better, although slower, way to calculate the
##### Integrated Autocorrelation Time (not plots though)
IAT <- function( info, par=0, from=1, to=info$Tr) {

	##-lag/log(GetAutoCorr( info, lag, par=par, from=from, to=to));
	
	
	#### we get the desired time series, the parameter and the from - to selection
	if (par>0) {
		if (info$n > 1)
			dt <-  info$output[from:to, par]
		else
			dt <-  info$output[from:to]
	}
	else
		dt <-  info$Us[from:to];

	n <- to-from;	
	mu <- mean(dt); ### with its mean and variance
	s2 <- var(dt);
	
	### The maximum lag is half the sample size
	maxlag <- max( 3, floor(n/2));
	
	#### The gammas are sums of two consecutive autocovariances
	Ga <- rep(0,2); ## two consecutive gammas
	
	lg <- 0;
	Ga[1] <- s2; #sum( (dt[1:(n-lg)]-mu)*(dt[(lg+1):n]-mu) )/n;
	lg <- 1;
	Ga[1] <- Ga[1] + sum( (dt[1:(n-lg)]-mu)*(dt[(lg+1):n]-mu) )/n;
	
	m <- 1;
	lg <- 2*m;
	Ga[2] <- sum( (dt[1:(n-lg)]-mu)*(dt[(lg+1):n]-mu) )/n;
	lg <- 2*m+1;
	Ga[2] <- Ga[2] + sum( (dt[1:(n-lg)]-mu)*(dt[(lg+1):n]-mu) )/n;

	IAT <- Ga[1]/s2; ### Add the autocorrelations
	
	
	### RULE: while Gamma stays positive and decreasing
	while  ((Ga[2] > 0.0) & (Ga[2] < Ga[1])) {
		m <- m+1;
		if (2*m+1 > maxlag) {
			cat("Not enough data, maxlag=", maxlag, "\n");
			break;
		}
		Ga[1] <- Ga[2];
		
		lg <- 2*m;
		Ga[2] <- sum( (dt[1:(n-lg)]-mu)*(dt[(lg+1):n]-mu) )/n;
		lg <- 2*m+1;
		Ga[2] <- Ga[2] + sum( (dt[1:(n-lg)]-mu)*(dt[(lg+1):n]-mu) )/n;
		
		IAT <- IAT + Ga[1]/s2;
	}
	 
	IAT <- -1 + 2*IAT;  ##Calculates the IAT from the gammas
	
	#cat("IAT: IAT=", IAT, ", last lag=", 2*m+1, ", last Gamma=", Ga[1], "\n");

	IAT;
	
}
	

PlotAutoCov <- function( dt, maxlag, main="", col="black", add=FALSE, ...)
{
	aut <- GetAutoCov( dt, maxlag);
	aut <- aut/aut[1];
	if (!add)
	 plot( 0:maxlag, aut, type="l",
		xlab="Lagg", ylab="Normalized Autocovariance",
		main=main, ylim=c(0,1), col=col, ...)
	else
	 lines( 0:maxlag, aut, col=col, ...);
	
	lines( c(0,maxlag), c(0,0));
	
	gams <- (aut[1:maxlag] + aut[2:(maxlag+1)])/2;
	
	lines((1:maxlag)/2, gams, lty="dashed", col=col, ...);	
	
	if (all(gams > 0))
		cutoff <- maxlag+1  ###no cutoff within maxlag
	else
		cutoff <- min(which(gams <= 0)); ## positive Gamas
	if (any(diff(gams) > 0))
		cutoff <- min( cutoff, which(diff(gams) > 0)); ##decreasing
	cutoff <- cutoff - 1;
	
		
	lines( c(0,maxlag), c(aut[cutoff],aut[cutoff]), lty=2, col=col);
	lines( c(cutoff-1,cutoff-1), c(0,aut[1]), col=col);
	
	axis( 3, cutoff-1, col=col);
		
	autv <- aut[1];
	if (cutoff > 1) 
		autv <- autv + 2*sum(aut[2:cutoff]);
	 	
	sqrt(autv);

}

PlotHistMCMC <- function( dt, maxlag, z=qnorm(0.95), main="", ...)
{
	aut <- GetAutoCov( dt, maxlag);
	
	gams <- (aut[1:maxlag] + aut[2:(maxlag+1)])/2;
		
	cutoff <- min(which(gams <= 0)); ## positive Gamas
	cutoff <- min( cutoff, which(diff(gams) > 0)); ##decreasing
	cutoff <- cutoff - 1;
				
	autsd <- aut[1];
	if (cutoff > 1) 
		autsd <- autsd + 2*sum(aut[2:cutoff]);
	 	
	autsd <- sqrt(autsd)/sqrt(length(dt));
	mu <- mean(dt);
	
	xl1 <- min( dt, mu-z*autsd);
	xl2 <- max( dt, mu+z*autsd);

	hist( dt, breaks=20, main=main,
		xlim=c( xl1, xl2), ...);

	axis( 3, c( mu-z*autsd, mu, mu+z*autsd));
	abline(v=c( mu-z*autsd, mu, mu+z*autsd), col="red")
	
	autsd;
}


