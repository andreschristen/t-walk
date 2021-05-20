

######### Functions to analyse the output of Runtwalk, saved in info.

######### To plot a time series of the log of the posterios
PlotUs <- function(info, from=0, to=length(info$Us))
{
	plot( from:(to-1), -info$Us[(from+1):to], type="l",
		ylab="Log of Objective", xlab="Iteration", main="");
}
######## log post of the x primes.
PlotUps <- function(info, from=0, to=length(info$Ups))
{
	plot( from:(to-1), -info$Ups[(from+1):to], type="l",
		ylab="Log of Objective", xlab="Iteration", main="");
}

######### To plot a histogram of any parameter
PlotHist <- function( info, par=1, from=0, xlab=paste("Parameter", par), main="", ...)
{
	hist( info$output[from:(info$Tr), par], xlab=xlab, main=main, ...);
}


############ To plot an outline of the output:

####  Plot autocorrelations and IAT times ... not calculated the best way
#### and other things
Ana <- function(info, from=1, to=info$Tr, maxlag=0, par=0, file="")
{
	sel <- from:to;
	
	layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))
	
	accrt <- rep(0, 4);
	for (h in 1:4)
	{
		selh <- which(info$recacc[sel,1] == h);
		accrt[h] <- sum(info$recacc[selh,2])/length(selh);
	}
	
	hist(info$recacc[sel,1], main="", xlab="", freq=FALSE);
	points( 1:4, accrt, col="red");
	##hist(info$recacc[selacc,1], col="red", add=TRUE);	
	
	if (n > 1)
	 hist(sqrt(as.matrix((info$output[sel,]-info$outputp[sel,])*(info$output[sel,]-info$outputp[sel,])) %*%   matrix(1, nrow=n)), main="", xlab="Step sizes", breaks=50)
	else
		hist(abs(info$output - info$outputp), main="", xlab="Step sizes", breaks=50);
	
	PlotUs(info, from, to);
	
	#PlotUps(info, from, to);
	
	Tint <- PlotAutoCorrTime( info, maxlag=maxlag, par=par);
	
	itmap = which(-info$Us == max(-info$Us))[1];
	
	
	cat( "Ratio of moved coodinates per it=\n",
		 accrt[1], accrt[2], accrt[3], accrt[4],
		 "\nn=", n, "AcceptanceRatio=", info$acc/info$Tr,
		"MAPlogPost=", -info$Us[itmap], "Tint/n=", Tint/n,"\n\n");

	if (file != "")
	 cat(file=file, n, 
		 accrt[1], accrt[2], accrt[3], accrt[4], info$acc/info$Tr,
		-info$Us[itmap], Tint/n,"\n");
	
}

####  Calculate IAT and acceptance ratios
Ana2 <- function(info, from=1, to=info$Tr, par=0, file="")
{
	sel <- from:to;
		
	accrt <- rep(0, 4);
	for (h in 1:4)
	{
		selh <- which(info$recacc[sel,1] == h);
		accrt[h] <- sum(info$recacc[selh,2])/length(selh);
	}

	
	#### No plots

	Tint <- IAT( info, par=par);
	
	itmap = which(-info$Us == max(-info$Us))[1];
	
	
	cat( "Ratio of moved coodinates per it=\n",
		 accrt[1], accrt[2], accrt[3], accrt[4],
		 "\nn=", info$n, "AcceptanceRatio=", info$acc/info$Tr,
		"MAPlogPost=", -info$Us[itmap], "IAT=", Tint, "IAT/n=", Tint/info$n,"\n\n");

	if (file != "")
	 cat(file=file, n, 
		 accrt[1], accrt[2], accrt[3], accrt[4], info$acc/info$Tr,
		-info$Us[itmap], Tint/info$n,"\n");
	
}
	

### Plot time series of parameters

TS <- function(info, pars=1:n, from=1, to=info$Tr, prime=FALSE)
{	
	sel <- from:to;
	
	if (length(pars) <= 10)
	{
		if (!prime)
			plot(as.ts(as.matrix(info$output[sel, pars])),  main="x")
		else
			plot(as.ts(as.matrix(info$outputp[sel, pars])), main="xp");
	}
	else
		cat("Cannot print time series for more than 10 parameters\n\n");
	
}




