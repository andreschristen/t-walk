
source("twalk.R");


### Para hacer los c·lculos del paper
### Shock absobrver data
dat <- read.table( "shock.dat", header=FALSE);
Scale <- 1000;
x <- dat[,1]/Scale; ##km x 1000
cen <- dat[,2];      ##1, actually observed, 0, random right censored

FailureUnits <- "km x 1000";


### NG prior pars
alpha <- 11.6;
beta <- 2.16;
NGm <- 22.8;
k <- 0.195;


##Utility function pars

ps <- 900;
c <- 700;
v <- ps - c;
te <- 10;
ta <- 20;
Ia <- 0.1; #pita/pite - 1
cr <- c;
q <- 0.05;

UtilityUnits <- "Net profit prop. to CGP"


################################## Program ########################################################


m <- length(x);  ##sample size

r <- sum(cen);  ## r= number of actually observed data points.
		   

##sum of logs of observed data points
logS <- sum(log(x)*cen);

WeibullU <- function(pars)
{
	nu <- pars[1];
	la <- pars[2];
		
	-r*log(la) + la*r*log(nu) - (la-1)*logS + sum((x/nu)^la) ;
}

WeibullSupp <- function(pars)
{
	nu <- pars[1];
	la <- pars[2];

	if ((nu > 0) && (la > 0)) TRUE else FALSE
}

Enu <- NGm;
Ela <- alpha/beta;
Sdnu <- sqrt(beta/(k*(alpha-1)));
Sdla <- sqrt(alpha)/beta;

 
PriorNGU <- function(pars)
{
	nu <- pars[1];
	la <- pars[2];
	
	-(alpha - 0.5)*log(la) + la*(k*0.5*(nu - NGm)^2 + beta);
}

ShockModelU <- function(pars) {  WeibullU(pars) + PriorNGU(pars); }


SimPriorNG <- function()
{
	la <- rgamma(1, shape=alpha, scale=1/beta);
	nu <- rnorm(1, mean=m, sd=1/sqrt(k*la));
	
	##Truncated NG
	while (nu < 0)
		nu <- rnorm(1, mean=m, sd=1/sqrt(k*la));
	
	c( nu, la);
}


SimWeibull <- function(pars)
{
	nu <- pars[1];
	la <- pars[2];
		
	st <- rweibull( 1, shape=la, scale=nu);

	#while (st > 100)
	#	st <- rweibull( 1, shape=la, scale=nu);

	st;
}

PriorPred <- function(ssize=20000)
{	
	output <- NULL;
	for (it in 1:ssize)
	{
		pars <- SimPriorNG(); 
		output <- append(output, SimWeibull(pars));
	}
	
	
	PriorPredOut <<- output;
	
	density( output, adjust=2, from=0);

}

PlotWeibull <- function(pars)
{
	nu <- pars[1];
	la <- pars[2];

	median <- nu*sqrt(log(2));
	xl <- c( max(0, median - nu), median + nu);
	
	plot( xlim=xl, function(t){ dweibull( t, shape=la, scale=nu); },
		xlab="Failure time", ylab=paste("Weibull(", nu, ",", la, ")"));
} 

TIMESSD <- 2.0;



RunMCMC <- function( it, usetwalk=FALSE)
{
	pars <- SimPriorNG(); ##initial values;
	if (!WeibullSupp(pars))
		cat("PRIOR SIMULATION OUT OF MODEL SUPPORT!!!!!!\n")
	CurrModelU <- WeibullU(pars);
	CurrPriorU <- PriorNGU(pars);
	CurrU <- CurrModelU + CurrPriorU;
	simt <- SimWeibull(pars);
	output <- c( pars, simt, -1*CurrU);
	
	#cat(paste("it: 0, ", pars, "\n"));
	
	if (usetwalk) {
		n <-2;
		
	  	info <- Runtwalk( Tr=it, Obj=ShockModelU, Supp=WeibullSupp,
	  		x0=pars, xp0=SimPriorNG(), PlotObj=FALSE);
	  		
	  	MCMCout <<- matrix( data=0, nrow=it+1, ncol=4, byrow=TRUE);
	  	MCMCout[,1:2] <<- info$output;
	  	MCMCout[,4] <<- -1*info$Us;
	  	
	  	for (i in 1:(it+1))
	  		MCMCout[i,3] <<- SimWeibull(info$output[i,]);

		# MCMC output selection, most be improved!!

		Sel <<- seq(from=1000, to=length(MCMCout[,1]), by=20);
	  	LenSel <<- length(Sel);

	}
	else {
	  for (it in (1:it))
	  {
		if (runif(1) < 0.1)
		{
			proppars <- SimPriorNG();  ##proposed jump from the prior
			
			if (WeibullSupp(proppars))
			{
				PropModelU <- WeibullU(proppars);
				PropPriorU <- PriorNGU(proppars);
				#The prior cancells
				A <- exp(CurrModelU -  PropModelU)
			}
			else
				A <- 0; ##out of support
			##cat(paste("it: ", it, "PropPrior: nu:", proppars[1], "la:", proppars[2], "A: ", A, "\n"));
			
		}
		else
		{	#proposed random walk jump
			proppars <- pars +
				c(rnorm( 1, 0, sd=TIMESSD*Sdnu), rnorm(1, 0, sd=TIMESSD*Sdla));

			if (WeibullSupp(proppars))
			{
				PropModelU <- WeibullU(proppars);
				PropPriorU <- PriorNGU(proppars);
				#Metropolis jump
				A <- exp(CurrU -  (PropModelU + PropPriorU))
			}
			else
				A <- 0; ##out of support
			##cat(paste("it: ", it, "PropRW: nu:", proppars[1], "la:", proppars[2], "A: ", A, "\n"));
		}
			 
		
		#print(paste( proppars[1], proppars[2], CurrU - PropU));
		
		if (runif(1) < A)
		{
			#accepted
			pars <- proppars;
			CurrModelU <- PropModelU;
			CurrPriorU <- PropPriorU;
			CurrU <- CurrModelU + CurrPriorU;
			simt <- SimWeibull(pars);
			
		}

		output <- append( output, c( pars, simt, -1*CurrU));
	  }

	  MCMCout <<- matrix( data=output, nrow=it+1, ncol=4, byrow=TRUE);
	  
   	# MCMC output selection, most be improved!!

	  Sel <<- seq(from=500, to=length(MCMCout[,1]), by=20);
	  LenSel <<- length(Sel);

	}	

	plot(ts(MCMCout[Sel,], names=c("nu", "la", "Pred", "-U")), main="");
 
}


PostPred <- function( sel=Sel)
{
	density(MCMCout[ sel, 3], adjust=2, from=0)
}


PlotAll <- function(xsup=-1, ysup=-1, main="", xlab=paste("Time to failure,", FailureUnits), ...)
{

	postpred <- PostPred();
	priorpred <- PriorPred();
	if (xsup <= 0) ##automatic axis
		xsup <- max(postpred$x, priorpred$x);
	if (ysup <= 0)
		ysup <- 1.15*max(postpred$y, priorpred$y);

	plot(postpred$x, postpred$y, type="l",
		main=main, xlab=xlab, ylab="Density",
		xlim=c(0,xsup), ylim=c(0,ysup), ...);

	hist(MCMCout[ Sel, 3], add=TRUE, freq=FALSE, breaks=50);
	
	points( x[which(cen == 0)], runif( m-r, max=0.06*ysup));
	points( x[which(cen == 1)], runif( r, min=0.07*ysup, max=0.12*ysup), pch=20);
	
	lines(priorpred$x, priorpred$y, lty=2);
}



#########  Utility

A1 <- log((Ia+1)/Ia)/te;
A3 <- 0;
A4 <- cr/v;
A5 <- q*ps/v;

DistFn <- function(t) 
{

	length(which(MCMCout[Sel,3] <= t))/LenSel;
  
}

Ind <- function(t,tw)
{
	if (t <= tw)
		1
	else
		0;
} 

Utility <- function(t,tw)
{
	(Ia+1)*(1 - exp(-A1*tw))*(1 -
	( A4*(1-(A3*t/tw)) + A5*(1-(t/tw)) )*Ind(t,tw));
}

ExpUtility <- function(twrange, sample=MCMCout[Sel,3], main="",
		xlab=paste("Warranty,", FailureUnits), ylab=paste("Expected utility,", UtilityUnits), ...)
{
	output <- NULL;

	for (ttw in twrange)
	{
		#Settw(ttw);

#tw <<- ttw;
#PostCDF <<- pnorm( tw, mean=22, sd=5);

		ex <- mean(sapply( sample, function(t) { Utility(t,ttw); }));

#ex <- mean(sapply( rnorm( 2000, mean=22, sd=5), Utility));

		output <- append(output,ex);
	}

	mx <- max(output);
	twstar <- twrange[which(mx == output)];

	plot( twrange, output, type="l",
		main=main, xlab=xlab, ylab=ylab, ...);
		
	list( tw=twrange, ExpUt=output, twop=twstar, MaxEx=mx);
}

#PlotAll(xsup=50, ysup=0.085, main="ComparaciÛn, roja=Mi Pred. Post., \nnegra=Pred. a priori" )
# ExpUtility(seq(from=1, to=20, by=0.1))

