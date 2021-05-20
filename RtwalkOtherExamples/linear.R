


#######  Fit linear models using the t-walk


### Model to be fitted, cubic
ln <- function(ts,x) { x[1] + x[2]*ts + x[3]*(ts)^2 + x[4]*(ts)^3; }
n  <- 5; ### Dimension: num of coeficients + la
X0 <- c( 5, -30, 5, 2, 1/(30^2));  ## True parameters
## Region
mn <- -6;
mx <- 6;
## sample size
m <- 20;

MakeDat <- function( tt=runif( m, min=mn, max=mx), la=X0[n])
{
	dat <<-
		matrix(
	  	append( tt, ln( tt, X0) + rnorm(m,sd=sqrt(1/la))),
		ncol=2); 



#### Good initial values

	#LMx0 <<- X0[1:(n-1)] + rnorm(n-1,sd=30);
	#LMx0 <<- append( LMx0, la*rgamma( 1, shape=1, rate=1));

	#LMxp0 <<- X0[1:(n-1)] + rnorm(n-1,sd=10);
	#LMxp0 <<- append( LMxp0, la*rgamma( 1, shape=1, rate=1));

## just any positive initial values

	#LMx0 <<- rgamma( n-1, shape=3, rate=1); ##rnorm(n-1, sd=3000);
	#LMx0 <<- append( LMx0, rgamma( 1, shape=3, rate=50));

	#LMxp0 <<- rgamma( n-1, shape=3, rate=10); ###rnorm(n-1, sd=3000);
	#LMxp0 <<- append( LMxp0, rgamma( 1, shape=3, rate=50));
	

## just any initial values

	#LMx0 <<- rnorm( n-1, mean=1000, sd=1000); 
	#LMx0 <<- append( LMx0, rgamma( 1, shape=3, rate=5));

	#LMxp0 <<- rnorm( n-1, mean=1000, sd=1000); 
	#LMxp0 <<- append( LMxp0, rgamma( 1, shape=3, rate=5));

### plot data
	
	PlotDatLM();
	PlotLn( X0, col="red");
	PlotLn(LMx0, col="green");
	PlotLn(LMxp0, col="green");

}


### Support
SuppLM <- function(x)
{
	rt <- TRUE;

	#if ((x[1] < 0)) ## && (x[1] > -2))
	#	rt <- TRUE
	#else 	
	#	rt <- FALSE;

	#rt <- (rt && (x[1] > 0)); ##Mean for the decay stuff

	#rt <- (rt && (x[2] > 0)); ##Std Dev for the decay stuff

	#rt <- (rt && (x[4] > 0)); ##Scale K

	rt <- (rt && (x[n] > 0));  ### precision

	rt;
	
	#all(x > 0);  ## for the dammped shock 

}




## x = (x1, x2, ... , x_{n-1}, la)
log2pi <- log(2*pi);
LMU <- function(x)
{
#cat("DU", x, SuppDecay(x), "\n");
	

	ll <- (m/2)*log2pi - (m/2)*log(x[n]) + (x[n]/2)*sum((dat[,2] - ln(dat[,1],x))^2);
	
	#CoefPr <- 0.5*(x[1] - 0.0002)^2/(0.001^2) + (1-5)*log(x[2]) + 0.0001*(x[2]) + (1-7)*log(x[4]) + 7*(x[4]);
	
	PresPr <- (1-5)*log(x[n]) + (5/0.001)*x[n];
	
	ll + PresPr;  # + CoefPr + PresPr;

} 



source("twalk.R");

PlotDatLM <- function(add = FALSE, type="p", pch=21, from=mn, to=mx)
{
	
	plot( dat[,1], dat[, 2], pch=pch,
		type=type, xlab="", ylab="", xlim=c(from,to));

}

PlotLn <- function( x, col="blue", from=mn, to=mx, add=TRUE)
{
	if (add)
		lines( seq( from, to, length=100), ln(seq( from, to, length=100), x), type="l", col=col)
	else
		plot( seq( from, to, length=100), ln(seq( from, to, length=100), x),
			type="l", col=col, xlab="", ylab="");
		
}

RunLM <- function(Tr=10000, x0=LMx0, xp0=LMxp0)
{
 	pl <<- Runtwalk( Tr=Tr, Obj=LMU, Supp=SuppLM, PlotObj=FALSE,
		x0=x0, xp0=xp0);

	Tr <<- Tr;	
	
	PlotFit();

	cat("x0:", x0, "\n", "xp0:", xp0, "\n");

}

PlotFit <- function(mp=-1)
{
	if (mp == -1)  ## MAP
	{
		locmp <- which(-1*(pl$Us) == max(-pl$Us))[1];
		mp <<- locmp;
	}
	PlotDatLM();
	PlotLn(X0, col="red");
	PlotLn(LMx0, col="green");
	PlotLn(LMxp0, col="green");
	PlotLn(pl$output[locmp,]);
	
	cat("True model = Red, MAP = blue, Strat = green\n")
	cat("MAP:", pl$output[locmp,], "LogPost=", -LMU(pl$output[locmp,]), "\n", "X0:", X0, "\n");

}

PlotUs <- function( lr=1, Tr=Tr, ...)
{
	 plot(lr:Tr, -pl$Us[lr:Tr], type="l",
		xlab="Iteration", ylab="LogPost", ...);

}


PlotHist <- function( lr, pr, xlab=paste("x",pr,sep=""), ...)
{
	hist(pl$output[lr:Tr,pr], xlab=xlab, ...);
}



### Other nonlinear models
## ln <- function(ts,x) { x[1] + x[2]*sin(x[3]*ts*(2*pi))*exp(-1*x[4]*ts); }
#n  <- 5; ### Dimension, four coeficients and la
#X0 <- c( 0, 10, 0.4, 0.34, 1/(0.1^2));  ## True parameters
## Region
#mn <- 0;
#mx <- 15;
#m <- 50;  ## sample points

#ln <- function(ts,x) { x[3] + x[4]*exp(-0.5*x[2]*(ts-x[1])^2); }
#n  <- 5; ### Dimension, coeficients and la
#X0 <- c( 0.0001, 86500, 0, 1, 1/(0.03^2));  ## True parameters
## Region
#mn <- 0;
#mx <- 0.02;
#m <- 500;
#LMx0  <- c( 0.001, 1/(0.005^2), 0, 1.1, 1/(0.001^2))
#LMxp0 <- c( 0.0001, 1/(0.004^2), 0.0001, 0.9, 1/(0.004^2))
 

