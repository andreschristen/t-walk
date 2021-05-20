
DecayDat <- function( file, mx=0)
{
	dat <<-  read.csv( file, skip=4);
	
	fnam <<- file;

	if (mx == 0) ### automatic maximum
	{
		mx <- length(dat[,1]);
		autpm <- TRUE;
	}
	else
		autpm <- FALSE;

	mx <<- mx;

	n <<- 5;
	m <<- mx;

#### initial values, ( mu, sig, a, y0, la)

	indmx <- which(dat[1:mx,2] == max(dat[1:mx,2]));
	indmn <- which(dat[1:mx,2] == min(dat[1:mx,2]));
	
	Mxt <<- dat[mx,1];

	Decayx0 <<- c(
		dat[indmx,1], ## mu = argmax(dat[,2])
		abs(dat[indmx,1]-dat[indmn,1])/2, ## range/2
		dat[indmx,2], ## a= max(dat[,2])
		dat[indmn,2], ## y0= min(dat[,2])
		0.001); #1/(mean(diff(dat[,2])^2)) ); ##c(1,2,3,-4,5);
		
	Decayx0 <<- Decayx0*rgamma(5, shape=80, rate=80);

	Decayxp0 <<- c( dat[indmn-100,1],## mu = argmax(dat[,2])
		abs(dat[indmx,1]-dat[indmn,1])/3, ## range/3
		dat[indmx,2]-dat[indmn,2], ## a= max(dat[,2])-min(dat[,2])
		dat[indmn,2]*2, ## y0= min(dat[,2])
		0.0011); #1/(mean(abs(diff(dat[,2])))) ); ##c(6,7,8,-9,10);
		
	Decayxp0 <<- Decayxp0*rgamma(5, shape=80, rate=80);
	

	PlotDecData();
	PlotDecLn(Decayx0, col="green");
	PlotDecLn(Decayxp0, col="green");

	cat("mx=", mx, "\n");
}


## x = ( mu, sig, a, y0, la)
SuppDecay <- function(x)
{
	if ((x[1] > 0) && (x[1] < Mxt))  ## mu
		rt <- TRUE
	else 	
		rt <- FALSE;

	rt <- (rt && (x[2] > 0)); ## sig
	rt <- (rt && (x[3] > 0)); ## a
	##rt <- (rt && (x[4] < 0)); ## y0
	rt <- (rt && (x[5] > 0)); ## la



	rt;
}


## y_0 + a e^{-frac{(t - \mu)^2}{2\sigma^2}}
Decln <- function(ts,x) { x[4] + x[3]*exp(-0.5*((ts - x[1])/x[2])^2); }

## x = ( mu, sig, a, y0, la)
log2pi <- log(2*pi);
DecayU <- function(x)
{
#cat("DU", x, SuppDecay(x), "\n");
	
	#(m/2)*log2pi - (m/2)*log(x[n]) +
	(x[n]/2)*sum((dat[1:mx,2] - Decln(dat[1:mx,1],x))^2);
	

	ll <- (m/2)*log2pi - (m/2)*log(x[n]) + (x[n]/2)*sum((dat[1:mx,2] - Decln(dat[1:mx,1],x))^2);
	
	CoefPr <- 0.5*(x[1] - 0.004)^2/(0.001^2);  ## mu
	CoefPr <- CoefPr + (1-5)*log(x[2]) + (5/120e-6)*(x[2]); ## sigma
	CoefPr <- CoefPr + (1-5)*log(x[3]) + (5/0.3)*(x[3]); ## a
	
	PresPr <- (1-5)*log(x[n]) + (5/0.001)*x[n];
	
	ll + CoefPr + PresPr;


} 

PlotGamma <- function( shape=5, rate)
{
	plot( function(x) { dgamma( x, shape=shape, rate=rate); }, xlim=c(0, shape/rate + 3*sqrt(shape)/rate),
	xlab="", ylab="Gamma density" );
}



source("twalk.R");

PlotDecData <- function(add = FALSE, type="l", main=fnam)
{
	
	plot( dat[,1], dat[,2],
		type=type, xlab="Time", ylab="Ampl.", main=main);

	lines( dat[1:mx,1], dat[1:mx,2], col="red");

}

PlotDecLn <- function( x, add=TRUE, col="blue", main="", ...)
{
	if (add)
		lines( dat[1:mx,1], Decln(dat[1:mx,1],x), col=col)
	else
		plot( dat[1:mx,1], Decln(dat[1:mx,1],x), col=col,
			type="l", xlab="Time", ylab="Ampl.", main=main, ...)
}

RunDec <- function(Tr=10000)
{
 	pl <<- Runtwalk( Tr=Tr, Obj=DecayU, Supp=SuppDecay, PlotObj=FALSE,
		x0=Decayx0, xp0=Decayxp0);

	Tr <<- Tr;	

	mp <<- which(-1*(pl$Us) == max(-pl$Us))[1];

	PlotDecData();
	PlotDecLn(pl$output[mp,]);
	PlotDecLn(Decayx0, col="green");
	PlotDecLn(Decayxp0, col="green");

	cat("Acc:", 100*(pl$acc/Tr), "%\n", "MAP:", pl$output[mp,], "LogPost:", -pl$Us[mp], "\n");
}


PlotUs <-  function(lr=0, main="")
{
	plot( lr:Tr, -pl$Us[(lr+1):(Tr+1)], type="l",
		xlab="Iteration", ylab="LogPost", main=main);
	lines( lr:Tr, -pl$Usp[(lr+1):(Tr+1)], col="blue")
}


###  taud = sig 
PlotHistTaud <- function(lr)
{
	hist( pl$output[lr:Tr,2], xlab="tau_d", main=fnam);
}



		