
####rats 

### Of course it takes far longer to program the a standard BUGS example in R
### and, at least for the very simple "rats" example, the twalk was rouhgly 10
### times less efficeint.  However, the twalk is completly flexible and permits the use
### of arbitrary distributions unlike BUGS.  For example, in the twalk
### we may easily constraint the support
### of some parameters.  However, this results in non standard full conditionals that severly
### undermines the use of BUGS.

### data here

source("twalk.R");
load("RatsData.Rdata");
al0jags <- read.table("RatsJagsAlpha0.out", header=FALSE);

xbar <- mean(Ratsx);
N <- dim(Rats)[1];
T <- dim(Rats)[2]; 

al <- 1.0E-3;
const2 <- 1.0E-4;



reg <- function( a, b, x) { a + b*(x-xbar); }


##        1         N  N+1       2N  2N+1   2N+2      2N+3       2N+4    2N+5
## x = ( a1, ... , aN, b1, ... , bN, tau.c, alpha.c, tau.alpha, beta.c, tau.beta)

RatslU <- function(par) {
	
	U <- 0;
	for (i in 1:N)
		U <- U + sum((Rats[i,] - reg( par[i], par[N+i], Ratsx))^2);
		
	-1*(N*T/2)*log(par[2*N+1]) + 0.5*par[2*N+1]*U;
}

RatsPriorU <- function(par) {
	
	Ua <- 0.0;
	Ub <- 0.0;
	
	for (i in 1:N) { 
		Ua <- Ua + (par[i] - par[2*N+2])^2;
		Ub <- Ub + (par[N+i] - par[2*N+4])^2;
	}
	
	U <- -1*(al-1)*log(par[2*N+1]) + al*par[2*N+1];
	U <- U + 0.5*(par[2*N+3]*(Ua + 2*al) + par[2*N+5]*(Ub + 2*al));
	U <- U + 0.5*const2*(par[2*N+2]^2 + par[2*N+4]^2);
	U <- U - (N/2 + al - 1)*(log(par[2*N + 3]) + log(par[2*N + 5]));
	
	U;
}

RatsU <- function(par) { RatslU(par) + RatsPriorU(par); }

n <- 2*N + 5;

##        1         N  N+1       2N  2N+1   2N+2      2N+3       2N+4    2N+5
## x = ( a1, ... , aN, b1, ... , bN, tau.c, alpha.c, tau.alpha, beta.c, tau.beta)

RatsInits <- function() {
	
	par <- rep(0,n);
	
	par[1:N] <- rnorm(N, mean=250,sd=3);
	par[(N+1):(2*N)] <- rnorm(N, mean=0,sd=0.5);
	par[(2*N+1)] <- rnorm(1, mean=1, sd=0.0001);
	par[(2*N+2)] <- rnorm(1, mean=250, sd=0.1);
	par[(2*N+3)] <- rnorm(1, mean=1, sd=0.0001);
	par[(2*N+4)] <- rnorm(1, mean=6, sd=0.1);
	par[(2*N+5)] <- rnorm(1, mean=1, sd=0.0001);

	par;
}

RatsSupp <- function(par) {
	
	all(par[c((2*N+1),(2*N+3),(2*N+5))] > 0);
}


PlotData <- function() {
	
	layout(matrix(1:6, ncol=2));
	for (b in 0:5) {
		plot( Ratsx-xbar, Rats[5*b+1,], type="b", ylim=c(130,380));
		for (i in (5*b+1):(5*(b+1)))
			lines(Ratsx-xbar, Rats[i,], type="b");
	}
}


#info <- Runtwalk(Tr=150000, Obj=RatsU, Supp=RatsSupp, x0=RatsInits(), xp0=RatsInits(),PlotObj=FALSE);  
#IAT <- PlotAutoCorrTime( list(n=1, Tr=10000, output=al0jags[,2]), par=1)


	
	

