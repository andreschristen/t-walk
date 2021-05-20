


#### This program is only for 4 resistors and is difficult to generalise

n <- 4; ### number of resistors and node ... dimension
  
### True resistances
truer <- c( 100, 150, 300, 350);


### Sets the admitance matrix
AdmMat <- function(r) {
	
	matrix( c( 1/r[1]+1/r[2],       -1/r[1],       -1/r[2],             0,
	                 -1/r[1], 1/r[1]+1/r[2],             0,       -1/r[3],
	                 -1/r[2],             0, 1/r[4]+1/r[2],       -1/r[4],
	                       0,       -1/r[3],       -1/r[4], 1/r[4]+1/r[3]),
	        ncol=4, nrow=4, byrow=TRUE);
}

SuppInvRes <- function(r) {
	
	if (all(0 < r)) {
		Chol <<- try(chol(AdmMat(truer)), silent=TRUE);
		
		is.matrix(Chol);
	}
	else
		FALSE;
}


### Sets the forward map
A <- function() {
	
	InvY <- chol2inv(Chol);
	
	cbind( (InvY %*% c(1,0,0,0))[-4,], (InvY %*% c(0,1,0,0))[-4,], (InvY %*% c(0,0,1,0))[-4,]);
}	

sigma <- 3;
##SuppInvRes(truer); Dat <- A() + matrix( rnorm( 9, sd=sigma), ncol=3, nrow=3);



UInvRes <- function(r) {
	
	0.5*sum((A() - Dat)^2)/(sigma^2); # + 0.5*sum((r - 150)^2)/(20^2);
}

InitVal <- function() {  rep( 400, 4) + rnorm( 4, sd=10); }
 
 
 
RunInvRes <- function(Tr=10000) {
	
	Runtwalk( Tr=Tr, Obj= UInvRes, Supp= SuppInvRes, xp0=InitVal(), x0=InitVal(), PlotObj=FALSE);
} 
	