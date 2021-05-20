InitGalaxy <- function() {
M<-m<-5;
n<-(m*3)-1+1;
alph<-2;
ikap<-0.0016;
xi<-21.73;
g<-0.2;
h<-0.016;
Delta<-rep(1,m);
Beta<-100;

x0<-c(rep(0.1,2),0.125,0.125,seq(from=0.1,to=48,length=5),rep(1.1,5),1)
xp0<-c(0.2,0.1,rep(0.09,2),seq(from=-0.01,to=62,length=5),rep(2,5),2)

}

InitGalaxy <- function() {

	M<<-m<<-5;
	n<<-(m*3)-1+1;
	alph<<-3;
	ikap<<-0.0016;
	xi<<-21.73;
	g<<-0.2;
	h<<-0.016;
	Delta<<-rep(1,m);
	Beta<<-20;
}

#### -log f(P|Delta); P~Di(Delta)
priorDIU<-function(P) {
  # sum(P)=1
  sum(log(gamma(Delta)))-log(gamma(sum(Delta)))+sum((1-Delta)*log(P))
}

priorNGGU<-function(Mu,Tau,beta=Beta,m=M) {
#m is the number of normals in the mix
#beta is a hyperparameter
#ikap = kapa^-1

  #-m*0.5*log(ikap/(2*pi)) -m*(alph+g-1)*log((beta)/(gamma(alph)))-m*log((h^g)/(gamma(g)))-
  #(alph-1)*sum(log(Tau)) + ikap*0.5*sum((Mu-xi)^2) + beta*sum(Tau+h)

  -m*0.5*log(ikap/(2*pi)) -m*(alph-1)*log((beta)/(gamma(alph)))-
  (alph-1)*sum(log(Tau)) + ikap*0.5*sum((Mu-xi)^2) + beta*sum(Tau)

}

priorU <- function(P, Mu, Tau, beta=Beta) {
  priorDIU(P) + priorNGGU(Mu, Tau, beta)
}

#### -log Likelihood(M5); M5<--five normals mix
#    sum^n {-log(sum^3{p_j*phi(y_i|mu_j,tao_j)})}
logLikhoodGalaxy<-function(P,Mu,Tau) {

  summ<-0;
  for(i in 1:length(y)) {
    summ = summ -log(sum(P*((Tau/(2*pi))^0.5)*exp( -(Tau/2)*(y[i]-Mu)^2)));
  }
  summ; 
}

### Ahora pars ya incluye a beta.
galaxyU<-function(pars,m=M) {
  
  P  <- pars[1:(m-1)];
  P2 <- c(P,1-sum(P));
  Mu <- pars[m:(2*m-1)];
  Tau<- pars[(2*m):(3*m-1)];
  beta<-pars[3*m];

# priorU(P2, Mu, Tau) + logLikhoodGalaxy(P2, Mu, Tau); 
  priorU(P2, Mu, Tau, beta) + logLikhoodGalaxy(P2, Mu, Tau);
}

suppGalaxy <- function(pars,m=M) {

  P  <- pars[1:(m-1)];
  Mu <- pars[(m):(2*m-1)];
  Tao<- pars[(2*m):(3*m-1)];
  beta<-pars[3*m];
  
  if( all(P>0) & all(Tao>0) & (sum(P)<1) & (!is.unsorted(Mu)) & (beta>0)) {
    TRUE;
  }
  else
    FALSE;
  
}


### Ploting functions
generaMuestra<-function(pars,n, m) {
	vect<-matrix(nrow=n ,ncol=m);
	u<-runif(n);
	p<-acumula(pars,m-1);

	for(i in 1:n){	

		v<-rep(0,m);
		if(u[i]<p[1]) {
			v[1]<-1;
			vect[i,]<-v;
		}

		for(j in 2:(m-1)) {
			if( (u[i]>p[j-1]) & (u[i]<p[j]) ) {
				v[j]<-1;
				vect[i,]<-v;
			}	
		}

		if(u[i]>p[m-1]) {
			v[m]<-1;
			vect[i,]<-v;
		}	
	}

	phi<-matrix(nrow=n,ncol=m);
	for(i in 1:m) {
		phi[,i]<-rnorm(n,mean=pars[i+m-1],sd=1/sqrt(pars[i+2*m-1]));
	}
	aux<-(phi*vect);
	
	y<-c(1:n);
	for(i in 1:n){
		y[i]<-sum(aux[i,])
	}
	
	y
}

plotDensity<-function(pars,m,space,a=1.2, add=FALSE, ...) {
  P  <- pars[1:(m-1)];
  P2 <- c(P,1-sum(P));
  Mu <- pars[m:(2*m-1)];
  Tao<- pars[(2*m):(3*m-1)];

  mTao<-max(1/sqrt(Tao));	
  lb<-min(Mu)-(3*mTao);
  ub<-max(Mu)+ (3*mTao);
  interval<-seq(lb,ub,space);
  yy<-interval;	

  for(i in 1:length(interval)){
	 yy[i]<-sum(P2*dnorm(interval[i],mean=Mu,sd=1/sqrt(Tao)));   
  } 

  if (add)
  	lines(interval,a*yy)
  else
  	plot(interval,a*yy,type='l', ...)
  #y;
}

evalDensity<-function(pars,m,x) {
  P  <- pars[1:(m-1)];
  P2 <- c(P,1-sum(P));
  Mu <- pars[m:(2*m-1)];
  Tao<- pars[(2*m):(3*m-1)];

  y<-x;	

  for(i in 1:length(x)){
	 y[i]<-sum(P2*dnorm(x[i],mean=Mu,sd=1/sqrt(Tao)));   
  } 
  y
}

#pars= a set of parameter vectors, ordered by row. i.e. one row is a parameter vector
#  nd= the number of densities to be ploted
#   m= the number of gaussians in the mix
#   x= a vector, where the densities are going to be evaluated
plotSevDens<-function(pars,nd,m,x) {
	n<-length(x);
	y<-matrix(ncol=nd,nrow=n);

	for(i in 1:nd) {
		y[,i]<-evalDensity(pars[i,],m,x);
	}
	
	matplot(x, y, pch = 1:nd, type = "l",  col = rainbow(nd), ylim=c(0,0.27),xlim=c(7,35),ylab="Density",xlab="Velocities")
}


Example4<-function(Tr=5000) {
	InitGalaxy();

	Tpars<-c(0.15,0.05,0.45,0.3,  #probabilidades
        	9, 16,19,23,33.5,     #medias
        	1/c(2,5,2,2.2,4), 1)   #varianzas, el ultimo es beta, no juega 

	x0 <- c(0.10,0.10,0.40,0.02,  #probabilidades
        	9.5, 14, 20, 25, 32.5,     #medias
        	1/c(2.1, 5.4, 2.1, 2.5, 4.1), 5)   #varianzas, el ultimo es beta, no juega 
	xp0 <- c(0.15,0.45,0.05,0.3,  #probabilidades
        	9, 16,19,23,33.5,     #medias
        	1/c(2,5,2,2.2,4), 1) 
		
	#x0<- Tpars + rnorm(n,mean=0,sd=0.02);
	#xp0<-Tpars + rnorm(n,mean=0,sd=0.02);
	
	#yy<<-read.table("galaxy.dat")
	#y<<-vector(length=82);
	#for(i in 1:82) {
	#	y[i]<<-yy[i,];
	#}
	
	y <<- scan("galaxy.dat");

	info<-Runtwalk(Tr=Tr,Obj=galaxyU,PlotObj=y,Supp=suppGalaxy,x0=x0,xp0=xp0)
	plot(-info$Us,type='l')
	for (i in 1:1000) exp(i); ##delay a bit
	itmap <- which(info$Us==min(info$Us))[1];
	map<-info$output[itmap,];
	plotDensity(map,m,0.1,a=1, ylim=c(0,0.27),xlim=c(7,35), ylab="Density",xlab="Velocities",col="red")
	hist(y,breaks=40,freq=FALSE,add=TRUE,lwd=30)
	
	info;

}

PlotMAP <- function(info)
{
	map<-info$output[info$Us==min(info$Us),][1,];
	plotDensity(map,m,0.1,a=1, ylim=c(0,0.27),xlim=c(7,35), ylab="Density",xlab="Velocities",col="red")
	hist(y,breaks=40,freq=FALSE,add=TRUE,lwd=30);
	lines( c(0,0), c(40,0));
}	

DiDens <- function(info)
{
	hist(y,breaks=40,freq=FALSE, ylim=c(0,0.27),xlim=c(7,35), ylab="Density",xlab="Velocities");
	
	for (tt in (100*(0:(info$Tr/100))+1))
	{
		hist(y,breaks=40,freq=FALSE, ylim=c(0,0.27),xlim=c(7,35), ylab="Density",xlab="Velocities");
		plotDensity( info$output[tt,],m,0.1,a=1, add=TRUE);
	}
}