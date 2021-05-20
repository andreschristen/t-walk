

##### Generic random walk


RunRW <- function(Tr=5000, sigma=c(1,1),
	Obj=RaU, Supp=SuppRa, PlotObj=PlotObjRa, pathcol="grey", add=FALSE, x0=c(13.5,-13.5))
{
	 x <- x0;
	 U <- Obj(x);
	 
	if (Supp(x) && Supp(xp))
	{	
		acc <- 0;
		rec <- matrix( 0, ncol=1+n, nrow=Tr+1);  ##save U and x 
		rec[1,] <- append( U, x);
		
		if (is.function(PlotObj))
		{
			PlotObj(add=add);
			Plpoints <- TRUE;
		}
		else
			Plpoints <- FALSE;
		propU <- U;

		y <- x;
	}
	else
	{
		cat(paste("Initial values out of support,\n  x=", x,"\n xp=", xp));
		Tr <- 0;
	}
	
	
	for (it in 1:Tr)
	{
		j <- ceiling(n*runif(1));  ### Choose a direction

		y <- x;
		y[j] <- y[j] + sigma[j]*rnorm(1);
		
		if (Supp(yp))
		{
			propU <- Obj(y);
		
			A <- exp(U - propU);
		}
		else
			A <- 0; ##out of support, not accepted
		
		if (runif(1) < A)
		{ ## accepted
			if (Plpoints)
			{
				#points( x[1], x[2], pch=20, col="blue");
				#points(xp[1],xp[2], pch=20, col="red");
			
				lines( c(x[1],y[1]), c(x[2],y[2]), col=pathcol);
				#lines( c(xp[1],yp[1]), c(xp[2],yp[2]), col="blue");
			}
			
			acc <-  acc + 1;
			
			x <- y;
			U <- propU;
			
		}
		
		rec[it+1,] <- append( U, x);
	}
	
		
	if (is.function(PlotObj))
		PlotObj(add=TRUE);
	
	list( n=n, Tr=Tr, acc=acc, Us=rec[,1], output=rec[,1+(1:n)]);
}	

	