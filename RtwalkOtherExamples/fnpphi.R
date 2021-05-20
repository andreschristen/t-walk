


#### Function to calculate the probability pphi of selecting a parameter to move
#### pphi <- fnpphi(n)/n;

n1phi <- 5; ### maximum expected number of 1's, fnpphi(n) <= 5.

### With laphi <- laphibound(2) (=0.2876821) we have fnpphi(2)=2
laphibound <- function(x) { (-1/(x-1))*log((n1phi-x)/(n1phi-1)) };
laphi <- laphibound(2);
)
#### fnpphi(n) rapidly tends to 5 (n1phi).  
fnpphi <- function(n) {  (n1phi-(n1phi-1)*exp(-laphi*(n-1))); }


#### The behaivour of fnpphi and the resulting pphi (=fnpphi(n)/n) may be
#### seen with this plot.
Plotnpphi <- function() {
	curve( fnpphi(x), xlim=c(1,n1phi+10), ylim=c(0,n1phi), axes=FALSE,
	ylab="Expected num of selected j's", xlab="n");
	
	curve( n1phi*fnpphi(x)/x, xlim=c(1,n1phi+10), col="red", add=TRUE);
	
	axis( 1, floor( seq(from=1, to=(n1phi+10), length=5)));
	axis( 2, 0:(n1phi));
	axis( 4, seq(from=0, to=n1phi, length=5),
		paste(seq( from=0, to=1, by=0.25)), col="red");
}

### Runtwalk runs this line to set pphi: pphi <- fnpphi(n)/n;

