"VSS.plot" <-
function(x,plottitle="Very Simple Structure",line=FALSE)
  {
  	n=dim(x)
  	symb=c(49,50,51,52)              #plotting sym
	plot(x$cfit.1,ylim=c(0,1),type="b",ylab="Very Simple Structure Fit",xlab="Number of Factors",pch=49)
    if (line) lines(x$fit)
	 

	title(main=plottitle)
	x$cfit.2[1]<-NA
	x$cfit.3[1]<-NA
	x$cfit.3[2]<-NA
	x$cfit.4[1]<-NA
	x$cfit.4[2]<-NA
	x$cfit.4[3]<-NA
	lines(x$cfit.2)
	points(x$cfit.2,pch=50)
	lines(x$cfit.3)
	points(x$cfit.3,pch=symb[3])
	lines(x$cfit.4)
	points(x$cfit.4,pch=symb[4])
}

