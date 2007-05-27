"fa.parallel" <-
function(x,ncases=0,main="Parallel Analysis Scree Plots")  {     
 nsub <- dim(x)[1]
 nvariables <- dim(x)[2]
 if (ncases >0) { nsub <- ncases 
  rx <- x} else {
  rx <- cor(x,use="pairwise")
  }
 if((ncases ==0) && (nsub==nvariables)) {warning("It seems as if you are using a correlation matrix, but have not specified the number of cases")}
  if(ncases == 0) {sampledata <- matrix(sample(unlist(x),size=nsub*nvariables,replace=TRUE),nrow=nsub,ncol=nvariables)
  values.samp <- eigen(cor(sampledata,use="pairwise"))$values}
  simdata=matrix(rnorm(nsub*nvariables),nrow=nsub,ncol=nvariables)    #make up simulated data
  values.sim <- eigen(cor(simdata))$values

 valuesx  <- eigen(rx)$values 

plot(valuesx,type="b", main = main ) 
points(values.sim,type ="b",lty="dotted")
if(ncases ==0) {points(values.samp,type ="b",lty="dashed")}
legend("topright", c("   Actual Data", "   Simulated Data", "  Resampled Data"), col = c(3,4,5),
       text.col = "green4", lty = c("solid","dotted", "dashed"),
       merge = TRUE, bg = 'gray90')
abline(h=1) } 

 