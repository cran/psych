"test.psych" <- 
function(first=1,last=5,short=TRUE) {  
s1 <- USArrests         #  Violent Crime Rates by US State  (4 variables)
s2 <- attitude          #The Chatterjee-Price Attitude Data
s3 <- Harman23.cor$cov     #   Harman Example 2.3 8 physical measurements
s4 <- Harman74.cor$cov     #   Harman Example 7.4  24 mental measurements
s5 <- ability.cov$cov       #  8 Ability and Intelligence Tests 

#convert covariance to correlation
d5 <- diag(1/sqrt(diag(s5)))
s5 <- d5 %*% s5 %*% d5

datasets <- list(s1,s2,s3,s4,s5)
out <- list()

for (i in first:last) {
   test.data <- datasets[[i]]
	pc <-   principal(test.data)
	pc2 <-    principal(test.data,2)
	fa2 <- factor.pa(test.data,2)
	fp <-    fa.parallel(test.data)
	ic <-   ICLUST(test.data)
	if(require(GPArotation)) {om <-  omega(test.data)} else {warning("Omega requires the GPArotation package to be loaded")
	  om <- NULL}
	fc <- factor.congruence(pc2,fa2)
	vss2 <- VSS(test.data)
	vsspc <- VSS(test.data,pc="pc")
	vss.scree <- VSS.scree(test.data)
	d <- describe(test.data)
	
	keys <- matrix(rep(0,dim(test.data)[2]*2),ncol=2)
	keys[,1] <- 1
	keys[1:3,2] <- 1
	if( dim(test.data)[1] != dim(test.data)[2]) {test.score <- score.items(keys,test.data)} else {test.score <- cluster.cor(keys,test.data)}
	
	out <- list(out,paste("test",i),pc,pc2,fa2,fp,ic,om,fc,vss2,vsspc,vss.scree,d,test.score)
  } #end loop
    #a few more tests
  
  simple <- item.sim(nvar=24)
  circ <-  circ.sim(nvar=24)
  simple.par <- fa.parallel(simple)
  fa.simple <- factor.pa(simple,2)
  f4 <- VSS.simulate()
 if (require(polycor)) { psych.d <- psycho.demo() } else {warning("psycho.demo requires the polycor package")
         psych.d <- NULL}
  cong <- congeneric.sim()
  pairs.panels(cong)
 if(require(Rgraphviz)) { fa.graph(factor.pa(item.sim(16),2)) } else {warning("fa.graph requires Rgraphviz") }
  out <- list(out,fa.simple,psych.d)
 if (!short) { return(out)}

}#end function
