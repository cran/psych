#Quality control function to run through hard problems 
"test.psych" <- 
function(first=1,last=5,short=TRUE,all=FALSE,fapc=FALSE) {  
s1 <- datasets::USArrests         #  Violent Crime Rates by US State  (4 variables)
s2 <- datasets::attitude          #The Chatterjee-Price Attitude Data
s3 <- datasets::Harman23.cor$cov     #   Harman Example 2.3 8 physical measurements
s4 <- datasets::Harman74.cor$cov     #   Harman Example 7.4  24 mental measurements
s5 <- datasets::ability.cov$cov      #  6 Ability and Intelligence Tests 


#convert covariance to correlation
d5 <- diag(1/sqrt(diag(s5)))
s5 <- d5 %*% s5 %*% d5

datasets <- list(s1,s2,s3,s4,s5)
out <- list()

for (i in first:last) {
   test.data <- datasets[[i]]
	pc <-   principal(test.data)
	pc2 <-    principal(test.data,2)
	if(i < 3) {
			fa2 <- fa(test.data,2)
			fp <-    fa.parallel(test.data)
			vss2 <- VSS(test.data)
			vsspc <- VSS(test.data,fm="pc")
		} else {
			fa2 <- fa(test.data,2,n.obs=200)
			cluster.plot(fa2)
			fp <-    fa.parallel(test.data,n.obs=200)
			vss2 <- VSS(test.data,n.obs=200)
			vsspc <- VSS(test.data,fm="pc",n.obs=200)
			

			}


	ic <-   ICLUST(test.data,plot=FALSE)
	if(requireNamespace('GPArotation')) {om <-  omega(test.data,plot=FALSE)} else {warning("Omega requires the GPArotation package to be loaded")
	  om <- NULL}
	fc <- factor.congruence(pc2,fa2)
	
	d <- describe(test.data)
	
	keys <- matrix(rep(0,dim(test.data)[2]*2),ncol=2)
	keys[,1] <- 1
	keys[1:3,2] <- 1
	if( dim(test.data)[1] != dim(test.data)[2]) {test.score <- scoreItems(keys,test.data)} else {test.score <- cluster.cor(keys,test.data)}
	
	out <- list(out,paste("test",i),pc,pc2,fa2,fp,ic,om,fc,vss2,vsspc,d,test.score)
  } #end loop
    #a few more tests
  
  cat("\n Testing cor plot and sim.item\n")
   set.seed(42)   #this way our simulaton will be consistent
  simple <- sim.item(nvar=24)
  circ <-  sim.circ(nvar=24)
  cor.plot(cor(circ),colors=TRUE,zlim=c(-1,1),main="24 variables in a circumplex")
  simple.par <- fa.parallel(simple)
  fa.simple <- fa(simple,2)
  cor.plot(fa.simple,TRUE,n=4)
  #fa.simple.keys <- ICLUST.sort(fa.simple,keys=TRUE) #why this way
 # simple.scores <-  scoreItems(fa.simple.keys$clusters,simple)
  fa.simple.keys <- factor2cluster(fa.simple)
  simple.scores <-  scoreItems(fa.simple.keys,simple)
 
 pairs.panels(simple.scores$scores)
 

 cat("\n Test of sim.VSS\n")
   f4 <- sim.VSS()  
    psych.d <- NULL 
  #the next test, phi.demo, throws multiple warnings that are from the polycor package and can not be found
  #if (!require(polycor)) { warning("psycho.demo requires the polycor package")  psych.d <- NULL  } else  {psych.d <- phi.demo() } 
  cong <- sim.congeneric()
  
if(all) { #test of factoring and scoring singular data  -- fails on some platforms
cat("\n Test of a singular matrix\n")
#a test example of a singular matrix
IRIS <- datasets::iris[,1:4]
IRIS[,5] <- datasets::iris[,1]+datasets::iris[,2]
f.iris <- fa(IRIS,5,scores=TRUE) #this is get around the failure of tenBerge for a singular matrix
p.iris <- principal(IRIS,5,scores=TRUE)
#this will fail if not using minres or pa
  
  }

   
cluster.plot(fa(sim.circ(nvar=24),nfactors=2),title="two circumplex factors")
 pairs.panels(cong) 
 #this section tests various functions that use Rgraphviz (if it is installed) 
 if(FALSE) { #{if(require(Rgraphviz) && !FALSE) {      
  fa.graph(fa(item.sim(16),2) ,title="Principal factor of a simple structure") 
  	ic.out <- ICLUST(s4,title="ICLUST of 24 Mental abilities")
  	v9 <-  omega(sim.hierarchical(),title="Omega with Schmid Leihman")
  	omega.graph(v9,sl=FALSE,title="Omega with hierarchical factors")
  	   #set up the parameters for the structure graph
		X6 <- matrix(c("a","b","c",rep(0,6),"d","e","f"),nrow=6)
		colnames(X6) <- c("L1","L2")
		rownames(X6) <- c("x1","x2","x3","x4","x5","x6")
		Y3 <- matrix(c("u","w","z"),ncol=1)
		colnames(Y3) <- "Y"
		rownames(Y3) <- c("y1","y2","y3")
		phi21 <- matrix(c(1,0,"r1",0,1,"r2",0,0,1),ncol=3)
		colnames(phi21) <- rownames(phi21) <-  c("L1","L2","Y")
	structure.graph(X6,phi21,Y3,title="Symbolic structural model")

 } else {warning("fa.graph, omega.graph, ICLUST.rgraph, structure.graph require Rgraphviz and were not tested") }
 
 fa.diagram(fa(item.sim(16),nfactors=2)) 
  	ic.out <- ICLUST(s4,title="ICLUST of 24 Mental abilities")
  	v9 <-  omega(sim.hierarchical(),title="Omega with Schmid Leihman")
  	omega.diagram(v9,sl=FALSE,main="Omega with hierarchical factors")
  	   #set up the parameters for the structure graph
		X6 <- matrix(c("a","b","c",rep(0,6),"d","e","f"),nrow=6)
		colnames(X6) <- c("L1","L2")
		rownames(X6) <- c("x1","x2","x3","x4","x5","x6")
		Y3 <- matrix(c("u","w","z"),ncol=1)
		colnames(Y3) <- "Y"
		rownames(Y3) <- c("y1","y2","y3")
		phi21 <- matrix(c(1,0,"r1",0,1,"r2",0,0,1),ncol=3)
		colnames(phi21) <- rownames(phi21) <-  c("L1","L2","Y")
	 example.model <- structure.diagram(X6,phi21,Y3,main="Symbolic structural model")

   R <- cor(sim.item(16))
    ss <- c(1,3,5,7,9,11,13,15)
   f <- fa(R[ss,ss],2)
   foe <- fa.extension(R[ss,-ss],f)
   fa.diagram(fa.results=f,fe.results=foe)
   
   #now test the iteration options  and the rotation options in fa
#not run by default for official testing
if(fapc) {
   data1 <- psych::bfi
 
   f3 <- fa(data1[1:15],3,n.iter=5)
   f3 <- fa(data1[1:15],3,n.iter=5,rotate="Varimax")
   f3 <- fa(data1[1:15],3,n.iter=5,rotate="varimax")
   f3 <- fa(data1[1:15],3,n.iter=5,rotate="quartimax")
   f3 <- fa(data1[1:15],3,n.iter=5,rotate="bentlerT")
   f3 <- fa(data1[1:15],3,n.iter=5,rotate="geominT")
    Target <- matrix(c(rep(1,5),rep(0,15),rep(1,5),rep(0,15),rep(1,5)),ncol=3)
   f3 <- fa(data1[1:15],3,n.iter=5,rotate="targetT",Target=Target)
   f3 <- fa(data1[1:15],3,n.iter=5,rotate="bifactor")
   f3 <- fa(data1[1:15],3,n.iter=5,rotate="TargetT",Target=Target)
   f3 <- fa(data1[1:15],3,n.iter=5,rotate="equamax")
   f3 <- fa(data1[1:15],3,n.iter=5,rotate="varimin")
   #f3 <- fa(data1[1:15],3,n.iter=5,rotate="specialQ")
   f3 <- fa(data1[1:15],3,n.iter=5,rotate="Promax")
   f3 <- fa(data1[1:15],3,n.iter=5,rotate="promax")
   f3 <- fa(data1[1:15],3,n.iter=5,rotate="cluster")
   f3 <- fa(data1[1:15],3,n.iter=5,rotate="biquartimin")
    Targ <- make.keys(15,list(f1=1:5,f2=6:10,f3=11:15)) 
 	Targ <- scrub(Targ,isvalue=1)  #fix the 0s, allow the NAs to be estimated
	Targ <- list(Targ)  #input must be a list
   f3 <- fa(data1[1:15],3,n.iter=5,rotate="TargetQ",Target=Targ)  #Michael Brown's 
   #f3 <- fa(data1[1:15],3,n.iter=5,rotate="specialQ")
   f3 <- fa(data1[1:15],3,n.iter=5,rotate="oblimin")
   f3 <- fa(data1[1:15],3,n.iter=5,rotate="quartimin")
   f3 <- fa(data1[1:15],3,n.iter=5,rotate="simplimax")
   f3 <- fa(data1[1:15],3,n.iter=5,rotate="geominQ")
   
   f3 <- fa(data1[1:15],3,n.iter=5,rotate="targetQ",Target=Target)
   f3 <- fa(data1[1:15],3,n.iter=5,rotate="bentlerQ")
   
     data2 <- psych::ability
    f1 <- fa(data2)
    fpoly <- fa(data2[1:10],2,n.iter=5,cor="poly")
    f1 <- fa(data2,n.iter=4)
    f1p <- fa(data2,n.iter=4,cor="tet")
    
   p3 <- principal(data1[1:15],3,n.iter=1)
   p3 <- principal(data1[1:15],3,rotate="Varimax")
   p3 <- principal(data1[1:15],3,rotate="varimax")
   p3 <- principal(data1[1:15],3,rotate="quartimax")
   p3 <- principal(data1[1:15],3,rotate="bentlerT")
    p3 <- principal(data1[1:15],3,rotate="geominT") 
  Target <- matrix(c(rep(1,5),rep(0,15),rep(1,5),rep(0,15),rep(1,5)),ncol=3)
    p3 <- principal(data1[1:15],3,rotate="targetT",Target=Target)
       p3 <- principal(data1[1:15],3,rotate="TargetT",Target=Target)
   p3 <- principal(data1[1:15],3,rotate="bifactor")
   p3 <- principal(data1[1:15],3,rotate="varimin")
   p3 <- principal(data1[1:15],3,rotate="bentlerT")
   p3 <- principal(data1[1:15],3,rotate="geominT")
   p3 <- principal(data1[1:15],3,rotate="equamax")
   p3 <- principal(data1[1:15],3,rotate="Promax")
      p3 <- principal(data1[1:15],3,rotate="promax")
   p3  <- principal(data1[1:15],3,rotate="cluster")
   p3  <- principal(data1[1:15],3,rotate="biquartimin")
   p3  <- principal(data1[1:15],3,rotate="equamax")
   Targ <- make.keys(15,list(f1=1:5,f2=6:10,f3=11:15)) 
 	Targ <- scrub(Targ,isvalue=1)  #fix the 0s, allow the NAs to be estimated
	Targ <- list(Targ)  #input must be a list
	p3  <- principal(data1[1:15],3,rotate="TargetQ",Target=Targ)
	 p3  <- principal(data1[1:15],3,rotate="oblimin")
   p3  <- principal(data1[1:15],3,rotate="quartimin")
   p3  <- principal(data1[1:15],3,rotate="simplimax")
    p3  <- principal(data1[1:15],3,rotate="geominQ")
   p3  <- principal(data1[1:15],3,rotate="biquartimin")
   p3  <- principal(data1[1:15],3,rotate="targetQ",Target=Target)
      p3  <- principal(data1[1:15],3,rotate="bentlerQ")
   
	
	
    fpoly <- principal(data1[1:10],2,cor="poly")
    
    f1 <- principal(data2)
    f1p <- principal(data2,cor="tet")
}
 


    out <- list(out,fa.simple,psych.d)
 if (!short) { return(out)}

}#end function
