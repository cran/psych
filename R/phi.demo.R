"phi.demo" <-
function(n=1000,r=.6 ,cuts =c(-2,-1,0,1,2)  ) {    
 #simulate correlation matrix with variable cut points -- psychometric demo
 #make up some correlations with different cut points
latent <-rnorm(n)  #make up some random normal theta scores
err<- rnorm(n)    #random normal error scores

observed <- latent*(r) + err*sqrt(1-r*r)  #observed = T + E

#convert to 0/1  with different values of cuts
trunc<- matrix(rep(observed,length(cuts)),ncol=length(cuts))  
for (i in 1:length(cuts)) {
   trunc[observed>cuts[i],i]<- 1
   trunc[observed< cuts[i],i]<- 0}
   
d.mat<- data.frame(latent,observed,trunc)  #combine into a data frame
pairs.panels(d.mat,main="Phi coefficients for extreme cut point normal data") 
trunc.cor<- cor(d.mat)                 #find the Pearson correlations
freq <- apply(d.mat,2,mean)                    #find the frequencies of scores


 #first demonstrate the tetrachoric function

 tetra <- tetrachoric(d.mat[,3:7],correct=FALSE)

for (i in 4:length(d.mat)) {
   for (j in 3:i) {
       trunc.cor[j,i]<- phi2tetra(trunc.cor[i,j],freq[i],freq[j]) 
       }}
  return(list(tetrachoric=tetra, phis=trunc.cor))
}

