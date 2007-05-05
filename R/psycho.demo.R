"psycho.demo" <-
function() {    
 #simulate correlation matrix with variable cut points -- psychometric demo
 #make up some correlations with different cut points
cuts <-c(-2,-1,0,1,2)
nsub<-1000    #how many subjects
r<-.6         #population correlation of observed with theta
theta <-rnorm(nsub)  #make up some random normal theta scores
err<- rnorm(nsub)    #random normal error scores

obser<- theta*(r) + err*sqrt(1-r*r)  #observed = T + E

#convert to 0/1  with different values of cuts
trunc<- matrix(rep(obser,length(cuts)),ncol=length(cuts))  
for (i in 1:length(cuts)) {
   trunc[obser>cuts[i],i]<- 1
   trunc[obser< cuts[i],i]<- 0}
   
d.mat<- data.frame(theta,obser,trunc)  #combine into a data frame
trunc.cor<- cor(d.mat)                 #find the correlations
freq <- mean(d.mat)                    #find the frequencies of scores

#now, convert the upper diagonal to polychorics using John Fox's polychor and my phi2poly

for (i in 4:length(d.mat)) {
   for (j in 3:i) {
       trunc.cor[j,i]<- phi2poly(trunc.cor[i,j],freq[i],freq[j]) 
       }}

}

