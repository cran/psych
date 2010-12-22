"comorbidity" <- 
function(d1,d2,com,labels=NULL) {
 cl <- match.call()
twobytwo <- matrix(c(com, d1-com,d2-com,1-d1-d2+com),ncol=2)
if(is.null(labels)) {
colnames(twobytwo) <- c("D1","d1")
rownames(twobytwo) <- c("D2","d2")} else { colnames(twobytwo) <- c(labels[1],paste("-",labels[1],sep=""))
                                          rownames(twobytwo) <- c(labels[2],paste("-",labels[2],sep=""))}
phi <- phi(twobytwo)
Yule <- Yule(twobytwo)
tetra<- tetrachoric(twobytwo)
answer <- list(twobytwo=twobytwo,phi=phi,Yule=Yule,tetra=tetra,Call=cl)
class(answer) <- c("psych","comorbid")
return(answer)
}